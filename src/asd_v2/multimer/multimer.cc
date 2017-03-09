//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multimer.cc
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/asd_v2/multimer/multimer.h>
#include <src/scf/hf/rhf.h>
#include <src/util/io/moldenout.h>
#include <src/integral/os/overlapbatch.h>
#include <src/mat1e/mixedbasis.h>

using namespace std;
using namespace bagel;

// construct Multimer class and perform rhf
Multimer::Multimer(shared_ptr<const PTree> input, shared_ptr<const Reference> ref) : prev_ref_(ref) {

  cout << " ===== Constructing Multimer Geometry ===== " << endl;
  const shared_ptr<const PTree> moldata = input->get_child("molecule");
  geom_ = make_shared<Geometry>(*ref->geom(), moldata);

  // direct rhf with larger basis
  auto HFinfo = input->get_child("hf") ? input->get_child("hf") : make_shared<PTree>();
  auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", HFinfo, geom_, nullptr));
  rhf->compute();
  rhf_ref_ = rhf->conv_to_ref();

#if 1  
  MoldenOut out("rhf_large.molden");
  out << geom_;
  out << rhf_ref_;
#endif
}


// preparation for ASD-DMRG
void Multimer::precompute(shared_ptr<const PTree> idata) {

  // reorder MO coeff to closed - active - virtual
  set_active(idata);
  // project active orbitals to fragments
  project_active(idata);
}


// pick active orbitals from reference actorbs
void Multimer::set_active(shared_ptr<const PTree> idata) {
  
  cout << "  o Forming multimer active orbitals from reference active list" << endl;

  auto isp = idata->get_child("multimer_active");
  set<int> RefActList;
  for (auto& s : *isp) { RefActList.insert(lexical_cast<int>(s->data())-1); };
  const int nactive = RefActList.size();

  // active orbitals with small basis
  auto prev_coeff = prev_ref_->coeff();
  auto prev_active = make_shared<Matrix>(prev_coeff->ndim(), nactive);
  int pos = 0;
  for (auto& iact : RefActList) {
    copy_n(prev_coeff->element_ptr(0, iact), prev_coeff->ndim(), prev_active->element_ptr(0, pos++));
  }

  // pick active orbitals with maximum overlap with small active orbitals from closed and virtual sets, respectively
  auto coeff = rhf_ref_->coeff();
  const MixedBasis<OverlapBatch> mix(geom_, prev_ref_->geom());

  const int multimerbasis = geom_->nbasis();
  const int nclosed_HF = rhf_ref_->nclosed();
  const int nclosed = nclosed_HF - nactive/2;
  const int nvirt = multimerbasis - nclosed - nactive;
  int closed_position = 0;
  int active_position = nclosed;
  int virt_position = nclosed + nactive;

  auto out_coeff = make_shared<Matrix>(multimerbasis, multimerbasis);

  // now deal with active orbitals
  // tuple information:
  // matrix --> reference active orbitals
  // pair   --> bounds to construct subspace
  // int    --> number of active orbitals to form
  // bool   --> closed(true) / virtual(false)
  vector<tuple<pair<int, int>, int, bool>> svd_info;

  svd_info.emplace_back(make_pair(0         ,    nclosed_HF), nactive/2,  true);
  svd_info.emplace_back(make_pair(nclosed_HF, multimerbasis), nactive/2, false);

  for (auto& subset : svd_info) {
    pair<int, int> bounds = get<0>(subset);
    const int nsubactive = get<1>(subset);
    const bool closed = get<2>(subset);

    auto subcoeff = coeff->slice_copy(bounds.first, bounds.second);
    auto overlap = make_shared<const Matrix>(*prev_active % mix * *subcoeff);
  
    multimap<double, int> norms;
    for (int i = 0; i != overlap->mdim(); ++i) {
      const double norm = blas::dot_product(overlap->element_ptr(0, i), overlap->ndim(), overlap->element_ptr(0, i));
      norms.emplace(norm, i);
    }
  
    set<int> subActList;
    active_thresh_ = idata->get<double>("active_thresh", 0.5);
    double max_overlap, min_overlap;
    {
      auto end = norms.rbegin(); advance(end, nsubactive);
      end = find_if(end, norms.rend(), [this] (const pair<double, int>& p) { return p.first < active_thresh_; });
      for_each(norms.rbegin(), end, [&subActList] (const pair<double, int>& p) { subActList.emplace(p.second); });
      auto mnmx = minmax_element(norms.rbegin(), end);
      tie(min_overlap, max_overlap) = make_tuple(mnmx.first->first, mnmx.second->first);
    }
    const int sublist_size = subActList.size();
    cout << "   - size of " << (closed ? "closed " : "virtual ") << "candidate space : " << sublist_size << endl;
    cout << "   - largest overlap with reference actorb : " << max_overlap << ", smallest : " << min_overlap << endl;

    if (sublist_size != nsubactive) {
      cout << "  o Performing SVD in candidate space" << endl;
      Matrix subspace(multimerbasis, sublist_size);

      int ip = 0;
      for (auto& imo : subActList)
        copy_n(rhf_ref_->coeff()->element_ptr(0, imo), multimerbasis, subspace.element_ptr(0, ip++));
  
      // by default active orbitals should be orthonormal, but here a general method is implemented
      Overlap SAO_prev(prev_ref_->geom());
      Matrix SMO_prev(*prev_active % SAO_prev * *prev_active);
      SMO_prev.inverse_symmetric();

      Matrix projector(SMO_prev * (*prev_active % mix * subspace));
      vector<double> singulars(sublist_size, 0.0);
      shared_ptr<Matrix> V, U;
      tie(ignore, V) = projector.svd(singulars.data());
      cout << "   - largest singular value : " << singulars[0] << ", smallest : " << singulars[nsubactive - 1] << endl;

      subspace = subspace ^ *V;

      // fill in closed / virtual orbitals
      for (int i = 0; i != subcoeff->mdim(); ++i) {
        if (subActList.find(i) == subActList.end())
          copy_n(subspace.element_ptr(0, i), multimerbasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }

      // fill in active orbitals
      copy_n(subspace.data(), multimerbasis * nsubactive, out_coeff->element_ptr(0, active_position));
      active_position += nsubactive;

      // non-active orbitals in subActList
      for (int i = nsubactive; i != sublist_size; ++i) {
        copy_n(subspace.element_ptr(0, i), multimerbasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }
    } else {
      for (int i = 0; i != subcoeff->mdim(); ++i) {
        if (subActList.find(i) != subActList.end()) copy_n(subcoeff->element_ptr(0, i), multimerbasis, out_coeff->element_ptr(0, active_position++));
        else copy_n(subcoeff->element_ptr(0, i), multimerbasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }
    }
  }

  active_ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}


void Multimer::project_active(shared_ptr<const PTree> idata) {

  // define fragments (regions) of the multimer
  vector<int> region_sizes = idata->get_vector<int>("region_sizes");
  const int nregions(region_sizes.size());
  vector<pair<int, int>> basis_bounds;
  {
    int nbasis = 0;
    int natoms = 0;
    for (int& rsize : region_sizes) {
      const int atombegin = natoms;
      const int basisbegin = nbasis;
      for (int iatom = atombegin; iatom != atombegin + rsize; ++iatom)
        nbasis += geom_->atoms()[iatom]->nbasis();

      natoms += rsize;
      if (basisbegin != nbasis)
        basis_bounds.emplace_back(basisbegin, nbasis);
    }
    if (natoms != count_if(geom_->atoms().begin(), geom_->atoms().end(), [](const shared_ptr<const Atom> a){return !a->dummy();}))
      throw logic_error("All atoms must be assigned to regions");
  }

  const int nclosed = active_ref_->nclosed();
  const int nact = active_ref_->nact();
  const int nvirt = active_ref_->nvirt();
  const int multimerbasis = active_ref_->geom()->nbasis();

  Overlap S(geom_);

  // project active orbitals to each fragement and do SVD
  vector<int> actsizes = idata->get_vector<int>("active_sizes");
  vector<int> fragbasis;
  int basisoffset = 0;
  int orboffset = 0;
  auto actcoeff = active_ref_->coeff()->get_submatrix(0, nclosed, multimerbasis, nact);
  auto tmp = actcoeff->clone();
  for (int i = 0; i != nregions; ++i) {
    const int nbasis = basis_bounds[i].second - basis_bounds[i].first;
    auto Sfrag = make_shared<Matrix>(nbasis, nbasis);
    auto Smix = make_shared<Matrix>(nbasis, multimerbasis);
    Sfrag->copy_block(0, 0, nbasis, nbasis, S.get_submatrix(basisoffset, basisoffset, nbasis, nbasis));
    Smix->copy_block(0, 0, nbasis, multimerbasis, S.get_submatrix(basisoffset, 0, nbasis, multimerbasis));

    auto S_inv = make_shared<Matrix>(*Sfrag);
    S_inv->inverse_symmetric();
 
    // forming projected orbitals of current fragment
    auto projected = make_shared<const Matrix>(*S_inv * *Smix * *actcoeff);

    // get rid of redundancy
    shared_ptr<Matrix> reduced;
    {
      auto CC = make_shared<Matrix>(*projected % *Sfrag *  *projected);
      VectorB eig(projected->mdim());
      CC->diagonalize(eig);
      auto P = CC->get_submatrix(0, (CC->mdim() - actsizes[i]), CC->ndim(), actsizes[i]);
      auto temp = make_shared<const Matrix>(*projected * *P);
      reduced = make_shared<Matrix>(multimerbasis, temp->mdim());
      reduced->copy_block(basisoffset, 0, nbasis, temp->mdim(), temp->data());
    }
    tmp->copy_block(0, orboffset, multimerbasis, actsizes[i], reduced->data());

    fragbasis.emplace_back(nbasis);
    basisoffset += nbasis;
    orboffset += actsizes[i];
  }
  assert(orboffset == nact);

  // normalization
  {
    auto csc = make_shared<Matrix>(*tmp % S * *tmp);
    for (int i = 0; i != nact; ++i)
      for_each(tmp->element_ptr(0, i), tmp->element_ptr(multimerbasis, i), [&i, &csc] (double& p) { p /= sqrt(*csc->element_ptr(i, i)); });
  }
  
  auto new_coeff = active_ref_->coeff()->copy();
  new_coeff->copy_block(0, nclosed, multimerbasis, nact, tmp->data());

  // lowdin orthogonalization
  cout << endl;
  cout << "    *** If linear dependency is detected, you shall try to find better initial active orbital guess ***   " << endl;
  auto tildex = make_shared<Matrix>(*new_coeff % S * *new_coeff);
  tildex->inverse_half();
  new_coeff = make_shared<Matrix>(*new_coeff * *tildex);

  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*new_coeff)), nclosed, nact, nvirt);

#if 1
  MoldenOut out("projected.molden");
  out << geom_;
  out << ref_;
#endif
}


