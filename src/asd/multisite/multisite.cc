//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/multisite/multisite.h>
#include <src/mat1e/overlap.h>
#include <src/util/muffle.h>
#include <src/scf/hf/rhf.h>
#include <src/wfn/construct_method.h>
#include <src/util/io/moldenout.h>
#include <src/mat1e/mixedbasis.h>
#include <src/integral/os/overlapbatch.h>
#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

extern "C" {
  void dgels_(const char* trans, const int* m, const int* n, const int* nrhs, double* a, const int* lda,
              double* b, const int* ldb, double* work, const int* lwork, int* info);
}

namespace {
  void dgels_(const char* trans, const int m, const int n, const int nrhs, double* a, const int lda, 
              double* b, const int ldb, double* work, const int lwork, int& info) {
    ::dgels_(trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info); }
}

MultiSite::MultiSite(shared_ptr<const PTree> input, vector<shared_ptr<const Reference>> refs) : input_(input), nsites_(refs.size()) {
  // build the super geometry
  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  for (auto& r : refs)
    geoms_.push_back(mdata ? make_shared<Geometry>(*r->geom(), mdata) : r->geom());

  auto envdata = input_->get_child_optional("environment");
  vector<shared_ptr<const Geometry>> geovec = geoms_;
  if (envdata) {
    Muffle hide_cout;
    geovec.push_back(make_shared<Geometry>(envdata));
  }

  auto sgeom = make_shared<Geometry>(geovec);

  // build combined references
  for (auto& r : refs)
    isolated_refs_.push_back(r->project_coeff(sgeom, false));
  const size_t norb = accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0ul, [] (size_t x, shared_ptr<const Reference> r) { return x + r->coeff()->mdim(); });
  Matrix coeff(sgeom->nbasis(), norb);

  // pull out all of the orbitals and put them in closed, active, virtual order
  int current = 0;
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,0),                      r->nclosed()*sgeom->nbasis(), coeff.element_ptr(0, current));
    closed_bounds_.push_back({current, current+r->nclosed()});
    current += r->nclosed();
  }
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,r->nclosed()+r->nact()), r->nvirt()*sgeom->nbasis(),   coeff.element_ptr(0, current));
    virt_bounds_.push_back({current, current+r->nvirt()});
    current += r->nvirt();
  }

  Overlap S(sgeom);
  Matrix transform(coeff % S * coeff);
  transform.inverse_half();
  auto scoeff = make_shared<Coeff>(coeff * transform);

  sref_ = make_shared<Reference>(sgeom, scoeff,
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nclosed(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nvirt(); }));
}


MultiSite::MultiSite(shared_ptr<const PTree> itree, shared_ptr<const Reference> ref) : input_(itree), prev_ref_(ref) {
  
  cout << " ===== Constructing MultiSite Geometry ===== " << endl;
  const shared_ptr<const PTree> moldata = itree->get_child("basis_info");
  geom_ = make_shared<Geometry>(*ref->geom(), moldata);

  // direct rhf with larger basis
  auto HFinfo = itree->get_child("hf_info") ? itree->get_child("hf_info") : make_shared<PTree>();
  auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", HFinfo, geom_, nullptr));
  rhf->compute();
  rhf_ref_ = rhf->conv_to_ref();
}


// preparation for ASD-DMRG
void MultiSite::precompute() {

  // reorder MO coeff to closed - active - virtual
  set_active_metal();
  
  // project active orbitals to fragments
  project_active();
  
  // canonicalize active orbitals in sub spaces
  canonicalize();
}


// pick active orbitals from reference actorbs
void MultiSite::set_active_metal() {
  
  cout << "    o Forming multisite active orbitals from reference active list" << endl << endl;

  auto isp = input_->get_child("multisite_active");
  set<int> RefActList;
  for (auto& s : *isp) { RefActList.insert(lexical_cast<int>(s->data())-1); };
  const int nactive = RefActList.size();

  // active orbitals with small basis
  auto prev_coeff = prev_ref_->coeff();
  auto prev_active = make_shared<Matrix>(prev_coeff->ndim(), nactive);
  const int prev_nclosed = prev_ref_->nclosed();
  int pos = 0; int prev_nactclo = 0;
  for (auto& iact : RefActList) {
    if (iact < prev_nclosed) { ++prev_nactclo; }
    copy_n(prev_coeff->element_ptr(0, iact), prev_coeff->ndim(), prev_active->element_ptr(0, pos++));
  }

  // pick active orbitals with maximum overlap with small active orbitals from closed and virtual sets, respectively
  auto coeff = rhf_ref_->coeff();
  const MixedBasis<OverlapBatch> mix(geom_, prev_ref_->geom());
  const int multisitebasis = geom_->nbasis();
  const int nclosed_HF = rhf_ref_->nclosed();
  const int nclosed = nclosed_HF - prev_nactclo;
  const int nvirt = multisitebasis - nclosed - nactive;

  int closed_position = 0;
  int active_position = nclosed;
  int virt_position = nclosed + nactive;

  auto out_coeff = make_shared<Matrix>(multisitebasis, multisitebasis);

  // now deal with active orbitals
  // tuple information:
  // matrix --> reference active orbitals
  // pair   --> bounds to construct subspace
  // int    --> number of active orbitals to form
  // bool   --> closed(true) / virtual(false)
  vector<tuple<pair<int, int>, int, bool>> svd_info;

  svd_info.emplace_back(make_pair(0         ,    nclosed_HF), prev_nactclo,  true);
  svd_info.emplace_back(make_pair(nclosed_HF, multisitebasis), nactive - prev_nactclo, false);

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
    active_thresh_ = input_->get<double>("active_thresh", 0.5);
    double max_overlap, min_overlap;
    {   
      auto end = norms.rbegin(); advance(end, nsubactive);
      end = find_if(end, norms.rend(), [this] (const pair<double, int>& p) { return p.first < active_thresh_; }); 
      for_each(norms.rbegin(), end, [&subActList] (const pair<double, int>& p) { subActList.emplace(p.second); }); 
      auto mnmx = minmax_element(norms.rbegin(), end);
      tie(min_overlap, max_overlap) = make_tuple(mnmx.first->first, mnmx.second->first);
    }   
    const int sublist_size = subActList.size();
    cout << "      - size of " << (closed ? "closed " : "virtual ") << "candidate space : " << sublist_size << endl;
    cout << "      - largest overlap with reference actorb : " << max_overlap << ", smallest : " << min_overlap << endl << endl;
    
    for (int i = 0; i != subcoeff->mdim(); ++i) {
      if (subActList.find(i) != subActList.end()) copy_n(subcoeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, active_position++));
      else copy_n(subcoeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
    }   

    if (sublist_size != nsubactive)
      throw runtime_error("currently SVD is disabled, try higher active_thresh");
/*
    // TODO currently disable SVD, may turn on when necessary
    if (sublist_size != nsubactive) {
      cout << "        o Performing SVD in " << (closed ? "closed" : "virtual") << " candidate space" << endl;
      Matrix subspace(multisitebasis, sublist_size);

      int ip = 0;
      for (auto& imo : subActList)
        copy_n(rhf_ref_->coeff()->element_ptr(0, imo), multisitebasis, subspace.element_ptr(0, ip++));

      // by default active orbitals should be orthonormal, but here a general method is implemented
      Overlap SAO_prev(prev_ref_->geom());
      Matrix SMO_prev_inv(*prev_active % SAO_prev * *prev_active);
      SMO_prev_inv.inverse_symmetric();

      Matrix projector(SMO_prev_inv * (*prev_active % mix * subspace));
      vector<double> singulars(sublist_size, 0.0);
      shared_ptr<Matrix> V;
      tie(ignore, V) = projector.svd(singulars.data());
      cout << "          - largest used singular value : " << singulars[0] << ", smallest : " << singulars[nsubactive - 1] << endl;
      cout << "          - largest ignored singular value : " << singulars[nsubactive] << endl << endl;

      subspace = subspace ^ *V;

      // fill in closed / virtual orbitals
      for (int i = 0; i != subcoeff->mdim(); ++i) {
        if (subActList.find(i) == subActList.end())
          copy_n(subspace.element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }   

      // fill in active orbitals
      copy_n(subspace.data(), multisitebasis * nsubactive, out_coeff->element_ptr(0, active_position));
      active_position += nsubactive;

      // non-active orbitals in subActList
      for (int i = nsubactive; i != sublist_size; ++i) {
        copy_n(subspace.element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }   
    } else {
      for (int i = 0; i != subcoeff->mdim(); ++i) {
        if (subActList.find(i) != subActList.end()) copy_n(subcoeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, active_position++));
        else copy_n(subcoeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed ? closed_position++ : virt_position++)));
      }   
    }
*/
  }

  active_ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}


void MultiSite::project_active() {

  // define fragments (regions) of the multisite
  active_electrons_ = input_->get_vector<int>("active_electrons");
  vector<int> region_sizes = input_->get_vector<int>("region_sizes");
  nsites_ = region_sizes.size();
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
  vector<int> actsizes = input_->get_vector<int>("active_sizes");
  int basisoffset = 0;
  int orboffset = 0;
  auto actcoeff = active_ref_->coeff()->get_submatrix(0, nclosed, multimerbasis, nact);
  auto tmp = actcoeff->clone();
  for (int i = 0; i != nsites_; ++i) {
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

    basisoffset += nbasis;
    orboffset += actsizes[i];
  }
  assert(orboffset == nact);
  
  // Now we want to solve min|| Act_coeff * x - tmp||, 
  // DGELS is used to obtain the transformation matrix, projected orbitals cannot have linear dependancy
  shared_ptr<Matrix> solution;
  {
    int N = actcoeff->ndim();
    int M = actcoeff->mdim();
    assert(N > M);
    assert(M == tmp->mdim());
    int lwork = 10 * N;
    unique_ptr<double[]> work(new double[lwork]);
    auto actcopy = actcoeff->copy();
    double* adata = actcopy->data();
    double* bdata = tmp->data();
    int info = 0;
    dgels_("N", N, M, M, adata, N, bdata, N, work.get(), lwork, info);
    if (info != 0) throw runtime_error("dgels failed in projecting coeff");
  
    solution = tmp->get_submatrix(0, 0, M, M);
  }

  // Lowdin orthogonalization for transformation matrix
  {
    auto tildeX = make_shared<Matrix>(*solution % *solution);
    tildeX->inverse_half();
    cout << "    *** If linear dependency is detected, you shall try to find better initial active orbital guess ***   " << endl;
    solution = make_shared<Matrix>(*solution * *tildeX);
  }

  actcoeff = make_shared<Matrix>(*actcoeff * *solution);

  auto new_coeff = active_ref_->coeff()->copy();
  new_coeff->copy_block(0, nclosed, multimerbasis, nact, actcoeff->data());

  // update multimer information
  active_sizes_ = actsizes;
  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*new_coeff)), nclosed, nact, nvirt);

}


void MultiSite::canonicalize() {

  auto out_coeff = ref_->coeff()->copy();
  
  const int multimerbasis = ref_->geom()->nbasis();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nactele = ref_->geom()->nele() - 2*nclosed;

  // projected active orbitals are partially occupied, modify density matrix accordingly
  auto density_coeff = ref_->coeff()->slice_copy(0, nclosed + nact);
  blas::scale_n(sqrt(0.5*nactele/nact), density_coeff->element_ptr(0, nclosed), multimerbasis * nact);
  auto density = density_coeff->form_density_rhf(nclosed + nact);
  auto fock = make_shared<const Fock<1>>(geom_, ref_->hcore(), density, density_coeff);

  // canonicalize active orbitals in each fragment
  {
    int offset = nclosed;
    for (int i = 0; i != active_sizes_.size(); ++i) {
      const int norb = active_sizes_[0];
      auto subspace = ref_->coeff()->slice_copy(offset, offset + norb);
      auto subfock = make_shared<Matrix>(*subspace % *fock * *subspace);
      VectorB eigs(norb);
      subfock->diagonalize(eigs);
      subspace = make_shared<Matrix>(*subspace * *subfock);
      copy_n(subspace->data(), multimerbasis * norb, out_coeff->element_ptr(0, offset));
      offset += norb;
    }
  }

  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*out_coeff)), nclosed, nact, ref_->nvirt());
  
#if 1
  MoldenOut out("canonicalized.molden");
  out << geom_;
  out << ref_;
#endif

}


shared_ptr<Reference> MultiSite::build_reference(const int site, const vector<bool> meanfield, bool metal) const {

  assert(meanfield.size()==nsites_ && site<nsites_ && site>=0);
  if (!metal) {
    vector<shared_ptr<const MatView>> closed_orbitals = {make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed()))};
    const int act_start = accumulate(active_refs_.begin(), active_refs_.begin() + site, sref_->nclosed(),
                                      [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); });
    const int act_fence = act_start + active_refs_[site]->nact();
    const MatView active_orbitals = sref_->coeff()->slice(act_start, act_fence);

    int current = sref_->nclosed();
    for (int i = 0; i < nsites_; ++i) {
      if (meanfield[i] && i!=site)
        closed_orbitals.push_back(make_shared<const MatView>(sref_->coeff()->slice(current, current+isolated_refs_[i]->nclosed()-active_refs_[i]->nclosed())));
      current += active_refs_[i]->nact();
    }

    const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [] (int x, shared_ptr<const MatView> m) { return x + m->mdim(); });
    const int nact = active_orbitals.mdim();

    auto out = make_shared<Matrix>(sref_->geom()->nbasis(), nclosed+nact);

    current = 0;
    closed_orbitals.push_back(make_shared<MatView>(active_orbitals));
    for (auto& orbitals : closed_orbitals) {
      copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
      current += orbitals->mdim();
    }

    return make_shared<Reference>(sref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);

  } else {// Raymond version

    int act_start = ref_->nclosed();
    for (int i = 0; i != site; ++i) act_start += active_sizes_[i];
    const int nact = active_sizes_[site];
    const Matrix active_orbitals = ref_->coeff()->slice(act_start, act_start + nact);

    vector<shared_ptr<const Matrix>> closed_orbitals = {make_shared<Matrix>(ref_->coeff()->slice(0, ref_->nclosed()))};
    int current = ref_->nclosed();
    for (int i = 0; i != nsites_; ++i) {
      if (meanfield[i] && i != site) {
        const int nele = active_electrons_[i];
        const int norb = active_sizes_[i];
        auto scale_coeff = ref_->coeff()->slice_copy(current, current + active_sizes_[i]);
        blas::scale_n(sqrt(0.5*nele/norb), scale_coeff->data(), scale_coeff->size());
        closed_orbitals.push_back(scale_coeff);
      }
      current += active_sizes_[i];
    }

    const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [] (int x, shared_ptr<const Matrix> m) { return x + m->mdim(); });

    auto out = make_shared<Matrix>(ref_->geom()->nbasis(), nclosed + nact);
    current = 0;
    closed_orbitals.push_back(make_shared<Matrix>(active_orbitals));
    for (auto& orbitals : closed_orbitals) {
      copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
      current += orbitals->mdim();
    }

    return make_shared<Reference>(ref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);
  }

}
