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

using namespace std;
using namespace bagel;

// construct Multimer class and perform rhf
Multimer::Multimer(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom) : geom_(geom) {
  
  // SCF
  auto HFinfo = input->get_child("hf") ? input->get_child("hf") : make_shared<PTree>();
  auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", HFinfo, geom_, rhf_ref_));
  rhf->compute();
  rhf_ref_ = rhf->conv_to_ref();
#if 1  
  MoldenOut out("out.molden");
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


void Multimer::set_active(shared_ptr<const PTree> idata) {
  
  const int multimerbasis = geom_->nbasis();
  const int nclosed_HF = rhf_ref_->nclosed();
  
  auto isp = idata->get_child("multimer_active");
  set<int> ActList;
  for (auto& s : *isp) { ActList.insert(lexical_cast<int>(s->data())-1); };

  int nactive = ActList.size();
  int nclosed = nclosed_HF;
  int nvirt = multimerbasis - nclosed;
  for (auto& amo : ActList) {
    if (amo < nclosed_HF) --nclosed;
    else --nvirt;
  }

  auto reorder_coeff = make_shared<Matrix>(multimerbasis, multimerbasis);
  { 
    auto coeff = rhf_ref_->coeff();

    int iclosed = 0;
    int iactive = nclosed;
    int ivirt = nclosed + nactive;

    auto cp = [&coeff, &reorder_coeff, &multimerbasis] (const int i, int& pos) { copy_n(coeff->element_ptr(0, i), multimerbasis, reorder_coeff->element_ptr(0, pos)); ++pos; };
  
    for (int i = 0; i != multimerbasis; ++i) {
      if (ActList.find(i) != ActList.end()) cp(i, iactive);
      else if (i < nclosed_HF) cp(i, iclosed);
      else cp(i, ivirt);
    }
  }

  active_ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*reorder_coeff)), nclosed, nactive, nvirt);
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
  auto tildex = make_shared<Matrix>(*new_coeff % S * *new_coeff);
  tildex->inverse_half();
  new_coeff = make_shared<Matrix>(*new_coeff * *tildex);

  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*new_coeff)), nclosed, nact, nvirt);

#if 0
  MoldenOut out("MoldenOut.molden");
  out << geom_;
  out << ref_;
#endif
}


