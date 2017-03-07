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

using namespace std;
using namespace bagel;

// construct Multimer class and perform rhf
Multimer::Multimer(shared_ptr<const PTree> input, shared_ptr<const Geometry> geom) : input_(input), geom_(geom) {
  
  // SCF
  auto HFinfo = input->get_child_optional("hf") ? input->get_child_optional("hf") : make_shared<PTree>();
  auto rhf = dynamic_pointer_cast<RHF>(construct_method("hf", HFinfo, geom_, rhf_ref_));
  rhf->compute();
  rhf_ref_ = rhf->conv_to_ref();

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
  vector<pair<int, int>> basis_bounds;
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
