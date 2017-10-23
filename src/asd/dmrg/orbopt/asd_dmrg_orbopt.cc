//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_orbopt.cc
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang
// Maintainer: Shiozaki Group
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

#include <src/asd/dmrg/orbopt/asd_dmrg_orbopt.h>
#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

ASD_DMRG_OrbOpt::ASD_DMRG_OrbOpt(shared_ptr<const PTree> idata, shared_ptr<const Reference> iref) : input_(idata) {
  
  print_header();
  
  // first construct multisite
  asd_info_ = idata->get_child("asd_info");
  multisite_ = make_shared<MultiSite>(asd_info_, iref);
  geom_ = multisite_->geom();
  ref_ = multisite_->ref();
  coeff_ = ref_->coeff();
  hcore_ = ref_->hcore();
  
  common_init();

}

void ASD_DMRG_OrbOpt::common_init() {

  nclosed_ = ref_->nclosed();
  nact_ = ref_->nact();
  nocc_ = nclosed_ + nact_;
  nvirt_ = ref_->nvirt();
  norb_ = ref_->coeff()->mdim();

  nstate_ = input_->get<int>("opt_nstate", 1);
  max_iter_ = input_->get<int>("opt_max_iter", 50);
  max_micro_iter_ = input_->get<int>("opt_max_micro_iter", 100);
  thresh_ = input_->get<double>("opt_thresh", 1.0e-8); // thresh for macro iteration
  thresh_micro_ = input_->get<double>("opt_thresh_micro", 5.0e-6); // thresh for micro iteration

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << nclosed_ << endl;
  cout << "    * nact     : " << setw(6) << nact_ << endl;
  cout << "    * nvirt    : " << setw(6) << nvirt_ << endl << endl;

  // DMRG with RHF orbitals
  asd_dmrg_ = make_shared<RASD>(asd_info_, multisite_);

  cout << "  ===== Orbital Optimization Iteration =====" << endl << endl;

}

void ASD_DMRG_OrbOpt::print_header() const {
  cout << "  --------------------------------------------------" << endl;
  cout << "     ASD-DMRG Second Order Orbital Optimization     " << endl;
  cout << "  --------------------------------------------------" << endl << endl;
}


shared_ptr<Matrix> ASD_DMRG_OrbOpt::compute_active_fock(const MatView acoeff, shared_ptr<const RDM<1>> rdm1) const {
  Matrix dkl(nact_, nact_);
  copy_n(rdm1->data(), dkl.size(), dkl.data());
  dkl.sqrt();
  dkl.scale(1.0/sqrt(2.0));
  return make_shared<Fock<1>>(geom_, hcore_->clone(), nullptr, acoeff*dkl, false, true);
}


shared_ptr<Matrix> ASD_DMRG_OrbOpt::compute_qvec(const MatView acoeff, shared_ptr<const RDM<2>> rdm2) const {
  
  auto half = geom_->df()->compute_half_transform(acoeff);

  // TODO MPI modification
  shared_ptr<const DFFullDist> full = half->apply_JJ()->compute_second_transform(coeff_->slice(nclosed_, nocc_));

  // [D|tu] = (D|xy) Gamma_{xy,tu}
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(*rdm2);

  // (r,u) = (rt|D) [D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  return make_shared<Matrix>(*coeff_ % *tmp);
}


