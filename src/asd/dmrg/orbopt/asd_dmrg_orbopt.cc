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

using namespace std;
using namespace bagel;

ASD_DMRG_Orbopt::ASD_DMRG_Orbopt(shared_ptr<const PTree> idata, shared_ptr<const Reference> iref) : input_(idata) {
  
  print_header();
  
  // first construct multisite
  asd_info_ = idata->get_child("asd_info");
  multisite_ = make_shared<MultiSite>(asd_info_, iref);
  ref_ = multisite_->ref();
  coeff_ = ref_->coeff();
  
  common_init();

}

void ASD_DMRG_Orbopt::common_init() {

  nclosed_ = ref_->nclosed();
  nact_ = ref_->nact();
  nocc_ = nclosed_ + nact_;
  nvirt_ = ref_->nvirt();
  nmo_ = ref_->coeff()->mdim();
  
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

void ASD_DMRG_Orbopt::print_header() const {
  cout << "  --------------------------------------------------" << endl;
  cout << "     ASD-DMRG Second Order Orbital Optimization     " << endl;
  cout << "  --------------------------------------------------" << endl << endl;
}

void ASD_DMRG_Orbopt::compute() {
  assert(nvirt_ && nact_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    
    // first obtain RDM from ASD_DMRG
    {
      if (iter) asd_dmrg_->update_multisite(coeff_);
      asd_dmrg_->compute();
      asd_dmrg_->compute_rdm12();
      energy_ = asd_dmrg_->energies();
    }


  
  }

}


