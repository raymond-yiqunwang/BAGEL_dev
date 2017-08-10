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

ASD_DMRG_Orbopt::ASD_DMRG_Orbopt(shared_ptr<const PTree> idata, shared_ptr<const Reference> iref) {
  // ASD-DMRG orbital optimization comes after RHF calculation
  
  // first construct multisite
  multisite_ = make_shared<MultiSite>(idata, iref);
  ref_ = multisite_->ref();

  common_init();

}

void ASD_DMRG_Orbopt::common_init() {

  print_header();

  nclosed_ = ref_->nclosed();
  nact_ = ref_->nact();
  nocc_ = nclosed_ + nact_;
  cout << " nclosed = " << nclosed_ << ", nact = " << nact_ << endl;

}

void ASD_DMRG_Orbopt::print_header() const {
  cout << "  --------------------------------------------------" << endl;
  cout << "     ASD-DMRG Second Order Orbital Optimization     " << endl;
  cout << "  --------------------------------------------------" << endl << endl;
}

void ASD_DMRG_Orbopt::compute() {

}


