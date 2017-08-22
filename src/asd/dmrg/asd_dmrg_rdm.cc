//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_rdm.cc
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang
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

#include <iostream>

#include <src/asd/dmrg/asd_dmrg.h>
#include <src/util/muffle.h>

using namespace std;
using namespace bagel;


void ASD_DMRG::compute_rdm12() {

  cout << endl << " " << string(10, '=') << " computing RDM12 " << string(10, '=') << endl << endl; 

  // compute ProdRASCivec after convergence
  {
    Muffle hide_cout("asd_dmrg_rdm.log", false);
    const int site = nsites_/2 - 1;
    auto left_block = (site==0) ? nullptr : left_blocks_[site-1];
    auto right_block = right_blocks_[nsites_-site-2];
    shared_ptr<const Reference> ref = multisite_->build_reference(site, vector<bool>(nsites_, false), metal_);
    // prepare input information
    shared_ptr<PTree> input = prepare_sweeping_input(site);
    {  
      input->put("nclosed", ref->nclosed());
      read_restricted(input, site);
      vector<int> actvec = multisite_->active_electrons();
      const int nactele = accumulate(actvec.begin(), actvec.end(), input->get<int>("charge"));
      input->put("nactele", nactele);
    }
    shared_ptr<const DMRG_Block> environment;
    if (!left_block) environment = right_block;
    else environment = make_shared<const DMRG_Block2>(left_block, right_block);
    auto prod_ras = make_shared<ProductRASCI>(input, ref, environment);
    prod_ras->compute();
    cc_ = prod_ras->civectors();
  }

  // for nstate==1, rdm1_av_ = rdm1_->at(0)
  // Needs initialization here because we use daxpy
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(nactorb_);
    rdm2_av_ = make_shared<RDM<2>>(nactorb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }

  for (int i = 0; i != nstate_; ++i)
    compute_rdm12(i, i);

  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weights_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weights_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0, 0);
    rdm2_av_ = rdm2_->at(0, 0);
  }

}


void ASD_DMRG::compute_rdm12(const int ist, const int jst) {
  shared_ptr<ProductRASCivec> ccbra = cc_[ist];
  shared_ptr<ProductRASCivec> ccket = cc_[jst];

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_prodcivec(ccbra, ccket);

  rdm1_->emplace(ist, jst, rdm1);
  rdm2_->emplace(ist, jst, rdm2);
}

tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> 
  ASD_DMRG::compute_rdm12_from_prodcivec(shared_ptr<const ProductRASCivec> cbra, shared_ptr<const ProductRASCivec> cket) const {
  auto rdm1 = make_shared<RDM<1>>(nactorb_);
  auto rdm2 = make_shared<RDM<2>>(nactorb_);
  
  return tie(rdm1, rdm2);
}
