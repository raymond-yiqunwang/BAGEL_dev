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

  // one additional sweeping after convergence to collect terms required to construct RDM
  shared_ptr<DMRG_Block1> left_block, right_block;
  for (int site = 0; site != nsites_; ++site) {
    left_block = (site==0) ? nullptr : left_blocks_[site-1];
    right_block = (site==nsites_-1) ? nullptr : right_blocks_[nsites_-site-2];

    // obtain ProductRASCivec 
    vector<shared_ptr<ProductRASCivec>> cc;
    {
      Muffle hide_cout("asd_dmrg_rdm.log", false);
      
      shared_ptr<const Reference> ref = multisite_->build_reference(site, vector<bool>(nsites_, false), metal_);
      shared_ptr<PTree> input = prepare_sweeping_input(site);
      {
        input->put("nclosed", ref->nclosed());
        input->put("metal", metal_);
        read_restricted(input, site);
        vector<int> actvec = multisite_->active_electrons();
        const int nactele = accumulate(actvec.begin(), actvec.end(), input->get<int>("charge"));
        input->put("nactele", nactele);
      }
      shared_ptr<ProductRASCI> prod_ras;
      // testing
      {
        cout << "site = " << site << endl;
        if (left_block)
          cout << "left block is not nullptr" << endl;
        if (right_block)
          cout << "right_block is not nullptr" << endl;
      }
      if ((left_block==nullptr) ^ (right_block==nullptr))
        prod_ras = make_shared<ProductRASCI>(input, ref, ((left_block==nullptr) ? right_block : left_block));
      else {
        assert((left_block!=nullptr) && (right_block!=nullptr));
        prod_ras = make_shared<ProductRASCI>(input, ref, make_shared<const DMRG_Block2>(left_block, right_block));
      }
      prod_ras->compute();
      cc = prod_ras->civectors();
    }

    // construct RDM by collecting terms during sweeping
    if ((site == 0) || (site == nsites_-1)) {
      cout << "  * special treatment for site[" << site << "]" << endl;
      vector<shared_ptr<Matrix>> ras_rdm1_mat;
      vector<shared_ptr<Matrix>> ras_rdm2_mat;
      tie(ras_rdm1_mat, ras_rdm2_mat) = compute_ras_rdm12(cc);
      // debugging
      {
        cout << "printing RAS RDM1 : " << endl;
        for (auto& mat1 : ras_rdm1_mat)
          mat1->print();
        cout << "printing RAS RDM2 : " << endl;
        for (auto& mat2 : ras_rdm2_mat)
          mat2->print();
      }
      // TODO implement -- copy this fragment of RDM into rdm_ with orbital range determined from site
      
    } else {
      // special treatment for first configuration as described by Garnet Chan, 2008
      if (site == 1) {
        cout << "  * special treatment for site[1]" << endl;
      }

      // general treatment
      {
        cout << "  * general treatment for site[" << site << "]" << endl;
      }

      // special treatment for final configuration
      if (site == nsites_-2) {
        cout << "  * special treatment for site[" << site << "]" << endl;
      }
    }

  }

/*
  // for nstate==1, rdm1_av_ = rdm1_->at(0)
  // Needs initialization here because we use daxpy
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(nactorb_);
    rdm2_av_ = make_shared<RDM<2>>(nactorb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }

  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weights_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weights_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }
*/
}


tuple<vector<shared_ptr<Matrix>>, vector<shared_ptr<Matrix>>> ASD_DMRG::compute_ras_rdm12(vector<shared_ptr<ProductRASCivec>> dvec) {
  vector<shared_ptr<Matrix>> rdm1_vec;
  vector<shared_ptr<Matrix>> rdm2_vec;
  const int norb = dvec.front()->space()->norb();
  const int nstate = dvec.size();
  shared_ptr<const RDM<2>> rdm2;
  shared_ptr<const RDM<1>> rdm1;
  for (int istate = 0; istate != nstate; ++istate) {
    auto rdm1_mat = make_shared<Matrix>(norb, norb);
    auto rdm2_mat = make_shared<Matrix>(norb*norb, norb*norb);
    for (auto& block : dvec[istate]->sectors()) {
      const int n_lr = block.second->mdim();
      for (int i = 0; i != n_lr; ++i) {
        auto rasvec = make_shared<RASCivec>(block.second->civec(i));
        tie(rdm1, rdm2) = rasvec->compute_rdm12_from_rasvec();
        assert(rdm1_mat->size() == rdm1->size());
        assert(rdm2_mat->size() == rdm2->size());
        blas::ax_plus_y_n(1.0, rdm1->data(), rdm1_mat->size(), rdm1_mat->data());
        blas::ax_plus_y_n(1.0, rdm2->data(), rdm2_mat->size(), rdm2_mat->data());
      }
    }

    rdm1_vec.push_back(rdm1_mat);
    rdm2_vec.push_back(rdm2_mat);
  }

  return tie(rdm1_vec, rdm2_vec);
}


