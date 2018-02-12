//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_rdm.cc
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu>
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

#include <src/asd/dmrg/asd_dmrg.h>
#include <src/asd/dmrg/product_rasci.h>
#include <src/util/muffle.h>

#include <src/ci/fci/knowles.h>

using namespace std;
using namespace bagel;


void ASD_DMRG::compute_rdm12() {

  cout << endl << "  * Computing ASD-DMRG Reduced Density Matrix.." << endl << endl;
  // initialize RDM12 with 0ull then do ax_plus_y
  const int nact = sref_->nact();
  auto rdm1 = make_shared<RDM<1>>(nact);
  auto rdm2 = make_shared<RDM<2>>(nact);
  for (int istate = 0; istate != nstate_; ++istate) {
    rdm1_->emplace(istate, istate, rdm1);
    rdm2_->emplace(istate, istate, rdm2);
  }

  // pick up one specific wavefunction and compute RDM from GammaForest
  auto right_block = right_blocks_[1];
  assert(right_block);
  const int site = 0;
  vector<shared_ptr<ProductRASCivec>> cc; // wavefunction
  {
    Muffle hide_cout("asd_dmrg_rdm.log", false);
    shared_ptr<const Reference> ref = build_reference(site, vector<bool>(nsites_, false));
    shared_ptr<PTree> input = prepare_sweeping_input(site);
    {
      input->put("nclosed", ref->nclosed());
      input->put("extern_nactele", true);
      const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), input->get<int>("charge"));
      input->put("nactele", nactele);
      read_restricted(input, site);
    }
    auto prod_ras = make_shared<ProductRASCI>(input, ref, right_block);
    prod_ras->compute();
    cc = prod_ras->civectors();
  }

  // compute RDM2
  {
    const list<list<GammaSQ>> gammalists = {
      {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
      {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
      {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
      {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}
    };

    for (int istate = 0; istate != nstate_; ++istate) {
      auto prod_vec = cc.at(istate);
      map<BlockKey, vector<shared_ptr<ProductRASCivec>>> states;
      BlockKey key(prod_vec->nelea(), prod_vec->neleb());
      states[key].push_back(prod_vec);

      GammaForestProdASD forest(states);
      forest.compute();

      for (auto& gammalist : gammalists) {
        auto mat = forest.get(key, key, gammalist);
        blas::ax_plus_y_n(1.0, mat->data(), mat->size(), rdm2_->at(istate)->data());
      }

      // compute RDM1 from RDM2
      const int nactelectrons = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
      for (int i = 0; i != nact; ++i) {
        for (int j = 0; j != nact; ++j) {
          for (int k = 0; k != nact; ++k) {
            *(rdm1_->at(istate)->element_ptr(j, i)) += 1.0/static_cast<double>(nactelectrons-1) * rdm2_->at(istate)->element(j, k, k, i);
          }
        }
      }
      
//      cout << "print RDM1" << endl;
//      rdm1_->at(istate)->print(0.01);
    } // end of loop over nstate
  } // end of compute RDM

/*
  auto fci_info = input_->get_child("fci");
  auto fci = make_shared<KnowlesHandy>(fci_info, sref_->geom(), sref_);
  fci->compute();
  auto fci_rdm1 = fci->rdm1_av();
  cout << "fci rdm1" << endl;
  fci_rdm1->print(0.01);
*/

  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(nact);
    rdm2_av_ = make_shared<RDM<2>>(nact);
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
}

