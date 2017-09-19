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
  // initialize RDM12 with 0.0 then do ax_plus_y
  auto rdm1 = make_shared<RDM<1>>(nactorb_);
  auto rdm2 = make_shared<RDM<2>>(nactorb_);
  for (int istate = 0; istate != nstate_; ++istate) {
    rdm1_->emplace(istate, istate, rdm1);
    rdm2_->emplace(istate, istate, rdm2);
  }

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
      // only collect on-site RASCI RDM from these two sites
      cout << "  * special treatment for site[" << site << "]" << endl;
      compute_rdm12_ras(cc, site);
    } else {
      // special treatment for first configuration as described by Garnet Chan, 2008
      if (site == 1) {
        cout << "  * special treatment for site[1]" << endl;
      }

      // general treatment
      {
        cout << "  * general treatment for site[" << site << "]" << endl;
        // compute RASCI RDM
        compute_rdm12_ras(cc, site);

        if (site == 1) compute_rdm12_130(cc, site);
      }

      // special treatment for final configuration
      if (site == nsites_-2) {
        cout << "  * special treatment for site[" << site << "]" << endl;
      }
    }

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

  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weights_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weights_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }

  // printing out results
#if 0
  for (int i = 0; i != nstate_; ++i) {
    cout << "rdm1_[" << i << "] : " << endl;
    rdm1_->at(i)->print(1e-3);
    cout << "rdm2_[" << i << "] : " << endl;
    rdm2_->at(i)->print(1e-3);
  }
#endif
}


void ASD_DMRG::compute_rdm12_ras(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  vector<shared_ptr<RDM<1>>> rdm1_vec;
  vector<shared_ptr<RDM<2>>> rdm2_vec;
  const int norb = dvec.front()->space()->norb();
  const int nstate = dvec.size();
  shared_ptr<const RDM<1>> rdm1_ptr;
  shared_ptr<const RDM<2>> rdm2_ptr;
  for (int istate = 0; istate != nstate; ++istate) {
    auto rdm1_tmp = make_shared<RDM<1>>(norb);
    auto rdm2_tmp = make_shared<RDM<2>>(norb);
    for (auto& block : dvec[istate]->sectors()) {
      const int n_lr = block.second->mdim();
      for (int i = 0; i != n_lr; ++i) {
        auto rasvec = make_shared<RASCivec>(block.second->civec(i));
        tie(rdm1_ptr, rdm2_ptr) = rasvec->compute_rdm12_from_rasvec();
        assert(rdm1_ptr->size() == rdm1_tmp->size());
        assert(rdm2_ptr->size() == rdm2_tmp->size());
        blas::ax_plus_y_n(1.0, rdm1_ptr->data(), rdm1_ptr->size(), rdm1_tmp->data());
        blas::ax_plus_y_n(1.0, rdm2_ptr->data(), rdm2_ptr->size(), rdm2_tmp->data());
      }
    }
    rdm1_vec.push_back(rdm1_tmp);
    rdm2_vec.push_back(rdm2_tmp);
  }

  // copy rdm12 elements into total rdm1_ and rdm2_
  vector<int> active_sizes = multisite_->active_sizes();
  const int orb_start = accumulate(active_sizes.begin(), active_sizes.begin()+site, 0);
  cout << "  * orb_start = " << orb_start << endl;
  for (int istate = 0; istate != nstate; ++istate) {
    auto rdm1 = rdm1_vec.at(istate);
    auto rdm2 = rdm2_vec.at(istate);
    const double* const rdm1_source = rdm1->data();
    auto rdm1_target = rdm1_->at(istate);
    auto rdm2_target = rdm2_->at(istate);
    for (int k = 0; k != norb; ++k) {
      copy_n(rdm1_source + k*norb, norb, rdm1_target->element_ptr(orb_start, orb_start + k));
      for (int j = 0; j != norb; ++j) {
        for (int i = 0; i != norb; ++i) {
          copy_n(rdm2->element_ptr(0,i,j,k), norb, rdm2_target->element_ptr(orb_start, (orb_start + i), (orb_start + j), (orb_start + k)));
        }
      }
    }
  }

}


void ASD_DMRG::compute_rdm12_130(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm12_130" << endl;
  const int nstate = dvec.size();
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    const int norb_site = multisite_->active_sizes().at(site);
    cout << "norb_site = " << norb_site << endl;
    const int norb_left = doubleblock->left_block()->norb();
    cout << "norb_left = " << norb_left << endl;
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site*norb_site, norb_left); // matrix to store RDM, use ax_plus_y...

    for (auto& ketblock : prod_civec->sectors()) {
      BlockKey ketkey = ketblock.first;
      BlockKey brakey(ketkey.nelea+1, ketkey.neleb); // <block_bra|a+|block_ket> <ras_bra|a+ a- a-|ras_ket>
      if (!prod_civec->contains_block(brakey)) continue;
      auto ketci_ptr = prod_civec->sector(ketkey);
      auto braci_ptr = prod_civec->sector(brakey);
      for (auto& ketpair : doubleblock->blockpairs(ketkey)) {
        const int ket_offset = ketpair.offset;
        BlockInfo ket_rblock = ketpair.right;
        BlockInfo ket_lblock = ketpair.left;
        // aaaa first
        for (auto& brapair : doubleblock->blockpairs(brakey)) {
          // only loop over blockpairs with the same right block
          if (!(brapair.right == ket_rblock)) continue;
          const int bra_offset = brapair.offset;
          for (int iright = 0; iright != ket_rblock.nstates; ++iright) {
            map<BlockKey, shared_ptr<const RASDvec>> states; // only store transition between left blocks with the same right block
            BlockInfo bra_lblock = brapair.left;
            vector<shared_ptr<RASCivec>> ket_tmpvec;
            for (int iket = 0; iket != ket_lblock.nstates; ++iket)
              ket_tmpvec.push_back(make_shared<RASCivec>(ketci_ptr->civec(ket_offset + iright*ket_lblock.nstates + iket)));
            vector<shared_ptr<RASCivec>> bra_tmpvec;
            for (int ibra = 0; ibra != bra_lblock.nstates; ++ibra)
              bra_tmpvec.push_back(make_shared<RASCivec>(braci_ptr->civec(bra_offset + iright*bra_lblock.nstates + ibra)));

            auto ket_rasdvec = make_shared<const RASDvec>(ket_tmpvec);
            auto bra_rasdvec = make_shared<const RASDvec>(bra_tmpvec);
            const int ket_ras_nelea = prod_civec->nelea() - ketkey.nelea;
            const int ket_ras_neleb = prod_civec->neleb() - ketkey.neleb;
            const int bra_ras_nelea = prod_civec->nelea() - brakey.nelea;
            const int bra_ras_neleb = prod_civec->neleb() - brakey.neleb;
            cout << "bra_ras_a = " << bra_ras_nelea << endl;
            cout << "bra_ras_b = " << bra_ras_neleb << endl;
            cout << "ket_ras_a = " << ket_ras_nelea << endl;
            cout << "ket_ras_b = " << ket_ras_neleb << endl;
            BlockInfo ket_rasblock(ket_ras_nelea, ket_ras_neleb, ket_lblock.nstates);
            BlockInfo bra_rasblock(bra_ras_nelea, bra_ras_neleb, bra_lblock.nstates);
            states[ket_rasblock] = ket_rasdvec;
            states[bra_rasblock] = bra_rasdvec;
            
            GammaForestASD<RASDvec> forest(states);
            forest.compute();

            cout << "RAS part" << endl;
            const size_t bratag = forest.block_tag(bra_rasblock);
            const size_t kettag = forest.block_tag(ket_rasblock);
            list<GammaSQ> gammalist_aaa = {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha};
#if 1
  // play with blas::transpose function

#endif
#if 0
            auto transition_mat_aaa = forest.template get<0>(kettag, bratag, gammalist_aaa); // conjugate
            btas::CRange<3> range(ket_rasblock.nstates, bra_rasblock.nstates, lrint(pow(norb_site, gammalist_aaa.size())));
            auto transition_tensor_aaa = make_shared<btas::Tensor3<double>>(range, move(transition_mat_aaa->storage()));
            // debug
            {
              cout << "printing original tensor" << endl;
              cout << "extent 0 = " << transition_tensor_aaa->extent(0) << endl;
              cout << "extent 1 = " << transition_tensor_aaa->extent(1) << endl;
              cout << "extent 2 = " << transition_tensor_aaa->extent(2) << endl;
              for (int i = 0; i != transition_tensor_aaa->extent(2); ++i) {
                for (int j = 0; j != transition_tensor_aaa->extent(1); ++j) {
                  for (int k = 0; k != transition_tensor_aaa->extent(0); ++k) {
                    if ((*transition_tensor_aaa)(k,j,i) > 1e-6) 
                      cout << (*transition_tensor_aaa)(k,j,i) << " " << "(" << k << ", " << j << ", " << i << ")" << endl;
                  }
                }
              }
            }
#endif
            // transpose the first two dimensions
            btas::CRange<3> range(bra_rasblock.nstates, ket_rasblock.nstates, lrint(pow(norb_site, gammalist_aaa.size())));
            auto transition_mat_aaa = forest.template get<0>(kettag, bratag, gammalist_aaa);
            auto transition_tensor_aaa = make_shared<btas::Tensor3<double>>(range, move(transition_mat_aaa->storage()));
            unique_ptr<double[]> buf(new double[ket_rasblock.nstates * bra_rasblock.nstates]);
            for (int i = 0; i != transition_tensor_aaa->extent(2); ++i) {
              copy_n(&(*transition_tensor_aaa)(0,0,i), ket_rasblock.nstates*bra_rasblock.nstates, buf.get());
              blas::transpose(buf.get(), ket_rasblock.nstates, bra_rasblock.nstates, &(*transition_tensor_aaa)(0,0,i));
            }
#if 0
            {
              cout << "printing new tensor" << endl;
              cout << "extent 0 = " << transition_tensor_aaa->extent(0) << endl;
              cout << "extent 1 = " << transition_tensor_aaa->extent(1) << endl;
              cout << "extent 2 = " << transition_tensor_aaa->extent(2) << endl;
              for (int i = 0; i != transition_tensor_aaa->extent(2); ++i) {
                for (int j = 0; j != transition_tensor_aaa->extent(1); ++j) {
                  for (int k = 0; k != transition_tensor_aaa->extent(0); ++k) {
                    if ((*transition_tensor_aaa)(k,j,i) > 1e-6) 
                      cout << (*transition_tensor_aaa)(k,j,i) << " " << "(" << k << ", " << j << ", " << i << ")" << endl;
                  }
                }
              }
            }
#endif

            cout << "block part" << endl;
            list<GammaSQ> gammalist_a = {GammaSQ::CreateAlpha};
            auto coupling_data_a = left_block->coupling(gammalist_a).at(make_pair(bra_lblock, ket_lblock)).data;

            auto target = rdm_mat->clone();
            contract(1.0, *transition_tensor_aaa, {2,3,0}, *coupling_data_a, {2,3,1}, 0.0, *target, {0,1});
            const int sign = 1 - (((bra_lblock.nelea + bra_lblock.neleb)%2) << 1);
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          } // end of looping over right blocks
        }
      }
    }
    cout << "printing rdm_mat" << endl;
    rdm_mat->print();
  }
}


