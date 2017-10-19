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
#include <src/asd/dmrg/gamma_forest_asd2.h>
#include <src/util/muffle.h>
#include <src/ci/fci/knowles.h>

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

    if (site != 0) {
      left_block->compute_left_index(site, multisite_->active_sizes());
    }

    // obtain ProductRASCivec 
    vector<shared_ptr<ProductRASCivec>> cc;
    {
      Muffle hide_cout("asd_dmrg_rdm.log", true);
      
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
      compute_rdm2_ras(cc, site);
    } else {
      // special treatment for first configuration as described by Garnet Chan, 2008
      if (site == 1) {
        cout << "  * special treatment for site[1]" << endl;
        // compute_310
        compute_rdm2_310(cc);

        // compute_301
        compute_rdm2_301(cc);
      }

      // general treatment
      {
        cout << "  * general treatment for site[" << site << "]" << endl;
        // compute RASCI RDM
        compute_rdm2_ras(cc, site);

        // compute_121
        compute_rdm2_121(cc, site);

        // compute_211
        compute_rdm2_211(cc, site);

        // compute_220
        compute_rdm2_220(cc, site);
        
        // compute_130
        compute_rdm2_130(cc, site);

        // compute_031
        compute_rdm2_031(cc, site);
      }

      // special treatment for final configuration
      if (site == nsites_-2) {
        cout << "  * special treatment for site[" << site << "]" << endl;
        // compute_013
        compute_rdm2_013(cc);

        // compute_103
        compute_rdm2_103(cc);

        // compute_022
        compute_rdm2_022(cc);

        // compute_202
        compute_rdm2_202(cc);

        // compute_112
        compute_rdm2_112(cc);
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

  // DEBUGGING -- compare results with FCI RDM
#if 1
  cout << string(12,'=') << endl;
  cout << "DEBUGGING..." << endl;
  cout << string(12,'=') << endl;
  auto fci_input = input_->get_child("fci");
  auto fci = make_shared<KnowlesHandy>(fci_input, multisite_->geom(), multisite_->ref());
  fci->compute();
  const int istate = 0;
  const int norb = fci->norb();
  auto dmrg_rdm2 = rdm2_->at(istate);
  auto fci_rdm2 = fci->rdm2()->at(istate);
  auto diff_tensor = make_shared<const RDM<2>>(*dmrg_rdm2 - *fci_rdm2);
  vector<int> active_size = multisite_->active_sizes();
  for (int site = 1; site != nsites_-1; ++site) {
    cout << "* site = " << site << endl;
    const int site_start = accumulate(active_size.begin(), active_size.begin()+site, 0);
    cout << "site_start = " << site_start << endl;
    const int norb_site = multisite_->active_sizes().at(site);
    const int right_start = site_start + norb_site;
    cout << "right_start = " << right_start << endl;
    
    pair<int, int> range1(0, site_start), range2(site_start, right_start), range3(right_start, norb);
    list<tuple<pair<int, int>, pair<int, int>, pair<int, int>, pair<int, int>>> list_tuplelist;
    
    // general list
    list<int> switchlist = {40/*RAS*/, 130, 31/*031*/, 220, 121, 211};
    if (site == 1) {
      switchlist.insert(switchlist.end(), 310);
      switchlist.insert(switchlist.end(), 301);
    } 
    if (site == nsites_-2) {
      switchlist.insert(switchlist.end(), 13/*013*/);
      switchlist.insert(switchlist.end(), 103);
      switchlist.insert(switchlist.end(), 22/*022*/);
      switchlist.insert(switchlist.end(), 202);
      switchlist.insert(switchlist.end(), 112);
    }

    for (int swch : switchlist) {
      vector<double> diff_vec;
      switch (swch) {
        // RAS
        case 40 : list_tuplelist = {
          {range2, range2, range2, range2}
        }; break;

        case 130 : list_tuplelist = {
          {range1, range2, range2, range2},
          {range2, range1, range2, range2},
          {range2, range2, range1, range2},
          {range2, range2, range2, range1}
        }; break;

        case 31 : list_tuplelist = {
          {range3, range2, range2, range2},
          {range2, range3, range2, range2},
          {range2, range2, range3, range2},
          {range2, range2, range2, range3}
        }; break;

        case 310 : list_tuplelist = {
          {range2, range1, range1, range1},
          {range1, range2, range1, range1},
          {range1, range1, range2, range1},
          {range1, range1, range1, range2}
        }; break;

        case 301 : list_tuplelist = {
          {range3, range1, range1, range1},
          {range1, range3, range1, range1},
          {range1, range1, range3, range1},
          {range1, range1, range1, range3}
        }; break;

        case 13 : list_tuplelist = {
          {range2, range3, range3, range3},
          {range3, range2, range3, range3},
          {range3, range3, range2, range3},
          {range3, range3, range3, range2}
        }; break;

        case 103 : list_tuplelist = {
          {range1, range3, range3, range3},
          {range3, range1, range3, range3},
          {range3, range3, range1, range3},
          {range3, range3, range3, range1}
        }; break;

        case 220 : list_tuplelist = {
          {range1, range1, range2, range2},
          {range2, range2, range1, range1},
          {range2, range1, range2, range1},
          {range1, range2, range1, range2},
          {range2, range1, range1, range2},
          {range1, range2, range2, range1}
        }; break;

        case 22 : list_tuplelist = {
          {range2, range2, range3, range3},
          {range3, range3, range2, range2},
          {range2, range3, range2, range3},
          {range3, range2, range3, range2},
          {range2, range3, range3, range2},
          {range3, range2, range2, range3}
        }; break;

        case 202 : list_tuplelist = {
          {range1, range1, range3, range3},
          {range3, range3, range1, range1},
          {range1, range3, range1, range3},
          {range3, range1, range3, range1},
          {range1, range3, range3, range1},
          {range3, range1, range1, range3}
        }; break;

        case 121 : list_tuplelist = {
          {range2, range2, range1, range3},
          {range1, range3, range2, range2},
          {range2, range2, range3, range1},
          {range3, range1, range2, range2},
          
          {range3, range2, range1, range2},
          {range1, range2, range3, range2},
          {range2, range3, range2, range1},
          {range2, range1, range2, range3},

          {range2, range1, range3, range2},
          {range3, range2, range2, range1},
          {range1, range2, range2, range3},
          {range2, range3, range1, range2}
        }; break;

        case 211 : list_tuplelist = {
          {range1, range1, range2, range3},
          {range2, range3, range1, range1},
          {range1, range1, range3, range2},
          {range3, range2, range1, range1},

          {range2, range1, range3, range1},
          {range3, range1, range2, range1},
          {range1, range2, range1, range3},
          {range1, range3, range1, range2},

          {range2, range1, range1, range3},
          {range1, range3, range2, range1},
          {range1, range2, range3, range1},
          {range3, range1, range1, range2}
        }; break;

        case 112 : list_tuplelist = {
          {range2, range1, range3, range3},
          {range3, range3, range2, range1},
          {range1, range2, range3, range3},
          {range3, range3, range1, range2},
          
          {range2, range3, range1, range3},
          {range1, range3, range2, range3},
          {range3, range2, range3, range1},
          {range3, range1, range3, range2},

          {range2, range3, range3, range1},
          {range3, range1, range2, range3},
          {range3, range2, range1, range3},
          {range1, range3, range3, range2}
        }; break;
      }

      cout << "switch : " << setfill('0') << setw(3) << swch << setfill(' ') << endl;
      for (auto& tuple_list : list_tuplelist) {
        for (int l = get<3>(tuple_list).first; l != get<3>(tuple_list).second; ++l) {
          for (int k = get<2>(tuple_list).first; k != get<2>(tuple_list).second; ++k) {
            for (int j = get<1>(tuple_list).first; j != get<1>(tuple_list).second; ++j) {
              for (int i = get<0>(tuple_list).first; i != get<0>(tuple_list).second; ++i) {
                diff_vec.push_back(fabs(diff_tensor->element(i,j,k,l)));
                if (diff_vec.back() > 1.0e-6) {
                  cout << "diff large (" << i << "," << j << "," << k << "," << l <<"), fci : " << fci_rdm2->element(i,j,k,l) << ", dmrg : " << dmrg_rdm2->element(i,j,k,l) << endl;
                }
              }
            }
          }
        }
      }
      
      double diff_sum = 0.0;
      std::for_each(diff_vec.begin(), diff_vec.end(), [&diff_sum] (const double v) { diff_sum += v*v; });
      const double rms = std::sqrt(diff_sum / static_cast<double>(diff_vec.size()));
      cout << "diff_rms = " << setw(16) << setprecision(12) << rms << ", RDM block size = " << diff_vec.size() << endl;
    } // end of looping over switchlist
  }

#endif
}


// TODO find a better algorithm
void ASD_DMRG::compute_rdm2_ras(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_ras" << endl;
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
        tie(rdm1_ptr, rdm2_ptr) = rasvec->compute_rdm2_from_rasvec();
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


void ASD_DMRG::compute_rdm2_130(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_130" << endl;
  const list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on site}, {ops on left}, {left nele change}}
    { {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {-1, 0} },
    { {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, {-1, 0} },
    { {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  { 0,-1} },
    { {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  { 0,-1} }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site*norb_site, norb_left); // matrix to store RDM, use ax_plus_y...
    auto unordered_rdm = rdm_mat->clone();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey ket_seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
          const int ketpairoffset = ketpair.offset;
          // left bra-ket pair
          BlockKey ket_leftkey = ketpair.left;
          const int ket_leftnstates = ketpair.left.nstates;
          BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
          if (!left_block->contains(bra_leftkey)) continue;
          const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
          // right block info
          BlockKey rightkey = ketpair.right;
          const int rightnstates = ketpair.right.nstates;
          // bra info
          BlockKey bra_seckey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
          if (!prod_civec->contains_block(bra_seckey)) continue;
          auto brapair = doubleblock->blockpairs(bra_seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
          const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
          const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
          const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;

          // site transition tensor
          shared_ptr<btas::Tensor3<double>> site_transition_tensor;
          {
            map<BlockKey, shared_ptr<const RASDvec>> ket_states;
            vector<shared_ptr<RASCivec>> ket_vecs;
            for (int iright = 0; iright != rightnstates; ++iright)
              for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + iright*ket_leftnstates)));
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);

            map<BlockKey, shared_ptr<const RASDvec>> bra_states;
            vector<shared_ptr<RASCivec>> bra_vecs;
            for (int iright = 0; iright != rightnstates; ++iright)
              for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + iright*bra_leftnstates)));
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
            bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);

            GammaForestASD2<RASDvec> forest(bra_states, ket_states);
            forest.compute();

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
            shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
            btas::CRange<3> range1(ket_leftnstates*bra_leftnstates, rightnstates, rightnstates*lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
            vector<double> buf1(ket_leftnstates*bra_leftnstates*rightnstates);
            for (int i = 0; i != tmp_tensor->extent(2); ++i) {
              copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
              blas::transpose(buf1.data(), bra_leftnstates*rightnstates, ket_leftnstates, &(*tmp_tensor)(0,0,i));
            }

            btas::CRange<3> site_range(ket_leftnstates*bra_leftnstates, rightnstates*rightnstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
          }

          // transposed left coupling tensor
          shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
          {
            btas::CRange<3> left_range(ket_leftnstates, bra_leftnstates, lrint(pow(norb_left, get<1>(gammalist_tuple).size())));
            shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_leftkey, bra_leftkey)).data;
            left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
          }

          // right coupling tensor
          shared_ptr<const Matrix> right_coupling_tensor;
          {
            auto eye = make_shared<Matrix>(rightnstates, rightnstates);
            eye->unit();
            right_coupling_tensor = eye;
          }

          // contraction
          auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
          btas::CRange<3> intermediate_range(site_transition_tensor->extent(1), site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          auto intermediate_tensor = make_shared<const btas::Tensor3<double>>(intermediate_range, move(contract_site_left->storage()));
          auto tmp_mat = make_shared<Matrix>(site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          auto gmat = group(*tmp_mat,0,2);
          contract(1.0, group(*right_coupling_tensor,0,2), {1}, group(*intermediate_tensor,1,3), {1,0}, 0.0, gmat, {0});
          const double sign = static_cast<double>(1 - (((ket_ras_nelea+ket_ras_neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), unordered_rdm->data());
        }
      }
    }

    // reorder left block orbitals
    for (int j = 0; j != rdm_mat->mdim(); ++j) {
      const int target = left_block->left_index(j);
      copy_n(unordered_rdm->element_ptr(0,j), rdm_mat->ndim(), rdm_mat->element_ptr(0,target));
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != norb_left; ++p) {
      for (int j = 0; j != norb_site; ++j) {
        for (int i = 0; i != norb_site; ++i) {
          for (int k = 0; k != norb_site; ++k) {
            const double value = *rdm_mat->element_ptr(k + i*norb_site + j*norb_site*norb_site, p);
            rdm2_target->element(i+norb_left, j+norb_left, k+norb_left, p) = value;
            rdm2_target->element(k+norb_left, p, i+norb_left, j+norb_left) = value;
            rdm2_target->element(j+norb_left, i+norb_left, p, k+norb_left) = value;
            rdm2_target->element(p, k+norb_left, j+norb_left, i+norb_left) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_130


void ASD_DMRG::compute_rdm2_031(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_031" << endl;
  const list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on site}, {ops on right}, {right nele change}}
    { {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {-1, 0} },
    { {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, {-1, 0} },
    { {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  { 0,-1} },
    { {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  { 0,-1} }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site*norb_site, norb_right); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey ket_seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
          const int ketpairoffset = ketpair.offset;
          // right bra-ket pair
          BlockKey ket_rightkey = ketpair.right;
          const int ket_rightnstates = ketpair.right.nstates;
          BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple).first, ket_rightkey.neleb + get<2>(gammalist_tuple).second);
          if (!right_block->contains(bra_rightkey)) continue;
          const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;

          // left block info
          BlockKey leftkey = ketpair.left;
          const int leftnstates = ketpair.left.nstates;
          // bra info
          BlockKey bra_seckey(leftkey.nelea + bra_rightkey.nelea, leftkey.neleb + bra_rightkey.neleb);
          if (!prod_civec->contains_block(bra_seckey)) continue;
          auto brapair = doubleblock->blockpairs(bra_seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
          const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
          const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
          const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;

          // site transition tensor
          shared_ptr<btas::Tensor3<double>> site_transition_tensor;
          {
            map<BlockKey, shared_ptr<const RASDvec>> ket_states;
            vector<shared_ptr<RASCivec>> ket_vecs;
            for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
              for (int ileft = 0; ileft != leftnstates; ++ileft)
                ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ileft + ket_ir*leftnstates)));
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);

            map<BlockKey, shared_ptr<const RASDvec>> bra_states;
            vector<shared_ptr<RASCivec>> bra_vecs;
            for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
              for (int ileft = 0; ileft != leftnstates; ++ileft)
                bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + ileft + bra_ir*leftnstates)));
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
            bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);

            GammaForestASD2<RASDvec> forest(bra_states, ket_states);
            forest.compute();

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
            shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
            btas::CRange<3> range1(leftnstates*leftnstates, bra_rightnstates, ket_rightnstates*lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
            vector<double> buf1(leftnstates*leftnstates*bra_rightnstates);
            for (int i = 0; i != tmp_tensor->extent(2); ++i) {
              copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
              blas::transpose(buf1.data(), leftnstates*bra_rightnstates, leftnstates, &(*tmp_tensor)(0,0,i));
            }

            btas::CRange<3> site_range(leftnstates*leftnstates, bra_rightnstates*ket_rightnstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
          }

          // left coupling tensor
          shared_ptr<const Matrix> left_coupling_tensor;
          {
            auto eye = make_shared<Matrix>(leftnstates, leftnstates);
            eye->unit();
            left_coupling_tensor = eye;
          }
          
          // right coupling tensor
          shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
          {
            btas::CRange<3> right_range(bra_rightnstates, ket_rightnstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
            shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;
            right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
            vector<double> buf(bra_rightnstates*ket_rightnstates);
            for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
              copy_n(&(*right_coupling_tensor)(0,0,i), buf.size(), buf.data());
              blas::transpose(buf.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,i));
            }
          }

          // contraction
          auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1), site_transition_tensor->extent(2));
          auto gmat = group(*contract_site_left,0,2);
          contract(1.0, group(*site_transition_tensor,1,3), {1,0}, group(*left_coupling_tensor,0,2), {1}, 0.0, gmat, {0});
          auto outmat = make_shared<Matrix>(site_transition_tensor->extent(2), right_coupling_tensor->extent(2));
          contract(1.0, *contract_site_left, {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *outmat, {0,1});
          const double sign = static_cast<double>(1 - (((ket_ras_nelea+ket_ras_neleb + leftkey.nelea + leftkey.neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, outmat->data(), outmat->size(), rdm_mat->data());
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_right; ++q) {
      for (int j = 0; j != norb_site; ++j) {
        for (int i = 0; i != norb_site; ++i) {
          for (int k = 0; k != norb_site; ++k) {
            const double value = *rdm_mat->element_ptr(k + i*norb_site + j*norb_site*norb_site, q);
            rdm2_target->element(i+norb_left, j+norb_left, k+norb_left, q+norb_left+norb_site) = value;
            rdm2_target->element(k+norb_left, q+norb_left+norb_site, i+norb_left, j+norb_left) = value;
            rdm2_target->element(j+norb_left, i+norb_left, q+rightoffset, k+norb_left) = value;
            rdm2_target->element(q+rightoffset, k+norb_left, j+norb_left, i+norb_left) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_031


void ASD_DMRG::compute_rdm2_310(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_310" << endl;
  const list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on site}, {ops on left}, {left nele change}}
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {-1, 0} },
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {-1, 0} },
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  { 0,-1} },
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, { 0,-1} }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(1); // site == 1
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site, lrint(pow(norb_left,3))); // matrix to store RDM, use ax_plus_y...
    auto unordered_rdm = rdm_mat->clone();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey ket_seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
          const int ketpairoffset = ketpair.offset;
          // left bra-ket pair
          BlockKey ket_leftkey = ketpair.left;
          const int ket_leftnstates = ketpair.left.nstates;
          BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
          if (!left_block->contains(bra_leftkey)) continue;
          const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
          // right block info
          BlockKey rightkey = ketpair.right;
          const int rightnstates = ketpair.right.nstates;
          // bra info
          BlockKey bra_seckey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
          if (!prod_civec->contains_block(bra_seckey)) continue;
          auto brapair = doubleblock->blockpairs(bra_seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
          const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
          const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
          const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;

          // site transition tensor
          shared_ptr<btas::Tensor3<double>> site_transition_tensor;
          {
            map<BlockKey, shared_ptr<const RASDvec>> ket_states;
            vector<shared_ptr<RASCivec>> ket_vecs;
            for (int iright = 0; iright != rightnstates; ++iright)
              for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + iright*ket_leftnstates)));
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);

            map<BlockKey, shared_ptr<const RASDvec>> bra_states;
            vector<shared_ptr<RASCivec>> bra_vecs;
            for (int iright = 0; iright != rightnstates; ++iright)
              for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + iright*bra_leftnstates)));
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
            bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);

            GammaForestASD2<RASDvec> forest(bra_states, ket_states);
            forest.compute();

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
            shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
            btas::CRange<3> range1(ket_leftnstates*bra_leftnstates, rightnstates, rightnstates*lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
            vector<double> buf1(ket_leftnstates*bra_leftnstates*rightnstates);
            for (int i = 0; i != tmp_tensor->extent(2); ++i) {
              copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
              blas::transpose(buf1.data(), bra_leftnstates*rightnstates, ket_leftnstates, &(*tmp_tensor)(0,0,i));
            }

            btas::CRange<3> site_range(ket_leftnstates*bra_leftnstates, rightnstates*rightnstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
          }

          // transposed left coupling tensor
          shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
          {
            btas::CRange<3> left_range(ket_leftnstates, bra_leftnstates, lrint(pow(norb_left, get<1>(gammalist_tuple).size())));
            shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_leftkey, bra_leftkey)).data;
            left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
          }

          // right coupling tensor
          shared_ptr<const Matrix> right_coupling_tensor;
          {
            auto eye = make_shared<Matrix>(rightnstates, rightnstates);
            eye->unit();
            right_coupling_tensor = eye;
          }

          // contraction
          auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
          btas::CRange<3> intermediate_range(site_transition_tensor->extent(1), site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          auto intermediate_tensor = make_shared<const btas::Tensor3<double>>(intermediate_range, move(contract_site_left->storage()));
          auto tmp_mat = make_shared<Matrix>(site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          auto gmat = group(*tmp_mat,0,2);
          contract(1.0, group(*right_coupling_tensor,0,2), {1}, group(*intermediate_tensor,1,3), {1,0}, 0.0, gmat, {0});
          const double sign = static_cast<double>(1 - (((ket_ras_nelea+ket_ras_neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), unordered_rdm->data());
        }
      }
    }

    // reorder left block orbitals
    for (int index3 = 0; index3 != norb_left; ++index3) {
      for (int index2 = 0; index2 != norb_left; ++index2) {
        for (int index1 = 0; index1 != norb_left; ++index1) {
          const int unordered_index = index1 + index2*norb_left + index3*norb_left*norb_left;
          const int target_idx1 = left_block->left_index(index1);
          const int target_idx2 = left_block->left_index(index2);
          const int target_idx3 = left_block->left_index(index3);
          const int target = target_idx1 + target_idx2*norb_left + target_idx3*norb_left*norb_left;
          copy_n(unordered_rdm->element_ptr(0,unordered_index), rdm_mat->ndim(), rdm_mat->element_ptr(0,target));
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_left; ++q) {
      for (int r = 0; r != norb_left; ++r) {
        for (int p = 0; p != norb_left; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            const double value = *rdm_mat->element_ptr(i, p + r*norb_left + q*norb_left*norb_left);
            rdm2_target->element(i+norb_left, p, q, r) = value;
            rdm2_target->element(q, r, i+norb_left, p) = value;
            rdm2_target->element(p, i+norb_left, r, q) = value;
            rdm2_target->element(r, q, p, i+norb_left) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_310



void ASD_DMRG::compute_rdm2_301(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_301" << endl;
  list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on left}, {ops on right}, {left nele change} }
    { {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {1,0} },
    { {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, {1,0} }, 
    { {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  {0,1} },
    { {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  {0,1} } 
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site =  multisite_->active_sizes().at(1); 
    const int rightoffset = norb_left + norb_site;
    auto rdm_mat = make_shared<Matrix>(norb_right, lrint(pow(norb_left, 3))); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = ketpair.offset;
          // left bra-ket pair
          BlockKey ket_leftkey = ketpair.left.key();
          const int ket_leftnstates = ketpair.left.nstates;
          BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
          if (!left_block->contains(bra_leftkey)) continue;
          const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockKey ket_rightkey = ketpair.right.key();
          const int ket_rightnstates = ketpair.right.nstates;
          BlockKey bra_rightkey(seckey.nelea - bra_leftkey.nelea, seckey.neleb - bra_leftkey.neleb);
          if (!right_block->contains(bra_rightkey)) continue;
          const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          // site transition matrix
          shared_ptr<Matrix> site_transition_mat;
          {
            site_transition_mat = make_shared<Matrix>(bra_leftnstates*ket_leftnstates, bra_rightnstates*ket_rightnstates);
            double* mat_ptr = site_transition_mat->data();
            auto ci_ptr = prod_civec->sector(seckey);
            for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                    *mat_ptr++ = ci_ptr->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates).dot_product(ci_ptr->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates));
          } // L-L'-R'-R

          // left coupling tensor L-L'
          shared_ptr<const btas::Tensor3<double>> left_coupling_tensor = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;
          
          // right coupling tensor R'-R
          shared_ptr<const btas::Tensor3<double>> right_coupling_tensor = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;

          // contraction>
          auto contract_site_left = make_shared<Matrix>(site_transition_mat->mdim(), left_coupling_tensor->extent(2));
          contract(1.0, *site_transition_mat, {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
          auto tmp = make_shared<Matrix>(right_coupling_tensor->extent(2), left_coupling_tensor->extent(2));
          contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *contract_site_left, {2,1}, 0.0, *tmp, {0,1});
          const double sign = static_cast<double>(1 - (((ket_leftkey.nelea+ket_leftkey.neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, tmp->data(), tmp->size(), rdm_mat->data());
        }
      }
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int r = 0; r != norb_left; ++r) {
      for (int q = 0; q != norb_left; ++q) {
        for (int p = 0; p != norb_left; ++p) {
          for (int t = 0; t != norb_right; ++t) {
            const double value = *rdm_mat->element_ptr(t, p+q*norb_left+r*norb_left*norb_left);
            rdm2_target->element(q, r, p, t+rightoffset) = value;
            rdm2_target->element(p, t+rightoffset, q, r) = value;
            rdm2_target->element(r, q, t+rightoffset, p) = value;
            rdm2_target->element(t+rightoffset, p, r, q) = value;
          }
        }
      }
    }

  } // end of looping over nstates
} // end of compute_301


void ASD_DMRG::compute_rdm2_013(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_013" << endl;
  const list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on site}, {ops on right}, {right nele change}}
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {-1, 0} },
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {-1, 0} },
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  { 0,-1} },
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, { 0,-1} }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_right*norb_right*norb_right, norb_site); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey ket_seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
          const int ketpairoffset = ketpair.offset;
          // right bra-ket pair
          BlockKey ket_rightkey = ketpair.right;
          const int ket_rightnstates = ketpair.right.nstates;
          BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple).first, ket_rightkey.neleb + get<2>(gammalist_tuple).second);
          if (!right_block->contains(bra_rightkey)) continue;
          const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;

          // left block info
          BlockKey leftkey = ketpair.left;
          const int leftnstates = ketpair.left.nstates;
          // bra info
          BlockKey bra_seckey(leftkey.nelea + bra_rightkey.nelea, leftkey.neleb + bra_rightkey.neleb);
          if (!prod_civec->contains_block(bra_seckey)) continue;
          auto brapair = doubleblock->blockpairs(bra_seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
          const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
          const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
          const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;

          // site transition tensor
          shared_ptr<btas::Tensor3<double>> site_transition_tensor;
          {
            map<BlockKey, shared_ptr<const RASDvec>> ket_states;
            vector<shared_ptr<RASCivec>> ket_vecs;
            for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
              for (int ileft = 0; ileft != leftnstates; ++ileft)
                ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ileft + ket_ir*leftnstates)));
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);

            map<BlockKey, shared_ptr<const RASDvec>> bra_states;
            vector<shared_ptr<RASCivec>> bra_vecs;
            for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
              for (int ileft = 0; ileft != leftnstates; ++ileft)
                bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + ileft + bra_ir*leftnstates)));
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
            bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);

            GammaForestASD2<RASDvec> forest(bra_states, ket_states);
            forest.compute();

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
            shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
            btas::CRange<3> range1(leftnstates*leftnstates, bra_rightnstates, ket_rightnstates*norb_site);
            auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
            vector<double> buf1(leftnstates*leftnstates*bra_rightnstates);
            for (int i = 0; i != tmp_tensor->extent(2); ++i) {
              copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
              blas::transpose(buf1.data(), leftnstates*bra_rightnstates, leftnstates, &(*tmp_tensor)(0,0,i));
            }

            btas::CRange<3> site_range(leftnstates*leftnstates, bra_rightnstates*ket_rightnstates, norb_site);
            site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
          }

          // left coupling tensor
          shared_ptr<const Matrix> left_coupling_tensor;
          {
            auto eye = make_shared<Matrix>(leftnstates, leftnstates);
            eye->unit();
            left_coupling_tensor = eye;
          }
          
          // right coupling tensor
          shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
          {
            btas::CRange<3> right_range(bra_rightnstates, ket_rightnstates, lrint(pow(norb_right, 3)));
            shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;
            right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
            vector<double> buf(bra_rightnstates*ket_rightnstates);
            for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
              copy_n(&(*right_coupling_tensor)(0,0,i), buf.size(), buf.data());
              blas::transpose(buf.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,i));
            }
          }

          // contraction
          auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1), site_transition_tensor->extent(2));
          auto gmat = group(*contract_site_left,0,2);
          contract(1.0, group(*site_transition_tensor,1,3), {1,0}, group(*left_coupling_tensor,0,2), {1}, 0.0, gmat, {0});
          auto outmat = make_shared<Matrix>(right_coupling_tensor->extent(2), site_transition_tensor->extent(2));
          contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *contract_site_left, {2,1}, 0.0, *outmat, {0,1});
          const double sign = static_cast<double>(1 - (((ket_ras_nelea+ket_ras_neleb + leftkey.nelea + leftkey.neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, outmat->data(), outmat->size(), rdm_mat->data());
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int i = 0; i != norb_site; ++i) {
      for (int q = 0; q != norb_right; ++q) {
        for (int r = 0; r != norb_right; ++r) {
          for (int p = 0; p != norb_right; ++p) {
            const double value = *rdm_mat->element_ptr(p + r*norb_right + q*norb_right*norb_right, i);
            rdm2_target->element(i+norb_left, p+rightoffset, q+rightoffset, r+rightoffset) = value;
            rdm2_target->element(q+rightoffset, r+rightoffset, i+norb_left, p+rightoffset) = value;
            rdm2_target->element(p+rightoffset, i+norb_left, r+rightoffset, q+rightoffset) = value;
            rdm2_target->element(r+rightoffset, q+rightoffset, p+rightoffset, i+norb_left) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_013


void ASD_DMRG::compute_rdm2_103(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_103" << endl;
  list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>>> gammalist_tuple_list = { 
    // { {ops on left}, {ops on right}, {left nele change} }
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {1,0} },
    { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {1,0} }, 
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {0,1} },
    { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {0,1} } 
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site =  multisite_->active_sizes().at(nsites_-2);
    const int rightoffset = norb_left + norb_site;
    auto rdm_mat = make_shared<Matrix>(lrint(pow(norb_right, 3)), norb_left); // matrix to store RDM, use ax_plus_y...
    auto unordered_rdm = rdm_mat->clone();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& ketpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = ketpair.offset;
          // left bra-ket pair
          BlockKey ket_leftkey = ketpair.left.key();
          const int ket_leftnstates = ketpair.left.nstates;
          BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
          if (!left_block->contains(bra_leftkey)) continue;
          const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockKey ket_rightkey = ketpair.right.key();
          const int ket_rightnstates = ketpair.right.nstates;
          BlockKey bra_rightkey(seckey.nelea - bra_leftkey.nelea, seckey.neleb - bra_leftkey.neleb);
          if (!right_block->contains(bra_rightkey)) continue;
          const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          // site transition matrix
          shared_ptr<Matrix> site_transition_mat;
          {
            site_transition_mat = make_shared<Matrix>(bra_leftnstates*ket_leftnstates, bra_rightnstates*ket_rightnstates);
            double* mat_ptr = site_transition_mat->data();
            auto ci_ptr = prod_civec->sector(seckey);
            for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                    *mat_ptr++ = ci_ptr->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates).dot_product(ci_ptr->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates));
          } // L-L'-R'-R

          // left coupling tensor L-L'
          shared_ptr<const btas::Tensor3<double>> left_coupling_tensor = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;
          
          // right coupling tensor R'-R
          shared_ptr<btas::Tensor3<double>> right_coupling_tensor = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;

          // contraction>
          auto contract_site_left = make_shared<Matrix>(site_transition_mat->mdim(), left_coupling_tensor->extent(2));
          contract(1.0, *site_transition_mat, {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
          auto tmp = make_shared<Matrix>(right_coupling_tensor->extent(2), left_coupling_tensor->extent(2));
          contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *contract_site_left, {2,1}, 0.0, *tmp, {0,1});
          const double sign = static_cast<double>(1 - (((ket_leftkey.nelea+ket_leftkey.neleb) & 1) << 1));
          blas::ax_plus_y_n(sign, tmp->data(), tmp->size(), unordered_rdm->data());
        }
      }
    }

    // reorder left block orbitals
    for (int j = 0; j != rdm_mat->mdim(); ++j) {
      const int target = left_block->left_index(j);
      copy_n(unordered_rdm->element_ptr(0,j), rdm_mat->ndim(), rdm_mat->element_ptr(0,target));
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int t = 0; t != norb_left; ++t) {
      for (int q = 0; q != norb_right; ++q) {
        for (int r = 0; r != norb_right; ++r) {
          for (int p = 0; p != norb_right; ++p) {
            const double value = *rdm_mat->element_ptr(p+r*norb_right+q*norb_right*norb_right, t);
            rdm2_target->element(t, p+rightoffset, q+rightoffset, r+rightoffset) = value;
            rdm2_target->element(q+rightoffset, r+rightoffset, t, p+rightoffset) = value;
            rdm2_target->element(p+rightoffset, t, r+rightoffset, q+rightoffset) = value;
            rdm2_target->element(r+rightoffset, q+rightoffset, p+rightoffset, t) = value;
          }
        }
      }
    }

  } // end of looping over nstates
} // end of compute_103


void ASD_DMRG::compute_rdm2_220(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_220" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>, bool, bool, bool>>> gammalist_tuple_list_list = { 
    // { {ops on site}, {ops on left}, {left nele change}, trans_site, trans_left, duplicate }
    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false },
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false },
    },

    {
      { {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha,  GammaSQ::AnnihilateAlpha},  { 2, 0}, false, true,  false },
      { {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha},  {GammaSQ::AnnihilateBeta,   GammaSQ::AnnihilateAlpha},  { 1, 1}, false, true,  true  },
      { {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta},   {GammaSQ::AnnihilateBeta,   GammaSQ::AnnihilateBeta},   { 0, 2}, false, true,  false }
    },

    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},       {-1, 1}, true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},       { 1,-1}, false, true,  false }
    }
    /// perhaps I have OCD...
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_left); // matrix to store RDM, use ax_plus_y...
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {
        const bool trans_site = get<3>(gammalist_tuple);
        const bool trans_left = get<4>(gammalist_tuple);
        const bool duplicate = get<5>(gammalist_tuple);
        for (auto& isec : prod_civec->sectors()) {
          BlockKey ket_seckey = isec.first;
          for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
            BlockKey ket_leftkey = ketpair.left;
            const int ket_leftnstates = ketpair.left.nstates;
            const int ketpairoffset = ketpair.offset;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
            if (!left_block->contains(bra_leftkey)) continue;
            const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
            // right block info
            BlockKey rightkey = ketpair.right;
            const int rightnstates = ketpair.right.nstates;
            // bra info
            BlockKey bra_seckey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!prod_civec->contains_block(bra_seckey)) continue;
            auto brapair = doubleblock->blockpairs(bra_seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &rightkey] (const DMRG::BlockPair& bp)
              { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;
  
            // site transition tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              map<BlockKey, shared_ptr<const RASDvec>> ket_states;
              vector<shared_ptr<RASCivec>> ket_vecs;
              for (int iright = 0; iright != rightnstates; ++iright)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + iright*ket_leftnstates)));
              const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
              const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;
              BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
              ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);
  
              map<BlockKey, shared_ptr<const RASDvec>> bra_states;
              vector<shared_ptr<RASCivec>> bra_vecs;
              for (int iright = 0; iright != rightnstates; ++iright)
                for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                  bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + iright*bra_leftnstates)));
              const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
              const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
              BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
              bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);
  
              if (trans_site) swap(bra_states, ket_states);
  
              GammaForestASD2<RASDvec> forest(bra_states, ket_states);
              forest.compute();
  
              const size_t bra_rastag = forest.block_tag(bra_raskey);
              const size_t ket_rastag = forest.block_tag(ket_raskey);
              const size_t bratag = trans_site ? ket_rastag : bra_rastag;
              const size_t kettag = trans_site ? bra_rastag : ket_rastag;
              assert(forest.template exist<0>(bratag, kettag, get<0>(gammalist_tuple)));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bratag, kettag, get<0>(gammalist_tuple));
              btas::CRange<3> range1(ket_leftnstates*bra_leftnstates, rightnstates, rightnstates*lrint(pow(norb_site, 2)));
              auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
              vector<double> buf1(ket_leftnstates*bra_leftnstates*rightnstates);
              const int dim1 = trans_site ? ket_leftnstates*rightnstates : bra_leftnstates*rightnstates;
              const int dim2 = trans_site ? bra_leftnstates : ket_leftnstates;
              for (int i = 0; i != tmp_tensor->extent(2); ++i) {
                copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
                blas::transpose(buf1.data(), dim1, dim2, &(*tmp_tensor)(0,0,i));
              } 
  
              btas::CRange<3> site_range(ket_leftnstates*bra_leftnstates, rightnstates*rightnstates, lrint(pow(norb_site, 2)));
              site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
            } // trans_site ? L-L'-R'-R : L'-L-R-R'
  
            // left coupling tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              const int dim1 = trans_site ? bra_leftnstates : ket_leftnstates;
              const int dim2 = trans_site ? ket_leftnstates : bra_leftnstates;
              btas::CRange<3> left_range(dim1, dim2, lrint(pow(norb_left, 2)));
              BlockKey brakey = trans_left ? ket_leftkey : bra_leftkey;
              BlockKey ketkey = trans_left ? bra_leftkey : ket_leftkey;
              shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
              if (!(trans_site^trans_left)) {
                vector<double> left_buf(bra_leftnstates*ket_leftnstates);
                for (int i = 0; i != left_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*left_coupling_tensor)(0,0,i), left_buf.size(), left_buf.data());
                  blas::transpose(left_buf.data(), left_coupling_tensor->extent(1), left_coupling_tensor->extent(0), &(*left_coupling_tensor)(0,0,i));
                }
              }
            } // trans_site ? L-L' : L'-L
  
            // right coupling tensor
            shared_ptr<const Matrix> right_coupling_tensor;
            {
              auto eye = make_shared<Matrix>(rightnstates, rightnstates);
              eye->unit();
              right_coupling_tensor = eye;
            }
  
            // contraction
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
            btas::CRange<3> intermediate_range(site_transition_tensor->extent(1), site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            auto intermediate_tensor = make_shared<const btas::Tensor3<double>>(intermediate_range, move(contract_site_left->storage()));
            auto tmp_mat = make_shared<Matrix>(site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            auto gmat = group(*tmp_mat,0,2);
            contract(1.0, group(*right_coupling_tensor,0,2), {1}, group(*intermediate_tensor,1,3), {1,0}, 0.0, gmat, {0});
            
            // swap indices
            if (trans_site) {
              vector<double> buf2(norb_site*norb_site);
              for (int i = 0; i != left_coupling_tensor->extent(2); ++i) {
                copy_n(tmp_mat->element_ptr(0,i), buf2.size(), buf2.data());
                blas::transpose(buf2.data(), norb_site, norb_site, tmp_mat->element_ptr(0,i));
              }
            }
            if (trans_left) {
              auto tmp2 = tmp_mat->clone();
              for (int i = 0; i != norb_left; ++i)
                for (int j = 0; j != norb_left; ++j)
                  copy_n(tmp_mat->element_ptr(0, j+i*norb_left), tmp_mat->ndim(), tmp2->element_ptr(0, i+j*norb_left));
              tmp_mat = tmp2;
            }
            if (duplicate) {
              auto dup_mat = tmp_mat->clone();
              for (int i = 0; i != norb_left; ++i)
                for (int j = 0; j != norb_left; ++j)
                  copy_n(tmp_mat->element_ptr(0,j+i*norb_left), tmp_mat->ndim(), dup_mat->element_ptr(0, i+j*norb_left));

              vector<double> buf3(norb_site*norb_site);
              for (int k = 0; k != dup_mat->mdim(); ++k) {
                copy_n(dup_mat->element_ptr(0,k), buf3.size(), buf3.data());
                blas::transpose(buf3.data(), norb_site, norb_site, dup_mat->element_ptr(0,k));
              }

              *tmp_mat += *dup_mat;
            }
            
            // copy this bra-ket pair density matrix to target
            const double sign = (scheme == 3) ? -1.0 : 1.0;
            blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), unordered_rdm->data());
          }
        }
      }
  
      // reorder left block orbitals
      for (int index2 = 0; index2 != norb_left; ++index2) {
        for (int index1 = 0; index1 != norb_left; ++index1) {
          const int unordered_index = index1 + index2*norb_left;
          const int target_idx1 = left_block->left_index(index1);
          const int target_idx2 = left_block->left_index(index2);
          const int target = target_idx1 + target_idx2*norb_left;
          copy_n(unordered_rdm->element_ptr(0,unordered_index), rdm_mat->ndim(), rdm_mat->element_ptr(0,target));
        }
      }
      
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_left);
                rdm2_target->element(i+norb_left, j+norb_left, p, q) = value;
                rdm2_target->element(p, q, i+norb_left, j+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_left);
                rdm2_target->element(i+norb_left, q, j+norb_left, p) = value;
                rdm2_target->element(q, i+norb_left, p, j+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_left);
                rdm2_target->element(i+norb_left, q, p, j+norb_left) = value;
                rdm2_target->element(q, i+norb_left, j+norb_left, p) = value;
              }
            }
          }
        }
      }

    } // end of operation list
  } // end of looping over nstates
} // end of compute_220


void ASD_DMRG::compute_rdm2_022(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_022" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>, bool, bool, bool>>> gammalist_tuple_list_list = { 
    // { {ops on site}, {ops on right}, {right nele change}, trans_site, trans_right, duplicate }
    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false },
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false }
    },

    {
      { {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha,  GammaSQ::AnnihilateAlpha},  { 2, 0}, false, true,  false },
      { {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha},  {GammaSQ::AnnihilateBeta,   GammaSQ::AnnihilateAlpha},  { 1, 1}, false, true,  true  },
      { {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta},   {GammaSQ::AnnihilateBeta,   GammaSQ::AnnihilateBeta},   { 0, 2}, false, true,  false }
    },

    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},       { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},       {-1, 1}, true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},        { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},       { 1,-1}, false, true,  false }
    }
  };

  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {
        const bool trans_site = get<3>(gammalist_tuple);
        const bool trans_right = get<4>(gammalist_tuple);
        const bool duplicate = get<5>(gammalist_tuple);
        for (auto& isec : prod_civec->sectors()) {
          BlockKey ket_seckey = isec.first;
          for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
            BlockKey ket_rightkey = ketpair.right;
            const int ket_rightnstates = ketpair.right.nstates;
            const int ketpairoffset = ketpair.offset;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple).first, ket_rightkey.neleb + get<2>(gammalist_tuple).second);
            if (!right_block->contains(bra_rightkey)) continue;
            const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
            // left block info
            BlockKey leftkey = ketpair.left;
            const int leftnstates = ketpair.left.nstates;
            // bra info
            BlockKey bra_seckey(leftkey.nelea + bra_rightkey.nelea, leftkey.neleb + bra_rightkey.neleb);
            if (!prod_civec->contains_block(bra_seckey)) continue;
            auto brapair = doubleblock->blockpairs(bra_seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
              { return make_pair(left_block->blockinfo(leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;

            // site transition tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              map<BlockKey, shared_ptr<const RASDvec>> ket_states;
              vector<shared_ptr<RASCivec>> ket_vecs;
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ileft = 0; ileft != leftnstates; ++ileft)
                  ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ileft + ket_ir*leftnstates)));
              const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
              const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;
              BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
              ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);
  
              map<BlockKey, shared_ptr<const RASDvec>> bra_states;
              vector<shared_ptr<RASCivec>> bra_vecs;
              for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
                for (int ileft = 0; ileft != leftnstates; ++ileft)
                  bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + ileft + bra_ir*leftnstates)));
              const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
              const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
              BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
              bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);
  
              if (trans_site) swap(bra_states, ket_states);
  
              GammaForestASD2<RASDvec> forest(bra_states, ket_states);
              forest.compute();
  
              const size_t bra_rastag = forest.block_tag(bra_raskey);
              const size_t ket_rastag = forest.block_tag(ket_raskey);
              const size_t bratag = trans_site ? ket_rastag : bra_rastag;
              const size_t kettag = trans_site ? bra_rastag : ket_rastag;
              assert(forest.template exist<0>(bratag, kettag, get<0>(gammalist_tuple)));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bratag, kettag, get<0>(gammalist_tuple));
              const int dim1 = trans_site ? ket_rightnstates : bra_rightnstates;
              const int dim2 = trans_site ? bra_rightnstates : ket_rightnstates;
              btas::CRange<3> range1(leftnstates*leftnstates, dim1, dim2*lrint(pow(norb_site, 2)));
              auto tmp_tensor = make_shared<btas::Tensor3<double>>(range1, transition_mat->storage());
              vector<double> buf1(leftnstates*leftnstates*dim1);
              for (int i = 0; i != tmp_tensor->extent(2); ++i) {
                copy_n(&(*tmp_tensor)(0,0,i), buf1.size(), buf1.data());
                blas::transpose(buf1.data(), leftnstates*dim1, leftnstates, &(*tmp_tensor)(0,0,i));
              } 
  
              btas::CRange<3> site_range(leftnstates*leftnstates, bra_rightnstates*ket_rightnstates, lrint(pow(norb_site, 2)));
              site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(tmp_tensor->storage()));
            } // trans_site ? L-L'-R'-R : L'-L-R-R'
  
            // left coupling tensor
            shared_ptr<Matrix> left_coupling_tensor;
            {
              auto eye = make_shared<Matrix>(leftnstates, leftnstates);
              eye->unit();
              left_coupling_tensor = eye;
            }
  
            // right coupling tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              const int dim1 = trans_site ? ket_rightnstates : bra_rightnstates;
              const int dim2 = trans_site ? bra_rightnstates : ket_rightnstates;
              btas::CRange<3> right_range(dim1, dim2, lrint(pow(norb_right, 2)));
              BlockKey brakey = trans_right ? ket_rightkey : bra_rightkey;
              BlockKey ketkey = trans_right ? bra_rightkey : ket_rightkey;
              shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
              if (trans_site^trans_right) {
                vector<double> right_buf(bra_rightnstates*ket_rightnstates);
                for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*right_coupling_tensor)(0,0,i), right_buf.size(), right_buf.data());
                  blas::transpose(right_buf.data(), right_coupling_tensor->extent(1), right_coupling_tensor->extent(0), &(*right_coupling_tensor)(0,0,i));
                }
              }
            } // trans_site ? R'-R : R-R'
  
            // contraction
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1), site_transition_tensor->extent(2));
            auto gmat = group(*contract_site_left,0,2);
            contract(1.0, group(*site_transition_tensor,1,3), {1,0}, group(*left_coupling_tensor,0,2), {1}, 0.0, gmat, {0});
            auto tmp_mat = make_shared<Matrix>(site_transition_tensor->extent(2), right_coupling_tensor->extent(2));
            contract(1.0, *contract_site_left, {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

            // swap indices
            if (trans_site) {
              vector<double> buf2(norb_site*norb_site);
              for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
                copy_n(tmp_mat->element_ptr(0,i), buf2.size(), buf2.data());
                blas::transpose(buf2.data(), norb_site, norb_site, tmp_mat->element_ptr(0,i));
              }
            }
            if (trans_right) {
              auto tmp2 = tmp_mat->clone();
              for (int i = 0; i != norb_right; ++i)
                for (int j = 0; j != norb_right; ++j)
                  copy_n(tmp_mat->element_ptr(0, j+i*norb_right), tmp_mat->ndim(), tmp2->element_ptr(0, i+j*norb_right));
              tmp_mat = tmp2;
            }
            if (duplicate) {
              auto dup_mat = tmp_mat->clone();
              for (int i = 0; i != norb_right; ++i)
                for (int j = 0; j != norb_right; ++j)
                  copy_n(tmp_mat->element_ptr(0, j+i*norb_right), tmp_mat->ndim(), dup_mat->element_ptr(0, i+j*norb_right));

              vector<double> buf3(norb_site*norb_site);
              for (int k = 0; k != dup_mat->mdim(); ++k) {
                copy_n(dup_mat->element_ptr(0,k), buf3.size(), buf3.data());
                blas::transpose(buf3.data(), norb_site, norb_site, dup_mat->element_ptr(0,k));
              }

              *tmp_mat += *dup_mat;
            }
            
            // copy this bra-ket pair density matrix to target
            const double sign = (scheme == 3) ? -1.0 : 1.0;
            blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), rdm_mat->data());
          }
        }
      }
  
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int q = 0; q != norb_right; ++q) {
          for (int p = 0; p != norb_right; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_right);
                rdm2_target->element(i+norb_left, j+norb_left, p+rightoffset, q+rightoffset) = value;
                rdm2_target->element(p+rightoffset, q+rightoffset, i+norb_left, j+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int q = 0; q != norb_right; ++q) {
          for (int p = 0; p != norb_right; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_right);
                rdm2_target->element(i+norb_left, q+rightoffset, j+norb_left, p+rightoffset) = value;
                rdm2_target->element(q+rightoffset, i+norb_left, p+rightoffset, j+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int q = 0; q != norb_right; ++q) {
          for (int p = 0; p != norb_right; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_right);
                rdm2_target->element(i+norb_left, q+rightoffset, p+rightoffset, j+norb_left) = value;
                rdm2_target->element(q+rightoffset, i+norb_left, j+norb_left, p+rightoffset) = value;
              }
            }
          }
        }
      }

    } // end of looping over operation list
  } // end of looping over nstates
} // end of compute_022


void ASD_DMRG::compute_rdm2_202(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_202" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, pair<int, int>, bool, bool, bool>>> gammalist_tuple_list_list = { 
    // { {ops on left}, {ops on right}, {left nele change}, trans_left, trans_right, duplicate }
    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     { 0, 0}, false, false, false },
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     { 0, 0}, false, false, false },
    },

    {
      { {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {-2, 0}, false, true,  false },
      { {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {-1,-1}, false, true,  true  },
      { {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  { 0,-2}, false, true,  false },
    },

    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     { 1,-1}, true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      { 0, 0}, false, false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {-1, 1}, false, true,  false }
    }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site =  multisite_->active_sizes().at(nsites_-2); 
    const int rightoffset = norb_left + norb_site;

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_right*norb_right, norb_left*norb_left); // matrix to store RDM, use ax_plus_y...
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {
        const bool trans_left = get<3>(gammalist_tuple);
        const bool trans_right = get<4>(gammalist_tuple);
        const bool duplicate = get<5>(gammalist_tuple);
        for (auto& isec : prod_civec->sectors()) {
          BlockKey seckey = isec.first;
          for (auto& ketpair : doubleblock->blockpairs(seckey)) {
            const int ketpairoffset = ketpair.offset;
            // left bra-ket pair
            BlockKey ket_leftkey = ketpair.left.key();
            const int ket_leftnstates = ketpair.left.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple).first, ket_leftkey.neleb + get<2>(gammalist_tuple).second);
            if (!left_block->contains(bra_leftkey)) continue;
            const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
            // right bra-ket pair
            BlockKey ket_rightkey = ketpair.right.key();
            const int ket_rightnstates = ketpair.right.nstates;
            BlockKey bra_rightkey(seckey.nelea - bra_leftkey.nelea, seckey.neleb - bra_leftkey.neleb);
            if (!right_block->contains(bra_rightkey)) continue;
            const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
            auto brapair = doubleblock->blockpairs(seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
              { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;
  
            // site transition matrix
            shared_ptr<Matrix> site_transition_tensor;
            {
              site_transition_tensor = make_shared<Matrix>(bra_leftnstates*ket_leftnstates, bra_rightnstates*ket_rightnstates);
              double* mat_ptr = site_transition_tensor->data();
              auto ci_ptr = prod_civec->sector(seckey);
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
                  for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                    for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                      *mat_ptr++ = ci_ptr->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates).dot_product(ci_ptr->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates));
            } // L-L'-R-R'
  
            // left coupling tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              BlockKey brakey = trans_left ? ket_leftkey : bra_leftkey;
              BlockKey ketkey = trans_left ? bra_leftkey : ket_leftkey;
              shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              btas::CRange<3> left_range(bra_leftnstates, ket_leftnstates, norb_left*norb_left);
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
              if (trans_left) {
                vector<double> buf(bra_leftnstates*ket_leftnstates);
                for (int i = 0; i != left_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*left_coupling_tensor)(0,0,i), bra_leftnstates*ket_leftnstates, buf.data());
                  blas::transpose(buf.data(), ket_leftnstates, bra_leftnstates, &(*left_coupling_tensor)(0,0,i));
                }
              }
            } // L-L'
            
            // right coupling tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              BlockKey brakey = trans_right ? ket_rightkey : bra_rightkey;
              BlockKey ketkey = trans_right ? bra_rightkey : ket_rightkey;
              shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              btas::CRange<3> range(bra_rightnstates, ket_rightnstates, norb_right*norb_right);
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(range, right_coupling->storage());
              if (trans_right) {
                vector<double> buf(bra_rightnstates*ket_rightnstates);
                for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*right_coupling_tensor)(0,0,i), bra_rightnstates*ket_rightnstates, buf.data());
                  blas::transpose(buf.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,i));
                }
              }
            } // R-R'
  
            // contraction
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->mdim(), left_coupling_tensor->extent(2));
            contract(1.0, *site_transition_tensor, {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
            auto tmp_mat = make_shared<Matrix>(right_coupling_tensor->extent(2), left_coupling_tensor->extent(2));
            contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *contract_site_left, {2,1}, 0.0, *tmp_mat, {0,1});

            // swap indices
            if (trans_left) {
              auto tmp = tmp_mat->clone();
              for (int i = 0; i != norb_left; ++i)
                for (int j = 0; j != norb_left; ++j)
                  copy_n(tmp_mat->element_ptr(0, j+i*norb_left), tmp_mat->ndim(), tmp->element_ptr(0, i+j*norb_left));
              tmp_mat = tmp;
            }
            if (trans_right) {
              vector<double> buf(norb_right*norb_right);
              for (int i = 0; i != tmp_mat->mdim(); ++i) {
                copy_n(tmp_mat->element_ptr(0, i), buf.size(), buf.data());
                blas::transpose(buf.data(), norb_right, norb_right, tmp_mat->element_ptr(0, i));
              }
            }
            if (duplicate){
              auto dup_mat = tmp_mat->clone();
              for (int i = 0; i != norb_left; ++i)
                for (int j = 0; j != norb_left; ++j)
                  copy_n(tmp_mat->element_ptr(0, j+i*norb_left), tmp_mat->ndim(), dup_mat->element_ptr(0, i+j*norb_left));
              
              vector<double> buf3(norb_right*norb_right);
              for (int i = 0; i != dup_mat->mdim(); ++i) {
                copy_n(dup_mat->element_ptr(0, i), buf3.size(), buf3.data());
                blas::transpose(buf3.data(), norb_right, norb_right, dup_mat->element_ptr(0, i));
              }
              
              *tmp_mat += *dup_mat;
            }

            // copy to target
            const double sign = (scheme == 3) ? -1.0 : 1.0;
            blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), unordered_rdm->data());
          }
        }
      }
      
      // redorder left orbitals
      for (int index2 = 0; index2 != norb_left; ++index2) {
        for (int index1 = 0; index1 != norb_left; ++index1) {
          const int unordered_index = index1 + index2*norb_left;
          const int target_idx1 = left_block->left_index(index1);
          const int target_idx2 = left_block->left_index(index2);
          const int target = target_idx1 + target_idx2*norb_left;
          copy_n(unordered_rdm->element_ptr(0,unordered_index), rdm_mat->ndim(), rdm_mat->element_ptr(0,target));
        }
      }
      
  
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int j = 0; j != norb_left; ++j) {
          for (int i = 0; i != norb_left; ++i) {
            for (int q = 0; q != norb_right; ++q) {
              for (int p = 0; p != norb_right; ++p) {
                const double value = *rdm_mat->element_ptr(p+q*norb_right, i+j*norb_left);
                rdm2_target->element(i, j, p+rightoffset, q+rightoffset) = value;
                rdm2_target->element(p+rightoffset, q+rightoffset, i, j) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int j = 0; j != norb_left; ++j) {
          for (int i = 0; i != norb_left; ++i) {
            for (int q = 0; q != norb_right; ++q) {
              for (int p = 0; p != norb_right; ++p) {
                const double value = *rdm_mat->element_ptr(p+q*norb_right, i+j*norb_left);
                rdm2_target->element(i, q+rightoffset, j, p+rightoffset) = value;
                rdm2_target->element(q+rightoffset, i, p+rightoffset, j) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int j = 0; j != norb_left; ++j) {
          for (int i = 0; i != norb_left; ++i) {
            for (int q = 0; q != norb_right; ++q) {
              for (int p = 0; p != norb_right; ++p) {
                const double value = *rdm_mat->element_ptr(p+q*norb_right, i+j*norb_left);
                rdm2_target->element(i, q+rightoffset, p+rightoffset, j) = value;
                rdm2_target->element(q+rightoffset, i, j, p+rightoffset) = value;
              }
            }
          }
        }
      }
  
    } // end of looping over operation list
  } // end of looping over nstates
} // end of compute_202


void ASD_DMRG::compute_rdm2_121(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_121" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, list<GammaSQ>, pair<int, int>, pair<int, int>, bool, bool, bool, bool>>> gammalist_tuple_list_list = {
    // { {ops on site}, {ops on left}, {ops on right}, {left nele change}, {right nele change}, trans_site, trans_left, trans_right, swap_site }
    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, { 1, 0}, {-1, 0}, false, false, true,  false },
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  { 0, 1}, { 0,-1}, false, false, true,  false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  { 0, 1}, { 0,-1}, false, false, true,  false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, { 1, 0}, {-1, 0}, false, false, true,  false }
    },

    {
      { {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, { 1, 0}, { 1, 0}, false, false, false, false },
      { {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha}, { 0, 1}, { 1, 0}, false, false, false, true  },
      { {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  { 0, 1}, { 0, 1}, false, false, false, false },
      { {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta},  { 1, 0}, { 0, 1}, false, false, false, false }
    },
    
    {
      { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {-1, 0}, { 1, 0}, false, true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta},  {-1, 0}, { 0, 1}, true,  true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  { 0,-1}, { 0, 1}, false, true,  false, false },
      { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha}, { 0,-1}, { 1, 0}, false, true,  false, false }
    }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site = multisite_->active_sizes().at(site);
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {   
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_right);
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {   
        const bool trans_site = get<5>(gammalist_tuple);
        const bool trans_left = get<6>(gammalist_tuple);
        const bool trans_right = get<7>(gammalist_tuple);
        const bool swap_site = get<8>(gammalist_tuple);
        // loop over product rasci sectors
        for (auto& isec : prod_civec->sectors()) {
          BlockKey ket_seckey = isec.first;
          for (auto& bpair : doubleblock->blockpairs(ket_seckey)) {
            const int ketpairoffset = bpair.offset;
            // left bra-ket pair
            BlockKey ket_leftkey = bpair.left.key();
            const int ket_leftnstates = bpair.left.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<3>(gammalist_tuple).first, ket_leftkey.neleb + get<3>(gammalist_tuple).second);
            if (!left_block->contains(bra_leftkey)) continue;
            const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
            // right bra-ket pair
            BlockKey ket_rightkey = bpair.right.key();
            const int ket_rightnstates = bpair.right.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<4>(gammalist_tuple).first, ket_rightkey.neleb + get<4>(gammalist_tuple).second);
            if (!right_block->contains(bra_rightkey)) continue;
            const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
            BlockKey bra_seckey(bra_leftkey.nelea+bra_rightkey.nelea, bra_leftkey.neleb+bra_rightkey.neleb);
            if (!prod_civec->contains_block(bra_seckey)) continue;
            auto brapair = doubleblock->blockpairs(bra_seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp) 
              { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;
  
            // site transition density tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              map<BlockKey, shared_ptr<const RASDvec>> ket_states;
              vector<shared_ptr<RASCivec>> ket_vecs;
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates)));
              const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
              const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;
              BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
              ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);

              map<BlockKey, shared_ptr<const RASDvec>> bra_states;
              vector<shared_ptr<RASCivec>> bra_vecs;
              for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
                for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                  bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates)));
              const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
              const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
              BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
              bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);

              if (trans_site) swap(bra_states, ket_states);
    
              // construct and compute GammaForest for site
              GammaForestASD2<RASDvec> forest(bra_states, ket_states);
              forest.compute();
    
              const size_t bra_rastag = forest.block_tag(bra_raskey);
              const size_t ket_rastag = forest.block_tag(ket_raskey);
              const size_t bratag = trans_site ? ket_rastag : bra_rastag;
              const size_t kettag = trans_site ? bra_rastag : ket_rastag;
              assert(forest.template exist<0>(bratag, kettag, get<0>(gammalist_tuple)));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bratag, kettag, get<0>(gammalist_tuple));
              const int dim1 = trans_site ? ket_leftnstates : bra_leftnstates;
              const int dim2 = trans_site ? bra_leftnstates : ket_leftnstates;
              const int dim3 = trans_site ? ket_rightnstates : bra_rightnstates;
              const int dim4 = trans_site ? bra_rightnstates : ket_rightnstates;
              btas::CRange<3> tmprange(dim1*dim2, dim3, dim4*norb_site*norb_site);
              auto tmp_transition_tensor = make_shared<btas::Tensor3<double>>(tmprange, transition_mat->storage()); 
   
              vector<double> buf1(dim1*dim2*dim3);
              for (int i = 0; i != tmp_transition_tensor->extent(2); ++i) {
                copy_n(&(*tmp_transition_tensor)(0,0,i), dim1*dim2*dim3, buf1.data());
                blas::transpose(buf1.data(), dim1*dim3, dim2, &(*tmp_transition_tensor)(0,0,i));
              } 
  
              btas::CRange<3> siterange(dim1*dim2, dim3*dim4, norb_site*norb_site);
              site_transition_tensor = make_shared<btas::Tensor3<double>>(siterange, move(tmp_transition_tensor->storage()));
            } // trans_site ? L-L'-R'-R : L'-L-R-R'
  
            // left coupling tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              const int dim1 = trans_site ? bra_leftnstates : ket_leftnstates;
              const int dim2 = trans_site ? ket_leftnstates : bra_leftnstates;
              btas::CRange<3> left_range(dim1, dim2, norb_left);
              BlockKey brakey = trans_left ? ket_leftkey : bra_leftkey;
              BlockKey ketkey = trans_left ? bra_leftkey : ket_leftkey;
              shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
              if (!(trans_site^trans_left)) {
                vector<double> buf2(dim1*dim2);
                for (int i = 0; i != norb_left; ++i) {
                  copy_n(&(*left_coupling_tensor)(0,0,i), dim1*dim2, buf2.data());
                  blas::transpose(buf2.data(), left_coupling_tensor->extent(1), left_coupling_tensor->extent(0), &(*left_coupling_tensor)(0,0,i));
                }
              }
            } // trans_site ? L-L' : L'-L
            
            // right coupling tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              const int dim1 = trans_site ? ket_rightnstates : bra_rightnstates;
              const int dim2 = trans_site ? bra_rightnstates : ket_rightnstates;
              btas::CRange<3> right_range(dim1, dim2, norb_right);
              BlockKey brakey = trans_right ? ket_rightkey : bra_rightkey;
              BlockKey ketkey = trans_right ? bra_rightkey : ket_rightkey;
              shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<2>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
              if (trans_site^trans_right) {
                vector<double> buf3(dim1*dim2);
                for (int j = 0; j != norb_right; ++j) {
                  copy_n(&(*right_coupling_tensor)(0,0,j), dim1*dim2, buf3.data());
                  blas::transpose(buf3.data(), right_coupling_tensor->extent(1), right_coupling_tensor->extent(0), &(*right_coupling_tensor)(0,0,j));
                }
              }
            } // trans_site ? R'-R : R-R'
            
            // contraction
            assert(rdm_mat->size() == site_transition_tensor->extent(2) * left_coupling_tensor->extent(2) * right_coupling_tensor->extent(2));
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
            btas::CRange<3> tmprange(site_transition_tensor->extent(1), site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            auto intermediate_tensor = make_shared<const btas::Tensor3<double>>(tmprange, move(contract_site_left->storage()));
            auto tmp_mat = make_shared<Matrix>(site_transition_tensor->extent(2)*left_coupling_tensor->extent(2), right_coupling_tensor->extent(2));
            contract(1.0, group(*intermediate_tensor,1,3), {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *tmp_mat, {0,1});
            auto tmp_rdm = rdm_mat->clone();
            copy_n(tmp_mat->data(), tmp_mat->size(), tmp_rdm->data());

            // swap indices
            if (trans_site^swap_site) {
              vector<double> buf4(tmp_rdm->ndim());
              for (int i = 0; i != tmp_rdm->mdim(); ++i) {
                copy_n(tmp_rdm->element_ptr(0,i), buf4.size(), buf4.data());
                blas::transpose(buf4.data(), norb_site, norb_site, tmp_rdm->element_ptr(0,i));
              }
            }
            
            // copy to target
            const double sign = static_cast<double>(1 - (((ket_leftkey.nelea + ket_leftkey.neleb) & 1) << 1)) * (1.0 - 2.0*static_cast<double>(swap_site));
            blas::ax_plus_y_n(sign, tmp_rdm->data(), tmp_rdm->size(), unordered_rdm->data());
          } 
        } 
      } 

      // reorder left orbitals
      for (int ir = 0; ir != norb_right; ++ir) {
        for (int index = 0; index != norb_left; ++index) {
          const int target = left_block->left_index(index);
          copy_n(unordered_rdm->element_ptr(0, index+ir*norb_left), rdm_mat->ndim(), rdm_mat->element_ptr(0, target+ir*norb_left));
        }
      }

      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int t = 0; t != norb_right; ++t) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+t*norb_left);
                rdm2_target->element(i+norb_left, j+norb_left, p, t+rightoffset) = value;
                rdm2_target->element(p, t+rightoffset, i+norb_left, j+norb_left) = value;
                rdm2_target->element(j+norb_left, i+norb_left, t+rightoffset, p) = value;
                rdm2_target->element(t+rightoffset, p, j+norb_left, i+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int t = 0; t != norb_right; ++t) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+t*norb_left);
                rdm2_target->element(t+rightoffset, i+norb_left, p, j+norb_left) = value;
                rdm2_target->element(p, j+norb_left, t+rightoffset, i+norb_left) = value;
                rdm2_target->element(i+norb_left, t+rightoffset, j+norb_left, p) = value;
                rdm2_target->element(j+norb_left, p, i+norb_left, t+rightoffset) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int t = 0; t != norb_right; ++t) {
          for (int p = 0; p != norb_left; ++p) {
            for (int j = 0; j != norb_site; ++j) {
              for (int i = 0; i != norb_site; ++i) {
                const double value = *rdm_mat->element_ptr(i+j*norb_site, p+t*norb_left);
                rdm2_target->element(i+norb_left, p, t+rightoffset, j+norb_left) = value;
                rdm2_target->element(t+rightoffset, j+norb_left, i+norb_left, p) = value;
                rdm2_target->element(p, i+norb_left, j+norb_left, t+rightoffset) = value;
                rdm2_target->element(j+norb_left, t+rightoffset, p, i+norb_left) = value;
              }
            }
          }
        }
      }

    } // end of looping over operation list
  } // end of looping over istate
} // end of compute_121


void ASD_DMRG::compute_rdm2_211(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_211" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, list<GammaSQ>, pair<int, int>, pair<int, int>, bool, bool, bool>>> gammalist_tuple_list_list = {
    // { {ops on site}, {ops on left}, {ops on right}, {left nele change}, {right nele change}, trans_left, trans_right, swap_left }
    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, { 0, 0}, {-1, 0}, false, true, false },
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateAlpha}, { 0, 0}, {-1, 0}, false, true, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta},  { 0, 0}, { 0,-1}, false, true, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta},  { 0, 0}, { 0,-1}, false, true, false }
    },

    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {-2, 0}, { 1, 0}, false, false, false },
      { {GammaSQ::CreateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  {-1,-1}, { 0, 1}, false, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  { 0,-2}, { 0, 1}, false, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {-1,-1}, { 1, 0}, false, false, true  }
    },
    
    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, { 0, 0}, {-1, 0}, false, true,  false },
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateBeta},  {-1, 1}, { 0,-1}, false, true,  false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {GammaSQ::CreateBeta},  { 0, 0}, { 0,-1}, false, true,  false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {GammaSQ::CreateAlpha}, { 1,-1}, {-1, 0}, true,  true,  false }
    }
  };
    
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site = multisite_->active_sizes().at(site);
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_right*norb_site, norb_left*norb_left);
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {
        const bool trans_left = get<5>(gammalist_tuple);
        const bool trans_right = get<6>(gammalist_tuple);
        const bool swap_left = get<7>(gammalist_tuple);
        // loop over product rasci sectors
        for (auto& isec : prod_civec->sectors()) {
          BlockKey ket_seckey = isec.first;
          for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
            const int ketpairoffset = ketpair.offset;
            // left bra-ket pair
            BlockKey ket_leftkey = ketpair.left.key();
            const int ket_leftnstates = ketpair.left.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<3>(gammalist_tuple).first, ket_leftkey.neleb + get<3>(gammalist_tuple).second);
            if (!left_block->contains(bra_leftkey)) continue;
            const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
            // right bra-ket pair
            BlockKey ket_rightkey = ketpair.right.key();
            const int ket_rightnstates = ketpair.right.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<4>(gammalist_tuple).first, ket_rightkey.neleb + get<4>(gammalist_tuple).second);
            if (!right_block->contains(bra_rightkey)) continue;
            const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
            BlockKey bra_seckey(bra_leftkey.nelea+bra_rightkey.nelea, bra_leftkey.neleb+bra_rightkey.neleb);
            if (!prod_civec->contains_block(bra_seckey)) continue;
            auto brapair = doubleblock->blockpairs(bra_seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
              { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;
            const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
            const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;
            
            // site transition density tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              map<BlockKey, shared_ptr<const RASDvec>> ket_states;
              vector<shared_ptr<RASCivec>> ket_vecs;
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates)));
              BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
              ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);
              map<BlockKey, shared_ptr<const RASDvec>> bra_states;
              vector<shared_ptr<RASCivec>> bra_vecs;
              for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
                for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                  bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates)));
              BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
              bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);
  
              // construct and compute GammaForest for site
              GammaForestASD2<RASDvec> forest(bra_states, ket_states);
              forest.compute();
    
              const size_t bra_rastag = forest.block_tag(bra_raskey);
              const size_t ket_rastag = forest.block_tag(ket_raskey);
              assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              btas::CRange<3> tmprange(ket_leftnstates*bra_leftnstates, bra_rightnstates, ket_rightnstates*norb_site);
              auto tmp_transition_tensor = make_shared<btas::Tensor3<double>>(tmprange, transition_mat->storage()); 
              vector<double> buf1(ket_leftnstates*bra_leftnstates*bra_rightnstates);
              for (int i = 0; i != tmp_transition_tensor->extent(2); ++i) {
                copy_n(&(*tmp_transition_tensor)(0,0,i),buf1.size(), buf1.data());
                blas::transpose(buf1.data(), bra_leftnstates*bra_rightnstates, ket_leftnstates, &(*tmp_transition_tensor)(0,0,i));
              }
  
              btas::CRange<3> siterange(ket_leftnstates*bra_leftnstates, bra_rightnstates*ket_rightnstates, norb_site);
              site_transition_tensor = make_shared<btas::Tensor3<double>>(siterange, move(tmp_transition_tensor->storage()));
            } // L'-L-R-R'
  
            // transposed left coupling tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              btas::CRange<3> left_range(ket_leftnstates, bra_leftnstates, norb_left*norb_left);
              BlockKey brakey = trans_left ? ket_leftkey : bra_leftkey;
              BlockKey ketkey = trans_left ? bra_leftkey : ket_leftkey;
              shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
              if (!trans_left) {
                vector<double> buf2(bra_leftnstates*ket_leftnstates);
                for (int j = 0; j != left_coupling->extent(2); ++j) {
                  copy_n(&(*left_coupling_tensor)(0,0,j), bra_leftnstates*ket_leftnstates, buf2.data());
                  blas::transpose(buf2.data(), bra_leftnstates, ket_leftnstates, &(*left_coupling_tensor)(0,0,j));
                }
              } // L'-L
            }
            
            // right coupling tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              btas::CRange<3> right_range(bra_rightnstates, ket_rightnstates, norb_right);
              BlockKey brakey = trans_right ? ket_rightkey : bra_rightkey;
              BlockKey ketkey = trans_right ? bra_rightkey : ket_rightkey;
              shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<2>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
              if (trans_right) {
                vector<double> buf3(bra_rightnstates*ket_rightnstates);
                for (int j = 0; j != right_coupling->extent(2); ++j) {
                  copy_n(&(*right_coupling_tensor)(0,0,j), buf3.size(), buf3.data());
                  blas::transpose(buf3.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,j));
                }
              }
            } // R-R'
            
            // contraction
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
            auto intermediate_mat = make_shared<Matrix>(site_transition_tensor->extent(1), site_transition_tensor->extent(2)*left_coupling_tensor->extent(2));
            copy_n(contract_site_left->data(), contract_site_left->size(), intermediate_mat->data());

            auto tmp_mat = make_shared<Matrix>(right_coupling_tensor->extent(2), intermediate_mat->mdim());
            contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *intermediate_mat, {2,1}, 0.0, *tmp_mat, {0,1});
            auto resize_mat = rdm_mat->clone();
            copy_n(tmp_mat->data(), tmp_mat->size(), resize_mat->data());

            // swap indices
            if (trans_left^swap_left) {
              auto tmp = resize_mat->clone();
              for (int i = 0; i != norb_left; ++i)
                for (int j = 0; j != norb_left; ++j)
                  copy_n(resize_mat->element_ptr(0, j+i*norb_left), resize_mat->ndim(), tmp->element_ptr(0, i+j*norb_left));
              resize_mat = tmp;
            }

            // copy to target
            const int sign_phase = ket_ras_nelea + ket_ras_neleb + ket_leftkey.nelea + ket_leftkey.neleb + ((scheme == 3) ? 1 : 0);
            const double sign = static_cast<double>(1 - ((sign_phase & 1) << 1)) * (1.0 - 2.0*static_cast<double>(swap_left));
            blas::ax_plus_y_n(sign, resize_mat->data(), resize_mat->size(), unordered_rdm->data());
          }
        }
      }

      // reorder left orbitals
      for (int index2 = 0; index2 != norb_left; ++index2) {
        for (int index1 = 0; index1 != norb_left; ++index1) {
          const int unordered_index = index1 + index2*norb_left;
          const int target_idx1 = left_block->left_index(index1);
          const int target_idx2 = left_block->left_index(index2);
          const int target = target_idx1 + target_idx2*norb_left;
          copy_n(unordered_rdm->element_ptr(0, unordered_index), unordered_rdm->ndim(), rdm_mat->element_ptr(0, target));
        }
      }
      
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int i = 0; i != norb_site; ++i) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+i*norb_right, p+q*norb_left);
                rdm2_target->element(i+norb_left, t+rightoffset, p, q) = value;
                rdm2_target->element(p, q, i+norb_left, t+rightoffset) = value;
                rdm2_target->element(t+rightoffset, i+norb_left, q, p) = value;
                rdm2_target->element(q, p, t+rightoffset, i+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int i = 0; i != norb_site; ++i) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+i*norb_right, p+q*norb_left);
                rdm2_target->element(i+norb_left, q, t+rightoffset, p) = value;
                rdm2_target->element(t+rightoffset, p, i+norb_left, q) = value;
                rdm2_target->element(q, i+norb_left, p, t+rightoffset) = value;
                rdm2_target->element(p, t+rightoffset, q, i+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int q = 0; q != norb_left; ++q) {
          for (int p = 0; p != norb_left; ++p) {
            for (int i = 0; i != norb_site; ++i) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+i*norb_right, p+q*norb_left);
                rdm2_target->element(i+norb_left, q, p, t+rightoffset) = value;
                rdm2_target->element(p, t+rightoffset, i+norb_left, q) = value;
                rdm2_target->element(q, i+norb_left, t+rightoffset, p) = value;
                rdm2_target->element(t+rightoffset, p, q, i+norb_left) = value;
              }
            }
          }
        }
      }
    
    } // end of looping over operation list
  } // end of looping over istate
} // end of compute_211


void ASD_DMRG::compute_rdm2_112(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_112" << endl;
  const list<list<tuple<list<GammaSQ>, list<GammaSQ>, list<GammaSQ>, pair<int, int>, pair<int, int>, bool, bool, bool>>> gammalist_tuple_list_list = {
    // { {ops on site}, {ops on left}, {ops on right}, {left nele change}, {right nele change}, trans_left, trans_right, swap_right }
    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {-1, 0}, { 0, 0}, true, false, false },
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      {-1, 0}, { 0, 0}, true, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      { 0,-1}, { 0, 0}, true, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     { 0,-1}, { 0, 0}, true, false, false }
    },

    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, { 1, 0}, {-2, 0}, false, false, false},
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, { 0, 1}, {-1,-1}, false, false, false},
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  { 0, 1}, { 0,-2}, false, false, false},
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, { 1, 0}, {-1,-1}, false, false, true }
    },
    
    {
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},     {-1, 0}, { 0, 0}, true, false, false },
      { {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     { 0,-1}, {-1, 1}, true, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},      { 0,-1}, { 0, 0}, true, false, false },
      { {GammaSQ::CreateBeta},  {GammaSQ::CreateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},     {-1, 0}, { 1,-1}, true, true,  false }
    }
  };

  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    int scheme = 0;
    for (auto& gammalist_tuple_list : gammalist_tuple_list_list) {
      ++scheme;
      auto rdm_mat = make_shared<Matrix>(norb_right*norb_right, norb_site*norb_left);
      auto unordered_rdm = rdm_mat->clone();
      for (auto& gammalist_tuple : gammalist_tuple_list) {
        const bool trans_left = get<5>(gammalist_tuple);
        const bool trans_right = get<6>(gammalist_tuple);
        const bool swap_right = get<7>(gammalist_tuple);
        // loop over product rasci sectors
        for (auto& isec : prod_civec->sectors()) {
          BlockKey ket_seckey = isec.first;
          for (auto& ketpair : doubleblock->blockpairs(ket_seckey)) {
            const int ketpairoffset = ketpair.offset;
            // left bra-ket pair
            BlockKey ket_leftkey = ketpair.left.key();
            const int ket_leftnstates = ketpair.left.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<3>(gammalist_tuple).first, ket_leftkey.neleb + get<3>(gammalist_tuple).second);
            if (!left_block->contains(bra_leftkey)) continue;
            const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
            // right bra-ket pair
            BlockKey ket_rightkey = ketpair.right.key();
            const int ket_rightnstates = ketpair.right.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<4>(gammalist_tuple).first, ket_rightkey.neleb + get<4>(gammalist_tuple).second);
            if (!right_block->contains(bra_rightkey)) continue;
            const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
            BlockKey bra_seckey(bra_leftkey.nelea+bra_rightkey.nelea, bra_leftkey.neleb+bra_rightkey.neleb);
            if (!prod_civec->contains_block(bra_seckey)) continue;
            auto brapair = doubleblock->blockpairs(bra_seckey);
            auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
              { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
            assert(braiter != brapair.end());
            const int brapairoffset = braiter->offset;
            const int bra_ras_nelea = tot_nelea - bra_seckey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_seckey.neleb;
            const int ket_ras_nelea = tot_nelea - ket_seckey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_seckey.neleb;
            
            // site transition density tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              map<BlockKey, shared_ptr<const RASDvec>> ket_states;
              vector<shared_ptr<RASCivec>> ket_vecs;
              for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
                for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                  ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(ket_seckey)->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates)));
              BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
              ket_states[ket_raskey] = make_shared<const RASDvec>(ket_vecs);
              map<BlockKey, shared_ptr<const RASDvec>> bra_states;
              vector<shared_ptr<RASCivec>> bra_vecs;
              for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
                for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                  bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(bra_seckey)->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates)));
              BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);
              bra_states[bra_raskey] = make_shared<const RASDvec>(bra_vecs);
  
              // construct and compute GammaForest for site
              GammaForestASD2<RASDvec> forest(bra_states, ket_states);
              forest.compute();
    
              const size_t bra_rastag = forest.block_tag(bra_raskey);
              const size_t ket_rastag = forest.block_tag(ket_raskey);
              assert(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              btas::CRange<3> tmprange(ket_leftnstates*bra_leftnstates, bra_rightnstates, ket_rightnstates*norb_site);
              auto tmp_transition_tensor = make_shared<btas::Tensor3<double>>(tmprange, transition_mat->storage()); 
              vector<double> buf1(ket_leftnstates*bra_leftnstates*bra_rightnstates);
              for (int i = 0; i != tmp_transition_tensor->extent(2); ++i) {
                copy_n(&(*tmp_transition_tensor)(0,0,i),buf1.size(), buf1.data());
                blas::transpose(buf1.data(), bra_leftnstates*bra_rightnstates, ket_leftnstates, &(*tmp_transition_tensor)(0,0,i));
              }
  
              btas::CRange<3> siterange(ket_leftnstates*bra_leftnstates, bra_rightnstates*ket_rightnstates, norb_site);
              site_transition_tensor = make_shared<btas::Tensor3<double>>(siterange, move(tmp_transition_tensor->storage()));
            } // L'-L-R-R'
  
            // transposed left coupling tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              btas::CRange<3> left_range(ket_leftnstates, bra_leftnstates, norb_left);
              BlockKey brakey = trans_left ? ket_leftkey : bra_leftkey;
              BlockKey ketkey = trans_left ? bra_leftkey : ket_leftkey;
              shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
              if (!trans_left) {
                vector<double> buf2(bra_leftnstates*ket_leftnstates);
                for (int j = 0; j != left_coupling->extent(2); ++j) {
                  copy_n(&(*left_coupling_tensor)(0,0,j), bra_leftnstates*ket_leftnstates, buf2.data());
                  blas::transpose(buf2.data(), bra_leftnstates, ket_leftnstates, &(*left_coupling_tensor)(0,0,j));
                }
              }
            } // L'-L
            
            // right coupling tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              btas::CRange<3> right_range(bra_rightnstates, ket_rightnstates, norb_right*norb_right);
              BlockKey brakey = trans_right ? ket_rightkey : bra_rightkey;
              BlockKey ketkey = trans_right ? bra_rightkey : ket_rightkey;
              shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<2>(gammalist_tuple)).at(make_pair(brakey, ketkey)).data;
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
              if (trans_right) {
                vector<double> buf3(bra_rightnstates*ket_rightnstates);
                for (int j = 0; j != right_coupling->extent(2); ++j) {
                  copy_n(&(*right_coupling_tensor)(0,0,j), buf3.size(), buf3.data());
                  blas::transpose(buf3.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,j));
                }
              }
            } // R-R'
            
            // contraction
            auto contract_site_left = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *contract_site_left, {0,1});
            auto intermediate_mat = make_shared<Matrix>(site_transition_tensor->extent(1), site_transition_tensor->extent(2)*left_coupling_tensor->extent(2));
            copy_n(contract_site_left->data(), contract_site_left->size(), intermediate_mat->data());

            auto tmp_mat = rdm_mat->clone();
            contract(1.0, group(*right_coupling_tensor,0,2), {2,0}, *intermediate_mat, {2,1}, 0.0, *tmp_mat, {0,1});

            // swap indices
            if (trans_right^swap_right) {
              auto tmp = tmp_mat->copy();
              vector<double> buf4(tmp_mat->ndim());
              for (int i = 0; i != tmp->mdim(); ++i) {
                copy_n(tmp->element_ptr(0,i), tmp->ndim(), buf4.data());
                blas::transpose(buf4.data(), norb_right, norb_right, tmp->element_ptr(0,i));
              }
              tmp_mat = tmp;
            }

            // copy to target
            const int sign_phase = ket_ras_nelea + ket_ras_neleb + ((scheme == 3) ? 1 : 0);
            const double sign = static_cast<double>(1 - ((sign_phase & 1) << 1)) * (1.0 - 2.0*static_cast<double>(swap_right));
            blas::ax_plus_y_n(sign, tmp_mat->data(), tmp_mat->size(), unordered_rdm->data());
          }
        }
      }

      // reorder left block orbitals
      for (int index = 0; index != norb_left; ++index) {
        const int target = left_block->left_index(index);
        for (int isite = 0; isite != norb_site; ++isite)
          copy_n(unordered_rdm->element_ptr(0, isite+index*norb_site), unordered_rdm->ndim(), rdm_mat->element_ptr(0,isite+target*norb_site));
      }
        
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      if (scheme == 1) {
        for (int p = 0; p != norb_left; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            for (int u = 0; u != norb_right; ++u) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+u*norb_right, i+p*norb_site);
                rdm2_target->element(i+norb_left, p, t+rightoffset, u+rightoffset) = value;
                rdm2_target->element(t+rightoffset, u+rightoffset, i+norb_left, p) = value;
                rdm2_target->element(p, i+norb_left, u+rightoffset, t+rightoffset) = value;
                rdm2_target->element(u+rightoffset, t+rightoffset, p, i+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 2) {
        for (int p = 0; p != norb_left; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            for (int u = 0; u != norb_right; ++u) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+u*norb_right, i+p*norb_site);
                rdm2_target->element(i+norb_left, u+rightoffset, p, t+rightoffset) = value;
                rdm2_target->element(p, t+rightoffset, i+norb_left, u+rightoffset) = value;
                rdm2_target->element(u+rightoffset, i+norb_left, t+rightoffset, p) = value;
                rdm2_target->element(t+rightoffset, p, u+rightoffset, i+norb_left) = value;
              }
            }
          }
        }
      } else if (scheme == 3) {
        for (int p = 0; p != norb_left; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            for (int u = 0; u != norb_right; ++u) {
              for (int t = 0; t != norb_right; ++t) {
                const double value = *rdm_mat->element_ptr(t+u*norb_right, i+p*norb_site);
                rdm2_target->element(i+norb_left, u+rightoffset, t+rightoffset, p) = value;
                rdm2_target->element(t+rightoffset, p, i+norb_left, u+rightoffset) = value;
                rdm2_target->element(u+rightoffset, i+norb_left, p, t+rightoffset) = value;
                rdm2_target->element(p, t+rightoffset, u+rightoffset, i+norb_left) = value;
              }
            }
          }
        }
      }
    
    } // end of looping over operation list
  } // end of looping over istate
} // end of compute_112


