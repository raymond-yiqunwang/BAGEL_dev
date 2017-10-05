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

    // obtain ProductRASCivec 
    vector<shared_ptr<ProductRASCivec>> cc;
    vector<double> energy;
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
      energy = prod_ras->energy();
      cc = prod_ras->civectors();
    }
    cout << "energy on site " << site << " : " << setw(16) << setprecision(12) << energy[0] << endl;

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
        // compute_rdm2_211(cc, site);

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
        //compute_rdm2_112(cc);
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
  auto fci_input = input_->get_child("fci");
  auto fci = make_shared<KnowlesHandy>(fci_input, multisite_->geom(), multisite_->ref());
  fci->compute();
  auto fcirdm2 = fci->rdm2();

  for (int istate = 0; istate != nstate_; ++istate) {
    cout << "rdm2_[" << istate << "] : " << endl;

    list<int> switchlist {111, 222, 333, 130, 310, 301, 31/*031*/, 13/*013*/, 103, 220, 22/*022*/, 202, 121};
    for (int swch : switchlist) {

      vector<double> tmpvec;
      pair<int, int> range1(0,2), range2(2,4), range3(4,6);
      list<tuple<pair<int, int>, pair<int, int>, pair<int, int>, pair<int, int>>> list_tuplelist;
      switch (swch) {
        // RAS0
        case 111: list_tuplelist = {
          {range1, range1, range1, range1}
        }; break;

        // RAS1
        case 222: list_tuplelist = {
          {range2, range2, range2, range2}
        }; break;

        // RAS2
        case 333: list_tuplelist = {
          {range3, range3, range3, range3}
        }; break;
        
        // 130
        case 130: list_tuplelist = {
          {range1, range2, range2, range2},
          {range2, range1, range2, range2},
          {range2, range2, range1, range2},
          {range2, range2, range2, range1}
        }; break;

        // 310
        case 310: list_tuplelist = {
          {range2, range1, range1, range1},
          {range1, range2, range1, range1},
          {range1, range1, range2, range1},
          {range1, range1, range1, range2}
        }; break;

        case 301: list_tuplelist = {
          {range3, range1, range1, range1},
          {range1, range3, range1, range1},
          {range1, range1, range3, range1},
          {range1, range1, range1, range3}
        }; break;

        case 31/*031*/: list_tuplelist = {
          {range3, range2, range2, range2},
          {range2, range3, range2, range2},
          {range2, range2, range3, range2},
          {range2, range2, range2, range3}
        }; break;

        case 13/*013*/: list_tuplelist = {
          {range2, range3, range3, range3},
          {range3, range2, range3, range3},
          {range3, range3, range2, range3},
          {range3, range3, range3, range2}
        }; break;

        case 103: list_tuplelist = {
          {range1, range3, range3, range3},
          {range3, range1, range3, range3},
          {range3, range3, range1, range3},
          {range3, range3, range3, range1}
        }; break;
        
        case 220: list_tuplelist = {
          {range2, range2, range1, range1}, 
          {range1, range1, range2, range2},
          {range2, range1, range2, range1},
          {range1, range2, range1, range2},
          {range2, range1, range1, range2},
          {range1, range2, range2, range1} 
        }; break;

        case 22/*022*/: list_tuplelist = {
          {range2, range2, range3, range3},
          {range3, range3, range2, range2},
          {range2, range3, range3, range2},
          {range3, range2, range2, range3},
          {range2, range3, range2, range3},
          {range3, range2, range3, range2} 
        }; break;

        case 202: list_tuplelist = {
          {range1, range1, range3, range3},
          {range3, range3, range1, range1},
          {range1, range3, range3, range1},
          {range3, range1, range1, range3},
          {range1, range3, range1, range3},
          {range3, range1, range3, range1}
        }; break;

        case 121: list_tuplelist = {
          {range2, range2, range1, range3},
          {range1, range3, range2, range2},
          {range3, range1, range2, range2},
          {range2, range2, range3, range1}
        }; break;

      }
 
      for (auto& tuple_list : list_tuplelist) {
        for (int i = get<3>(tuple_list).first; i != get<3>(tuple_list).second; ++i) {
          for (int j = get<2>(tuple_list).first; j != get<2>(tuple_list).second; ++j) {
            for (int k = get<1>(tuple_list).first; k != get<1>(tuple_list).second; ++k) {
              for (int l = get<0>(tuple_list).first; l != get<0>(tuple_list).second; ++l) {
                tmpvec.push_back(fcirdm2->at(istate)->element(l,k,j,i) - rdm2_->at(istate)->element(l,k,j,i));
              }
            }
          }
        }
      }
      double sum = 0.0;
      std::for_each(tmpvec.begin(), tmpvec.end(), [&sum] (const double v) { sum += v*v; });
      const double rms = std::sqrt(sum / static_cast<double>(tmpvec.size()));
      cout << "rms : with swich " << setfill('0') << setw(3) << swch << " = " << setfill(' ') << setw(16) << setprecision(12) << rms << ", tmpvec size = " << tmpvec.size() << endl;
    }
  }
#endif

}


// TODO find a better algorithm
void ASD_DMRG::compute_rdm2_ras(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
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
  const int nstate = dvec.size();
  // contains : site operator list, left_block operator list, change in alpha electrons at left_block, change in beta electrons at left_block
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    {{GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, 1, 0}, // <block_bra|a+|block_ket> <ras_bra|a+ a- a-|ras_ket>
    {{GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, 1, 0}, // <block_bra|a+|block_ket> <ras_bra|b+ b- a-|ras_ket>
    {{GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  0, 1}, // <block_bra|b+|block_ket> <ras_bra|a+ a- b-|ras_ket>
    {{GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  0, 1}  // <block_bra|b+|block_ket> <ras_bra|b+ b- b-|ras_ket>
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site*norb_site, norb_left); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // marco loop over right blocks since we have a delta_{r,r'}
      for (auto& rblock : right_block->blocks()) {
        BlockKey rightkey = rblock.key(); 
        const int rightnstates = rblock.nstates;
        // micro loop over one specific right block
        for (int ir = 0; ir != rightnstates; ++ir) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between left block states with the same right block state
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& lblock : left_block->blocks()) {
            BlockKey leftkey = lblock.key();
            const int leftnstates = lblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int ileft = 0; ileft != leftnstates; ++ileft)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + ir*leftnstates + ileft)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density matrices
          for (auto& lbinfo : left_block->blocks()) {
            BlockKey ket_leftkey = lbinfo.key();
            BlockKey ket_combinedkey(ket_leftkey.nelea + rightkey.nelea, ket_leftkey.neleb + rightkey.neleb);
            const int ket_nstates = lbinfo.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple), ket_leftkey.neleb + get<3>(gammalist_tuple));
            if (!(left_block->contains(bra_leftkey))) continue;
            BlockKey bra_combinedkey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = left_block->blockinfo(bra_leftkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple)))) continue;

            shared_ptr<btas::Tensor3<double>> transition_tensor;
            { // transpose the first two dimensions
              btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
              transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != transition_tensor->extent(2); ++i) {
                copy_n(&(*transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*transition_tensor)(0,0,i));
              }
            }
            
            shared_ptr<const btas::Tensor3<double>> coupling_data = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;
            
            auto target = rdm_mat->clone();
            assert (target->size() == transition_tensor->extent(2) * coupling_data->extent(2));
            contract(1.0, group(*transition_tensor,0,2), {2,0}, group(*coupling_data,0,2), {2,1}, 0.0, *target, {0,1});
            const double sign = static_cast<double>(1 - (((bra_ras_nelea + bra_ras_neleb)%2) << 1));
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int k = 0; k != norb_site; ++k) {
      for (int j = 0; j != norb_site; ++j) {
        for (int i = 0; i != norb_site; ++i) {
          for (int p = 0; p != norb_left; ++p) {
            const double value = *rdm_mat->element_ptr(i + j*norb_site + k*norb_site*norb_site, p);
            rdm2_target->element(p, i+norb_left, k+norb_left, j+norb_left) = value;
            rdm2_target->element(i+norb_left, p, j+norb_left, k+norb_left) = value;
            rdm2_target->element(k+norb_left, j+norb_left, p, i+norb_left) = value;
            rdm2_target->element(j+norb_left, k+norb_left, i+norb_left, p) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_130


void ASD_DMRG::compute_rdm2_310(vector<shared_ptr<ProductRASCivec>> dvec) {
  // site == 1
  cout << "  * compute_rdm2_310" << endl;
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = {
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 1, 0},
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  1, 0},
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 1},
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 1}
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(1/*site==1*/);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_left*norb_left*norb_left, norb_site);

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // macro loop over right blocks 
      for (auto& rblock : right_block->blocks()) {
        BlockKey rightkey = rblock.key();
        const int rightnstates = rblock.nstates;
        // micro loop over one specific right block
        for (int ir = 0; ir != rightnstates; ++ir) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between left blocks with the same right block state
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& lblock : left_block->blocks()) {
            BlockKey leftkey = lblock.key();
            const int leftnstates = lblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!prod_civec->contains_block(combinedkey)) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey to raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int ileft = 0; ileft != leftnstates; ++ileft)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + ir*leftnstates + ileft)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density amtrices
          for (auto& lbinfo : left_block->blocks()) {
            BlockKey ket_leftkey = lbinfo.key();
            BlockKey ket_combinedkey(ket_leftkey.nelea + rightkey.nelea, ket_leftkey.neleb + rightkey.neleb);
            const int ket_nstates = lbinfo.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple), ket_leftkey.neleb + get<3>(gammalist_tuple));
            if (!left_block->contains(bra_leftkey)) continue;
            BlockKey bra_combinedkey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!prod_civec->contains_block(bra_combinedkey)) continue;
            const int bra_nstates = left_block->blockinfo(bra_leftkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple)))) continue;

            // transpose the first two dimensions
            shared_ptr<btas::Tensor3<double>> transition_tensor;
            {
              btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
              transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != transition_tensor->extent(2); ++i) {
                copy_n(&(*transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*transition_tensor)(0,0,i));
              }
            }

            shared_ptr<const btas::Tensor3<double>> coupling_data = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;

            auto target = rdm_mat->clone();
            assert(target->size() == transition_tensor->extent(2) * coupling_data->extent(2));
            contract(1.0, group(*coupling_data,0,2), {2,0}, group(*transition_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            const double sign = static_cast<double>(1 - (((bra_ras_nelea + bra_ras_neleb)%2) << 1));
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          }
        
        }
      }
    }

    // copy rdm_mat into rdm2_ TODO make it more efficient
    auto rdm2_target = rdm2_->at(istate);
    for (int i = 0; i != norb_left; ++i) {
      for (int j = 0; j != norb_left; ++j) {
        for (int k = 0; k != norb_left; ++k) {
          for (int p = 0; p != norb_site; ++p) {
            const double value = *rdm_mat->element_ptr(k + j*norb_left + i*norb_left*norb_left, p);
            rdm2_target->element(p+norb_left, k, i, j) = value;
            rdm2_target->element(k, p+norb_left, j, i) = value;
            rdm2_target->element(i, j, p+norb_left, k) = value;
            rdm2_target->element(j, i, k, p+norb_left) = value;
          }
        }
      }
    }
  
  } // end of looping over istate
} // end of compute_310


void ASD_DMRG::compute_rdm2_301(vector<shared_ptr<ProductRASCivec>> dvec) {
  // site == 1
  cout << "  * compute_rdm2_301" << endl;
  const int nstate = dvec.size();
  // contains : site operator list, left_block operator list, change in alpha electrons at left_block, change in beta electrons at left_block
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    {{GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, 1, 0},
    {{GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, 1, 0}, 
    {{GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  0, 1}, 
    {{GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  0, 1} 
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int right_orboffset = norb_left + multisite_->active_sizes().at(1/*site*/);
    auto rdm_mat = make_shared<Matrix>(norb_left*norb_left*norb_left, norb_right); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = bpair.offset;
          // left bra-ket pair
          BlockInfo ket_leftinfo = bpair.left;
          const int lket_nstates = ket_leftinfo.nstates;
          BlockKey bra_leftkey(ket_leftinfo.nelea + get<2>(gammalist_tuple), ket_leftinfo.neleb + get<3>(gammalist_tuple));
          if (!left_block->contains(bra_leftkey)) continue;
          const int lbra_nstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockInfo ket_rightinfo = bpair.right;
          const int rket_nstates = ket_rightinfo.nstates;
          BlockKey bra_rightkey(ket_rightinfo.nelea - get<2>(gammalist_tuple), ket_rightinfo.neleb - get<3>(gammalist_tuple));
          if (!right_block->contains(bra_rightkey)) continue;
          const int rbra_nstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;
          
          shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftinfo.key())).data;
          
          shared_ptr<btas::Tensor3<double>> trans_right_coupling;
          { // transpose first two dimensions of right transition tensor
            shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightinfo.key(), bra_rightkey)).data;
            btas::CRange<3> range(rbra_nstates, rket_nstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
            trans_right_coupling = make_shared<btas::Tensor3<double>>(range, move(right_coupling->storage()));
            unique_ptr<double[]> buf(new double[rbra_nstates*rket_nstates]);
            for (int i = 0; i != trans_right_coupling->extent(2); ++i) {
              copy_n(&(*trans_right_coupling)(0,0,i), rbra_nstates*rket_nstates, buf.get());
              blas::transpose(buf.get(), rket_nstates, rbra_nstates, &(*trans_right_coupling)(0,0,i));
            }
          }

          // contract coefficient tensor
          auto contract_mat = make_shared<Matrix>(lbra_nstates*lket_nstates, rbra_nstates*rket_nstates);
          auto ciptr = prod_civec->sector(seckey);
          for (int irket = 0; irket != rket_nstates; ++irket)
            for (int irbra = 0; irbra != rbra_nstates; ++irbra)
              for (int ilket = 0; ilket != lket_nstates; ++ilket)
                for (int ilbra = 0; ilbra != lbra_nstates; ++ilbra)
                  *contract_mat->element_ptr(ilbra+ilket*lbra_nstates, irbra+irket*rbra_nstates) =
                     ciptr->civec(brapairoffset + ilbra + irbra*lbra_nstates).dot_product(ciptr->civec(ketpairoffset + ilket + irket*lket_nstates));

          auto tmp_mat = make_shared<Matrix>(rbra_nstates*rket_nstates, left_coupling->extent(2));
          contract(1.0, *contract_mat, {2,0}, group(*left_coupling,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

          auto target = rdm_mat->clone();
          contract(1.0, *tmp_mat, {2,0}, group(*trans_right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
          const double sign = static_cast<double>(1 - (((ket_leftinfo.nelea + ket_leftinfo.neleb) % 2) << 1));
          blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
        }
      }
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != norb_right; ++p)
      for (int i = 0; i != norb_left; ++i)
        for (int j = 0; j != norb_left; ++j)
          for (int k = 0; k != norb_left; ++k) {
            const double value = *rdm_mat->element_ptr(k+j*norb_left+i*norb_left*norb_left, p);
            rdm2_target->element(k,p+right_orboffset,j,i) = value;
            rdm2_target->element(p+right_orboffset,k,i,j) = value;
            rdm2_target->element(i,j,p+right_orboffset,k) = value;
            rdm2_target->element(j,i,k,p+right_orboffset) = value;
          }
  
  } // end of looping over nstates
} // end of compute_301


void ASD_DMRG::compute_rdm2_031(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_031" << endl;
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on right}, sitebra_deltaAlpha, sitebra_deltaBeta)
    {{GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, 1, 0},
    {{GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, 1, 0},
    {{GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  0, 1},
    {{GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  0, 1} 
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_right = right_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site*norb_site, norb_right); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // marco loop over left blocks since we have a delta_{l,l'}
      for (auto& lblock : left_block->blocks()) {
        BlockKey leftkey = lblock.key(); 
        const int leftnstates = lblock.nstates;
        // micro loop over one specific left block
        for (int il = 0; il != leftnstates; ++il) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between right block states with the same left block state
          // find all possible right blocks that can be coupled with the left state and compute GammaForestASD
          for (auto& rblock : right_block->blocks()) {
            BlockKey rightkey = rblock.key();
            const int rightnstates = rblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int iright = 0; iright != rightnstates; ++iright)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + il + iright*leftnstates)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over right blocks again to obtain all transition density matrices
          for (auto& rbinfo : right_block->blocks()) {
            BlockKey ket_rightkey = rbinfo.key();
            BlockKey ket_combinedkey(ket_rightkey.nelea + leftkey.nelea, ket_rightkey.neleb + leftkey.neleb);
            if (!(prod_civec->contains_block(ket_combinedkey))) continue;
            const int ket_nstates = rbinfo.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea - get<2>(gammalist_tuple), ket_rightkey.neleb - get<3>(gammalist_tuple));
            if (!(right_block->contains(bra_rightkey))) continue;
            BlockKey bra_combinedkey(bra_rightkey.nelea + leftkey.nelea, bra_rightkey.neleb + leftkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = right_block->blockinfo(bra_rightkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            if (!(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)))) continue;

            // build site transition density matrix
            shared_ptr<btas::Tensor3<double>> transition_tensor;
            {
              btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
            }

            // right block transition density matrix -- transpose the first two dimensions
            shared_ptr<btas::Tensor3<double>> right_coupling;
            {
              shared_ptr<const btas::Tensor3<double>> coupling_data = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;
              btas::CRange<3> right_range(bra_nstates, ket_nstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
              right_coupling = make_shared<btas::Tensor3<double>>(right_range, move(coupling_data->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != right_coupling->extent(2); ++i) {
                copy_n(&(*right_coupling)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*right_coupling)(0,0,i));
              }
            }
 
            auto target = rdm_mat->clone();
            assert (target->size() == transition_tensor->extent(2) * right_coupling->extent(2));
            contract(1.0, group(*transition_tensor,0,2), {2,0}, group(*right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
            const double sign = static_cast<double>(1 - (((ket_ras_nelea + ket_ras_neleb + leftkey.nelea + leftkey.neleb) % 2) << 1));
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != norb_right; ++p) {
      for (int i = 0; i != norb_site; ++i) {
        for (int j = 0; j != norb_site; ++j) {
          for (int k = 0; k != norb_site; ++k) {
            const double value = *rdm_mat->element_ptr(k + j*norb_site + i*norb_site*norb_site, p);
            rdm2_target->element(k+norb_left, p+norb_left+norb_site, j+norb_left, i+norb_left) = value;
            rdm2_target->element(p+norb_left+norb_site, k+norb_left, i+norb_left, j+norb_left) = value;
            rdm2_target->element(i+norb_left, j+norb_left, p+norb_left+norb_site, k+norb_left) = value;
            rdm2_target->element(j+norb_left, i+norb_left, k+norb_left, p+norb_left+norb_site) = value;
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_031


void ASD_DMRG::compute_rdm2_013(vector<shared_ptr<ProductRASCivec>> dvec) {
  // final configuration
  cout << "  * compute_rdm2_013" << endl;
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on right}, sitebra_deltaAlpha, sitebra_deltaBeta}
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 1, 0},
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  1, 0},
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 1},
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 1}
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_right = right_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_right*norb_right*norb_right, norb_site);

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // macro loop over right blocks
      for (auto& lblock : left_block->blocks()) {
        BlockKey leftkey = lblock.key();
        const int leftnstates = lblock.nstates;
        // micro loop over one specific left block
        for (int il = 0; il != leftnstates; ++il) {
          map<BlockKey, shared_ptr<const RASDvec>> states;
          // find all possible right blocks that can be coupled with the left state and compute GammaForestASD
          for (auto& rblock : right_block->blocks()) {
            BlockKey rightkey = rblock.key();
            const int rightnstates = rblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!prod_civec->contains_block(combinedkey)) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert (iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int iright = 0; iright != rightnstates; ++iright)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + il + iright*leftnstates)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density matrices
          for (auto& rbinfo : right_block->blocks()) {
            BlockKey ket_rightkey = rbinfo.key();
            BlockKey ket_combinedkey(ket_rightkey.nelea+leftkey.nelea, ket_rightkey.neleb+leftkey.neleb);
            if (!prod_civec->contains_block(ket_combinedkey)) continue;
            const int ket_nstates = rbinfo.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea - get<2>(gammalist_tuple), ket_rightkey.neleb - get<3>(gammalist_tuple));
            if (!right_block->contains(bra_rightkey)) continue;
            BlockKey bra_combinedkey(bra_rightkey.nelea+leftkey.nelea, bra_rightkey.neleb+leftkey.neleb);
            if (!prod_civec->contains_block(bra_combinedkey)) continue;
            const int bra_nstates = right_block->blockinfo(bra_rightkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t bra_rastag = forest.block_tag(bra_raskey);
            const size_t ket_rastag = forest.block_tag(ket_raskey);
            if (!forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple))) continue;

            // build site transition tensor
            shared_ptr<btas::Tensor3<double>> transition_tensor;
            {
              btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
            }

            // build right block transition density matrix
            shared_ptr<btas::Tensor3<double>> right_coupling;
            {
              shared_ptr<const btas::Tensor3<double>> coupling_data = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightkey, bra_rightkey)).data;
              btas::CRange<3> right_range(bra_nstates, ket_nstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
              right_coupling = make_shared<btas::Tensor3<double>>(right_range, move(coupling_data->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != right_coupling->extent(2); ++i) {
                copy_n(&(*right_coupling)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*right_coupling)(0,0,i));
              }
            }

            auto target = rdm_mat->clone();
            assert (target->size() == transition_tensor->extent(2) * right_coupling->extent(2));
            contract(1.0, group(*right_coupling,0,2), {2,0}, group(*transition_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            const double sign = static_cast<double>(1 - (((ket_ras_nelea + ket_ras_neleb + leftkey.nelea + leftkey.neleb) % 2) << 1));
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != norb_site; ++p) {
      for (int i = 0; i != norb_right; ++i) {
        for (int j = 0; j != norb_right; ++j) {
          for (int k = 0; k != norb_right; ++k) {
            const double value = *rdm_mat->element_ptr(k + j*norb_right + i*norb_right*norb_right, p);
            rdm2_target->element(p+norb_left, k+norb_left+norb_site, i+norb_left+norb_site, j+norb_left+norb_site) = value;
            rdm2_target->element(k+norb_left+norb_site, p+norb_left, j+norb_left+norb_site, i+norb_left+norb_site) = value;
            rdm2_target->element(i+norb_left+norb_site, j+norb_left+norb_site, p+norb_left, k+norb_left+norb_site) = value;
            rdm2_target->element(j+norb_left+norb_site, i+norb_left+norb_site, k+norb_left+norb_site, p+norb_left) = value;
          }
        }
      }
    }
  } // end of looping over istate
} // end of compute_013


void ASD_DMRG::compute_rdm2_103(vector<shared_ptr<ProductRASCivec>> dvec) {
  // site == 1
  cout << "  * compute_rdm2_103" << endl;
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    // ({operators on left}, {operators on right}, leftbra_deltaAlpha, leftbra_deltaBeta}
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 1, 0},
    {{GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  1, 0}, 
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 1}, 
    {{GammaSQ::CreateBeta},  {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 1} 
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int right_orboffset = norb_left + multisite_->active_sizes().at(1/*site*/);
    auto rdm_mat = make_shared<Matrix>(norb_right*norb_right*norb_right, norb_left); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = bpair.offset;
          // left bra-ket pair
          BlockInfo ket_leftinfo = bpair.left;
          const int lket_nstates = ket_leftinfo.nstates;
          BlockKey bra_leftkey(ket_leftinfo.nelea + get<2>(gammalist_tuple), ket_leftinfo.neleb + get<3>(gammalist_tuple));
          if (!left_block->contains(bra_leftkey)) continue;
          const int lbra_nstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockInfo ket_rightinfo = bpair.right;
          const int rket_nstates = ket_rightinfo.nstates;
          BlockKey bra_rightkey(ket_rightinfo.nelea - get<2>(gammalist_tuple), ket_rightinfo.neleb - get<3>(gammalist_tuple));
          if (!right_block->contains(bra_rightkey)) continue;
          const int rbra_nstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;
          
          shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftinfo.key())).data;
          
          shared_ptr<btas::Tensor3<double>> trans_right_coupling;
          {
            shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(ket_rightinfo.key(), bra_rightkey)).data;
            // transpose first two dimensions of right transition tensor
            btas::CRange<3> range(rbra_nstates, rket_nstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
            trans_right_coupling = make_shared<btas::Tensor3<double>>(range, move(right_coupling->storage()));
            unique_ptr<double[]> buf(new double[rbra_nstates*rket_nstates]);
            for (int i = 0; i != trans_right_coupling->extent(2); ++i) {
              copy_n(&(*trans_right_coupling)(0,0,i), rbra_nstates*rket_nstates, buf.get());
              blas::transpose(buf.get(), rket_nstates, rbra_nstates, &(*trans_right_coupling)(0,0,i));
            }
          }

          // contract coefficient tensor
          auto contract_mat = make_shared<Matrix>(lbra_nstates*lket_nstates, rbra_nstates*rket_nstates);
          auto ciptr = prod_civec->sector(seckey);
          for (int irket = 0; irket != rket_nstates; ++irket)
            for (int irbra = 0; irbra != rbra_nstates; ++irbra)
              for (int ilket = 0; ilket != lket_nstates; ++ilket)
                for (int ilbra = 0; ilbra != lbra_nstates; ++ilbra)
                  *contract_mat->element_ptr(ilbra+ilket*lbra_nstates, irbra+irket*rbra_nstates) =
                     ciptr->civec(brapairoffset + ilbra + irbra*lbra_nstates).dot_product(ciptr->civec(ketpairoffset + ilket + irket*lket_nstates));

          auto tmp_mat = make_shared<Matrix>(rbra_nstates*rket_nstates, left_coupling->extent(2));
          contract(1.0, *contract_mat, {2,0}, group(*left_coupling,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

          auto target = rdm_mat->clone()->transpose();
          contract(1.0, *tmp_mat, {2,0}, group(*trans_right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
          const double sign = static_cast<double>(1 - (((ket_leftinfo.nelea + ket_leftinfo.neleb) % 2) << 1));
          blas::ax_plus_y_n(sign, target->transpose()->data(), target->size(), rdm_mat->data());
        }
      }
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != norb_left; ++p) {
      for (int i = 0; i != norb_right; ++i) {
        for (int j = 0; j != norb_right; ++j) {
          for (int k = 0; k != norb_right; ++k) {
            const double value = *rdm_mat->element_ptr(k+j*norb_right+i*norb_right*norb_right, p);
            rdm2_target->element(p, k+right_orboffset, i+right_orboffset, j+right_orboffset) = value;
            rdm2_target->element(k+right_orboffset, p, j+right_orboffset, i+right_orboffset) = value;
            rdm2_target->element(i+right_orboffset, j+right_orboffset, p, k+right_orboffset) = value;
            rdm2_target->element(j+right_orboffset, i+right_orboffset, k+right_orboffset, p) = value;
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_103


void ASD_DMRG::compute_rdm2_220(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "  * compute_rdm2_220" << endl;
  
  // a^{\dagger}_{p \rho} a_{q \rho} on left
  compute_rdm2_220_part1(dvec, site);
  
  // a^{\dagger}_{p \rho} a_{q \sigma} on left
  compute_rdm2_220_part2(dvec, site);
  
  // a^{\dagger \rho} a^{\dagger \sigma} on left
  compute_rdm2_220_part3(dvec, site);
}


// orbital i,j on site, p,q on left block
// \Gamma_{ijpq} = \sum \left< A^{c}_{lr} | a^{\dagger}_{i \sigma} a_{j \sigma} | A^{c'}_{l'r} \right>
//                 \left< l \right| a^{\dagger}_{p \rho} a_{q \rho} \left| l' \right> 
void ASD_DMRG::compute_rdm2_220_part1(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on left}, leftbra_deltaAlpha, leftbra_deltaBeta}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0},
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0}
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_left); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // marco loop over right blocks since we have a delta_{r,r'}
      for (auto& rblock : right_block->blocks()) {
        BlockKey rightkey = rblock.key(); 
        const int rightnstates = rblock.nstates;
        // micro loop over one specific right block
        for (int ir = 0; ir != rightnstates; ++ir) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between left block states with the same right block state
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& lblock : left_block->blocks()) {
            BlockKey leftkey = lblock.key();
            const int leftnstates = lblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int ileft = 0; ileft != leftnstates; ++ileft)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + ir*leftnstates + ileft)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density matrices
          for (auto& lbinfo : left_block->blocks()) {
            BlockKey ket_leftkey = lbinfo.key();
            BlockKey ket_combinedkey(ket_leftkey.nelea + rightkey.nelea, ket_leftkey.neleb + rightkey.neleb);
            const int ket_nstates = lbinfo.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple), ket_leftkey.neleb + get<3>(gammalist_tuple));
            if (!(left_block->contains(bra_leftkey))) continue;
            BlockKey bra_combinedkey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = left_block->blockinfo(bra_leftkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)))) continue;

            shared_ptr<const btas::Tensor3<double>> site_transition_tensor;
            {
              btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              site_transition_tensor = make_shared<const btas::Tensor3<double>>(range, move(transition_mat->storage()));
            }

            shared_ptr<const btas::Tensor3<double>> left_coupling_data = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;
            
            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_tensor->extent(2) * left_coupling_data->extent(2));
            contract(1.0, group(*site_transition_tensor,0,2), {2,0}, group(*left_coupling_data,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_left; ++q) {
      for (int p = 0; p != norb_left; ++p) {
        for (int j = 0; j != norb_site; ++j) {
          for (int i = 0; i != norb_site; ++i) {
            const double value = *rdm_mat->element_ptr(i + j*norb_site, p + q*norb_left);
            rdm2_target->element(i+norb_left, j+norb_left, p, q) = value;
            rdm2_target->element(p, q, i+norb_left, j+norb_left) = value;
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_220_part1


// orbital i,j on site, p,q on left block
// \Gamma_{iqpj} = (-1) \sum \left< A^{c}_{lr} | a^{\dagger}_{i \rho} a_{j \sigma} | A^{c'}_{l'r} \right> 
//                 \left< l | a^{\dagger}_{p \sigma} a_{q \rho} | l' \right>
void ASD_DMRG::compute_rdm2_220_part2(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  const int nstate = dvec.size();
  bool TRANS_LEFT = false;
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on left}, leftbra_deltaAlpha, leftbra_deltaBeta, swap}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, 1,-1, false}, // transpose left transition tensor
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, 1,-1, true}   //  transpose left transition tensor
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      if (get<2>(gammalist_tuple)) TRANS_LEFT = true;
      const bool swap_idx = get<4>(gammalist_tuple);
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_left); // matrix to store RDM, use ax_plus_y...
      // marco loop over right blocks since we have a delta_{r,r'}
      for (auto& rblock : right_block->blocks()) {
        BlockKey rightkey = rblock.key(); 
        const int rightnstates = rblock.nstates;
        // micro loop over one specific right block
        for (int ir = 0; ir != rightnstates; ++ir) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between left block states with the same right block state
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& lblock : left_block->blocks()) {
            BlockKey leftkey = lblock.key();
            const int leftnstates = lblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int ileft = 0; ileft != leftnstates; ++ileft)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + ir*leftnstates + ileft)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density matrices
          for (auto& lbinfo : left_block->blocks()) {
            BlockKey ket_leftkey = lbinfo.key();
            BlockKey ket_combinedkey(ket_leftkey.nelea + rightkey.nelea, ket_leftkey.neleb + rightkey.neleb);
            const int ket_nstates = lbinfo.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple), ket_leftkey.neleb + get<3>(gammalist_tuple));
            if (!(left_block->contains(bra_leftkey))) continue;
            BlockKey bra_combinedkey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = left_block->blockinfo(bra_leftkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)))) continue;

            // site transition tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              btas::CRange<3> site_range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
              site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(transition_mat->storage()));
            }
 
            // left transition tensor
            shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
            {
              shared_ptr<const btas::Tensor3<double>> coupling_data = 
                  left_block->coupling(get<1>(gammalist_tuple)).at(make_pair((TRANS_LEFT ? ket_leftkey : bra_leftkey),
                                                                             (TRANS_LEFT ? bra_leftkey : ket_leftkey))).data;
              btas::CRange<3> left_range(bra_nstates, ket_nstates, lrint(pow(norb_left, get<1>(gammalist_tuple).size())));
              left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, move(coupling_data->storage()));
              if (TRANS_LEFT) {
                unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
                for (int i = 0; i != left_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*left_coupling_tensor)(0,0,i), bra_nstates*ket_nstates, buf.get());
                  blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*left_coupling_tensor)(0,0,i));
                }
              }
            }
            
            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_tensor->extent(2) * left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,0,2), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(-1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
      
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int j = 0; j != norb_site; ++j) {
        for (int i = 0; i != norb_site; ++i) {
          for (int q = 0; q != norb_left; ++q) {
            for (int p = 0; p != norb_left; ++p) {
              double value;
              if (swap_idx) value = *rdm_mat->element_ptr(j+i*norb_site, p+q*norb_left);
              else if (TRANS_LEFT) value = *rdm_mat->element_ptr(i+j*norb_site, q+p*norb_left);
              else value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_left);
              rdm2_target->element(i+norb_left, q, p, j+norb_left) += value;
              rdm2_target->element(q, i+norb_left, j+norb_left, p) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_220_part2


// orbital i,j on site, p,q on left block
// \Gamma_{iqjp} = \sum \left< A^{c}_{lr} | a^{\dagger}_{i \rho} a^{\dagger}_{j \sigma} | A^{c}_{l'r} \right> \left< l | a_{p \sigma} a_{q \rho} | l' \right>
void ASD_DMRG::compute_rdm2_220_part3(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on left}, leftbra_deltaAlpha, leftbra_deltaBeta, swap}
    {{GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},-2, 0, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  0,-2, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, true}
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(site);
    const int norb_left = left_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      const bool swap_idx = get<4>(gammalist_tuple);
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_left); // matrix to store RDM, use ax_plus_y...
      // marco loop over right blocks since we have a delta_{r,r'}
      for (auto& rblock : right_block->blocks()) {
        BlockKey rightkey = rblock.key(); 
        const int rightnstates = rblock.nstates;
        // micro loop over one specific right block
        for (int ir = 0; ir != rightnstates; ++ir) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices between left block states with the same right block state
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& lblock : left_block->blocks()) {
            BlockKey leftkey = lblock.key();
            const int leftnstates = lblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int ileft = 0; ileft != leftnstates; ++ileft)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + ir*leftnstates + ileft)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over left blocks again to obtain all transition density matrices
          for (auto& lbinfo : left_block->blocks()) {
            BlockKey ket_leftkey = lbinfo.key();
            BlockKey ket_combinedkey(ket_leftkey.nelea + rightkey.nelea, ket_leftkey.neleb + rightkey.neleb);
            const int ket_nstates = lbinfo.nstates;
            BlockKey bra_leftkey(ket_leftkey.nelea + get<2>(gammalist_tuple), ket_leftkey.neleb + get<3>(gammalist_tuple));
            if (!(left_block->contains(bra_leftkey))) continue;
            BlockKey bra_combinedkey(bra_leftkey.nelea + rightkey.nelea, bra_leftkey.neleb + rightkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = left_block->blockinfo(bra_leftkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple)))) continue;

            // site transition tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
              btas::CRange<3> site_range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(transition_mat->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != site_transition_tensor->extent(2); ++i) {
                copy_n(&(*site_transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*site_transition_tensor)(0,0,i));
              }
            }
            
            // left block transition tensor
            shared_ptr<const btas::Tensor3<double>> left_coupling_tensor = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;

            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_tensor->extent(2) * left_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,0,2), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
 
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int q = 0; q != norb_left; ++q) {
        for (int p = 0; p != norb_left; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            for (int j = 0; j != norb_site; ++j) {
              if (swap_idx) swap(i,j);
              if (swap_idx) swap(p,q);
              const double value = *rdm_mat->element_ptr(j+i*norb_site, p+q*norb_left);
              if (swap_idx) swap(i,j);
              if (swap_idx) swap(p,q);
              rdm2_target->element(i+norb_left, q, j+norb_left, p) += value;
              rdm2_target->element(q, i+norb_left, p, j+norb_left) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_220_part3


void ASD_DMRG::compute_rdm2_022(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_022" << endl;
  
  // a^{\dagger}_{p \rho} a_{q \rho} on right
  compute_rdm2_022_part1(dvec);
  
  // a^{\dagger}_{p \rho} a_{q \sigma} on right
  compute_rdm2_022_part2(dvec);

  // a^{\dagger}_{p \rho} a^{\dagger}_{q \sigma} on right
  compute_rdm2_022_part3(dvec);
}


// orbital i,j on site, p,q on right
// \Gamma_{ijpq} = \sum \left< A^{c}_{lr} | a^{\dagger}_{p \rho} a_{q \rho} | A^{c'}_{lr'}
//                 \right> \left< r | a^{\dagger}_{r \sigma} a_{s \sigma} | r' \right>
void ASD_DMRG::compute_rdm2_022_part1(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on right}, rightbra_deltaAlpha, rightbra_deltaBeta}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0},
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0}
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_left = left_block->norb();
    const int norb_rightoffset = norb_left + norb_site;
    const int norb_right = right_block->norb();
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // marco loop over left blocks since we have a delta_{l,l'}
      for (auto& lblock : left_block->blocks()) {
        BlockKey leftkey = lblock.key(); 
        const int leftnstates = lblock.nstates;
        // micro loop over one specific right block
        for (int il = 0; il != leftnstates; ++il) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store site transition density matrices
          // find all possible left blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& rblock : right_block->blocks()) {
            BlockKey rightkey = rblock.key();
            const int rightnstates = rblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int iright = 0; iright != rightnstates; ++iright)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + il + iright*leftnstates)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over right blocks again to obtain all transition density matrices
          for (auto& rbinfo : right_block->blocks()) {
            BlockKey ket_rightkey = rbinfo.key();
            BlockKey ket_combinedkey(ket_rightkey.nelea + leftkey.nelea, ket_rightkey.neleb + leftkey.neleb);
            BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple), ket_rightkey.neleb + get<3>(gammalist_tuple));
            if (!(right_block->contains(bra_rightkey))) continue;
            BlockKey bra_combinedkey(bra_rightkey.nelea + leftkey.nelea, bra_rightkey.neleb + leftkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)))) continue;

            // site transition tensor in matrix form
            shared_ptr<const Matrix> site_transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
            
            // right block transition tensor
            shared_ptr<const btas::Tensor3<double>> right_coupling_tensor = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_rightkey, ket_rightkey)).data;
            
            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_mat->extent(1) * right_coupling_tensor->extent(2));
            contract(1.0, *site_transition_mat, {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }
    
    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_right; ++q) {
      for (int p = 0; p != norb_right; ++p) {
        for (int j = 0; j != norb_site; ++j) {
          for (int i = 0; i != norb_site; ++i) {
            const double value = *rdm_mat->element_ptr(i + j*norb_site, p + q*norb_right);
            rdm2_target->element(i+norb_right, j+norb_right, p+norb_rightoffset, q+norb_rightoffset) = value;
            rdm2_target->element(p+norb_rightoffset, q+norb_rightoffset, i+norb_right, j+norb_right) = value;
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_022_part1


// orbital i,j on site, p,q on right
// \Gamma_{iqpj} = (-1) \sum \left< A^{c}_{lr} | a^{\dagger}_{i \rho} a_{j \sigma} | A^{c'}_{lr'} \right> 
//                 \left< r | a^{\dagger}_{\sigma} a_{q \rho} | r' \right>
void ASD_DMRG::compute_rdm2_022_part2(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  bool TRANS_RIGHT = false;
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on left}, leftbra_deltaAlpha, leftbra_deltaBeta, swap}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, 1,-1, false}, // transpose left transition tensor
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, 1,-1, true}   //  transpose left transition tensor and swap indices
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_right = right_block->norb();
    const int norb_rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      if (get<2>(gammalist_tuple)) TRANS_RIGHT = true;
      const bool swap_idx = get<4>(gammalist_tuple);
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...
      // marco loop over left blocks since we have a delta_{l,l'}
      for (auto& lblock : left_block->blocks()) {
        BlockKey leftkey = lblock.key(); 
        const int leftnstates = lblock.nstates;
        // micro loop over one specific left block
        for (int il = 0; il != leftnstates; ++il) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices
          // find all possible right blocks that can be coupled with the right state and compute GammaForestASD
          for (auto& rblock : right_block->blocks()) {
            BlockKey rightkey = rblock.key();
            const int rightnstates = rblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int iright = 0; iright != rightnstates; ++iright)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + il + iright*leftnstates)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over right blocks again to obtain all transition density matrices
          for (auto& rbinfo : right_block->blocks()) {
            BlockKey ket_rightkey = rbinfo.key();
            BlockKey ket_combinedkey(ket_rightkey.nelea + leftkey.nelea, ket_rightkey.neleb + leftkey.neleb);
            const int ket_nstates = rbinfo.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple), ket_rightkey.neleb + get<3>(gammalist_tuple));
            if (!(right_block->contains(bra_rightkey))) continue;
            BlockKey bra_combinedkey(bra_rightkey.nelea + leftkey.nelea, bra_rightkey.neleb + leftkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = right_block->blockinfo(bra_rightkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple)))) continue;

            // site transition tensor in matrix form
            auto site_transition_mat = forest.template get<0>(bra_rastag, ket_rastag, get<0>(gammalist_tuple));
 
            // right transition tensor
            shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
            {
              shared_ptr<const btas::Tensor3<double>> coupling_data =
                  right_block->coupling(get<1>(gammalist_tuple)).at(make_pair((TRANS_RIGHT ? ket_rightkey : bra_rightkey),
                                                                              (TRANS_RIGHT ? bra_rightkey : ket_rightkey))).data;
              btas::CRange<3> right_range(bra_nstates, ket_nstates, lrint(pow(norb_right, get<1>(gammalist_tuple).size())));
              right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, move(coupling_data->storage()));
              if (TRANS_RIGHT) {
                unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
                for (int i = 0; i != right_coupling_tensor->extent(2); ++i) {
                  copy_n(&(*right_coupling_tensor)(0,0,i), bra_nstates*ket_nstates, buf.get());
                  blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*right_coupling_tensor)(0,0,i));
                }
              }
            }
            
            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_mat->extent(1) * right_coupling_tensor->extent(2));
            contract(1.0, *site_transition_mat, {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(-1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
      
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int j = 0; j != norb_site; ++j) {
        for (int i = 0; i != norb_site; ++i) {
          for (int q = 0; q != norb_right; ++q) {
            for (int p = 0; p != norb_right; ++p) {
              double value;
              if (swap_idx) value = *rdm_mat->element_ptr(j+i*norb_site, p+q*norb_right);
              else if (TRANS_RIGHT) value = *rdm_mat->element_ptr(i+j*norb_site, q+p*norb_right);
              else value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_right);
              rdm2_target->element(i+norb_left, q+norb_rightoffset, p+norb_rightoffset, j+norb_left) += value;
              rdm2_target->element(q+norb_rightoffset, i+norb_left, j+norb_left, p+norb_rightoffset) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_220_part2


// orbital i,j on site, p,q on right
// \Gamma_{iqjp} = \left< A^{c}_{lr} | a^{\dagger}_{i \rho} a^{\dagger}_{j \sigma} | A^{c'}_{lr'} \right> 
//                 \left< r | a_{p \sigma} a_{q \rho} | r' \right>
void ASD_DMRG::compute_rdm2_022_part3(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = { 
    // ({operators on site}, {operators on left}, leftbra_deltaAlpha, leftbra_deltaBeta, swap}
    {{GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},-2, 0, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  0,-2, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, false},
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, true}
  };

  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_site = multisite_->active_sizes().at(nsites_-2);
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_rightoffset = norb_left + norb_site;
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      const bool swap_idx = get<4>(gammalist_tuple);
      auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...
      // marco loop over left blocks since we have a delta_{l,l'}
      for (auto& lblock : left_block->blocks()) {
        BlockKey leftkey = lblock.key(); 
        const int leftnstates = lblock.nstates;
        // micro loop over one specific left block
        for (int il = 0; il != leftnstates; ++il) {
          map<BlockKey, shared_ptr<const RASDvec>> states; // store transition density matrices
          // find all possible right blocks that can be coupled with the left state and compute GammaForestASD
          for (auto& rblock : right_block->blocks()) {
            BlockKey rightkey = rblock.key();
            const int rightnstates = rblock.nstates;
            BlockKey combinedkey(rightkey.nelea+leftkey.nelea, rightkey.neleb+leftkey.neleb);
            if (!(prod_civec->contains_block(combinedkey))) continue;
            auto bpair = doubleblock->blockpairs(combinedkey);
            auto iter = find_if(bpair.begin(), bpair.end(), [&lblock, &rblock] (const DMRG::BlockPair& bp)
              { return make_pair(lblock, rblock) == make_pair(bp.left, bp.right); });
            assert(iter != bpair.end());
            const int offset = iter->offset;
            // transform blockkey into raskey
            const int ras_nelea = tot_nelea - combinedkey.nelea;
            const int ras_neleb = tot_neleb - combinedkey.neleb;
            BlockKey ras_key(ras_nelea, ras_neleb);
            vector<shared_ptr<RASCivec>> tmpvec;
            for (int iright = 0; iright != rightnstates; ++iright)
              tmpvec.push_back(make_shared<RASCivec>(prod_civec->sector(combinedkey)->civec(offset + il + iright*leftnstates)));
            states[ras_key] = make_shared<const RASDvec>(tmpvec);
          }
          GammaForestASD<RASDvec> forest(states);
          forest.compute();

          // loop over right blocks again to obtain all transition density matrices
          for (auto& rbinfo : right_block->blocks()) {
            BlockKey ket_rightkey = rbinfo.key();
            BlockKey ket_combinedkey(ket_rightkey.nelea + leftkey.nelea, ket_rightkey.neleb + leftkey.neleb);
            const int ket_nstates = rbinfo.nstates;
            BlockKey bra_rightkey(ket_rightkey.nelea + get<2>(gammalist_tuple), ket_rightkey.neleb + get<3>(gammalist_tuple));
            if (!(right_block->contains(bra_rightkey))) continue;
            BlockKey bra_combinedkey(bra_rightkey.nelea + leftkey.nelea, bra_rightkey.neleb + leftkey.neleb);
            if (!(prod_civec->contains_block(bra_combinedkey))) continue;
            const int bra_nstates = right_block->blockinfo(bra_rightkey).nstates;
            // transform blockkey into raskey
            const int ket_ras_nelea = tot_nelea - ket_combinedkey.nelea;
            const int ket_ras_neleb = tot_neleb - ket_combinedkey.neleb;
            BlockKey ket_raskey(ket_ras_nelea, ket_ras_neleb);
            const int bra_ras_nelea = tot_nelea - bra_combinedkey.nelea;
            const int bra_ras_neleb = tot_neleb - bra_combinedkey.neleb;
            BlockKey bra_raskey(bra_ras_nelea, bra_ras_neleb);

            const size_t ket_rastag = forest.block_tag(ket_raskey);
            const size_t bra_rastag = forest.block_tag(bra_raskey);
            if (!(forest.template exist<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple)))) continue;

            // site transition tensor
            shared_ptr<btas::Tensor3<double>> site_transition_tensor;
            {
              shared_ptr<const Matrix> transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
              btas::CRange<3> site_range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
              site_transition_tensor = make_shared<btas::Tensor3<double>>(site_range, move(transition_mat->storage()));
              unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
              for (int i = 0; i != site_transition_tensor->extent(2); ++i) {
                copy_n(&(*site_transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
                blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*site_transition_tensor)(0,0,i));
              }
            }
            
            // right block transition tensor
            shared_ptr<const btas::Tensor3<double>> right_coupling_tensor = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_rightkey, ket_rightkey)).data;

            auto target = rdm_mat->clone();
            assert (target->size() == site_transition_tensor->extent(2) * right_coupling_tensor->extent(2));
            contract(1.0, group(*site_transition_tensor,0,2), {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *target, {0,1});
            blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
 
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int q = 0; q != norb_right; ++q) {
        for (int p = 0; p != norb_right; ++p) {
          for (int i = 0; i != norb_site; ++i) {
            for (int j = 0; j != norb_site; ++j) {
              if (swap_idx) { swap(i,j); swap(p,q); }
              const double value = *rdm_mat->element_ptr(j+i*norb_site, p+q*norb_right);
              if (swap_idx) { swap(i,j); swap(p,q); }
              rdm2_target->element(i+norb_left, q+norb_rightoffset, j+norb_left, p+norb_rightoffset) += value;
              rdm2_target->element(q+norb_rightoffset, i+norb_left, p+norb_rightoffset, j+norb_left) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_022_part3


void ASD_DMRG::compute_rdm2_202(vector<shared_ptr<ProductRASCivec>> dvec) {
  cout << "  * compute_rdm2_202" << endl;
  
  // a^{\dagger}_{p \rho} a_{q \rho} on right
  compute_rdm2_202_part1(dvec);
  
  // a^{\dagger}_{p \rho} a_{q \sigma} on right
  compute_rdm2_202_part2(dvec);
 
  // a^{\dagger}_{p \rho} a^{\dagger}_{q \sigma} on right
  compute_rdm2_202_part3(dvec);
}


// orbital i,j on left, p,q on right
// \Gamma_{ijpq} = \sum A^{c}_{lr} A^{c}_{l'r'} \left< l | a^{\dagger}_{i \rho} a_{j \rho} | l' \right> 
//                 \left< r | a^{\dagger}_{p \sigma} a_{q \sigma} | r' \right>
void ASD_DMRG::compute_rdm2_202_part1(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>>> gammalist_tuple_list = { 
    // ({operators on left}, {operators on right}}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}},
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}}
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_rightoffset = norb_left + multisite_->active_sizes().at(nsites_-2);
    auto rdm_mat = make_shared<Matrix>(norb_left*norb_left, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          const int pairoffset = bpair.offset;
          BlockInfo leftinfo = bpair.left;
          const int left_nstates = leftinfo.nstates;
          // right bra-ket pair
          BlockInfo rightinfo = bpair.right;
          const int right_nstates = rightinfo.nstates;
          
          shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(leftinfo, leftinfo)).data;
          
          shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(rightinfo, rightinfo)).data;

          // contract coefficient tensor
          auto contract_mat = make_shared<Matrix>(left_nstates*left_nstates, right_nstates*right_nstates);
          auto ciptr = prod_civec->sector(seckey);
          for (int irket = 0; irket != right_nstates; ++irket)
            for (int irbra = 0; irbra != right_nstates; ++irbra)
              for (int ilket = 0; ilket != left_nstates; ++ilket)
                for (int ilbra = 0; ilbra != left_nstates; ++ilbra)
                  *contract_mat->element_ptr(ilbra+ilket*left_nstates, irbra+irket*right_nstates) =
                     ciptr->civec(pairoffset + ilbra + irbra*left_nstates).dot_product(ciptr->civec(pairoffset + ilket + irket*left_nstates));

          auto tmp_mat = make_shared<Matrix>(right_nstates*right_nstates, left_coupling->extent(2));
          contract(1.0, *contract_mat, {2,0}, group(*left_coupling,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

          auto target = rdm_mat->clone();
          contract(1.0, *tmp_mat, {2,0}, group(*right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
          blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
        }
      }
    }

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_right; ++q) {
      for (int p = 0; p != norb_right; ++p) {
        for (int j = 0; j != norb_left; ++j) {
          for (int i = 0; i != norb_left; ++i) {
            const double value = *rdm_mat->element_ptr(i+j*norb_left, p+q*norb_right);
            rdm2_target->element(i, j, p+norb_rightoffset, q+norb_rightoffset) = value;
            rdm2_target->element(p+norb_rightoffset, q+norb_rightoffset, i, j) = value;
          }
        }
      }
    }
  
  } // end of looping over nstates
} // end of compute_202_part1


// orbital i,j on left, p,q on right
// \Gamma_{iqpj} = (-1) \sum A^{c}_{lr} A^{c}_{l'r'} \left< l | a^{\dagger}_{i \rho} a_{j \sigma} | l' \right> 
//                 \left< r | a^{\dagger}_{p \sigma} a_{q \rho} | r' \right>
void ASD_DMRG::compute_rdm2_202_part2(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = {
    // ({operators on left}, {operators on right}, rightbra_deltaAlpha, rightbra_deltaBeta, swap}
    {{GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, 0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  0, 0, false},
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},-1, 1, false}, // transpose left
    {{GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},-1, 1, true}   // transpose left and swap index
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_rightoffset = norb_left + multisite_->active_sizes().at(nsites_-2);

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      auto rdm_mat = make_shared<Matrix>(norb_left*norb_left, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...
      const bool TRANS_LEFT = get<3>(gammalist_tuple);
      const bool swap_idx = get<4>(gammalist_tuple);
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = bpair.offset;
          // left bra-ket pair
          BlockInfo ket_leftinfo = bpair.left;
          const int lket_nstates = ket_leftinfo.nstates;
          BlockKey bra_leftkey(ket_leftinfo.nelea - get<2>(gammalist_tuple), ket_leftinfo.neleb - get<3>(gammalist_tuple));
          if (!left_block->contains(bra_leftkey)) continue;
          const int lbra_nstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockInfo ket_rightinfo = bpair.right;
          const int rket_nstates = ket_rightinfo.nstates;
          BlockKey bra_rightkey(ket_rightinfo.nelea + get<2>(gammalist_tuple), ket_rightinfo.neleb + get<3>(gammalist_tuple));
          if (!right_block->contains(bra_rightkey)) continue;
          const int rbra_nstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          shared_ptr<btas::Tensor3<double>> left_transition_tensor;
          {
            shared_ptr<const btas::Tensor3<double>> left_coupling 
              = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair((TRANS_LEFT ? ket_leftinfo : bra_leftkey),
                                                                           (TRANS_LEFT ? bra_leftkey : ket_leftinfo))).data;
            btas::CRange<3> left_range(lbra_nstates, lket_nstates, lrint(pow(norb_left, get<0>(gammalist_tuple).size())));
            left_transition_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
            if (TRANS_LEFT) {
              unique_ptr<double[]> buf(new double[lbra_nstates*lket_nstates]);
              for (int i = 0; i != left_transition_tensor->extent(2); ++i) {
                copy_n(&(*left_transition_tensor)(0,0,i), lbra_nstates*lket_nstates, buf.get());
                blas::transpose(buf.get(), lket_nstates, lbra_nstates, &(*left_transition_tensor)(0,0,i));
              }
            }
          }

          shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_rightkey, ket_rightinfo)).data;

          // contract coefficient tensor
          auto contract_mat = make_shared<Matrix>(lbra_nstates*lket_nstates, rbra_nstates*rket_nstates);
          auto ciptr = prod_civec->sector(seckey);
          for (int irket = 0; irket != rket_nstates; ++irket)
            for (int irbra = 0; irbra != rbra_nstates; ++irbra)
              for (int ilket = 0; ilket != lket_nstates; ++ilket)
                for (int ilbra = 0; ilbra != lbra_nstates; ++ilbra)
                  *contract_mat->element_ptr(ilbra+ilket*lbra_nstates, irbra+irket*rbra_nstates) =
                     ciptr->civec(brapairoffset + ilbra + irbra*lbra_nstates).dot_product(ciptr->civec(ketpairoffset + ilket + irket*lket_nstates));

          auto tmp_mat = make_shared<Matrix>(rbra_nstates*rket_nstates, left_transition_tensor->extent(2));
          contract(1.0, *contract_mat, {2,0}, group(*left_transition_tensor,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

          auto target = rdm_mat->clone();
          contract(1.0, *tmp_mat, {2,0}, group(*right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
          blas::ax_plus_y_n(-1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
        }
      }
      
      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int q = 0; q != norb_right; ++q) {
        for (int p = 0; p != norb_right; ++p) {
          for (int j = 0; j != norb_left; ++j) {
            for (int i = 0; i != norb_left; ++i) {
              double value;
              if (swap_idx) value = *rdm_mat->element_ptr(i+j*norb_left, q+p*norb_right);
              else if (TRANS_LEFT) value = *rdm_mat->element_ptr(j+i*norb_left, p+q*norb_right);
              else value = *rdm_mat->element_ptr(i+j*norb_left, p+q*norb_right);
  
              rdm2_target->element(i, q+norb_rightoffset, p+norb_rightoffset, j) += value;
              rdm2_target->element(p+norb_rightoffset, j, i, q+norb_rightoffset) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_202_part2


// orbital i,j on left p,q on right
// \Gamma_{iqjp} = \sum A^{c}_{lr} A^{c}_{l'r'} \left< l | a^{\dagger}_{i \rho} a^{\dagger}_{j \sigma} | l' \right> 
//                 \left< r | a_{p \sigma} a_{q \rho} | r' \right>
void ASD_DMRG::compute_rdm2_202_part3(vector<shared_ptr<ProductRASCivec>> dvec) {
  const int nstate = dvec.size();
  list<tuple<list<GammaSQ>, list<GammaSQ>, int, int, bool>> gammalist_tuple_list = { 
    // ({operators on left}, {operators on right}, rightbra_deltaAlpha, rightbra_deltaBeta, swap}
    {{GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},-2, 0, false}, // transpose left
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  0,-2, false}, // transpose left
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, false}, // transpose left
    {{GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},-1,-1, true}  // transpose left and swap index
  };
  
  for (int istate = 0; istate != nstate; ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_rightoffset = norb_left + multisite_->active_sizes().at(nsites_-2/*site*/);

    for (auto& gammalist_tuple : gammalist_tuple_list) {
      const bool swap_idx = get<4>(gammalist_tuple);
      auto rdm_mat = make_shared<Matrix>(norb_left*norb_left, norb_right*norb_right); // matrix to store RDM, use ax_plus_y...
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          const int ketpairoffset = bpair.offset;
          // left bra-ket pair
          BlockInfo ket_leftinfo = bpair.left;
          const int lket_nstates = ket_leftinfo.nstates;
          BlockKey bra_leftkey(ket_leftinfo.nelea - get<2>(gammalist_tuple), ket_leftinfo.neleb - get<3>(gammalist_tuple));
          if (!left_block->contains(bra_leftkey)) continue;
          const int lbra_nstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockInfo ket_rightinfo = bpair.right;
          const int rket_nstates = ket_rightinfo.nstates;
          BlockKey bra_rightkey(ket_rightinfo.nelea + get<2>(gammalist_tuple), ket_rightinfo.neleb + get<3>(gammalist_tuple));
          if (!right_block->contains(bra_rightkey)) continue;
          const int rbra_nstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp)
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;
          
          shared_ptr<btas::Tensor3<double>> left_coupling;
          {
            shared_ptr<const btas::Tensor3<double>> lcoupling 
              = left_block->coupling(get<0>(gammalist_tuple)).at(make_pair(ket_leftinfo, bra_leftkey)).data;
            btas::CRange<3> range(lbra_nstates, lket_nstates, lrint(pow(norb_left, get<0>(gammalist_tuple).size())));
            left_coupling = make_shared<btas::Tensor3<double>>(range, move(lcoupling->storage()));
            unique_ptr<double[]> buf(new double[lbra_nstates*lket_nstates]);
            for (int i = 0; i != left_coupling->extent(2); ++i) {
              copy_n(&(*left_coupling)(0,0,i), lbra_nstates*lket_nstates, buf.get());
              blas::transpose(buf.get(), lket_nstates, lbra_nstates, &(*left_coupling)(0,0,i));
            }
          }
          
          shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_rightkey, ket_rightinfo)).data;
          
          // contract coefficient tensor
          auto contract_mat = make_shared<Matrix>(lbra_nstates*lket_nstates, rbra_nstates*rket_nstates);
          auto ciptr = prod_civec->sector(seckey);
          for (int irket = 0; irket != rket_nstates; ++irket)
            for (int irbra = 0; irbra != rbra_nstates; ++irbra)
              for (int ilket = 0; ilket != lket_nstates; ++ilket)
                for (int ilbra = 0; ilbra != lbra_nstates; ++ilbra)
                  *contract_mat->element_ptr(ilbra+ilket*lbra_nstates, irbra+irket*rbra_nstates) =
                     ciptr->civec(brapairoffset + ilbra + irbra*lbra_nstates).dot_product(ciptr->civec(ketpairoffset + ilket + irket*lket_nstates));

          auto tmp_mat = make_shared<Matrix>(rbra_nstates*rket_nstates, left_coupling->extent(2));
          contract(1.0, *contract_mat, {2,0}, group(*left_coupling,0,2), {2,1}, 0.0, *tmp_mat, {0,1});

          auto target = rdm_mat->clone();
          contract(1.0, *tmp_mat, {2,0}, group(*right_coupling,0,2), {2,1}, 0.0, *target, {0,1});
          blas::ax_plus_y_n(1.0/*sign*/, target->data(), target->size(), rdm_mat->data());
        }
      }

      // copy data into rdm2_
      auto rdm2_target = rdm2_->at(istate);
      for (int q = 0; q != norb_right; ++q) {
        for (int p = 0; p != norb_right; ++p) {
          for (int i = 0; i != norb_left; ++i) {
            for (int j = 0; j != norb_left; ++j) {
              if (swap_idx) { swap(i,j); swap(p,q); }
              const double value = *rdm_mat->element_ptr(j+i*norb_left, p+q*norb_right);
              if (swap_idx) { swap(i,j); swap(p,q); }
              rdm2_target->element(i, q+norb_rightoffset, j, p+norb_rightoffset) += value;
              rdm2_target->element(q+norb_rightoffset, i, p+norb_rightoffset, j) += value;
            }
          }
        }
      }
    }
  } // end of looping over nstates
} // end of compute_202_part3


void ASD_DMRG::compute_rdm2_121(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  cout << "compute_rdm2_121" << endl;

  compute_rdm2_121_part1(dvec, site);
  
}


// orbital i,j on site, p on left, q on right
void ASD_DMRG::compute_rdm2_121_part1(vector<shared_ptr<ProductRASCivec>> dvec, const int site) {
  const list<tuple<list<GammaSQ>, list<GammaSQ>, list<GammaSQ>, pair<int, int>, pair<int, int>>> gammalist_tuple_list = {
    // {ops on site}, {ops on left}, {ops on right}, {left nele change}, {right nele change}
    { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {1,0}, {-1,0} },
    { {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {0,1}, {0,-1} },
    { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateAlpha}, {GammaSQ::CreateAlpha}, {1,0}, {-1,0} },
    { {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  {GammaSQ::CreateBeta},  {GammaSQ::CreateBeta},  {0,1}, {0,-1} }
  };
  
  for (int istate = 0; istate != dvec.size(); ++istate) {
    auto prod_civec = dvec.at(istate);
    shared_ptr<const DMRG_Block2> doubleblock = dynamic_pointer_cast<const DMRG_Block2>(prod_civec->left());
    auto left_block = doubleblock->left_block();
    auto right_block = doubleblock->right_block();
    const int norb_left = left_block->norb();
    const int norb_right = right_block->norb();
    const int norb_site = multisite_->active_sizes().at(site);
    const int tot_nelea = prod_civec->nelea();
    const int tot_neleb = prod_civec->neleb();
    auto rdm_mat = make_shared<Matrix>(norb_site*norb_site, norb_left*norb_right);

    for (auto& gammalist_tuple : gammalist_tuple_list) {   
      // loop over product rasci sectors
      for (auto& isec : prod_civec->sectors()) {
        BlockKey seckey = isec.first;
        for (auto& bpair : doubleblock->blockpairs(seckey)) {
          cout << "start of a bpair " << endl;
          const int ketpairoffset = bpair.offset;
          // left bra-ket pair
          BlockInfo ket_leftinfo = bpair.left;
          const int ket_leftnstates = ket_leftinfo.nstates;
          BlockKey bra_leftkey(ket_leftinfo.nelea + get<3>(gammalist_tuple).first, ket_leftinfo.neleb + get<3>(gammalist_tuple).second);
          if (!left_block->contains(bra_leftkey)) continue;
          const int bra_leftnstates = left_block->blockinfo(bra_leftkey).nstates;
          // right bra-ket pair
          BlockInfo ket_rightinfo = bpair.right;
          const int ket_rightnstates = ket_rightinfo.nstates;
          BlockKey bra_rightkey(ket_rightinfo.nelea + get<4>(gammalist_tuple).first, ket_rightinfo.neleb + get<4>(gammalist_tuple).second);
          if (!right_block->contains(bra_rightkey)) continue;
          const int bra_rightnstates = right_block->blockinfo(bra_rightkey).nstates;
          auto brapair = doubleblock->blockpairs(seckey);
          auto braiter = find_if(brapair.begin(), brapair.end(), [&left_block, &right_block, &bra_leftkey, &bra_rightkey] (const DMRG::BlockPair& bp) 
            { return make_pair(left_block->blockinfo(bra_leftkey), right_block->blockinfo(bra_rightkey)) == make_pair(bp.left, bp.right); });
          assert(braiter != brapair.end());
          const int brapairoffset = braiter->offset;

          // site transition density tensor
          shared_ptr<btas::Tensor3<double>> site_transition_tensor;
          {
            const int ras_nelea = tot_nelea - seckey.nelea;
            const int ras_neleb = tot_neleb - seckey.neleb;
            BlockKey raskey(ras_nelea, ras_neleb);
            map<BlockKey, shared_ptr<const RASDvec>> ket_states;
            vector<shared_ptr<RASCivec>> ket_vecs;
            for (int ket_ir = 0; ket_ir != ket_rightnstates; ++ket_ir)
              for (int ket_il = 0; ket_il != ket_leftnstates; ++ket_il)
                ket_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(seckey)->civec(ketpairoffset + ket_il + ket_ir*ket_leftnstates)));
            ket_states[raskey] = make_shared<const RASDvec>(ket_vecs);
            map<BlockKey, shared_ptr<const RASDvec>> bra_states;
            vector<shared_ptr<RASCivec>> bra_vecs;
            for (int bra_ir = 0; bra_ir != bra_rightnstates; ++bra_ir)
              for (int bra_il = 0; bra_il != bra_leftnstates; ++bra_il)
                bra_vecs.push_back(make_shared<RASCivec>(prod_civec->sector(seckey)->civec(brapairoffset + bra_il + bra_ir*bra_leftnstates)));
            bra_states[raskey] = make_shared<const RASDvec>(bra_vecs);
  
            // construct and compute GammaForest for site
            GammaForestASD2<RASDvec> forest(bra_states, ket_states);
            forest.compute();
  
            const size_t rastag = forest.block_tag(raskey);
            assert(forest.template exist<0>(rastag, rastag, get<0>(gammalist_tuple)));
            shared_ptr<const Matrix> transition_mat = forest.template get<0>(rastag, rastag, get<0>(gammalist_tuple));
            //btas::CRange<3> siterange(bra_leftnstates*ket_leftnstates, bra_rightnstates*ket_rightnstates, norb_site*norb_site);
            btas::CRange<3> tmprange(ket_leftnstates*bra_leftnstates, bra_rightnstates, ket_rightnstates*norb_site*norb_site);
            // index : {bra_left, bra_right}, {ket_left}, {ket_right, i, j     }
            auto tmp_transition_tensor = make_shared<btas::Tensor3<double>>(tmprange, transition_mat->storage()); 
 
            vector<double> buf1(ket_leftnstates*bra_leftnstates*bra_rightnstates);
            for (int i = 0; i != tmp_transition_tensor->extent(2); ++i) {
              copy_n(&(*tmp_transition_tensor)(0,0,i),ket_leftnstates*bra_leftnstates*bra_rightnstates, buf1.data());
              blas::transpose(buf1.data(), bra_leftnstates*bra_rightnstates, ket_leftnstates, &(*tmp_transition_tensor)(0,0,i));
            } // index : ket_left, bra_left, bra_right, ket_right

            btas::CRange<3> siterange(ket_leftnstates*bra_leftnstates, bra_rightnstates*ket_rightnstates, norb_site*norb_site);
            site_transition_tensor = make_shared<btas::Tensor3<double>>(siterange, move(tmp_transition_tensor->storage()));
          }

          // transposed left coupling tensor
          shared_ptr<btas::Tensor3<double>> left_coupling_tensor;
          {
            btas::CRange<3> left_range(ket_leftnstates, bra_leftnstates, norb_left);
            shared_ptr<const btas::Tensor3<double>> left_coupling = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftinfo.key())).data;
            left_coupling_tensor = make_shared<btas::Tensor3<double>>(left_range, left_coupling->storage());
            vector<double> buf2(ket_leftnstates*bra_leftnstates);
            for (int i = 0; i != norb_left; ++i) {
              copy_n(&(*left_coupling_tensor)(0,0,i), ket_leftnstates*bra_leftnstates, buf2.data());
              blas::transpose(buf2.data(), bra_leftnstates, ket_leftnstates, &(*left_coupling_tensor)(0,0,i));
            }
          }
          
          // right coupling tensor
          shared_ptr<btas::Tensor3<double>> right_coupling_tensor;
          {
            btas::CRange<3> right_range(bra_rightnstates, ket_rightnstates, norb_right);
            shared_ptr<const btas::Tensor3<double>> right_coupling = right_block->coupling(get<2>(gammalist_tuple)).at(make_pair(ket_rightinfo, bra_rightkey)).data;
            right_coupling_tensor = make_shared<btas::Tensor3<double>>(right_range, right_coupling->storage());
            vector<double> buf3(bra_rightnstates*ket_rightnstates);
            for (int j = 0; j != norb_right; ++j) {
              copy_n(&(*right_coupling_tensor)(0,0,j), bra_rightnstates*ket_rightnstates, buf3.data());
              blas::transpose(buf3.data(), ket_rightnstates, bra_rightnstates, &(*right_coupling_tensor)(0,0,j));
            }
          }
          
          // contraction
          assert(rdm_mat->size() == site_transition_tensor->extent(2) * left_coupling_tensor->extent(2) * right_coupling_tensor->extent(2));
          auto tmpmat = make_shared<Matrix>(site_transition_tensor->extent(1)*site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          contract(1.0, group(*site_transition_tensor,1,3), {2,0}, group(*left_coupling_tensor,0,2), {2,1}, 0.0, *tmpmat, {0,1});
          btas::CRange<3> tmprange(site_transition_tensor->extent(1), site_transition_tensor->extent(2), left_coupling_tensor->extent(2));
          auto intermediate_tensor = make_shared<const btas::Tensor3<double>>(tmprange, (tmpmat->storage()));
          auto tmp_rdm_mat = make_shared<Matrix>(site_transition_tensor->extent(2)*left_coupling_tensor->extent(2), right_coupling_tensor->extent(2));
          contract(1.0, group(*intermediate_tensor,1,3), {2,0}, group(*right_coupling_tensor,0,2), {2,1}, 0.0, *tmp_rdm_mat, {0,1});
          const double sign = static_cast<double>(1 - (((ket_leftinfo.nelea + ket_leftinfo.neleb) % 2) << 1));
          blas::ax_plus_y_n(sign, tmp_rdm_mat->data(), tmp_rdm_mat->size(), rdm_mat->data());
        } // end of looping over one DMRG blockpair
      } // end of looping over sector with blockkey
    } // end of looping over gammalist_tuple

    // copy data into rdm2_
    auto rdm2_target = rdm2_->at(istate);
    for (int q = 0; q != norb_right; ++q) {
      for (int p = 0; p != norb_left; ++p) {
        for (int j = 0; j != norb_site; ++j) {
          for (int i = 0; i != norb_site; ++i) {
            const double value = *rdm_mat->element_ptr(i+j*norb_site, p+q*norb_left);
            rdm2_target->element(i+norb_left, j+norb_left, p, q+norb_left+norb_site) = value;
            rdm2_target->element(p, q+norb_left+norb_site, i+norb_left, j+norb_left) = value;
            rdm2_target->element(q+norb_left+norb_site, p, j+norb_left, i+norb_left) = value;
            rdm2_target->element(j+norb_left, i+norb_left, q+norb_left+norb_site, p) = value;
          }
        }
      }
    }

  } // end of looping over istate

}


