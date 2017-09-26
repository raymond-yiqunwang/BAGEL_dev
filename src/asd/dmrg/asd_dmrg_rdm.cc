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
      compute_rdm2_ras(cc, site);
    } else {
      // special treatment for first configuration as described by Garnet Chan, 2008
      if (site == 1) {
        cout << "  * special treatment for site[1]" << endl;
        // compute_310
        compute_rdm2_310(cc);
      }

      // general treatment
      {
        cout << "  * general treatment for site[" << site << "]" << endl;
        // compute RASCI RDM
        compute_rdm2_ras(cc, site);

        // compute_130
        compute_rdm2_130(cc, site);
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
#if 1
  for (int i = 0; i != nstate_; ++i) {
    cout << "rdm2_[" << i << "] : " << endl;
    rdm2_->at(i)->print(1e-6);
  }
#endif

}


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

            // transpose the first two dimensions
            btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            auto transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
            auto transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
            unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
            for (int i = 0; i != transition_tensor->extent(2); ++i) {
              copy_n(&(*transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
              blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*transition_tensor)(0,0,i));
            }
            
            auto coupling_data = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;
            
            auto target = rdm_mat->clone();
            assert (target->size() == transition_tensor->extent(2) * coupling_data->extent(2));
            contract(1.0, group(*transition_tensor,0,2), {2,0}, group(*coupling_data,0,2), {2,1}, 0.0, *target, {0,1});
            const double sign = static_cast<double>(1 - (((bra_ras_nelea + bra_ras_neleb)%2) << 1));
            blas::ax_plus_y_n(sign, target->data(), target->size(), rdm_mat->data());
          }
        }
      }
    }
    
    // copy data into rdm2_, looks stupid but leave it here for the moment 
    auto rdm2_target = rdm2_->at(istate);
    for (int p = 0; p != rdm_mat->mdim(); ++p)
      for (int i = 0; i != norb_site; ++i)
        for (int j = 0; j != norb_site; ++j)
          for (int k = 0; k != norb_site; ++k)
            rdm2_target->element(p,i+norb_left,k+norb_left,j+norb_left) = *rdm_mat->element_ptr(i + norb_site*j + norb_site*norb_site*k, p);

  } // end of looping over nstates
}


// \Gamma_{qprs}^{310} = \sum_{\{c\},\{lr\}} \sum_{\{c'\},\{l'r'\}} (-1)^{n_c} 
//     \left< A^{c}_{lr} | a_{q \rho} | A^{c'}_{l'r'} \right> \left< l \right| 
//     a^{\dagger}_{p \rho} a^{\dagger}_{r \sigma} a_{s \sigma} \left| l' \right>
void ASD_DMRG::compute_rdm2_310(vector<shared_ptr<ProductRASCivec>> dvec) {
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
            btas::CRange<3> range(bra_nstates, ket_nstates, lrint(pow(norb_site, get<0>(gammalist_tuple).size())));
            auto transition_mat = forest.template get<0>(ket_rastag, bra_rastag, get<0>(gammalist_tuple));
            auto transition_tensor = make_shared<btas::Tensor3<double>>(range, move(transition_mat->storage()));
            unique_ptr<double[]> buf(new double[bra_nstates*ket_nstates]);
            for (int i = 0; i != transition_tensor->extent(2); ++i) {
              copy_n(&(*transition_tensor)(0,0,i), ket_nstates*bra_nstates, buf.get());
              blas::transpose(buf.get(), ket_nstates, bra_nstates, &(*transition_tensor)(0,0,i));
            }

            auto coupling_data = left_block->coupling(get<1>(gammalist_tuple)).at(make_pair(bra_leftkey, ket_leftkey)).data;

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
    for (int l = 0; l != norb_left; ++l)
      for (int k = 0; k != norb_left; ++k)
        for (int i = 0; i != norb_left; ++i)
          for (int j = 0; j != norb_site; ++j)
            rdm2_target->element(i,j+norb_left,k,l) = *rdm_mat->element_ptr(i + norb_left*k + norb_left*norb_left*l, j);
  }

}


