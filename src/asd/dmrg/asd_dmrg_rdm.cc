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

  // compute ProdRASCivec after convergence TODO might need one more sweeping to the right for RDM<2>
  const int site = nsites_/2;
  auto left_block = left_blocks_[site-1];
  auto right_block = (nsites_==2) ? nullptr : right_blocks_[nsites_-site-2];
  {
    Muffle hide_cout("asd_dmrg_rdm.log", false);
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
    if (!right_block) environment = left_block;
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
    compute_rdm12(i, left_block, right_block);

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


void ASD_DMRG::compute_rdm12(const int ist, const shared_ptr<const DMRG_Block1> left, const shared_ptr<const DMRG_Block1> right) {
  shared_ptr<ProductRASCivec> ccbra = cc_[ist];

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_prodcivec(ccbra, left, right);

  rdm1_->emplace(ist, ist, rdm1);
  rdm2_->emplace(ist, ist, rdm2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> 
  ASD_DMRG::compute_rdm12_from_prodcivec(shared_ptr<const ProductRASCivec> cbra, shared_ptr<const DMRG_Block1> left, shared_ptr<const DMRG_Block1> right) const {

  auto rdm1 = make_shared<RDM<1>>(nactorb_);
  auto rdm2 = make_shared<RDM<2>>(nactorb_);

  const int site = nsites_/2;
  cout << "  * fragment " << site << " is used as site in RDM calculation..." << endl;

  vector<int> active_sizes = multisite_->active_sizes();
  const int site_offset = accumulate(active_sizes.begin(), active_sizes.begin() + site, 0);
  cout << "site_offset = " << site_offset << endl;
  // for now we assume we have right block, make it work first! TODO modify the code
  const int right_offset = site_offset + active_sizes[site];
  cout << "right_offset = " << right_offset << endl;

  // RDM1, calculate the lower half then fill_upper
  shared_ptr<const DMRG_Block2> block2 = dynamic_pointer_cast<const DMRG_Block2>(cbra->left());
  auto left_block = block2->left_block();
  auto right_block = block2->right_block();
  ///< both indices in left block
  for (auto sec : cbra->sectors()) {
    BlockKey bkey = sec.first;
    auto civec = sec.second;
    for (auto& sourcepair : block2->blockpairs(bkey)) {
      BlockInfo leftinfo = sourcepair.left;
      BlockInfo rightinfo = sourcepair.right;
      const int lstates = leftinfo.nstates;
      const int rstates = rightinfo.nstates;
      // transition density matrix
      auto rdm1_ll = make_shared<Matrix>(site_offset, site_offset);
      pair<BlockKey, BlockKey> left_cbkey = {leftinfo, leftinfo};
      list<GammaSQ> gammalist_alpha = {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha};
      list<GammaSQ> gammalist_beta =  {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta};
      shared_ptr<const btas::Tensor3<double>> tran_tensor_alpha = left_block->coupling(gammalist_alpha).at(left_cbkey).data;
      shared_ptr<const btas::Tensor3<double>> tran_tensor_beta = left_block->coupling(gammalist_beta).at(left_cbkey).data;
      auto transition_tensor = make_shared<const btas::Tensor3<double>>(*tran_tensor_alpha + *tran_tensor_beta);
      // contract DMRG coefficient tensor into {l l'}
      auto left_contract = make_shared<Matrix>(lstates, lstates);
      {
        auto contract_cr = [&civec, &lstates, &rstates] (const int j, const int i) {
          double result = 0.0;
          for (int iter = 0; iter != rstates; ++iter) {
            auto civec1 = civec->civec(lstates*iter+j);
            auto civec2 = civec->civec(lstates*iter+i);
            result += blas::dot_product(civec1.data(), civec->ndim(), civec2.data());
          }
          return result;
        };
        for (int ii = 0; ii != lstates; ++ii) {
          for (int jj = ii + 1; jj != lstates; ++jj) {
            *left_contract->element_ptr(jj, ii) = contract_cr(jj, ii);
          }
          *left_contract->element_ptr(ii, ii) = contract_cr(ii, ii);
        }
      }
      // fill in RDM<1>
      auto rdm1tmp = btas::group(*rdm1_ll,0,2);
      btas::contract(1.0, group(*transition_tensor, 0, 2), {0,1}, group(*left_contract, 0, 2), {0}, 0.0, rdm1tmp, {1});
      // TODO find a good way to add partial RDM to RDM<1>
    }
  }
  ///< one orbital on site, the other in left block
  const int nelea_tot = cbra->nelea();
  const int neleb_tot = cbra->neleb();
  cout << "nelea_tot = " << nelea_tot << endl;
  cout << "neleb_tot = " << nelea_tot << endl;
  for (auto sec : cbra->sectors()) {
    BlockKey bkey = sec.first;
    auto civec = sec.second;
    const int nelea_ci = nelea_tot - bkey.nelea;
    const int neleb_ci = neleb_tot - bkey.neleb;
    cout << "nelea_ci = " << nelea_ci << endl;
    cout << "neleb_ci = " << neleb_ci << endl;
    // left block wave function as ket
    for (auto& sourcepair : block2->blockpairs(bkey)) {
      BlockInfo leftinfo = sourcepair.left;
      BlockInfo rightinfo = sourcepair.right;
      const int lstates = leftinfo.nstates;
      const int rstates = rightinfo.nstates;
      cout << "lstates = " << lstates << endl;
      cout << "rstates = " << rstates << endl;
      cout << "left.nelea = " << leftinfo.nelea << endl;
      cout << "left.neleb = " << leftinfo.neleb << endl;
      cout << "right.nelea = " << rightinfo.nelea << endl;
      cout << "right.neleb = " << rightinfo.neleb << endl;
      BlockKey coupled_bkey_alpha = {leftinfo.nelea+1, leftinfo.neleb};
      BlockKey coupled_bkey_beta = {leftinfo.nelea, leftinfo.neleb+1};
      list<GammaSQ> gammalist_alpha = {GammaSQ::CreateAlpha};
      list<GammaSQ> gammalist_beta = {GammaSQ::CreateBeta};
      pair<BlockKey, BlockKey> left_cbkey_alpha = {coupled_bkey_alpha, leftinfo};
      pair<BlockKey, BlockKey> left_cbkey_beta = {coupled_bkey_beta, leftinfo};
      if (left_block->contains(coupled_bkey_alpha)) {
        shared_ptr<const btas::Tensor3<double>> tran_tensor_alpha = left_block->coupling(gammalist_alpha).at(left_cbkey_alpha).data;
        cout << "alpha ext0 = " << tran_tensor_alpha->extent(0) << ", ext1 = " << tran_tensor_alpha->extent(1) << ", ext2 = " << tran_tensor_alpha->extent(2) << endl;
      }
      if (left_block->contains(coupled_bkey_beta)) {
        shared_ptr<const btas::Tensor3<double>> tran_tensor_beta = left_block->coupling(gammalist_beta).at(left_cbkey_beta).data;
        cout << "beta ext0 = " << tran_tensor_beta->extent(0) << ", ext1 = " << tran_tensor_beta->extent(1) << ", ext2 = " << tran_tensor_beta->extent(2) << endl;
      }
    }
  }

  return tie(rdm1, rdm2);
}


