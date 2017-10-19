//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_rdm.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu>
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


#include <src/ci/ras/rasci.h>

using namespace std;
using namespace bagel;

void RASCI::compute_rdm12() {
  // Needs initialization here because we use daxpy.
  if (rdm1_av_ == nullptr && nstate_ > 1) {
    rdm1_av_ = make_shared<RDM<1>>(norb_);
    rdm2_av_ = make_shared<RDM<2>>(norb_);
  } else if (nstate_ > 1) {
    rdm1_av_->zero();
    rdm2_av_->zero();
  }

  for (int i = 0; i != nstate_; ++i)
    compute_rdm12(i, i);

  // calculate state averaged RDMs
  if (nstate_ != 1) {
    for (int ist = 0; ist != nstate_; ++ist) {
      rdm1_av_->ax_plus_y(weight_[ist], rdm1_->at(ist));
      rdm2_av_->ax_plus_y(weight_[ist], rdm2_->at(ist));
    }
  } else {
    rdm1_av_ = rdm1_->at(0,0);
    rdm2_av_ = rdm2_->at(0,0);
  }
}


void RASCI::compute_rdm12(const int ist, const int jst) {
  shared_ptr<RASCivec> rasdbra = cc_->data(ist);
  shared_ptr<RASCivec> rasdket = cc_->data(jst);

  shared_ptr<RDM<1>> rdm1;
  shared_ptr<RDM<2>> rdm2;
  tie(rdm1, rdm2) = compute_rdm12_from_rascivec(rasdbra, rasdket);

  rdm1_->emplace(ist, jst, rdm1);
  rdm2_->emplace(ist, jst, rdm2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  RASCI::compute_rdm12_from_rascivec(shared_ptr<const RASCivec> cbra, shared_ptr<const RASCivec> cket) const {

  auto dbra = make_shared<RASDvec>(cbra->det(), norb_*norb_);
  excite_alpha(cbra, dbra);
  excite_beta(cbra, dbra);

  shared_ptr<RASDvec> dket;
  if(cbra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb_*norb_);
    excite_alpha(cket, dket);
    excite_beta(cket, dket);
  } else {
    dket = dbra;
  }

  return compute_rdm12_last_step(dbra, dket, cbra);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  RASCI::compute_rdm12_last_step(shared_ptr<const RASDvec> dbra, shared_ptr<const RASDvec> dket, shared_ptr<const RASCivec> cibra) const {

  const int nri = cibra->det()->size();
  const int ij = norb_ * norb_;

  auto rdm1 = make_shared<RDM<1>>(norb_);
  auto rdm2 = make_shared<RDM<2>>(norb_);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, ij);
    for (int i = 0; i != norb_; ++i)
      for (int j = 0; j != norb_; ++j)
        copy_n(dket->data(j+i*norb_)->data(), nri, dket_data->element_ptr(0, i+j*norb_));
    auto rdm1tmp = btas::group(*rdm1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, rdm1tmp, {1});

    auto dbra_data = dket_data->clone();
    for (int k = 0; k != ij; ++k)
      copy_n(dbra->data(k)->data(), nri, dbra_data->element_ptr(0, k));
  }

  return tie(rdm1, rdm2);
}


void RASCI::excite_alpha(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  assert(cc->det() == d->det());
  auto det = cc->det();
  const double* const source_base = cc->data();

  for (auto& source_block : det->blockinfo()) {
    // source info
    auto source_astrings = source_block->stringsa();
    auto bstrings = source_block->stringsb();
    const size_t lenb = bstrings->size();
    const int nhb = bstrings->nholes();
    const int npb = bstrings->nparticles();
    const int source_offset = source_block->offset();

    for (auto& istring : *source_astrings) {
      // locate source
      const size_t source_lexoffset_a = source_astrings->lexical_offset(istring);
      const size_t source_lexzero_a = source_astrings->lexical_zero(istring);
      const double* const source_location = source_base + source_offset + source_lexzero_a*lenb;
      // excitation info
      for (auto& iterij : det->phia(source_lexoffset_a)) {
        // locate target
        const size_t ij = iterij.ij;
        const int sign = iterij.sign;
        const size_t target_lexoffset_a = iterij.source;
        const bitset<nbit__> target_astring = det->stringspacea()->strings(target_lexoffset_a);
        auto target_astrings = det->stringspacea()->find_string(target_astring);
        const size_t target_lexzero_a = target_astrings->lexical_zero(target_astring);
        const int target_nha = target_astrings->nholes();
        const int target_npa = target_astrings->nparticles();
        if (!det->allowed(target_nha, nhb, target_npa, npb)) continue; // ignore illegal excitations
        auto target_blockinfo = det->blockinfo(target_nha, nhb, target_npa, npb);
        const int target_offset = target_blockinfo->offset();
        double* const target_location = d->data(ij)->data() + target_offset + target_lexzero_a*lenb;
        
        blas::ax_plus_y_n(sign, source_location, lenb, target_location);
      }
    }
  }
}


void RASCI::excite_beta(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {

}
