//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/civector.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


#include <iomanip>

#include <src/util/taskqueue.h>
#include <src/ci/ras/civector.h>
#include <src/ci/ras/civec_spinop.h>
#include <src/ci/ras/apply_block.h>
#include <src/ci/fci/fci.h>

using namespace std;
using namespace bagel;


template<typename DataType>
RASCivector<DataType>::RASCivector(shared_ptr<const RASDeterminants> det) : RASCivector_impl<DataType>(det) {
  data_ = unique_ptr<DataType[]>(new DataType[size()]);
  fill_n(data_.get(), size(), 0.0);
  size_t sz = 0;
  for (auto& ipair : det->blockinfo())
    if (!ipair->empty()) {
      blocks_.push_back(make_shared<RBlock>(ipair->stringsa(), ipair->stringsb(), data_.get()+sz, sz));
      sz += blocks_.back()->size();
    } else {
      blocks_.push_back(nullptr);
    }
}


template<typename DataType>
RASCivector<DataType>& RASCivector<DataType>::operator=(RASCivector<DataType>&& o) {
  assert(*o.det() == *det());
  data_ = move(o.data_);
  blocks_ = move(o.blocks_);
  return *this;
}


template<typename DataType>
shared_ptr<RASCivector<DataType>> RASCivector<DataType>::apply(const int orbital, const bool action, const bool spin) const {
  // action: true -> create; false -> annihilate
  // spin: true -> alpha; false -> beta
  shared_ptr<const RASDeterminants> sdet = this->det();

  const int ras1 = sdet->ras(0);
  const int ras2 = sdet->ras(1);
  const int ras3 = sdet->ras(2);

  // 0 -> RASI, 1 -> RASII, 2 -> RASIII
  const int ras_space = ( orbital >= ras1 ) + (orbital >= ras1 + ras2);

  auto to_array = [] (shared_ptr<const RASBlock<DataType>> block) {
    auto sa = block->stringsa();
    auto sb = block->stringsb();
    return array<int, 6>({sa->nholes(), sb->nholes(), sa->nele2(), sb->nele2(), sa->nparticles(), sb->nparticles()});
  };

  auto op_on_array = [&ras_space, &action, &spin] ( array<int, 6> in ) {
    const int mod = ( action ? +1 : -1 ) * ( ras_space == 0 ? -1 : 1 );
    array<int, 6> out = in;
    out[2*ras_space] += (spin ? mod : 0);
    out[2*ras_space+1] += (spin ? 0 : mod);
    return out;
  };

  RAS::Apply_block apply_block(orbital, action, spin);

  const int mod = action ? +1 : -1;
  const int telea = sdet->nelea() + (spin ? mod : 0);
  const int teleb = sdet->neleb() + (spin ? 0 : mod);
  const int tholes = max(sdet->max_holes() - ((ras_space == 0) ? mod : 0), 0);
  const int tparts = max(sdet->max_particles() + ((ras_space == 2) ? mod : 0), 0);

  auto tdet = make_shared<const RASDeterminants>(ras1, ras2, ras3, telea, teleb, tholes, tparts, true);
  auto out = make_shared<RASCivector<DataType>>(tdet);

  for (shared_ptr<const RASBlock<double>> soblock : this->blocks()) {
    if (!soblock) continue;
    array<int, 6> tar_array = op_on_array(to_array(soblock));
    if (all_of(tar_array.begin(), tar_array.end(), [] (int i) { return i >= 0; })) {
      shared_ptr<RASBlock<double>> tarblock = out->block(tar_array[0], tar_array[1], tar_array[4], tar_array[5]);
      if (tarblock) apply_block(soblock, tarblock, false);
    }
  }
  return out;
}


template<typename DataType>
RASCivecView_<DataType>::RASCivecView_(shared_ptr<const RASDeterminants> det, double* const data)
 : RASCivector_impl<DataType>(det), data_ptr_(data), can_write_(true) {
  size_t sz = 0;
  for (auto& ipair : det->blockinfo()) {
    if (!ipair->empty()) {
      blocks_.push_back(make_shared<RBlock>(ipair->stringsa(), ipair->stringsb(), data+sz, sz));
      sz += blocks_.back()->size();
    } else {
      blocks_.push_back(nullptr);
    }
  }
}


template<>
shared_ptr<RASCivector<double>> RASCivector<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);
  RAS::spin_impl(*this, *out);
  return out;
}

template<>
shared_ptr<RASCivector<double>> RASCivecView_<double>::spin() const {
  auto out = make_shared<RASCivector<double>>(det_);
  RAS::spin_impl(*this, *out);
  return out;
}


template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert((tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_lower(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()-1, sdet->neleb()+1);
  assert((tdet->nelea() == sdet->nelea()-1) && (tdet->neleb() == sdet->neleb()+1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_lower_impl(*this, *out);
  return out;
}


template<> shared_ptr<RASCivector<double>> RASCivector<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert((tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}

template<> shared_ptr<RASCivector<double>> RASCivecView_<double>::spin_raise(shared_ptr<const RASDeterminants> tdet) const {
  shared_ptr<const RASDeterminants> sdet = det_;
  if (!tdet) tdet = sdet->clone(sdet->nelea()+1, sdet->neleb()-1);
  assert((tdet->nelea() == sdet->nelea()+1) && (tdet->neleb() == sdet->neleb()-1));
  auto out = make_shared<RASCivec>(tdet);
  RAS::spin_raise_impl(*this, *out);
  return out;
}


template<>
tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>>
  RASCivector<double>::compute_rdm12_from_rascivec(shared_ptr<const RASCivec> cket) const {
  
  auto cibra = shared_from_this();
  const int norb = cibra->det()->norb();
  const int nri = cibra->det()->size();
  const int ij = norb * norb;

  auto dbra = make_shared<RASDvec>(cibra->det(), norb*norb);
  excite_alpha(cibra, dbra);
  excite_beta(cibra, dbra);

  shared_ptr<RASDvec> dket;
  if(cibra != cket) {
    dket = make_shared<RASDvec>(cket->det(), norb*norb);
    excite_alpha(cket, dket);
    excite_beta(cket, dket);
  } else {
    dket = dbra;
  }

  auto rdm1 = make_shared<RDM<1>>(norb);
  auto rdm2_raw = make_shared<RDM<2>>(norb);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, ij);
    for (int i = 0; i != norb; ++i)
      for (int j = 0; j != norb; ++j)
        copy_n(dket->data(j+i*norb)->data(), nri, dket_data->element_ptr(0, i+j*norb));
    auto rdm1tmp = btas::group(*rdm1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, rdm1tmp, {1});

    auto dbra_data = dket_data->clone();
    for (int k = 0; k != ij; ++k)
      copy_n(dbra->data(k)->data(), nri, dbra_data->element_ptr(0, k));
    auto rdm2tmp = group(group(*rdm2_raw,2,4),0,2);
    btas::contract(1.0, *dbra_data, {2,0}, *dket_data, {2,1}, 0.0, rdm2tmp, {0,1});
  }

  auto rdm2 = make_shared<RDM<2>>(norb);
  // collect ij >= kl from rdm2_raw
  for (int k = 0; k != norb; ++k)
    for (int j = 0; j != norb; ++j)
      for (int l = 0; l != norb; ++l)
        for (int i = 0; i != norb; ++i) 
          if ((i-1)*norb+j >= (k-1)*norb+l)
            // tensor(i,j,k,l) == rdm(i,l,j,k)
            rdm2->element(i,l,j,k) = rdm2_raw->element(i,l,j,k);

  for (int k = 0; k != norb; ++k)
    for (int j = 0; j != norb; ++j)
      for (int l = 0; l != norb; ++l)
        for (int i = 0; i != norb; ++i)
          if ((i-1)*norb+j < (k-1)*norb+l)
            // tensor(i,j,k,l) == tensor(k,l,i,j)
            // rdm(i,l,j,k) == rdm(k,j,l,i)
            rdm2->element(i,l,j,k) = rdm2->element(k,j,l,i);

  for (int i = 0; i != norb; ++i)
    for (int k = 0; k != norb; ++k)
      for (int j = 0; j != norb; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}

template<>
void RASCivector<double>::excite_alpha(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
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


template<>
void RASCivector<double>::excite_beta(shared_ptr<const RASCivec> cc, shared_ptr<RASDvec> d) const {
  auto det = cc->det();
  const double* const source_base = cc->data();

  for (auto& source_block : det->blockinfo()) {
    // source info
    auto source_astrings = source_block->stringsa();
    auto source_bstrings = source_block->stringsb();
    const size_t source_lenb = source_bstrings->size();
    const int nha = source_astrings->nholes();
    const int npa = source_astrings->nparticles();
    const int source_offset = source_block->offset();

    for (auto& istringa : *source_astrings) {
      const size_t lexzero_a = source_astrings->lexical_zero(istringa);
      for (auto& istringb : *source_bstrings) {
        const size_t source_lexoffset_b = source_bstrings->lexical_offset(istringb);
        const size_t source_lexzero_b = source_bstrings->lexical_zero(istringb);
        // locate source
        const double* const source_location = source_base + source_offset + lexzero_a*source_lenb + source_lexzero_b;
        // excitation info
        for (auto& iterij : det->phib(source_lexoffset_b)) {
          // locate target
          const size_t ij = iterij.ij;
          const double sign = static_cast<double>(iterij.sign);
          const size_t target_lexoffset_b = iterij.source;
          const bitset<nbit__> target_bstring = det->stringspaceb()->strings(target_lexoffset_b);
          auto target_bstrings = det->stringspaceb()->find_string(target_bstring);
          const int target_lenb = target_bstrings->size();
          const size_t target_lexzero_b = target_bstrings->lexical_zero(target_bstring);
          const int target_nhb = target_bstrings->nholes();
          const int target_npb = target_bstrings->nparticles();
          if (!det->allowed(nha, target_nhb, npa, target_npb)) continue;
          auto target_blockinfo = det->blockinfo(nha, target_nhb, npa, target_npb);
          const int target_offset = target_blockinfo->offset();
          double* const target_location = d->data(ij)->data() + target_offset + lexzero_a*target_lenb + target_lexzero_b;
          
          *target_location += sign * *source_location;
        }
      }
    }
  }
}


template class bagel::RASCivector<double>;
template class bagel::RASCivecView_<double>;
