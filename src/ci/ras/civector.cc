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

template<> tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> RASCivector<double>::compute_rdm12_from_rasvec() const {
  const int norb = det_->norb();
  const int nelea = det_->nelea();
  const int neleb = det_->neleb();
  auto rdm2 = make_shared<RDM<2>>(norb);

  auto idet = make_shared<Determinants>(norb, nelea, neleb, false/*compress*/, true/*mute*/);
  auto dbra = make_shared<Dvec>(idet, norb*norb);
  auto fcivec = make_shared<Civec>(idet);
  
  // map RASCI vec into FCI vec
  for (auto& iblock : blocks_) {
    if (!iblock) continue;
    // beta_string runs first
    for (auto& abit : iblock->stringsa()->strings()) {
      const int block_index_a = iblock->stringsa()->lexical_zero(abit);
      for (auto & bbit : iblock->stringsb()->strings()) {
        // now I have a pair of legal bits, identify the source index in RASCivec and target index in Civec
        const int block_index_b = iblock->stringsb()->lexical_zero(bbit);
        double* const source_ptr = iblock->data() + block_index_b + (block_index_a * iblock->lenb());
        // FCI part
        const int fci_index_a = fcivec->det()->lexical<0>(abit);
        const int fci_index_b = fcivec->det()->lexical<1>(bbit);
        double* const target_ptr = fcivec->data() + fci_index_b + (fci_index_a * fcivec->det()->lenb());
        *target_ptr = *source_ptr;
      }
    }
  }
  // form RDM<2>
  { // sigma_2a1
    const int lb = dbra->lenb();
    const int ij = dbra->ij();
    const double* const source_base = fcivec->data();
    for (int ip = 0; ip != ij; ++ip) {
      double* const target_base = dbra->data(ip)->data();
      for (auto& iter : fcivec->det()->phia(ip)) {
        const double sign = static_cast<double>(iter.sign);
        double* const target_array = target_base + iter.source*lb;
        blas::ax_plus_y_n(sign, source_base + iter.target*lb, lb, target_array);
      }
    }
  }
  { // sigma_2a2
    const int la = dbra->lena();
    const int ij = dbra->ij();
    for (int i = 0; i < la; ++i) {
      const double* const source_array0 = fcivec->element_ptr(0, i);
      for (int ip = 0; ip != ij; ++ip) {
        double* const target_array0 = dbra->data(ip)->element_ptr(0, i);
        for (auto& iter : fcivec->det()->phib(ip)) {
          const double sign = static_cast<double>(iter.sign);
          target_array0[iter.source] += sign * source_array0[iter.target];
        }
      }
    }
  }

  return compute_rasrdm12_last_step(dbra, dbra, fcivec);
}

template<> tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> RASCivector<double>::compute_rasrdm12_last_step(shared_ptr<Dvec> dbra, shared_ptr<Dvec> dket,
                                                                                                         shared_ptr<Civec> cibra) const {
  const int nri = cibra->asize() * cibra->lenb();
  const int norb = det_->norb();
  const int ij = norb * norb;

  auto rdm1 = make_shared<RDM<1>>(norb);
  auto rdm2 = make_shared<RDM<2>>(norb);
  {
    auto cibra_data = make_shared<VectorB>(nri);
    copy_n(cibra->data(), nri, cibra_data->data());

    auto dket_data = make_shared<Matrix>(nri, ij);
    for (int i = 0; i != ij; ++i)
      copy_n(dket->data(i)->data(), nri, dket_data->element_ptr(0, i));
    auto rdm1t = btas::group(*rdm1,0,2);
    btas::contract(1.0, *dket_data, {0,1}, *cibra_data, {0}, 0.0, rdm1t, {1});

    auto dbra_data = dket_data;
    if (dbra != dket) {
      dbra_data = make_shared<Matrix>(nri, ij);
      for (int i = 0; i != ij; ++i)
        copy_n(dbra->data(i)->data(), nri, dbra_data->element_ptr(0, i));
    }
    auto rdm2t = group(group(*rdm2, 2,4), 0,2);
    btas::contract(1.0, *dbra_data, {1,0}, *dket_data, {1,2}, 0.0, rdm2t, {0,2});
  }

  // sorting... a bit stupid but cheap anyway
  // This is since we transpose operator pairs in dgemm - cheaper to do so after dgemm (usually Nconfig >> norb_**2).
  unique_ptr<double[]> buf(new double[norb*norb]);
  for (int i = 0; i != norb; ++i) {
    for (int k = 0; k != norb; ++k) {
      copy_n(&rdm2->element(0,0,k,i), norb*norb, buf.get());
      blas::transpose(buf.get(), norb, norb, rdm2->element_ptr(0,0,k,i));
    }
  }

  // put in diagonal into 2RDM
  // Gamma{i+ k+ l j} = Gamma{i+ j k+ l} - delta_jk Gamma{i+ l}
  for (int i = 0; i != norb; ++i)
    for (int k = 0; k != norb; ++k)
      for (int j = 0; j != norb; ++j)
        rdm2->element(j,k,k,i) -= rdm1->element(j,i);

  return tie(rdm1, rdm2);
}                                                                                                        


template class bagel::RASCivector<double>;
template class bagel::RASCivecView_<double>;
