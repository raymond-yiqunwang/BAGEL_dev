//
// BAGEL - Parallel electron correlation program.
// Filename: HZtasks.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef __BAGEL_FCI_HZTASKAA_H
#define __BAGEL_FCI_HZTASKAA_H

#include <src/util/f77.h>

namespace bagel {

class HZTaskAA {
  protected:
    std::shared_ptr<const Civec> cc_;
    const std::bitset<nbit__> targetstring_;
    double* const target_;
    const double* const h1_;
    const double* const h2_;

  public:
    HZTaskAA(std::shared_ptr<const Civec> cc, const std::bitset<nbit__> targetstring, double* const target, const double* const h1, const double* const h2) :
      cc_(cc), targetstring_(targetstring), target_(target), h1_(h1), h2_(h2) {}

    void compute() {
      std::shared_ptr<const Determinants> det = cc_->det();

      const int lb = cc_->lenb();
      const int norb = det->norb();

      // One-electron part
      for (int i = 0; i < norb; ++i) {
        if (!targetstring_[i]) continue;
        std::bitset<nbit__> ibs = targetstring_; ibs.reset(i);

        for (int j = 0; j < norb; ++j) {
          if (ibs[j]) continue;
          std::bitset <nbit__> sourcestring = ibs; sourcestring.set(j);

          const double hc = h1_[i + j*norb] * static_cast<double>(det->sign(sourcestring, i, j));
          daxpy_(lb, hc, cc_->element_ptr(0, det->lexical<0>(sourcestring)), 1, target_, 1);
        }
      }


      // Two-electron part
      for (int i = 0; i != norb; ++i) {
        if (!targetstring_[i]) continue;
        for (int j = 0; j < i; ++j) {
          if(!targetstring_[j]) continue;
          const int ij_phase = det->sign(targetstring_,i,j);
          std::bitset<nbit__> string_ij = targetstring_;
          string_ij.reset(i); string_ij.reset(j);
          for (int l = 0; l != norb; ++l) {
            if (string_ij[l]) continue;
            for (int k = 0; k < l; ++k) {
              if (string_ij[k]) continue;
              const int kl_phase = det->sign(string_ij,l,k);
              const double phase = -static_cast<double>(ij_phase*kl_phase);
              std::bitset<nbit__> string_ijkl = string_ij;
              string_ijkl.set(k); string_ijkl.set(l);
              const double temp = phase * h2_[l + k*norb + j*norb*norb + i*norb*norb*norb];
              const double* source = cc_->element_ptr(0, det->lexical<0>(string_ijkl));
              daxpy_(lb, temp, source, 1, target_, 1);
            }
          }
        }
      }
    }

};

}

#endif