//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_orbopt.h
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang
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

#ifndef __ASD_DMRG_ORBOPT_H
#define __ASD_DMRG_ORBOPT_H

#include <src/wfn/method.h>
#include <src/asd/multisite/multisite.h>

namespace bagel {

class ASD_DMRG_Orbopt : public std::enable_shared_from_this<ASD_DMRG_Orbopt> {

  protected:
  int nclosed_;
  int nact_;
  int nocc_; // sum of nclosed_ + nact_
  int nvirt_;
  int nmo_;
  int nstate_;

  int max_iter_;
  int max_micro_iter_;
  double thresh_;
  double thresh_micro_;
  
  std::shared_ptr<MultiSite> multisite_;
  std::shared_ptr<const Coeff> coeff_;
  std::vector<double> energy_;
  std::shared_ptr<const Reference> ref_;

  void print_header() const;
  void common_init();

  public:
  ASD_DMRG_Orbopt(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> iref);

  void compute();

  // return functions
  std::shared_ptr<const Reference> conv_to_ref() const { std::cout << "111" << std::endl; return ref_; }

  int nclosed() const { return nclosed_; }
  int nact() const { return nact_; }
  int nocc() const { return nocc_; }
  int nvirt() const { return nvirt_; }
  int nmo() const { return nmo_; }
  int nstate() const { return nstate_; }
  int max_iter() const { return max_iter_; }
  int max_micro_iter() const { return max_micro_iter_; }
  double thresh() const { return thresh_; }
  double thresh_micro() const { return thresh_micro_; }
  
  double energy(const int i) const { return energy_[i]; }
  double energy_av() const { return blas::average(energy_); }
  const std::vector<double>& energy() const { return energy_; }
  std::shared_ptr<const Coeff> coeff() const { return coeff_; }

};

}

#endif
