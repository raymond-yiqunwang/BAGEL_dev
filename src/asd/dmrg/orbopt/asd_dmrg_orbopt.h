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

namespace bagel {

class ASD_DMRG_Orbopt : public Method/*, public : std::enable_shared_from_this<ASD_DMRG_Orbopt> */{

  protected:
  int nclosed_;
  int nact_;
  int nocc_;
  int nvirt_;
  int nbasis_;
  int nstate_;

  public:
  ASD_DMRG_Orbopt(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> input_ref);

  void compute() override;

  // return functions
  std::shared_ptr<const Reference> conv_to_ref() const override { std::cout << "111" << std::endl; return ref_; }

};

}

#endif
