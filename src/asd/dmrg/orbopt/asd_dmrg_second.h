//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_second.h
// Copyright (C) 2017 Raymond Wang
//
// Author: Raymond Wang <raymondwang@u.northwestern.edu>
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

#ifndef __ASD_DMRG_SECOND_H
#define __ASD_DMRG_SECOND_H

#include <src/asd/dmrg/orbopt/asd_dmrg_orbopt.h>

namespace bagel {

class ASD_DMRG_Second : public ASD_DMRG_OrbOpt {
  protected:
    // convergence threshold for micro iteration relative to stepsize
    double thresh_microstep_;

    // compute orbital gradient
    std::shared_ptr<ASD_DMRG_RotFile> compute_gradient(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock,
                                                       std::shared_ptr<const Matrix> qxr) const;
  
    // compute exact diagonal Hessian
    std::shared_ptr<ASD_DMRG_RotFile> compute_denom(std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> half_1j,
                                                    std::shared_ptr<const DFHalfDist> halfa, std::shared_ptr<const DFHalfDist> halfa_JJ,
                                                    std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock) const;

    // apply denominator in micro-iterations
    std::shared_ptr<ASD_DMRG_RotFile> apply_denom(std::shared_ptr<const ASD_DMRG_RotFile> grad, std::shared_ptr<const ASD_DMRG_RotFile> denom, 
                                                  const double shift, const double scale) const;
    
    // compute Hessian times trial vector
    std::shared_ptr<ASD_DMRG_RotFile> compute_hess_trial(std::shared_ptr<const ASD_DMRG_RotFile> trot, std::shared_ptr<const DFHalfDist> half, 
                                                         std::shared_ptr<const DFHalfDist> halfa, std::shared_ptr<const DFHalfDist> halfa_JJ, 
                                                         std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock,
                                                         std::shared_ptr<const Matrix> qxr) const;

  public:
    ASD_DMRG_Second(std::shared_ptr<const PTree> idata, std::shared_ptr<const Reference> iref)
      : ASD_DMRG_OrbOpt(idata, iref) {
      std::cout << "  * Using the second-order algorithm for orbital optimization" << std::endl << std::endl;  
      // overwriting thresh_micro
      thresh_micro_ = idata->get<double>("opt_thresh_micro", thresh_*0.5);
      thresh_microstep_ = idata->get<double>("opt_thresh_microstep", 1.0e-4);
    }

    void compute() override;

    void trans_natorb_block();
    void semi_canonicalize_block();
};

}

#endif


