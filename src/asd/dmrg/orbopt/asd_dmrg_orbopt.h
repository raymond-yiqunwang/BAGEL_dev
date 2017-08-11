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

#include <src/asd/dmrg/rasd.h>
#include <src/asd/dmrg/orbopt/asd_dmrg_rotfile.h>

namespace bagel {

class ASD_DMRG_Orbopt : public std::enable_shared_from_this<ASD_DMRG_Orbopt> {

  protected:
    int nclosed_;
    int nact_;
    int nocc_; // sum of nclosed_ + nact_
    int nvirt_;
    int nmo_;
    int nstate_;
  
    // parameters for iteration
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;
   
    std::shared_ptr<RASD> asd_dmrg_; // should have DMRG member
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const PTree> asd_info_; // information for doing ASD-DMRG
    std::shared_ptr<MultiSite> multisite_;
    std::shared_ptr<const Coeff> coeff_;
    std::vector<double> energy_;
    std::shared_ptr<const Reference> ref_;
  
    // util functions
    void print_header() const;
    void common_init();
  
    // second-order optimization functions
    // compute orbital gradient
    std::shared_ptr<ASD_DMRG_RotFile> compute_gradient(std::shared_ptr<const Matrix> cfock, std::shared_ptr<const Matrix> afock,
                                                       std::shared_ptr<const Matrix> qxr) const;
  
    // compute exact diagonal Hessian
    std::shared_ptr<ASD_DMRG_RotFile> compute_denom(std::shared_ptr<const DFHalfDist> half, std::shared_ptr<const DFHalfDist> half_1j,
                                                    std::shared_ptr<const DFHalfDist> halfa, std::shared_ptr<const Matrix> cfock,
                                                    std::shared_ptr<const Matrix> afock) const;
  
    // compute Hessian times trial vector
    std::shared_ptr<ASD_DMRG_RotFile> compute_hess_trial(std::shared_ptr<const ASD_DMRG_RotFile> trot, std::shared_ptr<const DFHalfDist> hafl,
                                                         std::shared_ptr<const DFHalfDist> halfa, std::shared_ptr<const Matrix> cfock,
                                                         std::shared_ptr<const Matrix> afock, std::shared_ptr<const Matrix> qxr) const;
  
    // apply denominator in micro-iterations
    std::shared_ptr<ASD_DMRG_RotFile> apply_denom(std::shared_ptr<const ASD_DMRG_RotFile> grad, std::shared_ptr<const ASD_DMRG_RotFile> denom, 
                                                  const double shift, const double scale) const;

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
