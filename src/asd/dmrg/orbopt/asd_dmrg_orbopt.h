//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_orbopt.h
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

#ifndef __ASD_DMRG_ORBOPT_H
#define __ASD_DMRG_ORBOPT_H

#include <src/asd/dmrg/rasd.h>
#include <src/asd/dmrg/orbopt/asd_dmrg_rotfile.h>

namespace bagel {

class ASD_DMRG_OrbOpt : public std::enable_shared_from_this<ASD_DMRG_OrbOpt> {

  protected:
    int naa_; // size of active-active part
    // parameters for iteration
    int max_iter_;
    int max_micro_iter_;
    double thresh_;
    double thresh_micro_;
   
    std::shared_ptr<RASD> asd_dmrg_;
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Coeff> coeff_;

#ifdef AAROT
    // active-active rotation parameters
    std::vector<ASD_ActRotBlock> act_rotblocks_;
#endif

    std::shared_ptr<const Coeff> semi_canonical_orb() const;

    // util functions
    void print_header() const;
    void print_iteration(const int iter, const std::vector<double>& energy, const double error) const;

  public:
    ASD_DMRG_OrbOpt(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> iref);
  
    virtual void compute() = 0;
  
    // return functions
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }

    std::shared_ptr<Matrix> compute_active_fock(const MatView acoeff, std::shared_ptr<const RDM<1>> rdm1) const;
    std::shared_ptr<Matrix> compute_qvec(const MatView acoeff, std::shared_ptr<const RDM<2>> rdm2) const;
};

}

#endif
