//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multimer.h 
// Copyright (C) 2017 Raymond Wang 
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

#ifndef __SRC_ASD_V2_MULTIMER_MULTIMER_H
#define __SRC_ASD_V2_MULTIMER_MULTIMER_H

#include <vector>
#include <src/wfn/construct_method.h>

namespace bagel {

// preparation step for ASD calculation, computes the MO coefficient matrix and cispace

class Multimer : public std::enable_shared_from_this<Multimer> {
  protected:
    std::shared_ptr<const Geometry> geom_; // only super geometry is provided for generality
    
    const std::shared_ptr<const Reference> prev_ref_;
    std::shared_ptr<const Reference> ref_;
    std::shared_ptr<const Reference> rhf_ref_;
    std::shared_ptr<const Reference> active_ref_;

    double active_thresh_;
//  double region_thresh_;
//  bool print_orbitals_;

    // reorder MO coeff to closed - active - virtual
    void set_active(std::shared_ptr<const PTree> idata);
    // project active orbitals to fragments
    void project_active(std::shared_ptr<const PTree> idata);
    
    // compute CI space
//    template <class VecType>
//    std::shared_ptr<MultimerCISpace<<VecType>> compute_cispace(std::shared_ptr<const PTree> idata);

  public:
    // constructor
    Multimer(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> ref);

    // utility functions
    void precompute(std::shared_ptr<const PTree> idata);

    std::shared_ptr<const Reference> ref() const { return ref_; } 

};

} // end of namespace TODO --delete

#endif
