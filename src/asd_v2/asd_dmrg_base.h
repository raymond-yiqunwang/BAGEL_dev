//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_base.h
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

#ifndef __SRC_ASD_V2_ASD_DMRG_BASE_H
#define __SRC_ASD_V2_ASD_DMRG_BASE_H

#include <src/asd_v2/multimer/multimer.h>
#include <src/asd_v2/asd_dmrg_block.h>

namespace bagel {

// Base class for ASD_DMRG_METAL, not to be confused with previous version of ASD_DMRG
class ASD_DMRG_base {
  protected:
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Multimer> multimer_;

    // ASD_DMRG_Block representing L block containing l sites is stored in left_blocks_[l-1]
    std::vector<std::shared_ptr<ASD_DMRG_Block1>> left_blocks_;
    // ASD_DMRG_Block representing R block containing r sites is stored in right_blocks_[r-1]
    std::vector<std::shared_ptr<ASD_DMRG_Block1>> right_blocks_;

    std::vector<double> weights_; // weights to use when building RDM

    std::vector<std::vector<double>> sweep_energies_; // stores the energies of each states for each step of the sweep
    std::vector<double> energies_; // final energies

    int nsites_;  // number of sites(fragments) in DMRG model
    int nstate_;  // number of target states
    int maxiter_; // maximum number of full sweeps to perform
    int ntrunc_;  // number of states to keep in each DMRG block

    double thresh_; // TODO add description
    double perturb_;
    double perturb_thresh_;
    double perturb_min_;

    // print graphical depiction of sweep process
    std::string print_progress(const int position, const std::string left_symbol, const std::string right_symbol) const;

/*    // implemented in CAS/RAS models, respectively
    virtual std::shared_ptr<ASD_DMRG_Block1> compute_first_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref) = 0;
    // add one site to the block
    virtual std::shared_ptr<ASD_DMRG_Block1> grow_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref,
                                                          std::shared_ptr<ASD_DMRG_Block1> left, const int site) = 0;
    virtual std::shared_ptr<ASD_DMRG_Block1> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref,
                                                              std::shared_ptr<ASD_DMRG_Block1> system, std::shared_ptr<ASD_DMRG_Block1> environment, const int site) = 0;
*/
  private:
    // prepare several input files used for growing the chain
    std::vector<std::shared_ptr<PTree>> prepare_growing_input(const int site) const;
    // prepare one input to be used during the sweep
    std::shared_ptr<PTree> prepare_sweeping_input(const int site) const;
  
  public:
    ASD_DMRG_base(const std::shared_ptr<const PTree> input, std::shared_ptr<const Multimer> multimer);

    void compute();
    
    const std::vector<double>& energies() const { return energies_; }
    double energies(const int i) const { return energies_.at(i); }

    // utility functions
    std::shared_ptr<const Multimer> multimer() const { return multimer_; }
};

} 

#endif
