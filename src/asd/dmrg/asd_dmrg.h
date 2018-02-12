//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg.h
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __ASD_DMRG_ASD_DMRG_H
#define __ASD_DMRG_ASD_DMRG_H

#include <src/asd/dmrg/dmrg_block.h>

namespace bagel {

/// Base class for DMRG using ASD
class ASD_DMRG {
  protected:
    std::shared_ptr<const PTree> input_;
    std::shared_ptr<const Reference> sref_;

    std::vector<std::shared_ptr<DMRG_Block1>> left_blocks_;
    std::vector<std::shared_ptr<DMRG_Block1>> right_blocks_;

    int nsites_;
    int nstate_;
    int charge_;
    int nspin_;
    int ntrunc_;
    int maxiter_;
    
    std::vector<int> active_electrons_;
    std::vector<int> active_sizes_;
    std::vector<int> region_sizes_;
    std::vector<double> energies_;
    std::vector<std::vector<double>> sweep_energies_;

    double thresh_; ///< convergence threshold for initial portion of calculation
    double perturb_; ///< magnitude of perturbation added to RDM
    double perturb_thresh_; ///< threshold at which to decrease the perturbation
    double perturb_min_; ///< minimum value of perturbation (below this, it is just set to zero)

    double down_thresh_; ///< convergence threshold for sweeping downwards with smaller M. Should probably be tighter than thresh_
    bool down_sweep_; ///< controls whether to sweep with decreasing values of ntrunc_ after the main calculation
    std::vector<int> down_sweep_truncs_; ///< descending list of values to use for ntrunc_

    std::string print_progress(const int position, const std::string left_symbol, const std::string right_symbol) const;

    virtual std::shared_ptr<DMRG_Block1> compute_first_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref) = 0;
    virtual std::shared_ptr<DMRG_Block1> grow_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> left, const int site) = 0;
    virtual std::shared_ptr<DMRG_Block1> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref, std::shared_ptr<DMRG_Block1> system, std::shared_ptr<DMRG_Block1> environment, const int site) = 0;

    // RDM
    std::shared_ptr<VecRDM<1>> rdm1_;
    std::shared_ptr<VecRDM<2>> rdm2_;
    // state averaged RDM
    std::vector<double>  weights_;
    std::shared_ptr<RDM<1>> rdm1_av_;
    std::shared_ptr<RDM<2>> rdm2_av_;

  public:
    ASD_DMRG(const std::shared_ptr<const PTree> input, std::shared_ptr<const Reference> ref);

    void sweep(const bool restart = true);
    void project_active();
    void down_sweep();

    // return functions
    int nstate() const { return nstate_; }
    int nsites() const { return nsites_; }
    const std::vector<double>& energies() const { return energies_; }
    double energies(const int i) const { return energies_.at(i); }
    std::shared_ptr<const Reference> sref() const { return sref_; }
    std::vector<int> active_sizes() const { return active_sizes_; }

    void read_restricted(std::shared_ptr<PTree> input, const int site) const;
    void update_coeff(std::shared_ptr<const Coeff> coeff) { sref_ = std::make_shared<const Reference>(*sref_, coeff); }

    // compute RDM
    void compute_rdm12();

    std::shared_ptr<VecRDM<1>> rdm1() { return rdm1_; }
    std::shared_ptr<VecRDM<2>> rdm2() { return rdm2_; }
    std::shared_ptr<RDM<1>> rdm1(const int i, const int j) { return rdm1_->at(i, j); }
    std::shared_ptr<RDM<2>> rdm2(const int i, const int j) { return rdm2_->at(i, j); }
    std::shared_ptr<RDM<1>> rdm1(const int i) { return rdm1(i, i); }
    std::shared_ptr<RDM<2>> rdm2(const int i) { return rdm2(i, i); }
    std::shared_ptr<RDM<1>> rdm1(const int i, const int j) const { return rdm1_->at(i, j); }
    std::shared_ptr<RDM<2>> rdm2(const int i, const int j) const { return rdm2_->at(i, j); }
    std::shared_ptr<RDM<1>> rdm1(const int i) const { return rdm1(i, i); }
    std::shared_ptr<RDM<2>> rdm2(const int i) const { return rdm2(i, i); }
    std::shared_ptr<RDM<1>> rdm1_av() { return rdm1_av_; }
    std::shared_ptr<RDM<2>> rdm2_av() { return rdm2_av_; }
    std::shared_ptr<const RDM<1>> rdm1_av() const { return rdm1_av_; }
    std::shared_ptr<const RDM<2>> rdm2_av() const { return rdm2_av_; }

  private:
    void rearrange_orbitals(std::shared_ptr<const Reference> iref);
    std::shared_ptr<Reference> build_reference(const int site, const std::vector<bool> meanfield) const;
    std::vector<std::shared_ptr<PTree>> prepare_growing_input(const int site) const;
    std::shared_ptr<PTree> prepare_sweeping_input(const int site) const;
};

}

#endif
