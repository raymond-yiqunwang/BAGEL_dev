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

#include <src/util/muffle.h>
#include <src/asd/dmrg/rasd.h>
#include <src/asd/dmrg/orbopt/asd_dmrg_rotfile.h>
#include <src/asd/dmrg/orbopt/rotblock.h>

namespace bagel {

class ASD_DMRG_OrbOpt : public std::enable_shared_from_this<ASD_DMRG_OrbOpt> {

  protected:
    int nclosed_;
    int nact_;
    int nocc_; // sum of nclosed_ + nact_
    int nvirt_;
    int naa_; // size of active-active part
    int norb_;
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
    std::shared_ptr<const Matrix> hcore_;
    std::shared_ptr<const Geometry> geom_;
    
    // active-active rotation parameters
    int nsites_;
    std::vector<ASD_RotBlock> rotblocks_;

    std::shared_ptr<const Coeff> update_coeff(const std::shared_ptr<const Matrix> cold, std::shared_ptr<const Matrix> natorb) const;
    std::shared_ptr<const Coeff> semi_canonical_orb() const;

    // util functions
    void print_header() const;
    void print_iteration(const int iter, const std::vector<double>& energy, const double error) const;
    void common_init();

    // mask some of the output
    mutable std::shared_ptr<Muffle> muffle_;
  
  public:
    ASD_DMRG_OrbOpt(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> iref);
  
    virtual void compute() = 0;
  
    // return functions
    int nclosed() const { return nclosed_; }
    int nact() const { return nact_; }
    int nocc() const { return nocc_; }
    int nvirt() const { return nvirt_; }
    int norb() const { return norb_; }
    int nstate() const { return nstate_; }
    int max_iter() const { return max_iter_; }
    int max_micro_iter() const { return max_micro_iter_; }
    double thresh() const { return thresh_; }
    double thresh_micro() const { return thresh_micro_; }
    
    double energy(const int i) const { return energy_[i]; }
    double energy_av() const { return blas::average(energy_); }
    const std::vector<double>& energy() const { return energy_; }
    std::shared_ptr<const Coeff> coeff() const { return coeff_; }

    std::shared_ptr<Matrix> compute_active_fock(const MatView acoeff, std::shared_ptr<const RDM<1>> rdm1) const;
    std::shared_ptr<Matrix> compute_qvec(const MatView acoeff, std::shared_ptr<const RDM<2>> rdm2) const;

};

}

#endif
