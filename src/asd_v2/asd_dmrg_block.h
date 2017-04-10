//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_block.h
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

#ifndef __SRC_ASD_V2_ASD_DMRG_BLOCK_H
#define __SRC_ASD_V2_ASD_DMRG_BLOCK_H

namespace bagel {

class ASD_DMRG_Block {
  protected:
    std::set<BlockInfo> blocks_;
    std::shared_ptr<const Matrix> coeff_; // coefficients for orbitals stored in DMRG block

  public:
    ASD_DMRG_Block() {}
    ASD_DMRG_Block(std::shared_ptr<const Matrix> c) : coeff_(c) {}

//    virtual std::shared_ptr<Matrix> spin(const BlockKey k) const = 0;
//    virtual std::shared_ptr<Matrix> spin_lower(const BlockKey k) const = 0;
//    virtual std::shared_ptr<Matrix> spin_raise(const BlockKey k) const = 0;

//    virtual std::shared_ptr<const BlockOperator> compute_block_ops(std::shared_ptr<MultimerJop> jop) const = 0; TODO build MultimerJop class

};

// Store matrix representations fo all of the operators needed to apply the Hamiltonian to a tensor-product state
class ASD_DMRG_Block1 : public std::enable_shared_from_this<ASD_DMRG_Block1>, public ASD_DMRG_Block {
  protected:
    using SparseMap = std::map<std::list<GammaSQ>, std::map<std::pair<BlockKey, BlockKey>, CouplingBlock>>;

    SparseMap sparse_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> H2e_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> spin_;

  public:
    ASD_DMRG_Block1() {}

    ASD_DMRG_Block1(GammaForestASD<>)
}

} // end of namespace TODO --delete

#endif
