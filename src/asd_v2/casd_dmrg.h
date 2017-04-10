//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_metal_base.h
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

#ifndef __SRC_ASD_V2_CASD_DMRG_H
#define __SRC_ASD_V2_CASD_DMRG_H

#include <src/asd_v2/asd_dmrg_base.h>

namespace bagel {

class CASD_DMRG : public ASD_DMRG_base {
  protected:
/*    std::shared_ptr<ASD_DMRG_Block1> compute_first_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref);
    std::shared_ptr<ASD_DMRG_Block1> grow_block(std::vector<std::shared_ptr<PTree>> inputs, std::shared_ptr<const Reference> ref,
                                                  std::shared_ptr<ASD_DMRG_Block1> left, const int site);
    std::shared_ptr<ASD_DMRG_Block1> decimate_block(std::shared_ptr<PTree> input, std::shared_ptr<const Reference> ref,
                                                      std::shared_ptr<ASD_DMRG_Block1> system, std::shared_ptr<ASD_DMRG_Block1> environment, const int site);
*/
  public:
    // constructor
    CASD_DMRG(const std::shared_ptr<const PTree> itree, std::shared_ptr<const Multimer> multimer) : ASD_DMRG_base(itree, multimer) { };

  private:

};

}

#endif
