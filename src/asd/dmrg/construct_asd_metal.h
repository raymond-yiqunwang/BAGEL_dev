//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_asd_metal.h 
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

#ifndef __SRC_ASD_DMRG_CONSTRUCT_ASD_METAL_H
#define __SRC_ASD_DMRG_CONSTRUCT_ASD_METAL_H

#include <src/asd/dmrg/rasd.h>

namespace bagel {
  extern std::shared_ptr<RASD> construct_ASD_METAL(std::shared_ptr<const PTree> itree, std::shared_ptr<const Reference> ref);
}

#endif
