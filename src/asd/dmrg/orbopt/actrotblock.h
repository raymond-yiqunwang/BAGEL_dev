//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: actrotblock.h
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

#ifndef BAGEL_ASD_DMRG_ORBOPT_ACTROTBLOCK_H
#define BAGEL_ASD_DMRG_ORBOPT_ACTROTBLOCK_H

namespace bagel {

struct ASD_ActRotBlock {
    int iorbstart;
    int norb_i;
    int jorbstart;
    int norb_j;
    int offset;
    int size;
  
    ASD_ActRotBlock(std::vector<int> active_sizes, const int sj, int& off) : offset(off) {
      norb_j = active_sizes.at(sj);
      jorbstart = accumulate(active_sizes.begin(), active_sizes.begin()+sj, 0);
      iorbstart = jorbstart + norb_j;
      norb_i = accumulate(active_sizes.begin()+sj+1, active_sizes.end(), 0);
      size = norb_i * norb_j;
      off += size;
    }
};

}


#endif