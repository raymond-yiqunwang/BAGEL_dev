//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_asd_metal.cc
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

#include <src/asd/multisite/multisite.h>
#include <src/asd/dmrg/rasd.h>

using namespace std;
using namespace bagel;

namespace bagel {

shared_ptr<RASD> construct_ASD_METAL(shared_ptr<const PTree> itree, shared_ptr<const Reference> ref) {
  
  shared_ptr<const PTree> multisite_info = itree->get_child("multisite");
  auto multisite = make_shared<MultiSite>(multisite_info, ref);

  auto out = make_shared<RASD>(itree, multisite);

  return out;

}

}
