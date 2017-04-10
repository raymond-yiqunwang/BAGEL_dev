//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: construct_asd_matal.cc
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

#include <src/asd_v2/casd_dmrg.h>

using namespace std;
using namespace bagel;

namespace bagel {
  
// Since we project the total active space to fragments, FCI should be a better starting point than RASCI TODO --delete
shared_ptr<ASD_DMRG_base> construct_ASD_Metal(shared_ptr<const PTree> itree, shared_ptr<const Reference> ref) {

  // construct multimer first
  shared_ptr<const PTree> multimer_info = itree->get_child("multimer");
  auto multimer = make_shared<Multimer>(multimer_info, ref);
  multimer->precompute(multimer_info);

  shared_ptr<ASD_DMRG_base> out;
  string method = itree->get<string>("method", "cas");

  if (method == "cas" || method == "fci") {
    out = make_shared<CASD_DMRG>(itree, multimer);
  } else if (method == "ras") {
    // TODO RAS version of ASD_DMRG_METAL
  } else {
    throw logic_error("Unrecognized method for ASD_METAL");
  }

  return out;
}
  
}


