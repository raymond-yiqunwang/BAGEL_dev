//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_base.cc
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

#include <src/asd_v2/asd_dmrg_base.h>

using namespace std;
using namespace bagel;

ASD_DMRG_base::ASD_DMRG_base(const std::shared_ptr<const PTree> input, std::shared_ptr<const Multimer> multimer) : input_(input), multimer_(multimer) {
  nsites_ = multimer->nsites();
  nstate_ = input_->get<int>("nstate", 1);
  ntrunc_ = input_->get<int>("ntrunc");
  thresh_ = input_->get<int>("thresh", 1.0e-6);
  maxiter_ = input_->get<int>("maxiter", 50);

  // TODO add perturbation
  // perturb = input_->get<double>("perturb", 0.001);
  // perturb_thresh_ = input_->get<double>("perturb_thresh", 0.0001);
  // perturb_min = input_->get<double>("perturb_min", 1.0e-5);

  auto weighinput = input_->get_child_optional("weights");
  if (weighinput)
    weights_ = input_->get_vector<double>("weights", nstate_);
  else
    weights_.resize(nstate_, 1.0/static_cast<double>(nstate_));

  energies_.resize(nstate_);
  sweep_energies_.resize(nstate_);
}
