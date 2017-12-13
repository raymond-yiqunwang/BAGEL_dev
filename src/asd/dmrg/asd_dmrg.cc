//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_base.cc
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

#include <iostream>

#include <src/asd/dmrg/asd_dmrg.h>

using namespace std;
using namespace bagel;

ASD_DMRG::ASD_DMRG(shared_ptr<const PTree> input, shared_ptr<const MultiSite> multisite) : input_(input), multisite_(multisite) {
  nsites_ = multisite->nsites();
  nstate_ = input_->get<int>("nstate", 1);
  if (metal_) nactorb_ = multisite_->ref()->nact();
  thresh_ = input_->get<double>("thresh", 1.0e-6);
  maxiter_ = input_->get<int>("maxiter", 50);

  perturb_ = input_->get<double>("perturb", 0.001);
  perturb_thresh_ = input_->get<double>("perturb_thresh", 0.0001);
  perturb_min_ = input_->get<double>("perturb_min", 1.0e-5);

  down_thresh_ = input_->get<double>("down_thresh", 1.0e-8);
  auto down = input_->get_child_optional("down_sweep_truncs");
  down_sweep_ = static_cast<bool>(down);
  if (down_sweep_)
    down_sweep_truncs_ = input_->get_vector<int>("down_sweep_truncs");

  auto winput = input_->get_child_optional("weights");
  if (winput)
    weights_ = input_->get_vector<double>("weights", nstate_);
  else
    weights_.resize(nstate_, 1.0/static_cast<double>(nstate_));

  energies_.resize(nstate_);
  sweep_energies_.resize(nstate_);

  // initialize VecRDM
  rdm1_ = make_shared<VecRDM<1>>();
  rdm2_ = make_shared<VecRDM<2>>();
}


string ASD_DMRG::print_progress(const int position, const string left_symbol, const string right_symbol) const {
  stringstream out;
  for (int i = 0; i < position; ++i) out << left_symbol << " ";
  out << "** ";
  for (int i = position+1; i < nsites_; ++i) out << right_symbol << " ";

  return out.str();
}


vector<shared_ptr<PTree>> ASD_DMRG::prepare_growing_input(const int site) const {
  vector<shared_ptr<PTree>> out;

  shared_ptr<PTree> base_input = input_->get_child_optional("ras");
  if (!base_input) base_input = make_shared<PTree>();

  base_input->erase("charge");
  base_input->erase("nspin");
  base_input->erase("nstate");

  auto space = input_->get_child_optional("spaces");
  if (!space) throw runtime_error("spaces must be specified");
  else if ( !(space->size()==1 || space->size()==nsites_) )
    throw runtime_error("Must specify either one \"space\" object or one per site");

  auto space_iter = space->begin();
  if (space->size() == nsites_)
    advance(space_iter, site);

  // Quirk of PTree class requires this intermediate step
  auto space_input = make_shared<PTree>(**space_iter);

  for (auto& siter : *space_input) {
    auto tmp = make_shared<PTree>(*base_input);
    tmp->put("charge", siter->get<string>("charge"));
    tmp->put("nspin", siter->get<string>("nspin"));
    tmp->put("nstate", siter->get<string>("nstate"));
    out.push_back(tmp);
  }

  return out;
}


shared_ptr<PTree> ASD_DMRG::prepare_sweeping_input(const int site) const {
  shared_ptr<PTree> out = input_->get_child_optional("ras");
  if (!out) out = make_shared<PTree>();

  out->erase("charge"); out->put("charge", multisite_->charge());
  out->erase("nspin");  out->put("nspin", multisite_->nspin());
  out->erase("nstate"); out->put("nstate", input_->get<string>("nstate", "1"));

  return out;
}


shared_ptr<const Reference> ASD_DMRG::conv_to_ref() const {
  // TODO this function is to be modified
  auto mref = multisite_->ref();
  return make_shared<Reference>(mref->geom(), mref->coeff(), mref->nclosed(), mref->nact(), mref->nvirt(), energies_, make_shared<VecRDM<1>>(), make_shared<VecRDM<2>>());
}


void ASD_DMRG::rotate_rdms(shared_ptr<const Matrix> trans) {
  for (auto& i : *rdm1_)
    i.second->transform(trans);
  for (auto& i : *rdm2_)
    i.second->transform(trans);

  // Only when #state > 1
  if (rdm1_->size() > 1) rdm1_av_->transform(trans);
  if (rdm2_->size() > 1) rdm2_av_->transform(trans);
  assert(rdm1_->size() > 1 || rdm1_->at(0) == rdm1_av_);
  assert(rdm2_->size() > 1 || rdm2_->at(0) == rdm2_av_);
}


