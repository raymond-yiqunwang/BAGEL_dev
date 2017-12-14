//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_orbopt.cc
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

#include <src/asd/dmrg/orbopt/asd_dmrg_orbopt.h>
#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

ASD_DMRG_OrbOpt::ASD_DMRG_OrbOpt(shared_ptr<const PTree> idata, shared_ptr<const Reference> iref) : input_(idata) {
  
  print_header();
  
  // collect ASD-DMRG info
  auto asd_dmrg_info = input_->get_child_optional("asd_dmrg_info");
  if (!asd_dmrg_info) throw runtime_error("ASD-DMRG info has to be provided for orbital optimization");

  // collect MultiSite info
  shared_ptr<const MultiSite> multisite;
  {
    auto multisite_info = input_->get_child_optional("multisite");
    if (!multisite_info) throw runtime_error("MultiSite info has to be provided for ASD-DMRG orbital optimization");
    const int nsites = multisite_info->get<int>("nsites");
    auto ms = make_shared<MultiSite>(multisite_info, iref, nsites);
    ms->compute();
    multisite = ms;
  }

  auto mref = multisite->sref();
  coeff_ = mref->coeff();

#ifdef AAROT
  cout << " *** Active-active rotation turned on!" << endl;
  // initialize active-active rotation parameters
  int offset = 0;
  for (int sj = 0; sj != multisite->nsites()-1; ++sj)
    act_rotblocks_.emplace_back(multisite->active_sizes(), sj, offset);
  naa_ = offset;
#else
  naa_ = 0;
#endif

  nstate_ = input_->get<int>("opt_nstate", 1);
  max_iter_ = input_->get<int>("opt_maxiter", 50);
  max_micro_iter_ = input_->get<int>("opt_max_micro_iter", 100);
  thresh_ = input_->get<double>("opt_thresh", 1.0e-8); // thresh for macro iteration
  thresh_micro_ = input_->get<double>("opt_thresh_micro", 5.0e-6); // thresh for micro iteration

  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << mref->nclosed() << endl;
  cout << "    * nact     : " << setw(6) << mref->nact() << endl;
  cout << "    * nvirt    : " << setw(6) << mref->nvirt() << endl << endl;
  assert(mref->nact() && mref->nvirt());

  muffle_ = make_shared<Muffle>("asd_dmrg_orbopt.log");
  muffle_->unmute();

  // DMRG with RHF orbitals
  asd_dmrg_ = make_shared<RASD>(asd_dmrg_info, multisite);
  
  cout << "  ===== Orbital Optimization Iteration =====" << endl << endl;
}


void ASD_DMRG_OrbOpt::print_header() const {
  cout << "  --------------------------------------------------" << endl;
  cout << "     ASD-DMRG Second Order Orbital Optimization     " << endl;
  cout << "  --------------------------------------------------" << endl << endl;
}


void ASD_DMRG_OrbOpt::print_iteration(const int iter, const vector<double>& energy, const double error) const {
  if (energy.size() != 1 && iter) cout << endl;

  int i = 0;
  for (auto& e : energy) {
    cout << "  " << setw(5) << iter << setw(3) << i << setw(19) << fixed << setprecision(12) << e << "   "
                 << setw(10) << scientific << setprecision(2) << (i==0 ? error : 0.0) << endl;
    ++i;
  }
}


shared_ptr<Matrix> ASD_DMRG_OrbOpt::compute_active_fock(const MatView acoeff, shared_ptr<const RDM<1>> rdm1) const {
  auto mref = asd_dmrg_->multisite()->sref();
  const int nact = asd_dmrg_->multisite()->sref()->nact();
  Matrix dkl(nact, nact);
  copy_n(rdm1->data(), dkl.size(), dkl.data());
  dkl.sqrt();
  dkl.scale(1.0/sqrt(2.0));
  return make_shared<Fock<1>>(mref->geom(), mref->hcore()->clone(), nullptr, acoeff*dkl, false, true);
}


shared_ptr<Matrix> ASD_DMRG_OrbOpt::compute_qvec(const MatView acoeff, shared_ptr<const RDM<2>> rdm2) const {
  
  auto mref = asd_dmrg_->multisite()->sref();
  const int nocc = mref->nclosed() + mref->nact();
  auto half = mref->geom()->df()->compute_half_transform(acoeff);

  // TODO MPI modification
  shared_ptr<const DFFullDist> full = half->apply_JJ()->compute_second_transform(coeff_->slice(mref->nclosed(), nocc));

  // [D|tu] = (D|xy) Gamma_{xy,tu}
  shared_ptr<const DFFullDist> prdm = full->apply_2rdm(*rdm2);

  // (r,u) = (rt|D) [D|tu]
  shared_ptr<const Matrix> tmp = half->form_2index(prdm, 1.0);

  return make_shared<Matrix>(*coeff_ % *tmp);
}


shared_ptr<const Coeff> ASD_DMRG_OrbOpt::semi_canonical_orb() const {
  auto mref = asd_dmrg_->multisite()->sref();
  const int nclosed = mref->nclosed();
  const int nact = mref->nact();
  const int nocc = nclosed + nact;
  const int nvirt = mref->nvirt();
  const int norb = coeff_->mdim();
  
  auto rdm1_mat = make_shared<Matrix>(nact, nact);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1_mat->size(), rdm1_mat->data());
  rdm1_mat->sqrt();
  rdm1_mat->scale(1.0/sqrt(2.0));

  const MatView ccoeff = coeff_->slice(0, nclosed);
  const MatView acoeff = coeff_->slice(nclosed, nocc);
  const MatView vcoeff = coeff_->slice(nocc, norb);

  VectorB eig(coeff_->mdim());
  auto core_fock = nclosed ? make_shared<Fock<1>>(mref->geom(), mref->hcore(), nullptr, coeff_->slice(0, nclosed), false/*store*/, true/*rhf*/)
                            : make_shared<Matrix>(*mref->hcore());
  Fock<1> fock(mref->geom(), core_fock, nullptr, acoeff * *rdm1_mat, false, true);

  Matrix trans(norb, norb);
  trans.unit();
  if (nclosed) {
    Matrix ofock = ccoeff % fock * ccoeff;
    ofock.diagonalize(eig);
    trans.copy_block(0, 0, nclosed, nclosed, ofock);
  }
  Matrix vfock = vcoeff % fock * vcoeff;
  vfock.diagonalize(eig);
  trans.copy_block(nocc, nocc, nvirt, nvirt, vfock);
  
  return make_shared<Coeff>(*coeff_ * trans);
}


