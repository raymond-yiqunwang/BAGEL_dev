//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/caspt2/MSCASPT2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;


MSCASPT2::MSCASPT2::MSCASPT2(const CASPT2::CASPT2& cas) {
  info_    = cas.info_;
  virt_    = cas.virt_;
  active_  = cas.active_;
  closed_  = cas.closed_;
  rvirt_   = cas.rvirt_;
  ractive_ = cas.ractive_;
  rclosed_ = cas.rclosed_;
  heff_    = cas.heff_;
  fockact_ = cas.fockact_;
  e0all_   = cas.e0all_;
  xmsmat_  = cas.xmsmat_;

  t2all_ = cas.t2all_;
  lall_  = cas.lall_;
  rall_  = cas.rall_;
  h1_ = cas.h1_;
  f1_ = cas.f1_;
  v2_ = cas.v2_;
  den1 = cas.h1_->clone();
  den2 = cas.h1_->clone();
  Den1 = cas.init_residual();

  rdm0all_ = cas.rdm0all_;
  rdm1all_ = cas.rdm1all_;
  rdm2all_ = cas.rdm2all_;
  rdm3all_ = cas.rdm3all_;
  rdm4all_ = cas.rdm4all_;

  den0ci = cas.rdm0_->clone();
  den1ci = cas.rdm1_->clone();
  den2ci = cas.rdm2_->clone();
  den3ci = cas.rdm3_->clone();
  den4ci = cas.rdm3_->clone();

  // Total tensor (that is summed up)
  den0cit = cas.rdm0_->clone();
  den1cit = cas.rdm1_->clone();
  den2cit = cas.rdm2_->clone();
  den3cit = cas.rdm3_->clone();
  den4cit = cas.rdm3_->clone();

  den0ciall = make_shared<Vec<Tensor>>();
  den1ciall = make_shared<Vec<Tensor>>();
  den2ciall = make_shared<Vec<Tensor>>();
  den3ciall = make_shared<Vec<Tensor>>();
  den4ciall = make_shared<Vec<Tensor>>();
}

void MSCASPT2::MSCASPT2::solve_dm() {
  {
    const int nstates = info_->ciwfn()->nstates();
    const int targetJ = info_->target();
    const int targetI = info_->target2();

    shared_ptr<Tensor> resultv = den1->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      const double jheffJ = (*heff_)(jst, targetJ);
      const double jheffI = (*heff_)(jst, targetI);
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          if (info_->sssr() && ist != istate)
            continue;
          const double isheffJ = (*heff_)(istate, targetJ);
          const double isheffI = (*heff_)(istate, targetI);
          const double ijvJI = (jheffJ*isheffI - jheffI*isheffJ) * 0.5;
          l2 = t2all_[istate]->at(ist); // careful

          shared_ptr<Queue> queue = make_density1q(true, ist == jst);
          while (!queue->done())
            queue->next_compute();
          resultv->ax_plus_y(ijvJI, den1);
        }
      }
    }
    vden1_ = resultv->matrix();
  }
}

void MSCASPT2::MSCASPT2::add_total(double factor) {
  den0cit->ax_plus_y(factor, den0ci);
  den1cit->ax_plus_y(factor, den1ci);
  den2cit->ax_plus_y(factor, den2ci);
  den3cit->ax_plus_y(factor, den3ci);
  den4cit->ax_plus_y(factor, den4ci);
}

void MSCASPT2::MSCASPT2::zero_total() {
  den0cit->zero();
  den1cit->zero();
  den2cit->zero();
  den3cit->zero();
  den4cit->zero();
}

void MSCASPT2::MSCASPT2::do_rdm_deriv(double factor) {
  const int nstates = info_->ciwfn()->nstates();
  Timer timer(1);
  for (int nst = 0; nst != nstates; ++nst) {
    const size_t ndet = ci_deriv_->data(nst)->size();
    const size_t nact  = info_->nact();
    const size_t norb2 = nact*nact;
    const size_t ijmax = 317520001;
    const size_t ijnum = ndet * norb2 * norb2;
    const size_t npass = (ijnum-1) / ijmax + 1;
    const size_t nsize = (ndet-1) / npass + 1;
    if (npass > 1)
      cout << "       - CI derivative contraction (state " << setw(2) << nst + 1 << ") will be done with " << npass << " passes" << endl;

    for (int ipass = 0; ipass != npass; ++ipass) {
      const size_t ioffset = ipass * nsize;
      const size_t isize = (ipass != (npass - 1)) ? nsize : ndet - ioffset;

      // If it is the first time for this state, compute fock-weighted 2RDM derivative and save it
      tie(ci_, rci_, rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3fderiv_, rdm2fderiv_)
        = SpinFreeMethod<double>::feed_rdm_deriv_3(info_, active_, fockact_, nst, ioffset, isize, /*reset=*/(ipass==0), rdm2fderiv_);
      for (int mst = 0; mst != nstates; ++mst) {
        den0cit = den0ciall->at(nst, mst);
        den1cit = den1ciall->at(nst, mst);
        den2cit = den2ciall->at(nst, mst);
        den3cit = den3ciall->at(nst, mst);
        den4cit = den4ciall->at(nst, mst);
        mpi__->barrier();

        deci = make_shared<Tensor>(vector<IndexRange>{ci_});
        deci->allocate();
        auto bdata = make_shared<VectorB>(ndet);
        shared_ptr<Queue> queue = contract_rdm_deriv(/*ciwfn=*/info_->ciwfn(), bdata, /*offset=*/ioffset, /*size=*/isize, /*reset=*/true);
        while (!queue->done())
          queue->next_compute();

        blas::ax_plus_y_n(factor, deci->vectorb()->data(), isize, ci_deriv_->data(mst)->data()+ioffset);
        blas::ax_plus_y_n(factor, bdata->data(), ndet, ci_deriv_->data(mst)->data());
      }
      if (npass > 1) {
        stringstream ss; ss << "Multipassing (" << setw(2) << ipass + 1 << " / " << npass << ")";
        timer.tick_print(ss.str());
      }
    }
  }
}


void MSCASPT2::MSCASPT2::solve_deriv() {
  Timer timer;
  const int nstates = info_->ciwfn()->nstates();
  const int target = info_->target();

  // first-order energy from the energy expression
  // TODO these can be computed more efficiently using mixed states
  {
    shared_ptr<Tensor> result = den1->clone();
    shared_ptr<Tensor> result2 = Den1->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      const double jheff = (*heff_)(jst, target);
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          if (info_->sssr() && ist != istate)
            continue;
          const double isheff = (*heff_)(istate, target);
          l2 = t2all_[istate]->at(ist); // careful
          shared_ptr<Queue> queue = make_density1q(true, ist == jst);
          while (!queue->done())
            queue->next_compute();
          result->ax_plus_y(isheff*jheff, den1);

          shared_ptr<Queue> queue2 = make_density2q(true, ist == jst);
          while (!queue2->done())
            queue2->next_compute();
          result2->ax_plus_y(isheff*jheff, Den1);
        }
      }
    }
    den1_ = result->matrix();
    Den1_ = result2;
  }
  // second-order contribution from the lambda terms
  {
    den2->zero();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          l2 = lall_[istate]->at(ist);
          t2 = t2all_[istate]->at(jst);
          shared_ptr<Queue> queue = make_densityq(false, ist == jst);
          while (!queue->done())
            queue->next_compute();
        }
      }
    }
    den2_ = den2->matrix();
  }
  // first-order contribution from the lambda terms
  {
    den1->zero();
    Den1->zero();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        if (info_->sssr() && jst != ist)
          continue;
        set_rdm(jst, ist);

        l2 = lall_[jst]->at(ist);
        shared_ptr<Queue> queue = make_density1q(false, ist == jst);
        while (!queue->done())
          queue->next_compute();

        shared_ptr<Queue> queue2 = make_density2q(false, ist == jst);
        while (!queue2->done())
          queue2->next_compute();
      }
    }
    den1_->ax_plus_y(1.0, den1->matrix());
    Den1_->ax_plus_y(1.0, Den1);
  }
  // because of the convention...
  den1_->scale(0.5);
  Den1_->scale(0.5);
  timer.tick_print("Correlated density matrix evaluation");

  // CI derivative..
  ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), nstates);

  for (int mst = 0; mst != nstates; ++mst) {
    const double mheff = (*heff_)(mst, target);
    shared_ptr<Queue> dec;

    for (int nst = 0; nst != nstates; ++nst) {
      const double nheff = (*heff_)(nst, target);
      zero_total();

      for (int lst = 0; lst != nstates; ++lst) {
        const double lheff = (*heff_)(lst, target);
        if (!info_->sssr() || nst == lst) {
          l2 = t2all_[lst]->at(nst);
          dec = make_deci3q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(lheff*mheff);
        }

        if (!info_->sssr() || mst == lst) {
          l2 = t2all_[lst]->at(mst);
          dec = make_deci4q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(lheff*nheff);
        }

        if (!info_->sssr() || (mst == lst && nst == lst)) {
          e0_ = 2.0*info_->shift();
          l2 = t2all_[lst]->at(nst);
          t2 = t2all_[lst]->at(mst);
          dec = make_deci2q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(lheff*lheff);
        }
      }

      if (!info_->sssr() || nst == mst) {
        l2 = lall_[mst]->at(nst);
        dec = make_deci3q(/*zero*/true);
        while (!dec->done())
          dec->next_compute();

        l2 = lall_[nst]->at(mst);
        dec = make_deci4q(false);
        while (!dec->done())
          dec->next_compute();

        for (int lst = 0; lst != nstates; ++lst) {
          if (info_->sssr() && (nst != lst || mst != lst))
            continue;

          e0_ = e0all_[lst] - info_->shift();
          l2 = lall_[lst]->at(nst);
          t2 = t2all_[lst]->at(mst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();
          dec = make_deci2q(false);
          while (!dec->done())
            dec->next_compute();

          l2 = t2all_[lst]->at(nst);
          t2 = lall_[lst]->at(mst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();
          dec = make_deci2q(false);
          while (!dec->done())
            dec->next_compute();
        }
        add_total(1.0);
      }

      // when active is divided into the blocks, den4cit is evaluated (activeblock)**2 times
      double den4factor = 1.0 / static_cast<double>(active_.nblock() * active_.nblock());
      den4cit->scale(den4factor);

      den0ciall->emplace(nst, mst, den0cit->copy());
      den1ciall->emplace(nst, mst, den1cit->copy());
      den2ciall->emplace(nst, mst, den2cit->copy());
      den3ciall->emplace(nst, mst, den3cit->copy());
      den4ciall->emplace(nst, mst, den4cit->copy());
    }

    stringstream ss; ss << "CI derivative evaluation   (" << setw(2) << mst+1 << " /" << setw(2) << nstates << ")";
    timer.tick_print(ss.str());
  }

  do_rdm_deriv(1.0);
  timer.tick_print("CI derivative contraction");
}

void MSCASPT2::MSCASPT2::solve_nacme() {
  Timer timer;
  const int nstates = info_->ciwfn()->nstates();
  const int targetJ = info_->target();
  const int targetI = info_->target2();

  // first-order energy from the energy expression
  {
    shared_ptr<Tensor> result = den1->clone();
    shared_ptr<Tensor> resultv = den1->clone();
    shared_ptr<Tensor> result2 = Den1->clone();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      const double jheffJ = (*heff_)(jst, targetJ);
      const double jheffI = (*heff_)(jst, targetI);
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          if (info_->sssr() && ist != istate)
            continue;
          const double isheffJ = (*heff_)(istate, targetJ);
          const double isheffI = (*heff_)(istate, targetI);
          const double ijhJI = (jheffJ*isheffI + jheffI*isheffJ) * 0.5;
          const double ijvJI = (jheffJ*isheffI - jheffI*isheffJ) * 0.5;
          l2 = t2all_[istate]->at(ist); // careful

          shared_ptr<Queue> queue = make_density1q(true, ist == jst);
          while (!queue->done())
            queue->next_compute();
          result->ax_plus_y(ijhJI, den1);
          resultv->ax_plus_y(ijvJI, den1);

          shared_ptr<Queue> queue2 = make_density2q(true, ist == jst);
          while (!queue2->done())
            queue2->next_compute();
          result2->ax_plus_y(ijhJI, Den1);
        }
      }
    }
    den1_ = result->matrix();
    vden1_ = resultv->matrix();
    Den1_ = result2;
  }

  // second-order contribution from the lambda terms
  {
    den2->zero();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        set_rdm(jst, ist);
        for (int istate = 0; istate != nstates; ++istate) { // state of T
          if (info_->sssr() && (jst != istate || ist != istate))
            continue;
          l2 = lall_[istate]->at(ist);
          t2 = t2all_[istate]->at(jst);
          shared_ptr<Queue> queue = make_densityq(false, ist == jst);
          while (!queue->done())
            queue->next_compute();
        }
      }
    }
    den2_ = den2->matrix();
  }
  // first-order contribution from the lambda terms
  {
    den1->zero();
    Den1->zero();
    for (int jst = 0; jst != nstates; ++jst) { // bra
      for (int ist = 0; ist != nstates; ++ist) { // ket
        if (info_->sssr() && jst != ist)
          continue;
        set_rdm(jst, ist);

        l2 = lall_[jst]->at(ist);
        shared_ptr<Queue> queue = make_density1q(false, ist == jst);
        while (!queue->done())
          queue->next_compute();

        shared_ptr<Queue> queue2 = make_density2q(false, ist == jst);
        while (!queue2->done())
          queue2->next_compute();
      }
    }
    den1_->ax_plus_y(1.0, den1->matrix());
    Den1_->ax_plus_y(1.0, Den1);
  }
  // because of the convention...
  den1_->scale(0.5);
  Den1_->scale(0.5);
  timer.tick_print("Correlated density matrix evaluation");

  // CI derivative..
  ci_deriv_ = make_shared<Dvec>(info_->ref()->ciwfn()->det(), nstates);

  for (int mst = 0; mst != nstates; ++mst) {
    const double mheffJ = (*heff_)(mst, targetJ);
    const double mheffI = (*heff_)(mst, targetI);
    shared_ptr<Queue> dec;

    for (int nst = 0; nst != nstates; ++nst) {
      const double nheffJ = (*heff_)(nst, targetJ);
      const double nheffI = (*heff_)(nst, targetI);
      zero_total();

      for (int lst = 0; lst != nstates; ++lst) {
        const double lheffJ = (*heff_)(lst, targetJ);
        const double lheffI = (*heff_)(lst, targetI);
        const double lnhJI  = (lheffJ * nheffI + lheffI * nheffJ) * 0.5;
        const double llhJI  = (lheffJ * lheffI + lheffI * lheffJ) * 0.5;
        const double lmhJI  = (lheffJ * mheffI + lheffI * mheffJ) * 0.5;

        if (!info_->sssr() || nst == lst) {
          l2 = t2all_[lst]->at(nst);
          dec = make_deci3q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(lmhJI);
        }

        if (!info_->sssr() || mst == lst) {
          l2 = t2all_[lst]->at(mst);
          dec = make_deci4q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(lnhJI);
        }

        if (!info_->sssr() || (mst == lst && nst == lst)) {
          e0_ = 2.0*info_->shift();
          l2 = t2all_[lst]->at(nst);
          t2 = t2all_[lst]->at(mst);
          dec = make_deci2q(/*zero=*/true);
          while (!dec->done())
            dec->next_compute();
          add_total(llhJI);
        }
      }

      if (!info_->sssr() || nst == mst) {
        l2 = lall_[mst]->at(nst);
        dec = make_deci3q(/*zero*/true);
        while (!dec->done())
          dec->next_compute();

        l2 = lall_[nst]->at(mst);
        dec = make_deci4q(false);
        while (!dec->done())
          dec->next_compute();

        for (int lst = 0; lst != nstates; ++lst) {
          if (info_->sssr() && (nst != lst || mst != lst))
            continue;

          e0_ = e0all_[lst] - info_->shift();
          l2 = lall_[lst]->at(nst);
          t2 = t2all_[lst]->at(mst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();
          dec = make_deci2q(false);
          while (!dec->done())
            dec->next_compute();

          l2 = t2all_[lst]->at(nst);
          t2 = lall_[lst]->at(mst);
          dec = make_deciq(false);
          while (!dec->done())
            dec->next_compute();
          dec = make_deci2q(false);
          while (!dec->done())
            dec->next_compute();
        }
        add_total(1.0);
      }

      // when active is divided into the blocks, den4cit is evaluated (activeblock)**2 times
      double den4factor = 1.0 / static_cast<double>(active_.nblock() * active_.nblock());
      den4cit->scale(den4factor);

      den0ciall->emplace(nst, mst, den0cit->copy());
      den1ciall->emplace(nst, mst, den1cit->copy());
      den2ciall->emplace(nst, mst, den2cit->copy());
      den3ciall->emplace(nst, mst, den3cit->copy());
      den4ciall->emplace(nst, mst, den4cit->copy());
    }

    stringstream ss; ss << "CI derivative evaluation   (" << setw(2) << mst+1 << " /" << setw(2) << nstates << ")";
    timer.tick_print(ss.str());
  }

  do_rdm_deriv(1.0);
  timer.tick_print("CI derivative contraction");
}

#endif
