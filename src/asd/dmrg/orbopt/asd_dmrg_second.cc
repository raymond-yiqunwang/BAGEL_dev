//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_second.cc
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

#include <src/asd/dmrg/orbopt/asd_dmrg_second.h>
#include <src/scf/hf/fock.h>
#include <src/util/math/aughess.h>

#define DEBUG

#ifdef DEBUG
#include <src/multi/casscf/cassecond.h>
#endif

using namespace std;
using namespace bagel;

void ASD_DMRG_Second::compute() {
  assert(nvirt_ && nact_);

#ifdef DEBUG
  cout << string(12,'=') << endl;
  auto casscf_input = input_->get_child("casscf");
  auto casscf = make_shared<CASSecond>(casscf_input, geom_, ref_); // same orbital ordering with ASD
  casscf->compute();
#endif

  muffle_->mute();
  for (int iter = 0; iter != max_iter_; ++iter) {
    
    // first obtain RDM from ASD_DMRG
    {
      if (iter) asd_dmrg_->update_multisite(coeff_);
      asd_dmrg_->compute(!iter);
      asd_dmrg_->compute_rdm12();
      // convert to natrual orbitals
      auto natorb = asd_dmrg_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
      // debugging
      if (iter < 6) {  
        auto rdm1 = make_shared<Matrix>(nact_, nact_);
        copy_n(asd_dmrg_->rdm1_av()->data(), rdm1->size(), rdm1->data());
        auto diff_rdm = make_shared<Matrix>(*(casscf->cfock0().at(iter)) - *rdm1);
        cout << "rdm1 diff rms = " << setw(16) << setprecision(12) << diff_rdm->rms() << endl;
        auto diff_coeff = make_shared<Matrix>(*(casscf->coeff0().at(iter)) - *coeff_);
        cout << "coeff diff rms = " << setw(16) << setprecision(12) << diff_coeff->rms() << endl;
      }
      energy_ = asd_dmrg_->energies();
    }
    
    shared_ptr<const Matrix> cfockao = nclosed_ ? make_shared<Fock<1>>(geom_, hcore_, nullptr, coeff_->slice(0, nclosed_), true/*store*/, true/*rhf*/) : hcore_;
    shared_ptr<const Matrix> afockao = compute_active_fock(coeff_->slice(nclosed_, nocc_), asd_dmrg_->rdm1_av());
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const Matrix> qxr = compute_qvec(coeff_->slice(nclosed_, nocc_), asd_dmrg_->rdm2_av());

    shared_ptr<const ASD_DMRG_RotFile> grad = compute_gradient(cfock, afock, qxr);

    // check gradient and break if converged
    const double gradient = grad->rms();
    print_iteration(iter, energy_, gradient);
    if (gradient < thresh_) {
      muffle_->unmute();
      cout << endl << "    * Second-Order Optimization Converged. *" << endl << endl;
      break;
    }

    // half-transformed integrals (with JJ)
    shared_ptr<const DFHalfDist> half_1j = nclosed_ ? dynamic_pointer_cast<const Fock<1>>(cfockao)->half() : nullptr;
    shared_ptr<const DFHalfDist> half = nclosed_ ? half_1j->apply_J() : nullptr;
    shared_ptr<const DFHalfDist> halfa = geom_->df()->compute_half_transform(coeff_->slice(nclosed_, nocc_));
    shared_ptr<const DFHalfDist> halfa_JJ = halfa->apply_JJ();
    
    // compute_denominator
    shared_ptr<const ASD_DMRG_RotFile> denom = compute_denom(half, half_1j, halfa, halfa_JJ, cfock, afock);

    AugHess<ASD_DMRG_RotFile> solver(max_micro_iter_, grad);

    // initial trial vector
    shared_ptr<ASD_DMRG_RotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    trot->normalize();

    // debugging
    if (iter < 6) {
      auto casscf_trot = trot->clone();
      copy_n(casscf->trot0().at(iter)->data(), casscf_trot->size(), casscf_trot->data());
      auto diff_trot = make_shared<ASD_DMRG_RotFile>(*casscf_trot - *trot);
      cout << "  * trot diff rms = " << setw(16) << setprecision(12) << diff_trot->rms() << endl;
    }

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      
      shared_ptr<const ASD_DMRG_RotFile> sigma = compute_hess_trial(trot, half, halfa_JJ, cfock, afock, qxr);
      shared_ptr<const ASD_DMRG_RotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma);
      const double err = residual->norm() / lambda;
      muffle_->unmute();
      if (!miter) cout << endl;
      cout << "         res : " << setw(8) << setprecision(2) << scientific << err
           <<       "   lamb: " << setw(8) << setprecision(2) << scientific << lambda
           <<       "   eps : " << setw(8) << setprecision(2) << scientific << epsilon
           <<       "   step: " << setw(8) << setprecision(2) << scientific << stepsize << endl;
      muffle_->mute();
      if (err < max(thresh_micro_, stepsize*thresh_microstep_))
        break;

      trot = apply_denom(residual, denom, -epsilon, 1.0/lambda);
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    } // end of micro iter

    shared_ptr<const ASD_DMRG_RotFile> sol = solver.civec();
    shared_ptr<const Matrix> a = sol->unpack();
    Matrix w(*a * *a);
    VectorB eig(a->ndim());
    w.diagonalize(eig);
    Matrix wc(w);
    Matrix ws(w);
    for (int i = 0; i != a->mdim(); ++i) {
      const double tau = sqrt(fabs(eig(i)));
      blas::scale_n(cos(tau), wc.element_ptr(0,i), wc.ndim());
      blas::scale_n(tau > 1.0e-15 ? sin(tau)/tau : 1.0, ws.element_ptr(0,i), ws.ndim());
    }
    const Matrix R = (wc ^ w) + (ws ^ w) * *a;

    // debugging
    if (iter < 6) {
      //auto diff_coeff = make_shared<Matrix>(*(casscf->coeff0().at(iter)) - *coeff_);
      //cout << "coeff diff rms = " << setw(16) << setprecision(12) << diff_coeff->rms() << endl;

      //auto diff_cfock = make_shared<Matrix>(*(casscf->cfock0().at(iter)) - *cfock);
      //cout << "cfock diff rms = " << setw(16) << setprecision(12) << diff_cfock->rms() << endl;

      auto diff_afock = make_shared<Matrix>(*(casscf->afock0().at(iter)) - *afock);
      cout << "afock diff rms = " << setw(16) << setprecision(12) << diff_afock->rms() << endl;

      auto diff_qxr = make_shared<Matrix>(*(casscf->qxr0().at(iter)) - *qxr);
      cout << "qxr diff rms = " << setw(16) << setprecision(12) << diff_qxr->rms() << endl;

      auto casscf_grad = grad->clone();
      copy_n(casscf->grad0().at(iter)->data(), casscf_grad->size(), casscf_grad->data());
      auto diff_grad = make_shared<ASD_DMRG_RotFile>(*casscf_grad - *grad);
      cout << "grad diff rms = " << setw(16) << setprecision(12) << diff_grad->rms() << endl;

      auto casscf_denom = denom->clone();
      copy_n(casscf->denom0().at(iter)->data(), casscf_denom->size(), casscf_denom->data());
      auto diff_denom = make_shared<ASD_DMRG_RotFile>(*casscf_denom - *denom);
      cout << "denom diff rms = " << setw(16) << setprecision(12) << diff_denom->rms() << endl;

      auto diff_R = make_shared<Matrix>(*(casscf->R0().at(iter)) - R);
      cout << "R diff rms = " << setw(16) << setprecision(12) << diff_R->rms() << endl;
    }

    coeff_ = make_shared<Coeff>(*coeff_ * R);

    if (iter == max_iter_-1) {
      muffle_->unmute();
      cout << endl << "    * Max iteration reached during the second-order optimization.  Convergence not reached! *   " << endl << endl;
    }
  
  } // end of macro iter
  muffle_->unmute();

  // block diagonalize coeff_ in nclosed and nvirt
  coeff_ = semi_canonical_orb();
  
  // TODO maybe one more ASD iteration
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_gradient(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto grad = make_shared<ASD_DMRG_RotFile>(nclosed_, nact_, nvirt_, naa_);
  shared_ptr<const RDM<1>> rdm1 = asd_dmrg_->rdm1_av();
  
  // closed-active section, closed runs first
  if (nclosed_) {
    double* target = grad->ptr_ca();
    for (int t = 0; t != nact_; ++t, target += nclosed_) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(0, nclosed_+t), nclosed_, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(0, nclosed_+t), nclosed_, target);
      blas::ax_plus_y_n(-2.0, qxr->element_ptr(0, t), nclosed_, target);
      for (int u = 0; u != nact_; ++u)
        blas::ax_plus_y_n(-2.0*rdm1->element(u, t), cfock->element_ptr(0, nclosed_+u), nclosed_, target);
    }
  }
  // virtual-active section, virtual runs first
  {
    double* target = grad->ptr_va();
    for (int t = 0; t != nact_; ++t, target += nvirt_) {
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nocc_, t), nvirt_, target);
      for (int u = 0; u != nact_; ++u)
        blas::ax_plus_y_n(2.0*rdm1->element(u, t), cfock->element_ptr(nocc_, nclosed_+u), nvirt_, target);
    }
  }
  // virtual-closed asection, virtual runs firsgt
  if (nclosed_){
    double* target = grad->ptr_vc();
    for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(nocc_, i), nvirt_, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(nocc_, i), nvirt_, target);
    }
  }
  // active-active part
  for (auto& block : act_rotblocks_) {
    const int istart = block.iorbstart;
    const int jstart = block.jorbstart;
    const int inorb = block.norb_i;
    const int jnorb = block.norb_j;
    const int offset = block.offset;
    double* target = grad->ptr_aa_offset(offset);

    for (int j = 0; j != jnorb; ++j, target += inorb) {
      for (int v = 0; v != nact_; ++v) {
        blas::ax_plus_y_n(2.0*rdm1->element(v, jstart+j), cfock->element_ptr(nclosed_+istart, nclosed_+v), inorb, target);
        blas::ax_plus_y_n(-2.0*cfock->element(v, jstart+j), rdm1->element_ptr(istart, v), inorb, target);
      }
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nclosed_+istart, jstart+j), inorb, target);
      blas::ax_plus_y_n(-2.0, qxr->transpose()->element_ptr(istart, nclosed_+jstart+j), inorb, target);
    }
  }

  return grad;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_denom(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> half_1j, shared_ptr<const DFHalfDist> halfa,
    shared_ptr<const DFHalfDist> halfa_JJ, shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock) const {

  auto denom = make_shared<ASD_DMRG_RotFile>(nclosed_, nact_, nvirt_, naa_);
  const MatView ccoeff = coeff_->slice(0, nclosed_);
  const MatView acoeff = coeff_->slice(nclosed_, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  Matrix rdm1(nact_, nact_);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1.size(), rdm1.data());

  // Fock related part
  const Matrix fcd = *cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * rdm1;
  const Matrix fock = *cfock + *afock;
  {
    // closed-active
    for (int i = 0; i != nact_; ++i) 
      for (int j = 0; j != nclosed_; ++j)
        denom->ele_ca(j, i) = 4.0 * fock(i+nclosed_, i+nclosed_) - 4.0 * fock(j, j) - 2.0 * fcd(i, i) + 2.0 * (*cfock)(j, j) * rdm1(i, i);
  
    // virtual-active
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_va(j, i) = 2.0 * (*cfock)(j+nocc_, j+nocc_) * rdm1(i, i) - 2.0 * fcd(i, i);
  
    // virtual-closed
    for (int i = 0; i != nclosed_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_vc(j, i) = 4.0 * fock(j+nocc_, j+nocc_) - 4.0 * fock(i, i);
  }

  const int nao = coeff_->ndim();
  
  // rdm-integral part  
  // [tt|pq] = \Gamma_{vw,tt}(vw|pq)
  shared_ptr<const DFFullDist> vaa = halfa_JJ->compute_second_transform(acoeff);
  const int nri = vaa->block(0)->asize();
  shared_ptr<const DFFullDist> vgaa = vaa->apply_2rdm(*asd_dmrg_->rdm2_av());
  {
    Matrix tmp_ao(nao, nao);
    for (int i = 0; i != nact_; ++i) {
      dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact_*i), 1, 0.0, tmp_ao.data(), 1);
      // tmp_ao.allreduce();
      Matrix tmp_virt = vcoeff % tmp_ao * vcoeff;
      blas::ax_plus_y_n(2.0, tmp_virt.diag().get(), nvirt_, denom->ptr_va()+nvirt_*i);
      if (nclosed_) {
        Matrix tmp_clo = ccoeff % tmp_ao * ccoeff;
        blas::ax_plus_y_n(2.0, tmp_clo.diag().get(), nclosed_, denom->ptr_ca()+nclosed_*i);
      }
    }
  }
  
  // [t,t] = \Gamma_{vw,xt}(vw|xt)
  shared_ptr<const DFFullDist> vaa_exc = halfa->compute_second_transform(acoeff);
  shared_ptr<const Matrix> mo2e = vaa->form_4index(vaa_exc, 1.0);
  {
    for (int i = 0; i != nact_; ++i) {
      const double e2 = -2.0 * blas::dot_product(mo2e->element_ptr(0, i*nact_), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,i));
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_va(j, i) += e2;
      for (int k = 0; k != nclosed_; ++k)
        denom->ele_ca(k, i) += e2;
    }
  }

  // mixed rdm2
  Matrix rdmk(nact_*nact_, nact_); // stores \Gamma_{k,i,j,i} + \Gamma_{k,i,i,j}
  for (int i = 0; i != nact_; ++i)
    for (int j = 0; j != nact_; ++j)
      for (int k = 0; k != nact_; ++k)
        rdmk(k+nact_*j, i) = asd_dmrg_->rdm2_av()->element(k, i, j, i) + asd_dmrg_->rdm2_av()->element(k, i, i, j);
  {
    shared_ptr<const DFFullDist> vav = halfa->compute_second_transform(vcoeff)->apply_J();
    denom->ax_plus_y_va(2.0, *(rdmk % *vav->form_4index_diagonal_part()).transpose());
    if (nclosed_) {
      shared_ptr<const DFFullDist> vac = halfa->compute_second_transform(ccoeff)->apply_J();
      shared_ptr<const Matrix> mcaa = vac->form_4index_diagonal_part()->transpose();
      denom->ax_plus_y_ca(2.0, *mcaa * rdmk);
      shared_ptr<Matrix> mcaad = mcaa->copy();
      dgemm_("N", "N", nclosed_*nact_, nact_, nact_, -1.0, mcaa->data(), nclosed_*nact_, rdm1.data(), nact_, 1.0, mcaad->data(), nclosed_*nact_);
      for (int i = 0; i != nact_; ++i)
        blas::ax_plus_y_n(12.0, mcaad->element_ptr(0, i+nact_*i), nclosed_, denom->ptr_ca()+i*nclosed_);
      
      Matrix tmp(nao, nao);
      shared_ptr<DFFullDist> vgaa = vaa->copy();
      vgaa->rotate_occ1(make_shared<Matrix>(rdm1));
      vgaa->ax_plus_y(-1.0, vaa);
      for (int i = 0; i != nact_; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact_*i), 1, 0.0, tmp.data(), 1);
        // tmp.allreduce();
        Matrix tmp0 = ccoeff % tmp * ccoeff;
        blas::ax_plus_y_n(4.0, tmp0.diag().get(), nclosed_, denom->ptr_ca()+nclosed_*i);
      }
    }
  }

  // 4-index integral part
  {
    // virtual-closed
    if (nclosed_) {
      auto vvc = half_1j->compute_second_transform(vcoeff)->form_4index_diagonal()->transpose();
      denom->ax_plus_y_vc(12.0, *vvc);

      shared_ptr<const DFFullDist> vgcc = half->compute_second_transform(ccoeff);
      const int nri = vgcc->block(0)->asize();
      Matrix tmp_ao(nao, nao);
      for (int i = 0; i != nclosed_; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgcc->block(0)->data()+nri*(i+nclosed_*i), 1, 0.0, tmp_ao.data(), 1);
        // tmp_ao.allreduce();
        Matrix tmp_virt = vcoeff % tmp_ao * vcoeff;
        blas::ax_plus_y_n(-4.0, tmp_virt.diag().get(), nvirt_, denom->ptr_vc()+nvirt_*i);
      }
    }
  }
  
  // active-active part
  shared_ptr<const Matrix> maa = vaa_exc->apply_J()->form_4index_diagonal_part();
  Matrix mgaa(nact_, nact_);
  {
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nact_; ++j)
        mgaa.element(j, i) = blas::dot_product(vaa_exc->block(0)->data()+nri*(j+nact_*i), nri, vgaa->block(0)->data()+nri*(j+nact_*i));
  }
  for (auto& block : act_rotblocks_) {
    const int istart = block.iorbstart;
    const int jstart = block.jorbstart;
    const int inorb = block.norb_i;
    const int jnorb = block.norb_j;
    const int offset = block.offset;
    for (int j = 0; j != jnorb; ++j) {
      
      // [t,t] = \Gamma_{vw,xt}(vw|xt)
      const double e2j = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact_*(jstart+j)), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,jstart+j));

      // Fock related part
      for (int i = 0; i != inorb; ++i) {
        const double e2i = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact_*(istart+i)), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,istart+i));
        denom->ele_aa_offset(i, inorb, j, offset) = 2.0 * (*cfock)(nclosed_+istart+i, nclosed_+istart+i) * rdm1(jstart+j, jstart+j)
                                                    + 2.0 * (*cfock)(nclosed_+jstart+j, nclosed_+jstart+j) * rdm1(istart+i, istart+i)
                                                    - 4.0 * (*cfock)(nclosed_+istart+i, nclosed_+jstart+j) * rdm1(istart+i, jstart+j)
                                                    - 2.0 * fcd(istart+i, istart+i) - 2.0 * fcd(jstart+j, jstart+j)
                                                    + e2j + e2i;
      }
      
      // [tt|pq] = \Gamma_{vw,tt}(vw|pq)
      Matrix tmp1_act(nact_, nact_);
      dgemv_("T", nri, nact_*nact_, 1.0, vaa_exc->block(0)->data(), nri, vgaa->block(0)->data()+nri*(jstart+j+nact_*(jstart+j)), 1, 0.0, tmp1_act.data(), 1);
      dgemv_("T", nri, nact_*nact_, 1.0, vgaa->block(0)->data(), nri, vaa_exc->block(0)->data()+nri*(jstart+j+nact_*(jstart+j)), 1, 1.0, tmp1_act.data(), 1);
      blas::ax_plus_y_n(2.0, tmp1_act.diag().get()+istart, inorb, denom->ptr_aa_offset(offset)+j*inorb);
      blas::ax_plus_y_n(-4.0, mgaa.data()+istart+nact_*(jstart+j), inorb, denom->ptr_aa_offset(offset)+j*inorb);

      // mixed rdm2
      blas::ax_plus_y_n(2.0, ((rdmk % *maa) + (*maa % rdmk)).data()+istart+nact_*(jstart+j), inorb, denom->ptr_aa_offset(offset)+j*inorb);
      {
        auto rdm_mat = make_shared<Matrix>(nact_*nact_, nact_*nact_);
        {
          btas::CRange<3> range1(nact_, nact_, nact_*nact_);
          auto tmptensor1 = make_shared<btas::Tensor3<double>>(range1, asd_dmrg_->rdm2_av()->storage());
          vector<double> buf1(nact_*nact_);
          for (int i = 0; i != range1.extent(2); ++i) {
            copy_n(&(*tmptensor1)(0,0,i), buf1.size(), buf1.data());
            blas::transpose(buf1.data(), nact_, nact_, &(*tmptensor1)(0,0,i));
          }
          copy_n(asd_dmrg_->rdm2_av()->data(), rdm_mat->size(), rdm_mat->data());
          blas::ax_plus_y_n(1.0, tmptensor1->data(), rdm_mat->size(), rdm_mat->data());
          btas::CRange<3> range2(nact_*nact_, nact_, nact_);
          auto tmptensor2 = make_shared<btas::Tensor3<double>>(range2, move(tmptensor1->storage()));
          vector<double> buf2(nact_*nact_*nact_);
          for (int i = 0; i != range2.extent(2); ++i) {
            copy_n(&(*tmptensor2)(0,0,i), buf2.size(), buf2.data());
            blas::transpose(buf2.data(), nact_*nact_, nact_, &(*tmptensor2)(0,0,i));
          } 
          btas::CRange<3> range3(nact_, nact_, nact_*nact_);
          auto tmptensor3 = make_shared<btas::Tensor3<double>>(range3, move(tmptensor2->storage()));
          vector<double> buf3(nact_*nact_);
          for (int i = 0; i != range3.extent(2); ++i) {
            copy_n(&(*tmptensor3)(0,0,i), buf3.size(), buf3.data());
            blas::transpose(buf3.data(), nact_, nact_, &(*tmptensor3)(0,0,i));
          }
          copy_n(tmptensor3->data(), rdm_mat->size(), rdm_mat->data());
        } // now we have \Gamma_{txuy} + \Gamma_{xtuy} stored in rdm_mat(tu,xy)
        auto mo2ep = mo2e->clone();
        {
          btas::CRange<3> range4(nact_*nact_, nact_, nact_);
          auto tmptensor4 = make_shared<btas::Tensor3<double>>(range4, mo2e->storage());
          vector<double> buf4(nact_*nact_*nact_);
          for (int i = 0; i != range4.extent(2); ++i) {
            copy_n(&(*tmptensor4)(0,0,i), buf4.size(), buf4.data());
            blas::transpose(buf4.data(), nact_*nact_, nact_, &(*tmptensor4)(0,0,i));
          }
          copy_n(tmptensor4->data(), mo2ep->size(), mo2ep->data());
        } // (ux|ty) stored in mo2ep(tu,xy)
        auto outmat = mo2ep->clone();
        contract(1.0, *rdm_mat, {0,2}, *mo2ep, {1,2}, 0.0, *outmat, {0,1});
        blas::ax_plus_y_n(-4.0, outmat->diag().get()+istart+nact_*(jstart+j), inorb, denom->ptr_aa_offset(offset)+j*inorb);
      }
    } // end of looping over second active index
  } // end of looping over blocks

  return denom;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::apply_denom(shared_ptr<const ASD_DMRG_RotFile> grad, shared_ptr<const ASD_DMRG_RotFile> denom, 
                                                          const double shift, const double scale) const {
  shared_ptr<ASD_DMRG_RotFile> out = grad->copy();
  for (int i = 0; i != out->size(); ++i)
    if (fabs(denom->data(i)*scale + shift) > 1.0e-12)
      out->data(i) /= denom->data(i)*scale + shift;
  return out;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_hess_trial(shared_ptr<const ASD_DMRG_RotFile> trot, shared_ptr<const DFHalfDist> half,
    shared_ptr<const DFHalfDist> halfa_JJ, shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {

  shared_ptr<ASD_DMRG_RotFile> sigma = trot->clone();

  shared_ptr<const Matrix> va = trot->va_mat();
  shared_ptr<const Matrix> ca = nclosed_ ? trot->ca_mat() : nullptr;
  shared_ptr<const Matrix> vc = nclosed_ ? trot->vc_mat() : nullptr;

  shared_ptr<const Matrix> fcaa = cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
  shared_ptr<const Matrix> faaa = afock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
  shared_ptr<const Matrix> fcva = cfock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
  shared_ptr<const Matrix> fava = afock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
  shared_ptr<const Matrix> fcvv = cfock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
  shared_ptr<const Matrix> favv = afock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
  shared_ptr<const Matrix> fccc = nclosed_ ? cfock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
  shared_ptr<const Matrix> facc = nclosed_ ? afock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
  shared_ptr<const Matrix> fcca = nclosed_ ? cfock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
  shared_ptr<const Matrix> faca = nclosed_ ? afock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
  shared_ptr<const Matrix> fcvc = nclosed_ ? cfock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;
  shared_ptr<const Matrix> favc = nclosed_ ? afock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;

  const MatView ccoeff = coeff_->slice(0, nclosed_);
  const MatView acoeff = coeff_->slice(nclosed_, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  Matrix rdm1(nact_, nact_);
  copy_n(asd_dmrg_->rdm1_av()->data(), nact_*nact_, rdm1.data());

  // lambda for computing g(D)
  auto compute_gd = [&, this] (shared_ptr<const DFHalfDist> halft, shared_ptr<const DFHalfDist> halfjj, const MatView pcoeff) {
    shared_ptr<const Matrix> pcoefft = make_shared<Matrix>(pcoeff)->transpose();
    shared_ptr<Matrix> gd = geom_->df()->compute_Jop(halft, pcoefft);
    shared_ptr<Matrix> ex0 = halfjj->form_2index(halft, 1.0);
    ex0->symmetrize();
    gd->ax_plus_y(-0.5, ex0);
    return gd;
  };

  // g(t_vc) operator and g(t_ca) operator
  if (nclosed_) {
    const Matrix tcoeff = vcoeff * *vc + acoeff * *ca->transpose();
    auto halft = geom_->df()->compute_half_transform(tcoeff);
    const Matrix gt = *compute_gd(halft, half, ccoeff);
    sigma->ax_plus_y_va(16.0, vcoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_ca(32.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_ca(-16.0, ccoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_vc(32.0, vcoeff % gt * ccoeff);
  }

  const Matrix tcoeff = nclosed_ ? (vcoeff * *va - ccoeff * *ca) : (vcoeff *  *va);
  shared_ptr<const DFHalfDist> halfta = geom_->df()->compute_half_transform(tcoeff);
  
  // g(t_va - t_ca)
  if (nclosed_) {
    shared_ptr<DFHalfDist> halftad = halfta->copy();
    halftad->rotate_occ(make_shared<Matrix>(rdm1));
    const Matrix gt = *compute_gd(halftad, halfa_JJ, acoeff);
    sigma->ax_plus_y_ca(16.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_vc(16.0, vcoeff % gt * ccoeff);
  }

  // terms with Qvec
  {
    shared_ptr<const Matrix> qaa = qxr->cut(nclosed_, nocc_);
    sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
    sigma->ax_plus_y_va(-2.0, *va * *qaa);
    if (nclosed_) {
      shared_ptr<const Matrix> qva = qxr->cut(nocc_, nocc_+nvirt_);
      shared_ptr<const Matrix> qca = qxr->cut(0, nclosed_);
      sigma->ax_plus_y_va(-2.0, *vc * *qca);
      sigma->ax_plus_y_ca(-2.0, *ca ^ *qaa);
      sigma->ax_plus_y_ca(-2.0, *ca * *qaa);
      sigma->ax_plus_y_ca(-2.0, *vc % *qva);
      sigma->ax_plus_y_vc(-2.0, *va ^ *qca);
      sigma->ax_plus_y_vc(-2.0, *qva ^ *ca);
    }
  }

  // Q' and Q'' part
  {
    shared_ptr<const DFFullDist> fullaa = halfa_JJ->compute_second_transform(acoeff);
    shared_ptr<DFFullDist> fullta = halfta->compute_second_transform(acoeff);
    shared_ptr<const DFFullDist> fulltas = fullta->swap();
    fullta->ax_plus_y(1.0, fulltas);
    shared_ptr<const DFFullDist> fullaaD = fullaa->apply_2rdm(*asd_dmrg_->rdm2_av());
    shared_ptr<const DFFullDist> fulltaD = fullta->apply_2rdm(*asd_dmrg_->rdm2_av());
    shared_ptr<const Matrix> qp = halfa_JJ->form_2index(fulltaD, 1.0);
    shared_ptr<const Matrix> qpp = halfta->form_2index(fullaaD, 1.0);

    sigma->ax_plus_y_va(4.0, vcoeff % (*qp + *qpp));
    if (nclosed_)
      sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp + *qpp));
  }

  // Fock related terms
  {
    sigma->ax_plus_y_va( 4.0, *fcvv * *va * rdm1);
    sigma->ax_plus_y_va(-2.0, *va * (rdm1 * *fcaa + *fcaa * rdm1));
    if (nclosed_) {
      sigma->ax_plus_y_vc( 8.0, (*fcvv + *favv) * *vc);
      sigma->ax_plus_y_vc(-8.0, *vc * (*fccc + *facc));
      sigma->ax_plus_y_vc( 8.0, (*fcva + *fava) ^ *ca);
      sigma->ax_plus_y_vc(-2.0, *fcva * rdm1 ^ *ca);
      sigma->ax_plus_y_vc(-4.0, *va ^ (*fcca + *faca));
      sigma->ax_plus_y_vc(-2.0, *va * rdm1 ^ *fcca);
      sigma->ax_plus_y_ca( 8.0, *vc % (*fcva + *fava));
      sigma->ax_plus_y_ca(-2.0, *vc % *fcva * rdm1);
      sigma->ax_plus_y_ca( 4.0, (*fcvc + *favc) % *va);
      sigma->ax_plus_y_ca(-4.0, *fcvc % *va * rdm1);
      sigma->ax_plus_y_ca( 8.0, *ca * (*fcaa + *faaa));
      sigma->ax_plus_y_ca(-2.0, *ca * (rdm1 * *fcaa + *fcaa * rdm1));
      sigma->ax_plus_y_ca(-8.0, (*fccc + *facc) * *ca);
      sigma->ax_plus_y_ca( 4.0, *fccc * *ca * rdm1);
      sigma->ax_plus_y_va(-4.0, *vc * (*fcca + *faca));
      sigma->ax_plus_y_va(-2.0, *vc * *fcca * rdm1);
      sigma->ax_plus_y_va( 4.0, (*fcvc + *favc) * *ca);
      sigma->ax_plus_y_va(-4.0, *fcvc * *ca * rdm1);
    }
  }

  // active-active part (P.S. use x to represent arbitrary active orbitals)
  for (auto& block : act_rotblocks_) {
    const int istart = block.iorbstart;
    const int jstart = block.jorbstart;
    const int inorb = block.norb_i;
    const int jnorb = block.norb_j;
    const int offset = block.offset;
    const int bsize = block.size;

    // cut the trial rotation vector
    auto rotblock_aa = make_shared<Matrix>(inorb, jnorb);
    copy_n(trot->ptr_aa_offset(offset), bsize, rotblock_aa->data());

    // Fock related part
    shared_ptr<const Matrix> rdmxj = rdm1.slice_copy(jstart, jstart+jnorb);
    shared_ptr<const Matrix> rdmxi = rdm1.slice_copy(istart, istart+inorb);
    { // (at, uv)
      shared_ptr<const Matrix> fcvai = fcva->slice_copy(istart, istart+inorb);
      shared_ptr<const Matrix> fcvaj = fcva->slice_copy(jstart, jstart+jnorb);
      shared_ptr<const Matrix> vai = va->slice_copy(istart, istart+inorb);
      shared_ptr<const Matrix> vaj = va->slice_copy(jstart, jstart+jnorb);
  
      sigma->ax_plus_y_va_offset(2.0, *fcva * *rdmxj ^ *rotblock_aa, istart);

      sigma->ax_plus_y_va(4.0, *fcvai * *rotblock_aa ^ (*rdmxj));

      sigma->ax_plus_y_va_offset(-2.0, *fcva * *rdmxi * *rotblock_aa, jstart);

      sigma->ax_plus_y_va(-4.0, *fcvaj ^ (*rdmxi * *rotblock_aa));

      sigma->ax_plus_y_aa_offset(2.0, *vai % *fcva * *rdmxj, offset);

      sigma->ax_plus_y_aa_offset(4.0, *fcvai % *va * *rdmxj, offset);

      sigma->ax_plus_y_aa_offset(-2.0, (*fcva * *rdmxi) % *vaj, offset);

      sigma->ax_plus_y_aa_offset(-4.0, (*va * *rdmxi) % *fcvaj, offset);
    }

    { // (ti, uv)
      shared_ptr<const Matrix> fccai = fcca->get_submatrix(0, istart, nclosed_, inorb);
      shared_ptr<const Matrix> fccaj = fcca->get_submatrix(0, jstart, nclosed_, jnorb);
      shared_ptr<const Matrix> cai = ca->get_submatrix(0, istart, nclosed_, inorb);
      shared_ptr<const Matrix> caj = ca->get_submatrix(0, jstart, nclosed_, jnorb); // TODO compare performance with slice() etc.
      
      sigma->ax_plus_y_ca_offset(4.0, *fccai * *rotblock_aa, jstart);

      sigma->ax_plus_y_ca_offset(-4.0, *fccaj ^ *rotblock_aa, istart);

      sigma->ax_plus_y_ca_offset(2.0, *fcca * *rdmxi * *rotblock_aa, istart);

      sigma->ax_plus_y_ca_offset(-2.0, *fcca * *rdmxj ^ *rotblock_aa, istart);

      sigma->ax_plus_y_ca(4.0, *fccaj ^ (*rdmxi * *rotblock_aa));

      sigma->ax_plus_y_ca(-4.0, *fccai * *rotblock_aa ^ *rdmxj);

      
      sigma->ax_plus_y_aa_offset(4.0, *fccai % *caj, offset);

      sigma->ax_plus_y_aa_offset(-4.0, *cai % *fccaj, offset);

      sigma->ax_plus_y_aa_offset(2.0, (*fcca * *rdmxi) % *caj, offset);

      sigma->ax_plus_y_aa_offset(-2.0, *cai % *fcca * *rdmxj, offset);

      sigma->ax_plus_y_aa_offset(4.0, (*ca * *rdmxi) % *fccaj, offset);

      sigma->ax_plus_y_aa_offset(-4.0, *fccai % *ca * *rdmxj, offset);
    }

    { // (tu, vw)
      shared_ptr<const Matrix> fcaai = fcaa->get_submatrix(0, istart, nact_, inorb);
      shared_ptr<const Matrix> fcaaj = fcaa->get_submatrix(0, jstart, nact_, jnorb);

      // need another loop over (vw) active pairs
      for (auto& block2 : act_rotblocks_) {
        const int istart2 = block2.iorbstart;
        const int jstart2 = block2.jorbstart;
        const int inorb2 = block2.norb_i;
        const int jnorb2 = block2.norb_j;
        const int offset2 = block2.offset;
        const int bsize2 = block2.size;
        
        // used capitalized char for block2 orbitals
        shared_ptr<const Matrix> fcaaiI = fcaa->get_submatrix(istart, istart2, inorb, inorb2);
        shared_ptr<const Matrix> fcaaiJ = fcaa->get_submatrix(istart, jstart2, inorb, jnorb2);
        shared_ptr<const Matrix> fcaaJ = fcaa->get_submatrix(0, jstart2, nact_, jnorb2);
        shared_ptr<const Matrix> fcaajI = fcaa->get_submatrix(jstart, istart2, jnorb, inorb2);
        shared_ptr<const Matrix> fcaajJ = fcaa->get_submatrix(jstart, jstart2, jnorb, jnorb2);
        auto rotblock_aa2 = rotblock_aa->clone();
        copy_n(trot->ptr_aa_offset(offset2), bsize2, rotblock_aa2->data());
        shared_ptr<const Matrix> rdmjJ = rdm1.get_submatrix(jstart, jstart2, jnorb, jnorb2);
        shared_ptr<const Matrix> rdmiJ = rdm1.get_submatrix(istart, jstart2, inorb, jnorb2);
        shared_ptr<const Matrix> rdmiI = rdm1.get_submatrix(istart, istart2, inorb, inorb2);
        shared_ptr<const Matrix> rdmxJ = rdm1.get_submatrix(0, jstart2, nact_, jnorb2);
        shared_ptr<const Matrix> rdmjI = rdm1.get_submatrix(jstart, istart2, jnorb, inorb2);

        sigma->ax_plus_y_aa_offset(4.0, *fcaaiI * *rotblock_aa2 ^ *rdmjJ, offset);
        sigma->ax_plus_y_aa_offset(-4.0, *fcaaiJ ^ (*rdmjI * *rotblock_aa2), offset);
        sigma->ax_plus_y_aa_offset(-4.0, *rdmiJ ^ (*fcaajI * *rotblock_aa2), offset);
        sigma->ax_plus_y_aa_offset(4.0, *rdmiI * *rotblock_aa2 ^ *fcaajJ, offset);

        // \delta_{uv} part
        if (jstart >= istart2) {
          auto tmpmat = make_shared<const Matrix>(*fcaai % *rdmxJ ^ *rotblock_aa2);
          auto tmpmat2 = make_shared<const Matrix>(*rdmxi % *fcaaJ ^ *rotblock_aa2);
          for (int j = 0; j != jnorb; ++j) {
            blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j+jstart-istart2), inorb, sigma->ptr_aa_offset(offset)+j*inorb);
            blas::ax_plus_y_n(2.0, tmpmat2->element_ptr(0, j+jstart-istart2), inorb, sigma->ptr_aa_offset(offset)+j*inorb);
          }
        }
        
        // \delta_{uw} part, two blocks are identical
        if (jstart == jstart2) {
          sigma->ax_plus_y_aa_offset(-2.0, *fcaai % *rdmxi * *rotblock_aa, offset);
          sigma->ax_plus_y_aa_offset(-2.0, *rdmxi % *fcaai * *rotblock_aa, offset);
        }

        // \delta_{tw} part
        if (jstart2 >= istart) {
          auto tmpmat = make_shared<const Matrix>((*fcaai * *rotblock_aa2) % *rdmxj);
          auto tmpmat2 = make_shared<const Matrix>((*rdmxj * *rotblock_aa2) % *fcaaj);
          for (int j = 0; j != jnorb; ++j) {
            blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j), jnorb2, sigma->ptr_aa_offset(offset)+jstart2-istart+j*inorb);
            blas::ax_plus_y_n(2.0, tmpmat2->element_ptr(0, j), jnorb2, sigma->ptr_aa_offset(offset)+jstart2-istart+j*inorb);
          }
        }

        // \delta_{tv} part
        {
          auto tmpmat = make_shared<const Matrix>(*rotblock_aa2 ^ *fcaaJ * *rdmxj);
          auto tmpmat2 = make_shared<const Matrix>(*rotblock_aa2 ^ *rdmxJ * *fcaaj);
          if (istart < istart2) {
            for (int j = 0; j != jnorb; ++j) {
              blas::ax_plus_y_n(-2.0, tmpmat->element_ptr(0, j), inorb2, sigma->ptr_aa_offset(offset)+j*inorb+istart2-istart);
              blas::ax_plus_y_n(-2.0, tmpmat2->element_ptr(0, j), inorb2, sigma->ptr_aa_offset(offset)+j*inorb+istart2-istart);
            }
          } else {
            auto tmp = tmpmat->get_submatrix(istart-istart2, 0, inorb, jnorb);
            auto tmp2 = tmpmat2->get_submatrix(istart-istart2, 0, inorb, jnorb);
            sigma->ax_plus_y_aa_offset(-2.0, *tmp, offset);
            sigma->ax_plus_y_aa_offset(-2.0, *tmp2, offset);
          }
        }
      }
    }
    
  } // end of looping over blocks
  
  sigma->scale(0.5);

  return sigma;
}
