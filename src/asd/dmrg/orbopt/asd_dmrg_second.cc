//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_second.cc
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

  for (int iter = 0; iter != max_iter_; ++iter) {
    
    // first obtain RDM from ASD_DMRG
    {
      if (iter) asd_dmrg_->update_multisite(coeff_);
      asd_dmrg_->compute();
      asd_dmrg_->compute_rdm12();
      // convert to natrual orbitals
      auto natorb = asd_dmrg_->natorb_convert();
      coeff_ = update_coeff(coeff_, natorb.first);
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
    if (gradient < thresh_) {
      // muffle->unmute();
      cout << endl << "    * Second-Order Optimization Converged. *" << endl << endl;
      break;
    }

    // half-transformed integrals (with JJ)
    shared_ptr<const DFHalfDist> half_1j = nclosed_ ? dynamic_pointer_cast<const Fock<1>>(cfockao)->half() : nullptr;
    shared_ptr<const DFHalfDist> half = nclosed_ ? half_1j->apply_J() : nullptr;
    // TODO is this the right way to do it?
    shared_ptr<const DFHalfDist> halfa = geom_->df()->compute_half_transform(coeff_->slice(nclosed_, nocc_));
    shared_ptr<const DFHalfDist> halfa_JJ = halfa->apply_JJ();
    
    // compute_denominator
    shared_ptr<const ASD_DMRG_RotFile> denom = compute_denom(half, half_1j, halfa, halfa_JJ, cfock, afock);

    AugHess<ASD_DMRG_RotFile> solver(max_micro_iter_, grad);

    // initial trial vector
    shared_ptr<ASD_DMRG_RotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    trot->normalize();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      
      shared_ptr<const ASD_DMRG_RotFile> sigma = compute_hess_trial(trot, half, halfa_JJ, cfock, afock, qxr);
      shared_ptr<const ASD_DMRG_RotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma);
      const double err = residual->norm() / lambda;
      if (!miter) cout << endl;
      cout << "         res : " << setw(8) << setprecision(2) << scientific << err
           <<       "   lamb: " << setw(8) << setprecision(2) << scientific << lambda
           <<       "   eps : " << setw(8) << setprecision(2) << scientific << epsilon
           <<       "   step: " << setw(8) << setprecision(2) << scientific << stepsize << endl;
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

    coeff_ = make_shared<Coeff>(*coeff_ * R);

    if (iter == max_iter_-1) {
      cout << endl << "    * Max iteration reached during the second-order optimization.  Convergence not reached! *   " << endl << endl;
    }

    // block diagonalize coeff_ in nclosed and nvirt
    coeff_ = semi_canonical_orb();

    // TODO maybe one more ASD iteration

  } // end of macro iter
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_gradient(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto grad = make_shared<ASD_DMRG_RotFile>(nclosed_, nact_, nvirt_);
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

  return grad;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_denom(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> half_1j, shared_ptr<const DFHalfDist> halfa,
    shared_ptr<const DFHalfDist> halfa_JJ, shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock) const {

  auto denom = make_shared<ASD_DMRG_RotFile>(nclosed_, nact_, nvirt_);
  const MatView ccoeff = coeff_->slice(0, nclosed_);
  const MatView acoeff = coeff_->slice(nclosed_, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  Matrix rdm1(nact_, nact_);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1.size(), rdm1.data());

  // Fock related part
  {
    // closed-active
    const Matrix fcd = *cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * rdm1;
    const Matrix fock = *cfock + *afock;
    for (int i = 0; i != nact_; ++i) 
      for (int j = 0; j != nclosed_; ++j)
      // TODO += should also work, check later
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
  {
    // [tt|pq] = \Gamma_{vw,tt}(vw|pq)
    shared_ptr<const DFFullDist> vaa = halfa_JJ->compute_second_transform(acoeff);
    const int nri = vaa->block(0)->asize();
    shared_ptr<const DFFullDist> vgaa = vaa->apply_2rdm(*asd_dmrg_->rdm2_av());
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
    
    // [t,t] = \Gamma_{vw,xt}(vw|xt)
    shared_ptr<const DFFullDist> vaa_exc = halfa->compute_second_transform(acoeff);
    shared_ptr<const Matrix> mo2e = vaa->form_4index(vaa_exc, 1.0);
    for (int i = 0; i != nact_; ++i) {
      const double e2 = -2.0 * blas::dot_product(mo2e->element_ptr(0, i*nact_), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,i));
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_va(j, i) += e2;
      for (int k = 0; k != nclosed_; ++k)
        denom->ele_ca(k, i) += e2;
    }

    // mixed rdm2
    Matrix rdmk(nact_*nact_, nact_); // stores \Gamma_{k,i,j,i} + \Gamma_{k,i,i,j}
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nact_; ++j)
        for (int k = 0; k != nact_; ++k)
          rdmk(k+nact_*j, i) = asd_dmrg_->rdm2_av()->element(k, i, j, i) + asd_dmrg_->rdm2_av()->element(k, i, i, j);
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
  sigma->scale(0.5);

  return sigma;
}
