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

  //muffle_->mute();
  for (int iter = 0; iter != max_iter_; ++iter) {
    
    // first obtain RDM from ASD_DMRG
    {
      muffle_->mute();
      if (iter) asd_dmrg_->update_multisite(coeff_);
      asd_dmrg_->compute(!iter);
      asd_dmrg_->compute_rdm12();
      trans_natorb();
      energy_ = asd_dmrg_->energies();
      muffle_->unmute();
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

#ifdef DEBUG
  // build Hessian matrix
  const int rotsize = denom->size();
  cout << "rotsize = " << rotsize << endl;
  Matrix rdm1(nact_, nact_);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1.size(), rdm1.data());
  auto rdm2 = *asd_dmrg_->rdm2_av();
  Matrix mo2e(norb_, norb_);
  {
    const MatView coeff(*coeff_);
    shared_ptr<const DFHalfDist> halfx = geom_->df()->compute_half_transform(coeff);
    shared_ptr<const DFFullDist> fullx_1j = halfx->compute_second_transform(coeff)->apply_J();
    mo2e = *fullx_1j->form_4index(fullx_1j, 1.0);
  }
  const int va_offset = nclosed_ * nact_;
  const int vc_offset = va_offset + nvirt_ * nact_;
  cout << "va_offset = " << va_offset << endl;
  cout << "vc_offset = " << vc_offset << endl;
#ifdef AAROT
  const int aa_offset = vc_offset + nvirt_ * nclosed_;
  cout << "aa_offset = " << aa_offset << endl;
#endif
  const Matrix fock = *cfock + *afock;
  shared_ptr<const Matrix> cfk_xt = cfock->get_submatrix(0, nclosed_, norb_, nact_);
  const Matrix cfkd = *cfk_xt * rdm1;
  
  VectorB grad_check(rotsize);
  {
    // closed-active
    for (int t = 0; t != nact_; ++t) {
      for (int i = 0; i != nclosed_; ++i) {
        double value = 4.0 * fock(i, nclosed_+t) - 2.0 * cfkd(i, t);
        for (int u = 0; u != nact_; ++u) {
          for (int v = 0; v != nact_; ++v) {
            for (int w = 0; w != nact_; ++w) {
              value -= 2.0 * rdm2(t, u, v, w) * mo2e(nclosed_+u+i*norb_, nclosed_+v+(nclosed_+w)*norb_);
            }
          }
        }
        grad_check(i+t*nclosed_) = value;
      }
    }
    // virtual-active
    for (int t = 0; t != nact_; ++t) {
      for (int a = 0; a != nvirt_; ++a) {
        double value = 2.0 * cfkd(nocc_+a, t);
        for (int u = 0; u != nact_; ++u) {
          for (int v = 0; v != nact_; ++v) {
            for (int w = 0; w != nact_; ++w) {
              value += 2.0 * rdm2(t, u, v, w) * mo2e(nocc_+a+(nclosed_+u)*norb_, nclosed_+v+(nclosed_+w)*norb_);
            }
          }
        }
        grad_check(va_offset+a+t*nvirt_) = value;
      }
    }
    // virtual-closed
    for (int i = 0; i != nclosed_; ++i) {
      for (int a = 0; a != nvirt_; ++a) {
        grad_check(vc_offset+a+i*nvirt_) = 4.0 * fock(nocc_+a, i);
      }
    }
#ifdef AAROT
    // active-active
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norb_i;
      const int jnorb = block.norb_j;
      const int offset = block.offset;

      for (int t = 0; t != inorb; ++t) {
        for (int u = 0; u != jnorb; ++u) {
          double value = 2.0 * (cfkd(nclosed_+istart+t, jstart+u) - cfkd(nclosed_+jstart+u, istart+t));
          for (int v = 0; v != nact_; ++v) {
            for (int w = 0; w != nact_; ++w) {
              for (int x = 0; x != nact_; ++x) {
                value += 2.0 * rdm2(jstart+u, v, w, x) * mo2e(nclosed_+istart+t+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_)
                       - 2.0 * rdm2(istart+t, v, w, x) * mo2e(nclosed_+jstart+u+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_);
              }
            }
          }
          grad_check(aa_offset+offset+t+u*inorb) = value;
        }
      }
    }
#endif
  }

  auto hessian = make_shared<Matrix>(rotsize, rotsize);
  {
    // Fock part
    {
      // (at, bu)
      for (int u = 0; u != nact_; ++u) {
        for (int b = 0; b != nvirt_; ++b) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 2.0 * rdm1(t, u) * cfock->element(nocc_+a, nocc_+b);
              if (a == b) value -= cfkd(nclosed_+u, t) + cfkd(nclosed_+t, u);
              hessian->element(va_offset+a+t*nvirt_, va_offset+b+u*nvirt_) = value;
            }
          }
        }
      }
      // (at, ui)
      for (int u = 0; u != nact_; ++u) {
        for (int i = 0; i != nclosed_; ++i) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = -2.0 * rdm1(t, u) * cfock->element(nocc_+a, i);
              if (t == u) value += 2.0 * fock(nocc_+a, i);
              hessian->element(va_offset+a+t*nvirt_, i+u*nclosed_) = value;
              hessian->element(i+u*nclosed_, va_offset+a+t*nvirt_) = value;
            }
          }
        }
      }
      // (bi, at)
      for (int t = 0; t != nact_; ++t) {
        for (int a = 0; a != nvirt_; ++a) {
          for (int i = 0; i != nclosed_; ++i) {
            for (int b = 0; b != nvirt_; ++b) {
              if (a == b) {
                double value = -2.0 * fock(i, nclosed_+t) - cfkd(i, t);
                hessian->element(vc_offset+b+i*nvirt_, va_offset+a+t*nvirt_) = value;
                hessian->element(va_offset+a+t*nvirt_, vc_offset+b+i*nvirt_) = value;
              }
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact_; ++u) {
        for (int j = 0; j != nclosed_; ++j) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              double value = 2.0 * rdm1(t, u) * cfock->element(i, j);
              if (i == j) {
                value += 4.0 * fock(nclosed_+t, nclosed_+u) - cfkd(nclosed_+t, u) - cfkd(nclosed_+u, t);
              }
              if (t == u) {
                value -= 4.0 * fock(i, j);
              }
              hessian->element(i+t*nclosed_, j+u*nclosed_) = value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed_; ++j) {
        for (int a = 0; a != nvirt_; ++a) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              if (i == j) {
                double value = 4.0 * fock(nocc_+a, nclosed_+t) - cfkd(nocc_+a, t);
                hessian->element(i+t*nclosed_, vc_offset+a+j*nvirt_) = value;
                hessian->element(vc_offset+a+j*nvirt_, i+t*nclosed_) = value;
              }
            }
          }
        }
      }
      // (ai, bj)
      for (int j = 0; j != nclosed_; ++j) {
        for (int b = 0; b != nvirt_; ++b) {
          for (int i = 0; i != nclosed_; ++i) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 0;
              if (a == b) {
                value -= 4.0 * fock(i, j);
              }
              if (i == j) {
                value += 4.0 * fock(nocc_+a, nocc_+b);
              }
              hessian->element(vc_offset+a+i*nvirt_, vc_offset+b+j*nvirt_) = value;
            }
          }
        }
      }
#ifdef AAROT
      for (auto& block : act_rotblocks_) {
        const int istart = block.iorbstart;
        const int jstart = block.jorbstart;
        const int inorb = block.norb_i;
        const int jnorb = block.norb_j;
        const int offset = block.offset;

        // (ti, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact_; ++t) {
              for (int i = 0; i != nclosed_; ++i) {
                double value = 0.0;
                if (t == jstart+v) {
                  value += 2.0 * cfock->element(i, nclosed_+istart+u) + cfkd(i, istart+u);
                }
                if (t == istart+u) {
                  value -= 2.0 * cfock->element(i, nclosed_+jstart+v) + cfkd(i, jstart+v);
                }
                value += 2.0 * rdm1(t, istart+u) * cfock->element(i, nclosed_+jstart+v) - 2.0 * rdm1(t, jstart+v) * cfock->element(i, nclosed_+istart+u);
                hessian->element(i+t*nclosed_, aa_offset+offset+u+v*inorb) = value;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed_) = value;
              }
            }
          }
        }
        // (at, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact_; ++t) {
              for (int a = 0; a != nvirt_; ++a) {
                double value = 2.0 * rdm1(t, jstart+v) * cfock->element(nocc_+a, nclosed_+istart+u)
                             - 2.0 * rdm1(t, istart+u) * cfock->element(nocc_+a, nclosed_+jstart+v);
                if (t == istart+u) {
                  value += cfkd(nocc_+a, jstart+v);
                }
                if (t == jstart+v) {
                  value -= cfkd(nocc_+a, istart+u);
                }
                hessian->element(va_offset+a+t*nvirt_, aa_offset+offset+u+v*inorb) = value;
                hessian->element(aa_offset+offset+u+v*inorb, va_offset+a+t*nvirt_) = value;
              }
            }
          }
        }
        // (tu, vw)
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norb_i;
          const int jnorb2 = block2.norb_j;
          const int offset2 = block2.offset;

          for (int w = 0; w != jnorb2; ++w) {
            for (int v = 0; v != inorb2; ++v) {
              for (int u = 0; u != jnorb; ++u) {
                for (int t = 0; t != inorb; ++t) {
                  double value = 2.0 * rdm1(jstart2+w, jstart+u) * cfock->element(nclosed_+istart+t, nclosed_+istart2+v)
                               - 2.0 * rdm1(istart2+v, jstart+u) * cfock->element(nclosed_+istart+t, nclosed_+jstart2+w)
                               - 2.0 * rdm1(jstart2+w, istart+t) * cfock->element(nclosed_+jstart+u, nclosed_+istart2+v)
                               + 2.0 * rdm1(istart2+v, istart+t) * cfock->element(nclosed_+jstart+u, nclosed_+jstart2+w);
                  if (jstart+u == istart2+v) {
                    value += cfkd(nclosed_+istart+t, jstart2+w) + cfkd(nclosed_+jstart2+w, istart+t);
                  }
                  if (jstart+u == jstart2+w) {
                    value -= cfkd(nclosed_+istart+t, istart2+v) + cfkd(nclosed_+istart2+v, istart+t);
                  }
                  if (istart+t == istart2+v) {
                    value -= cfkd(nclosed_+jstart+u, jstart2+w) + cfkd(nclosed_+jstart2+w, jstart+u);
                  }
                  if (istart+t == jstart2+w) {
                    value += cfkd(nclosed_+jstart+u, istart2+v) + cfkd(nclosed_+istart2+v, jstart+u);
                  }
                  hessian->element(aa_offset+offset+t+u*inorb, aa_offset+offset2+v+w*inorb2) = value;
                }
              }
            }
          }
        }
      } // end end looping over actrotblocks
#endif
    } // end of Fock part of Hessian matrix

    // 2-electron integral part
    {
      // (at, ui)
      for (int u = 0; u != nact_; ++u) {
        for (int i = 0; i != nclosed_; ++i) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 0;
              for (int v = 0; v != nact_; ++v) {
                assert(mo2e(nocc_+a+(nclosed_+v)*norb_, nclosed_+u+i*norb_) - mo2e(nclosed_+v+(nocc_+a)*norb_, nclosed_+u+i*norb_) < 1.0e-15);
                value += 2.0 * rdm1(t, v) * (4.0 * mo2e(nocc_+a+(nclosed_+v)*norb_, nclosed_+u+i*norb_) - mo2e(nocc_+a+(nclosed_+u)*norb_, nclosed_+v+i*norb_)
                                                                                                        - mo2e(nocc_+a+i*norb_, nclosed_+u+(nclosed_+v)*norb_));
              }
              hessian->element(va_offset+a+t*nvirt_, i+u*nclosed_) += value;
              hessian->element(i+u*nclosed_, va_offset+a+t*nvirt_) += value;
            }
          }
        }
      }
      // (at, bi)
      for (int i = 0; i != nclosed_; ++i) {
        for (int b = 0; b != nvirt_; ++b) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 0;
              for (int u = 0; u != nact_; ++u) {
                value += 2.0 * rdm1(t, u) * (4.0 * mo2e(nocc_+a+(nclosed_+u)*norb_, nocc_+b+i*norb_) - mo2e(nocc_+a+(nocc_+b)*norb_, nclosed_+u+i*norb_)
                                                                                                     - mo2e(nocc_+a+i*norb_, nclosed_+u+(nocc_+b)*norb_));
              }
              hessian->element(va_offset+a+t*nvirt_, vc_offset+b+i*nvirt_) += value;
              hessian->element(vc_offset+b+i*nvirt_, va_offset+a+t*nvirt_) += value;
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact_; ++u) {
        for (int j = 0; j != nclosed_; ++j) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              double value = 16.0 * mo2e(nclosed_+t+i*norb_, nclosed_+u+j*norb_) - 4.0 * mo2e(nclosed_+t+j*norb_, i+(nclosed_+u)*norb_)
                                                                                 - 4.0 * mo2e(nclosed_+t+(nclosed_+u)*norb_, i+j*norb_);
              for (int v = 0; v != nact_; ++v) {
                value -= 2.0 * rdm1(u, v) * (4.0 * mo2e(nclosed_+t+i*norb_, nclosed_+v+j*norb_) - mo2e(nclosed_+t+j*norb_, i+(nclosed_+v)*norb_)
                                                                                                - mo2e(nclosed_+t+(nclosed_+v)*norb_, i+j*norb_));
                value -= 2.0 * rdm1(t, v) * (4.0 * mo2e(nclosed_+v+i*norb_, nclosed_+u+j*norb_) - mo2e(nclosed_+v+j*norb_, nclosed_+u+i*norb_)
                                                                                                - mo2e(nclosed_+v+(nclosed_+u)*norb_, i+j*norb_));
              }
              hessian->element(i+t*nclosed_, j+u*nclosed_) += value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed_; ++j) {
        for (int a = 0; a != nvirt_; ++a) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              double value = 16.0 * mo2e(nclosed_+t+i*norb_, nocc_+a+j*norb_) - 4.0 * mo2e(nclosed_+t+j*norb_, nocc_+a+i*norb_)
                                                                              - 4.0 * mo2e(nclosed_+t+(nocc_+a)*norb_, i+j*norb_);
              for (int u = 0; u != nact_; ++u) {
                value -= 2.0 * rdm1(t, u) * (4.0 * mo2e(nclosed_+u+i*norb_, nocc_+a+j*norb_) - mo2e(nclosed_+u+j*norb_, nocc_+a+i*norb_)
                                                                                             - mo2e(nclosed_+u+(nocc_+a)*norb_, i+j*norb_));
              }
              hessian->element(i+t*nclosed_, vc_offset+a+j*nvirt_) += value;
              hessian->element(vc_offset+a+j*nvirt_, i+t*nclosed_) += value;
            }
          }
        }
      }
      // (ai, bj)
      for (int j = 0; j != nclosed_; ++j) {
        for (int b = 0; b != nvirt_; ++b) {
          for (int i = 0; i != nclosed_; ++i) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 16.0 * mo2e(nocc_+a+i*norb_, nocc_+b+j*norb_) - 4.0 * mo2e(nocc_+a+(nocc_+b)*norb_, i+j*norb_)
                                                                           - 4.0 * mo2e(nocc_+a+j*norb_, nocc_+b+i*norb_);
              hessian->element(vc_offset+a+i*nvirt_, vc_offset+b+j*nvirt_) += value;
            }
          }
        }
      }
#ifdef AAROT
      for (auto& block : act_rotblocks_) {
        const int istart = block.iorbstart;
        const int jstart = block.jorbstart;
        const int inorb = block.norb_i;
        const int jnorb = block.norb_j;
        const int offset = block.offset;

        // (ti, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact_; ++t) {
              for (int i = 0; i != nclosed_; ++i) {
                double value = 0.0;
                for (int w = 0; w != nact_; ++w) {
                  value += rdm1(w, jstart+v) * (8.0 * mo2e(nclosed_+istart+u+(nclosed_+w)*norb_, nclosed_+t+i*norb_)
                                              - 2.0 * mo2e(nclosed_+istart+u+(nclosed_+t)*norb_, nclosed_+w+i*norb_)
                                              - 2.0 * mo2e(nclosed_+istart+u+i*norb_, nclosed_+w+(nclosed_+t)*norb_));
                  value -= rdm1(w, istart+u) * (8.0 * mo2e(nclosed_+jstart+v+(nclosed_+w)*norb_, nclosed_+t+i*norb_)
                                              - 2.0 * mo2e(nclosed_+jstart+v+(nclosed_+t)*norb_, nclosed_+w+i*norb_)
                                              - 2.0 * mo2e(nclosed_+jstart+v+i*norb_, nclosed_+w+(nclosed_+t)*norb_));
                }
                if (t == jstart+v) {
                  for (int x = 0; x != nact_; ++x) {
                    for (int y = 0; y != nact_; ++y) {
                      value += rdm1(x, y) * (2.0 * mo2e(nclosed_+x+(nclosed_+y)*norb_, nclosed_+istart+u+i*norb_)
                                                 - mo2e(nclosed_+istart+u+(nclosed_+y)*norb_, nclosed_+x+i*norb_));
                    }
                  }
                }
                if (t == istart+u) {
                  for (int x = 0; x != nact_; ++x) {
                    for (int y = 0; y != nact_; ++y) {
                      value -= rdm1(x, y) * (2.0 * mo2e(nclosed_+x+(nclosed_+y)*norb_, nclosed_+jstart+v+i*norb_)
                                                 - mo2e(nclosed_+jstart+v+(nclosed_+y)*norb_, nclosed_+x+i*norb_));
                    }
                  }
                }
                hessian->element(i+t*nclosed_, aa_offset+offset+u+v*inorb) += value;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed_) += value;
              }
            }
          }
        }
        // (ai, tu)
        for (int u = 0; u != jnorb; ++u) {
          for (int t = 0; t != inorb; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              for (int a = 0; a != nvirt_; ++a) {
                double value = 0.0;
                for (int v = 0; v != nact_; ++v) {
                  value += rdm1(jstart+u, v) * (8.0 * mo2e(nclosed_+istart+t+(nclosed_+v)*norb_, nocc_+a+i*norb_)
                                              - 2.0 * mo2e(nclosed_+istart+t+(nocc_+a)*norb_, nclosed_+v+i*norb_)
                                              - 2.0 * mo2e(nclosed_+istart+t+i*norb_, nclosed_+v+(nocc_+a)*norb_));
                  value -= rdm1(istart+t, v) * (8.0 * mo2e(nclosed_+jstart+u+(nclosed_+v)*norb_, nocc_+a+i*norb_)
                                              - 2.0 * mo2e(nclosed_+jstart+u+(nocc_+a)*norb_, nclosed_+v+i*norb_)
                                              - 2.0 * mo2e(nclosed_+jstart+u+i*norb_, nclosed_+v+(nocc_+a)*norb_));
                }
                hessian->element(vc_offset+a+i*nvirt_, aa_offset+offset+t+u*inorb) += value;
                hessian->element(aa_offset+offset+t+u*inorb, vc_offset+a+i*nvirt_) += value;
              }
            }
          }
        }
      }
#endif
    } // end of 2-electron integral part

    // rdm2 and integral part
    {
      // (at, bu)
      for (int u = 0; u != nact_; ++u) {
        for (int b = 0; b != nvirt_; ++b) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 0.0;
              for (int v = 0; v != nact_; ++v) {
                for (int w = 0; w != nact_; ++w) {
                  value += 2.0 * rdm2(t, u, v, w) * mo2e(nocc_+a+(nocc_+b)*norb_, nclosed_+v+(nclosed_+w)*norb_)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nocc_+b+(nclosed_+v)*norb_, nocc_+a+(nclosed_+w)*norb_);
                  if (a == b) {
                    for (int x = 0; x != nact_; ++x) {
                      value -= rdm2(t, v, w, x) * mo2e(nclosed_+u+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_)
                             + rdm2(u, v, w, x) * mo2e(nclosed_+t+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_);
                    }
                  }
                }
              }
              hessian->element(va_offset+a+t*nvirt_, va_offset+b+u*nvirt_) += value;
            }
          }
        }
      }
      // (at, ui)
      for (int u = 0; u != nact_; ++u) {
        for (int i = 0; i != nclosed_; ++i) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              double value = 0.0;
              for (int v = 0; v != nact_; ++v) {
                for (int w = 0; w != nact_; ++w) {
                  value -= 2.0 * rdm2(t, u, v, w) * mo2e(nocc_+a+i*norb_, nclosed_+v+(nclosed_+w)*norb_)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nclosed_+v+i*norb_, nclosed_+w+(nocc_+a)*norb_);
                }
              }
              hessian->element(va_offset+a+t*nvirt_, i+u*nclosed_) += value;
              hessian->element(i+u*nclosed_, va_offset+a+t*nvirt_) += value;
            }
          }
        }
      }
      // (at, bi)
      for (int i = 0; i != nclosed_; ++i) {
        for (int b = 0; b != nclosed_; ++b) {
          for (int t = 0; t != nact_; ++t) {
            for (int a = 0; a != nvirt_; ++a) {
              if (a == b) {
                double value = 0.0;
                for (int u = 0; u != nact_; ++u) {
                  for (int v = 0; v != nact_; ++v) {
                    for (int w = 0; w != nact_; ++w) {
                      value -= rdm2(t, u, v, w) * mo2e(nclosed_+u+i*norb_, nclosed_+v+(nclosed_+w)*norb_);
                    }
                  }
                }
                hessian->element(va_offset+a+t*nvirt_, vc_offset+b+i*nvirt_) += value;
                hessian->element(vc_offset+b+i*nvirt_, va_offset+a+t*nvirt_) += value;
              }
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact_; ++u) {
        for (int j = 0; j != nclosed_; ++j) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              double value = 0.0;
              for (int v = 0; v != nact_; ++v) {
                for (int w = 0; w != nact_; ++w) {
                  value += 2.0 * rdm2(t, u, v, w) * mo2e(i+j*norb_, nclosed_+v+(nclosed_+w)*norb_)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nclosed_+v+j*norb_, nclosed_+w+i*norb_);
                  if (i == j) {
                    for (int x = 0; x != nact_; ++x) {
                      value -= rdm2(u, v, w, x) * mo2e(nclosed_+t+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_)
                             + rdm2(t, v, w, x) * mo2e(nclosed_+u+(nclosed_+v)*norb_, nclosed_+w+(nclosed_+x)*norb_);
                    }
                  }
                }
              }
              hessian->element(i+t*nclosed_, j+u*nclosed_) += value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed_; ++j) {
        for (int a = 0; a != nvirt_; ++a) {
          for (int t = 0; t != nact_; ++t) {
            for (int i = 0; i != nclosed_; ++i) {
              double value = 0.0;
              if (i == j) {
                for (int u = 0; u != nact_; ++u) {
                  for (int v = 0; v != nact_; ++v) {
                    for (int w = 0; w != nact_; ++w) {
                      value -= rdm2(t, u, v, w) * mo2e(nocc_+a+(nclosed_+u)*norb_, nclosed_+v+(nclosed_+w)*norb_);
                    }
                  }
                }
              }
              hessian->element(i+t*nclosed_, vc_offset+a+j*nvirt_) += value;
              hessian->element(vc_offset+a+j*nvirt_, i+t*nclosed_) += value;
            }
          }
        }
      }
    } // end of rdm2 and integral part

    // check for compute_gradient
    {
      cout << " * checking compute_grad" << endl;
      VectorB grad_vec(rotsize);
      copy_n(grad->data(), rotsize, grad_vec.data());
      VectorB grad_diff(rotsize);
      VectorB tmp(grad_vec - grad_check);
      for (int i = 0; i != rotsize; ++i) {
        if (fabs(grad_vec(i)) >= 1.0e-15) {
          grad_diff(i) = tmp(i) / grad_vec(i);
          if (grad_diff(i) > 1.0e-8) cout << i << " : " << grad_diff(i) << ", grad_vec = " << grad_vec(i) << endl;
        }
      }
      cout << "diff grad rms = " << grad_diff.rms() << endl;
    }
    // check for compute_denom
    {
      cout << " * checking compute_denom" << endl;
      VectorB denom_vec(rotsize);
      copy_n(denom->data(), rotsize, denom_vec.data());
      VectorB hessian_diag(rotsize);
      copy_n(hessian->diag().get(), rotsize, hessian_diag.data());
      VectorB denom_diff(rotsize);
      VectorB tmp(denom_vec - hessian_diag);
      for (int i = 0; i != rotsize; ++i) {
        if (fabs(denom_vec(i) >= 1.0e-15)) {
          denom_diff(i) = tmp(i) / denom_vec(i);
          if (denom_diff(i) > 1.0e-8) cout << i << " : " << denom_diff(i) << ", denom_vec = " << denom_vec(i) << endl;
        }
      }
      cout << "diff denom rms = " << denom_diff.rms() << endl;
    }
    // check for hessian_trial
    {
      cout << " * checking compute_hess_trial" << endl;
      auto rot = trot->clone();
      std::iota(rot->begin(), rot->begin()+va_offset, 1.0);
      std::iota(rot->begin()+va_offset, rot->begin()+vc_offset, 1.0);
#ifdef AAROT
      std::iota(rot->begin()+vc_offset, rot->begin()+aa_offset, 1.0);
      std::iota(rot->begin()+aa_offset, rot->end(), 1.0);
#else
      std::iota(rot->begin()+vc_offset, rot->end(), 1.0);
#endif
      auto hess_try = compute_hess_trial(rot, half, halfa, halfa_JJ, cfock, afock, qxr);
      VectorB hess_try_vec(rotsize);
      copy_n(hess_try->data(), rotsize, hess_try_vec.data());
      VectorB rot_vec(rot->size());
      copy_n(rot->data(), rotsize, rot_vec.data());
      VectorB hess_debug_vec(rotsize);
      dgemv_("N", hessian->ndim(), hessian->mdim(), 1.0, hessian->data(), rotsize, rot_vec.data(), 1, 0.0, hess_debug_vec.data(), 1);
      VectorB hess_t_diff(rotsize);
      VectorB tmp(hess_try_vec - hess_debug_vec);
      for (int i = 0; i != rotsize; ++i) {
        if (fabs(hess_try_vec(i) >= 1.0e-15)) {
          hess_t_diff(i) = tmp(i) / hess_try_vec(i);
          if (hess_t_diff(i) > 1.0e-8) cout << i << " : " << hess_t_diff(i) << ", hess_t_vec = " << hess_try_vec(i) << endl;
        }
      }
      cout << "diff hess_trial rms : " << hess_t_diff.rms() << endl;
    }
  }
  
#endif

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      
      shared_ptr<const ASD_DMRG_RotFile> sigma = compute_hess_trial(trot, half, halfa, halfa_JJ, cfock, afock, qxr);
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
      //muffle_->mute();
      if (err < max(thresh_micro_, stepsize*thresh_microstep_))
        break;

      trot = apply_denom(residual, denom, -epsilon, 1.0/lambda);
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    } // end of micro iter

    shared_ptr<const ASD_DMRG_RotFile> sol = solver.civec();
#ifdef AAROT
    shared_ptr<const Matrix> a = sol->unpack(act_rotblocks_);
#else
    shared_ptr<const Matrix> a = sol->unpack();
#endif
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
#ifdef AAROT
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
        blas::ax_plus_y_n(-2.0*cfock->element(nclosed_+v, nclosed_+jstart+j), rdm1->element_ptr(istart, v), inorb, target);
      }
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nclosed_+istart, jstart+j), inorb, target);
      blas::ax_plus_y_n(-2.0, qxr->transpose()->element_ptr(istart, nclosed_+jstart+j), inorb, target);
    }
  }
#endif

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
      vgaa = vgaa->transform_occ1(make_shared<Matrix>(rdm1));
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

#ifdef AAROT
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
//      const double e2j = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact_*(jstart+j)), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,jstart+j));

      // Fock related part
      for (int i = 0; i != inorb; ++i) {
//        const double e2i = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact_*(istart+i)), nact_*nact_*nact_, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,istart+i));
/*
        denom->ele_aa_offset(i, inorb, j, offset) = 2.0 * (*cfock)(nclosed_+istart+i, nclosed_+istart+i) * rdm1(jstart+j, jstart+j)
                                                    + 2.0 * (*cfock)(nclosed_+jstart+j, nclosed_+jstart+j) * rdm1(istart+i, istart+i)
                                                    - 4.0 * (*cfock)(nclosed_+istart+i, nclosed_+jstart+j) * rdm1(istart+i, jstart+j)
                                                    - 2.0 * fcd(istart+i, istart+i) - 2.0 * fcd(jstart+j, jstart+j)
                                                    + e2j + e2i;
*/
        denom->ele_aa_offset(i, inorb, j, offset) = 2.0 * (*cfock)(nclosed_+istart+i, nclosed_+istart+i) * rdm1(jstart+j, jstart+j)
                                                    + 2.0 * (*cfock)(nclosed_+jstart+j, nclosed_+jstart+j) * rdm1(istart+i, istart+i)
                                                    - 4.0 * (*cfock)(nclosed_+istart+i, nclosed_+jstart+j) * rdm1(istart+i, jstart+j)
                                                    - 2.0 * fcd(istart+i, istart+i) - 2.0 * fcd(jstart+j, jstart+j);
      }

/*
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
*/
    } // end of looping over second active index
  } // end of looping over blocks
#endif

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


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_hess_trial(shared_ptr<const ASD_DMRG_RotFile> trot, shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> halfa,
    shared_ptr<const DFHalfDist> halfa_JJ, shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {

  shared_ptr<ASD_DMRG_RotFile> sigma = trot->clone();

  shared_ptr<const Matrix> va = trot->va_mat();
  shared_ptr<const Matrix> ca = nclosed_ ? trot->ca_mat() : nullptr;
  shared_ptr<const Matrix> vc = nclosed_ ? trot->vc_mat() : nullptr;

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
    const Matrix tcoeff2c = vcoeff * *vc + acoeff * *ca->transpose();
    auto halft = geom_->df()->compute_half_transform(tcoeff2c);
    const Matrix gt = *compute_gd(halft, half, ccoeff);
    sigma->ax_plus_y_va(16.0, vcoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_ca(32.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_ca(-16.0, ccoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_vc(32.0, vcoeff % gt * ccoeff);
  }

  const Matrix tcoeff2a = nclosed_ ? (vcoeff * *va - ccoeff * *ca) : (vcoeff *  *va);
  shared_ptr<const DFHalfDist> halfta = geom_->df()->compute_half_transform(tcoeff2a);
  
  // g(t_va - t_ca)
  if (nclosed_) {
    shared_ptr<DFHalfDist> halftad = halfta->copy();
    halftad = halftad->transform_occ(make_shared<Matrix>(rdm1));
    const Matrix gt = *compute_gd(halftad, halfa_JJ, acoeff);
    sigma->ax_plus_y_ca(16.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_vc(16.0, vcoeff % gt * ccoeff);
  }

#ifdef AAROT
  // active-active 2-electron integral part
  if (nclosed_){
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norb_i;
      const int jnorb = block.norb_j;
      const int offset = block.offset;
      const int bsize = block.size;

      auto rotblock_aa = make_shared<Matrix>(inorb, jnorb);
      copy_n(trot->ptr_aa_offset(offset), bsize, rotblock_aa->data());

      const MatView acoeffi = coeff_->slice(nclosed_+istart, nclosed_+istart+inorb);
      const MatView acoeffj = coeff_->slice(nclosed_+jstart, nclosed_+jstart+jnorb);

      const MatView rdmxi = rdm1.slice(istart, istart+inorb);
      const MatView rdmxj = rdm1.slice(jstart, jstart+jnorb);

      const MatView cai = ca->slice(istart, istart+inorb);
      const MatView caj = ca->slice(jstart, jstart+jnorb);

      { // ai->tu
        const Matrix tcoeffv2c = vcoeff * *vc;
        auto halftv2c = geom_->df()->compute_half_transform(tcoeffv2c);
        const Matrix gtv2c = *compute_gd(halftv2c, half, ccoeff);
        sigma->ax_plus_y_aa_offset(16.0, acoeffi % gtv2c * acoeff * rdmxj , offset);
        sigma->ax_plus_y_aa_offset(-16.0, rdmxi % (acoeff % gtv2c * acoeffj), offset);
      }

      shared_ptr<DFHalfDist> halfxy = halfa->copy();
      halfxy = halfxy->transform_occ(make_shared<Matrix>(rdm1));
      const Matrix gtaa = *compute_gd(halfxy, halfa_JJ, acoeff);

      { // ti->uv
        const Matrix tcoeffc2a = ccoeff * *ca;
        auto halftc2a = geom_->df()->compute_half_transform(tcoeffc2a);
        const Matrix gtc2a = *compute_gd(halftc2a, halfa_JJ, acoeff);
        sigma->ax_plus_y_aa_offset(16.0, acoeffi % gtc2a * acoeff * rdmxj, offset);
        sigma->ax_plus_y_aa_offset(-16.0, rdmxi % (acoeff % gtc2a * acoeffj), offset);

        // \delta_{tv} part
        sigma->ax_plus_y_aa_offset(4.0, (ccoeff % gtaa * acoeffi) % caj, offset);
        // \delta_{tu} part
        sigma->ax_plus_y_aa_offset(-4.0, (ccoeff * cai) % gtaa * acoeffj, offset);
      }

      { // tu->ai and uv->it
        // \gamma_{uv} and \gamma_{vw} part
        const Matrix tcoeffi2x = acoeffi * *rotblock_aa ^ rdmxj;
        auto halftix = geom_->df()->compute_half_transform(tcoeffi2x);
        const Matrix gt1 = *compute_gd(halftix, halfa_JJ, acoeff);
        sigma->ax_plus_y_vc(16.0, vcoeff % gt1 * ccoeff);
        sigma->ax_plus_y_ca(16.0, ccoeff % gt1 * acoeff);
        // \gamma_{tu} and \gamma_{uw} part
        const Matrix tcoeffj2x = acoeffj ^ (rdmxi * *rotblock_aa);
        auto halftjx = geom_->df()->compute_half_transform(tcoeffj2x);
        const Matrix gt2 = *compute_gd(halftjx, halfa_JJ, acoeff);
        sigma->ax_plus_y_vc(-16.0, vcoeff % gt2 * ccoeff);
        sigma->ax_plus_y_ca(-16.0, ccoeff % gt2 * acoeff);
        // \delta_{tv} part
        sigma->ax_plus_y_ca_offset(4.0, ccoeff % gtaa * acoeffi * *rotblock_aa, jstart);
        sigma->ax_plus_y_ca_offset(-4.0, ccoeff % gtaa * acoeffj ^ *rotblock_aa, istart);
      }
    }
  }
#endif

  // terms with Qvec
  {
    shared_ptr<const Matrix> qaa = qxr->cut(nclosed_, nocc_);
    shared_ptr<const Matrix> qva = qxr->cut(nocc_, nocc_+nvirt_);
    shared_ptr<const Matrix> qca = qxr->cut(0, nclosed_);
    sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
    sigma->ax_plus_y_va(-2.0, *va * *qaa);
    if (nclosed_) {
      sigma->ax_plus_y_va(-2.0, *vc * *qca);
      sigma->ax_plus_y_ca(-2.0, *ca ^ *qaa);
      sigma->ax_plus_y_ca(-2.0, *ca * *qaa);
      sigma->ax_plus_y_ca(-2.0, *vc % *qva);
      sigma->ax_plus_y_vc(-2.0, *va ^ *qca);
      sigma->ax_plus_y_vc(-2.0, *qva ^ *ca);
    }

//#ifdef AAROT
#if 0
    // active-active part
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norb_i;
      const int jnorb = block.norb_j;
      const int offset = block.offset;
      const int bsize = block.size;

      auto rotblock_aa = make_shared<Matrix>(inorb, jnorb);
      copy_n(trot->ptr_aa_offset(offset), bsize, rotblock_aa->data());

      const MatView qvaj = qva->slice(jstart, jstart+jnorb);
      const MatView vai = va->slice(istart, istart+inorb);
      sigma->ax_plus_y_va_offset(2.0, qvaj ^ *rotblock_aa, istart);
      sigma->ax_plus_y_aa_offset(2.0, vai % qvaj, offset);

      if (nclosed_) {
        const MatView qcai = qca->slice(istart, istart+inorb);
        const MatView qcaj = qca->slice(jstart, jstart+jnorb);
        const MatView caj = ca->slice(jstart, jstart+jnorb);
        const MatView cai = ca->slice(istart, istart+inorb);
        sigma->ax_plus_y_ca_offset(2.0, qcai * *rotblock_aa, jstart);
        sigma->ax_plus_y_ca_offset(-2.0, qcaj ^ *rotblock_aa, istart);
        sigma->ax_plus_y_aa_offset(2.0, qcai % caj, offset);
        sigma->ax_plus_y_aa_offset(-2.0, cai % qcaj, offset);
      }

      for (auto& block2 : act_rotblocks_) {
        const int istart2 = block2.iorbstart;
        const int jstart2 = block2.jorbstart;
        const int inorb2 = block2.norb_i;
        const int jnorb2 = block2.norb_j;
        const int offset2 = block2.offset;
        const int bsize2 = block2.size;

        auto rotblock2_aa = make_shared<Matrix>(inorb2, jnorb2);
        copy_n(trot->ptr_aa_offset(offset2), bsize2, rotblock2_aa->data());

        // \delta_{vu} part
        if (jstart >= istart2) {
          shared_ptr<const Matrix> qaaiJ = qaa->get_submatrix(istart, jstart2, inorb, jnorb2);
          shared_ptr<const Matrix> qaaJi = qaa->get_submatrix(jstart2, istart, jnorb2, inorb);
          auto tmpmat = make_shared<const Matrix>((*qaaiJ ^ *rotblock2_aa) + *((*rotblock2_aa * *qaaJi).transpose()));
          for (int j = 0; j != jnorb; ++j)
            blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j+jstart-istart2), inorb, sigma->ptr_aa_offset(offset)+j*inorb);
        }

        // \delta_{tw} part
        if (jstart2 >= istart) {
          shared_ptr<const Matrix> qaaIj = qaa->get_submatrix(istart2, jstart, inorb2, jnorb);
          shared_ptr<const Matrix> qaajI = qaa->get_submatrix(jstart, istart2, jnorb, inorb2);
          auto tmpmat = make_shared<const Matrix>((*rotblock2_aa % *qaaIj) + *((*qaajI * *rotblock2_aa).transpose()));
          for (int j = 0; j != jnorb; ++j)
            blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j), jnorb2, sigma->ptr_aa_offset(offset)+j*inorb+jstart2-istart);
        }

        // \delta_{wu} part, two blocks are identical
        if (jstart2 == jstart) {
          assert(istart == istart2 && inorb == inorb2);
          shared_ptr<const Matrix> qaaiI = qaa->get_submatrix(istart, istart2, inorb, inorb2);
          shared_ptr<const Matrix> qaaIi = qaa->get_submatrix(istart2, istart, inorb2, inorb);
          auto tmpmat = make_shared<const Matrix>((*qaaiI * *rotblock2_aa) + (*qaaIi % *rotblock2_aa));
          sigma->ax_plus_y_aa_offset(-2.0, *tmpmat, offset);
        }

        // \delta_{tv} part
        {
          shared_ptr<const Matrix> qaaJj = qaa->get_submatrix(jstart2, jstart, jnorb2, jnorb);
          shared_ptr<const Matrix> qaajJ = qaa->get_submatrix(jstart, jstart2, jnorb, jnorb2);
          auto tmpmat = make_shared<const Matrix>((*rotblock2_aa * *qaaJj) + (*rotblock2_aa ^ *qaajJ));
          if (istart < istart2) {
            for (int j = 0; j != jnorb; ++j)
              blas::ax_plus_y_n(-2.0, tmpmat->element_ptr(0, j), inorb2, sigma->ptr_aa_offset(offset)+j*inorb+istart2-istart);
          } else {
            auto tmpsub = tmpmat->get_submatrix(istart-istart2, 0, inorb, jnorb);
            sigma->ax_plus_y_aa_offset(-2.0, *tmpsub, offset);
          }
        }
      } // end of looping over block2
    } // end of looping over block
#endif
  } // end of Qvec part

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

//#ifdef AAROT
#if 0
    // active-active part
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norb_i;
      const int jnorb = block.norb_j;
      const int offset = block.offset;
      const int bsize = block.size;

      auto rotblock_aa = make_shared<Matrix>(inorb, jnorb);
      copy_n(trot->ptr_aa_offset(offset), bsize, rotblock_aa->data());

      const MatView acoeffi = coeff_->slice(nclosed_+istart, nclosed_+istart+inorb);
      const MatView acoeffj = coeff_->slice(nclosed_+jstart, nclosed_+jstart+jnorb);

      // (at, uv) and (ti, uv) terms
      shared_ptr<Matrix> rotmat2i;
      shared_ptr<Matrix> rotmat2j;
      {
        rotmat2i = make_shared<Matrix>(nact_, inorb);
        for (int i = 0; i != inorb; ++i)
          *rotmat2i->element_ptr(istart+i, i) = 1.0;
        rotmat2j = make_shared<Matrix>(nact_, jnorb);
        for (int j = 0; j != jnorb; ++j)
          *rotmat2j->element_ptr(jstart+j, j) = 1.0;
        const Matrix tcoeff2i = acoeffj ^ *rotblock_aa;
        const Matrix tcoeff2j = acoeffi * *rotblock_aa;
        shared_ptr<const DFHalfDist> halfta_2i = geom_->df()->compute_half_transform(tcoeff2i);
        shared_ptr<const DFHalfDist> halfta_2j = geom_->df()->compute_half_transform(tcoeff2j);
        shared_ptr<DFFullDist> fullta_2i = halfta_2i->compute_second_transform(acoeff);
        shared_ptr<DFFullDist> fullta_2j = halfta_2j->compute_second_transform(acoeff);
        fullta_2i = fullta_2i->transform_occ1(rotmat2i->transpose());
        fullta_2j = fullta_2j->transform_occ1(rotmat2j->transpose());
        shared_ptr<const DFFullDist> fulltas_2i = fullta_2i->swap();
        shared_ptr<const DFFullDist> fulltas_2j = fullta_2j->swap();
        fullta_2i->ax_plus_y(1.0, fulltas_2i);
        fullta_2j->ax_plus_y(1.0, fulltas_2j);
        shared_ptr<const DFFullDist> fulltaD_2i = fullta_2i->apply_2rdm(*asd_dmrg_->rdm2_av());
        shared_ptr<const DFFullDist> fulltaD_2j = fullta_2j->apply_2rdm(*asd_dmrg_->rdm2_av());
        shared_ptr<DFFullDist> fulltaD_ij = fulltaD_2j->copy();
        fulltaD_ij->ax_plus_y(-1.0, fulltaD_2i); // TODO check here later
        shared_ptr<const Matrix> qp = halfa_JJ->form_2index(fulltaD_ij, 1.0);
        
        shared_ptr<DFFullDist> fullaaD_2i = fullaaD->copy();
        shared_ptr<DFFullDist> fullaaD_2j = fullaaD->copy();
        fullaaD_2i = fullaaD_2i->transform_occ1(rotmat2i);
        fullaaD_2j = fullaaD_2j->transform_occ1(rotmat2j);
        shared_ptr<const Matrix> qpp2i = halfta_2i->form_2index(fullaaD_2i, 1.0);
        shared_ptr<const Matrix> qpp2j = halfta_2j->form_2index(fullaaD_2j, 1.0);

        sigma->ax_plus_y_va(4.0, vcoeff % (*qp + *qpp2j - *qpp2i));
        if (nclosed_)
          sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp + *qpp2j - *qpp2i));

        shared_ptr<const DFFullDist> fullaaDs_2i = fullaaD_2i->swap();
        shared_ptr<const DFFullDist> fullaaDs_2j = fullaaD_2j->swap();
        shared_ptr<const DFFullDist> fulltta_2i = halfta->compute_second_transform(acoeffi);
        shared_ptr<const DFFullDist> fulltta_2j = halfta->compute_second_transform(acoeffj);
        shared_ptr<const Matrix> Qpp1 = fulltta_2i->form_2index(fullaaDs_2j, 1.0);
        shared_ptr<const Matrix> Qpp2 = fullaaDs_2i->form_2index(fulltta_2j, 1.0);
        shared_ptr<const Matrix> Qpp = make_shared<Matrix>(*Qpp1 - *Qpp2);

        auto tmp_2i = fulltaD->copy();
        tmp_2i = tmp_2i->transform_occ1(rotmat2i);
        shared_ptr<const DFFullDist> fullttaD_2i = tmp_2i->swap();
        auto tmp_2j = fulltaD->copy();
        tmp_2j = tmp_2j->transform_occ1(rotmat2j);
        shared_ptr<const DFFullDist> fullttaD_2j = tmp_2j->swap();
        shared_ptr<const DFFullDist> fullaa_2i = halfa_JJ->compute_second_transform(acoeffi);
        shared_ptr<const DFFullDist> fullaa_2j = halfa_JJ->compute_second_transform(acoeffj);
        shared_ptr<const Matrix> Qp1 = fullaa_2i->form_2index(fullttaD_2j, 1.0);
        shared_ptr<const Matrix> Qp2 = fullttaD_2i->form_2index(fullaa_2j, 1.0);
        shared_ptr<const Matrix> Qp = make_shared<Matrix>(*Qp1 - *Qp2);

        sigma->ax_plus_y_aa_offset(4.0, (*Qpp + *Qp), offset);
      }

      // (tu, vw) part
      {
        auto qp_aa = make_shared<Matrix>(inorb, jnorb);
        auto qpp_aa = qp_aa->clone();
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norb_i;
          const int jnorb2 = block2.norb_j;
          const int offset2 = block2.offset;
          const int bsize2 = block2.size;

          auto rotblock2_aa = make_shared<Matrix>(inorb2, jnorb2);
          copy_n(trot->ptr_aa_offset(offset2), bsize2, rotblock2_aa->data());

          const MatView acoeffi2 = coeff_->slice(nclosed_+istart2, nclosed_+istart2+inorb2);
          const MatView acoeffj2 = coeff_->slice(nclosed_+jstart2, nclosed_+jstart2+jnorb2);

          auto rotmat2i2 = make_shared<Matrix>(nact_, inorb2);
          for (int i = 0; i != inorb2; ++i)
            *rotmat2i2->element_ptr(istart2+i, i) = 1.0;
          auto rotmat2j2 = make_shared<Matrix>(nact_, jnorb2);
          for (int j = 0; j != jnorb2; ++j)
            *rotmat2j2->element_ptr(jstart2+j, j) = 1.0; 
          
          const Matrix matv2w = acoeffi2 * *rotblock2_aa;
          shared_ptr<const DFHalfDist> halfav2w = geom_->df()->compute_half_transform(matv2w);
          shared_ptr<const Matrix> tmpout;
          {
            shared_ptr<DFFullDist> full_wu;
            shared_ptr<DFFullDist> full_wt;
            auto fullaaD1 = fullaaD->copy();
            fullaaD1 = fullaaD1->transform_occ1(rotmat2j);
            fullaaD1 = fullaaD1->swap();
            full_wu = fullaaD1->copy();
            full_wu = full_wu->transform_occ1(rotmat2j2);
            full_wt = halfav2w->compute_second_transform(acoeffi);
            tmpout = full_wt->form_2index(full_wu, 1.0);
            assert(tmpout->size() == tmpout->size());
            copy_n(tmpout->data(), tmpout->size(), qpp_aa->data());

            full_wu = fullaaD1->copy();
            full_wu = full_wu->transform_occ1(make_shared<Matrix>(*rotmat2i2 * *rotblock2_aa));
            shared_ptr<const DFHalfDist> halfa2 = geom_->df()->compute_half_transform(acoeffj2);
            full_wt = halfa2->compute_second_transform(acoeffi);
            tmpout = full_wt->form_2index(full_wu, 1.0);
            assert(tmpout->size() == qpp_aa->size());
            blas::ax_plus_y_n(-1.0, tmpout->data(), tmpout->size(), qpp_aa->data());

            auto fullaaD2 = fullaaD->copy();
            fullaaD2 = fullaaD2->transform_occ1(rotmat2i);
            fullaaD2 = fullaaD2->swap();
            full_wt = fullaaD2->copy();
            full_wt = full_wt->transform_occ1(rotmat2j2);
            full_wu = halfav2w->compute_second_transform(acoeffj);
            tmpout = full_wt->form_2index(full_wu, 1.0);
            assert(tmpout->size() == qpp_aa->size());
            blas::ax_plus_y_n(-1.0, tmpout->data(), tmpout->size(), qpp_aa->data());

            full_wt = fullaaD2->copy();
            full_wt = full_wt->transform_occ1(make_shared<Matrix>(*rotmat2i2 * *rotblock2_aa));
            full_wu = halfa2->compute_second_transform(acoeffj);
            tmpout = full_wt->form_2index(full_wu, 1.0);
            assert(tmpout->size() == qpp_aa->size());
            blas::ax_plus_y_n(1.0, tmpout->data(), tmpout->size(), qpp_aa->data());
          }

          {
            shared_ptr<DFFullDist> full_yu;
            shared_ptr<DFFullDist> full_yt;
            shared_ptr<DFFullDist> full_wx = halfav2w->compute_second_transform(acoeff);
            full_wx = full_wx->transform_occ1(rotmat2j2->transpose());
            shared_ptr<const DFFullDist> full_wxs = full_wx->swap();
            full_wx->ax_plus_y(1.0, full_wxs);
            shared_ptr<const DFFullDist> full_dwx = full_wx->apply_2rdm(*asd_dmrg_->rdm2_av());
            
            auto full_uy = full_dwx->copy();
            full_uy = full_uy->transform_occ1(rotmat2j);
            full_yu = full_uy->swap();
            auto full_ty = fullaa->copy();
            full_ty = full_ty->transform_occ1(rotmat2i);
            full_yt = full_ty->swap();
            tmpout = full_yt->form_2index(full_yu, 1.0);
            assert(tmpout->size() == qp_aa->size());
            copy_n(tmpout->data(), tmpout->size(), qp_aa->data());

            full_ty = full_dwx->copy();
            full_ty = full_ty->transform_occ1(rotmat2i);
            full_yt = full_ty->swap();
            full_uy = fullaa->copy();
            full_uy = full_uy->transform_occ1(rotmat2j);
            full_yu = full_uy->swap();
            tmpout = full_yt->form_2index(full_yu, 1.0);
            assert(tmpout->size() == qp_aa->size());
            blas::ax_plus_y_n(-1.0, tmpout->data(), tmpout->size(), qp_aa->data());

            const Matrix w2v = acoeffj2 ^ *rotblock2_aa;
            auto full_xv = halfa->compute_second_transform(w2v);
            auto full_vx = full_xv->swap();
            full_vx = full_vx->transform_occ1(rotmat2i2->transpose());
            full_xv = full_vx->swap();
            full_vx->ax_plus_y(1.0, full_xv);
            shared_ptr<const DFFullDist> full_dvx = full_vx->apply_2rdm(*asd_dmrg_->rdm2_av());

            full_uy = full_dvx->copy();
            full_uy = full_uy->transform_occ1(rotmat2j);
            full_yu = full_uy->swap();
            full_ty = fullaa->copy();
            full_ty = full_ty->transform_occ1(rotmat2i);
            full_yt = full_ty->swap();
            tmpout = full_yt->form_2index(full_yu, 1.0);
            assert(tmpout->size() == qp_aa->size());
            blas::ax_plus_y_n(-1.0, tmpout->data(), tmpout->size(), qp_aa->data());

            full_ty = full_dvx->copy();
            full_ty = full_ty->transform_occ1(rotmat2i);
            full_yt = full_ty->swap();
            full_uy = fullaa->copy();
            full_uy = full_uy->transform_occ1(rotmat2j);
            full_yu = full_uy->swap();
            tmpout = full_yt->form_2index(full_yu, 1.0);
            assert(tmpout->size() == qp_aa->size());
            blas::ax_plus_y_n(1.0, tmpout->data(), tmpout->size(), qp_aa->data());
          }

          sigma->ax_plus_y_aa_offset(4.0, (*qp_aa + *qpp_aa), offset);
          qp_aa->zero();
          qpp_aa->zero(); // TODO not necessary, just in case, to be deleted..
        }
      }
    } // end of looping over block
#endif
  } // end of Q' and Q'' part

  // Fock related terms
  {
    // construct submatrices
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

#ifdef AAROT
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
      const MatView rdmxi = rdm1.slice(istart, istart+inorb);
      const MatView rdmxj = rdm1.slice(jstart, jstart+jnorb);

      { // (at, uv)
        const MatView fcvai = fcva->slice(istart, istart+inorb);
        const MatView fcvaj = fcva->slice(jstart, jstart+jnorb);
        const MatView vai = va->slice(istart, istart+inorb);
        const MatView vaj = va->slice(jstart, jstart+jnorb);
    
        sigma->ax_plus_y_va_offset(2.0, *fcva * rdmxj ^ *rotblock_aa, istart);
        sigma->ax_plus_y_va(4.0, fcvai * *rotblock_aa ^ rdmxj);
        sigma->ax_plus_y_va_offset(-2.0, *fcva * rdmxi * *rotblock_aa, jstart);
        sigma->ax_plus_y_va(-4.0, fcvaj ^ (rdmxi * *rotblock_aa));
        sigma->ax_plus_y_aa_offset(2.0, vai % *fcva * rdmxj, offset);
        sigma->ax_plus_y_aa_offset(4.0, fcvai % *va * rdmxj, offset);
        sigma->ax_plus_y_aa_offset(-2.0, (*fcva * rdmxi) % vaj, offset);
        sigma->ax_plus_y_aa_offset(-4.0, (*va * rdmxi) % fcvaj, offset);
      }
  
      if (nclosed_) { // (ti, uv)
        const MatView fccai = fcca->slice(istart, istart+inorb);
        const MatView fccaj = fcca->slice(jstart, jstart+jnorb);
        const MatView cai = ca->slice(istart, istart+inorb);
        const MatView caj = ca->slice(jstart, jstart+jnorb);
        
        sigma->ax_plus_y_ca_offset(4.0, fccai * *rotblock_aa, jstart);
        sigma->ax_plus_y_ca_offset(-4.0, fccaj ^ *rotblock_aa, istart);
        sigma->ax_plus_y_ca_offset(2.0, *fcca * rdmxi * *rotblock_aa, jstart);
        sigma->ax_plus_y_ca_offset(-2.0, *fcca * rdmxj ^ *rotblock_aa, istart);
        sigma->ax_plus_y_ca(4.0, fccaj ^ (rdmxi * *rotblock_aa));
        sigma->ax_plus_y_ca(-4.0, fccai * *rotblock_aa ^ rdmxj);
        sigma->ax_plus_y_aa_offset(4.0, fccai % caj, offset);
        sigma->ax_plus_y_aa_offset(-4.0, cai % fccaj, offset);
        sigma->ax_plus_y_aa_offset(2.0, (*fcca * rdmxi) % caj, offset);
        sigma->ax_plus_y_aa_offset(-2.0, (cai % *fcca) * rdmxj, offset);
        sigma->ax_plus_y_aa_offset(4.0, (*ca * rdmxi) % fccaj, offset);
        sigma->ax_plus_y_aa_offset(-4.0, fccai % *ca * rdmxj, offset);
      }
  
      { // (tu, vw)
        shared_ptr<const Matrix> fcaai = fcaa->slice_copy(istart, istart+inorb);
        shared_ptr<const Matrix> fcaaj = fcaa->slice_copy(jstart, jstart+jnorb);
  
        // need another loop over (vw) active pairs
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norb_i;
          const int jnorb2 = block2.norb_j;
          const int offset2 = block2.offset;
          const int bsize2 = block2.size;
          
          auto rotblock_aa2 = make_shared<Matrix>(inorb2, jnorb2);
          copy_n(trot->ptr_aa_offset(offset2), bsize2, rotblock_aa2->data());
  
          // used capitalized char for block2 orbitals
          const MatView fcaaI = fcaa->slice(istart2, istart2+inorb2);
          const MatView fcaaJ = fcaa->slice(jstart2, jstart2+jnorb2);
          shared_ptr<const Matrix> fcaaiI = fcaa->get_submatrix(istart, istart2, inorb, inorb2);
          shared_ptr<const Matrix> fcaaiJ = fcaa->get_submatrix(istart, jstart2, inorb, jnorb2);
          shared_ptr<const Matrix> fcaajI = fcaa->get_submatrix(jstart, istart2, jnorb, inorb2);
          shared_ptr<const Matrix> fcaajJ = fcaa->get_submatrix(jstart, jstart2, jnorb, jnorb2);
          
          shared_ptr<const Matrix> rdmiI = rdm1.get_submatrix(istart, istart2, inorb, inorb2);
          shared_ptr<const Matrix> rdmiJ = rdm1.get_submatrix(istart, jstart2, inorb, jnorb2);
          shared_ptr<const Matrix> rdmjI = rdm1.get_submatrix(jstart, istart2, jnorb, inorb2);
          shared_ptr<const Matrix> rdmjJ = rdm1.get_submatrix(jstart, jstart2, jnorb, jnorb2);
          const MatView rdmxI = rdm1.slice(istart2, istart2+inorb2);
          const MatView rdmxJ = rdm1.slice(jstart2, jstart2+jnorb2);
  
          sigma->ax_plus_y_aa_offset(4.0, *fcaaiI * *rotblock_aa2 ^ *rdmjJ, offset);
          sigma->ax_plus_y_aa_offset(-4.0, *fcaaiJ ^ (*rdmjI * *rotblock_aa2), offset);
          sigma->ax_plus_y_aa_offset(-4.0, *rdmiJ ^ (*fcaajI * *rotblock_aa2), offset);
          sigma->ax_plus_y_aa_offset(4.0, *rdmiI * *rotblock_aa2 ^ *fcaajJ, offset);
  
          // \delta_{uv} part
          if (jstart >= istart2) {
            auto tmpmat = make_shared<const Matrix>((*fcaai % rdmxJ ^ *rotblock_aa2) + (rdmxi % fcaaJ ^ *rotblock_aa2));
            for (int j = 0; j != jnorb; ++j)
              blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j+jstart-istart2), inorb, sigma->ptr_aa_offset(offset)+j*inorb);
          }
          
          // \delta_{uw} part, two blocks are identical
          if (jstart == jstart2) {
            auto tmpmat = make_shared<const Matrix>(*fcaai % rdmxI * *rotblock_aa + rdmxi % fcaaI * *rotblock_aa);
            sigma->ax_plus_y_aa_offset(-2.0, *tmpmat, offset);
          }
  
          // \delta_{tw} part
          if (jstart2 >= istart) {
            auto tmpmat = make_shared<const Matrix>((fcaaI * *rotblock_aa2) % rdmxj + (rdmxI * *rotblock_aa2) % *fcaaj);
            for (int j = 0; j != jnorb; ++j)
              blas::ax_plus_y_n(2.0, tmpmat->element_ptr(0, j), jnorb2, sigma->ptr_aa_offset(offset)+j*inorb+jstart2-istart);
          }
  
          // \delta_{tv} part
          {
            auto tmpmat = make_shared<const Matrix>(((*rotblock_aa2 ^ fcaaJ) * rdmxj) + ((*rotblock_aa2 ^ rdmxJ) * *fcaaj));
            if (istart < istart2) {
              for (int j = 0; j != jnorb; ++j)
                blas::ax_plus_y_n(-2.0, tmpmat->element_ptr(0, j), inorb2, sigma->ptr_aa_offset(offset)+j*inorb+istart2-istart);
            } else {
              auto tmpsub = tmpmat->get_submatrix(istart-istart2, 0, inorb, jnorb);
              sigma->ax_plus_y_aa_offset(-2.0, *tmpsub, offset);
            }
          }
        } // end of looping over block2
      }
    } // end of looping over blocks
#endif
  } // end of Fock part
  
  sigma->scale(0.5);

  return sigma;
}


void ASD_DMRG_Second::trans_natorb() {
  auto trans = make_shared<Matrix>(nact_, nact_);
  trans->add_diag(2.0);
  blas::ax_plus_y_n(-1.0, asd_dmrg_->rdm1_av()->data(), nact_*nact_, trans->data());

  VectorB occup(nact_);
  trans->diagonalize(occup);

  asd_dmrg_->rotate_rdms(trans);

  auto cnew = make_shared<Coeff>(*coeff_);
  cnew->copy_block(0, nclosed_, cnew->ndim(), nact_, coeff_->slice(nclosed_, nocc_) * *trans);
  coeff_ = cnew;
}


