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

#ifdef DEBUG_Hess
#include <src/multi/casscf/cassecond.h>
#endif

using namespace std;
using namespace bagel;

void ASD_DMRG_Second::compute() {

#ifdef DEBUG_Hess
  auto casscf_input = input_->get_child("casscf");
  auto casscf = make_shared<CASSecond>(casscf_input, asd_dmrg_->multisite()->sref()->geom(), asd_dmrg_->multisite()->sref()); // same orbital ordering with ASD
  casscf->compute();
#endif

  for (int iter = 0; iter != max_iter_; ++iter) {
    
    // first obtain RDM from ASD_DMRG
    {
      if (iter) asd_dmrg_->update_multisite(coeff_);
      asd_dmrg_->compute(!iter);
      asd_dmrg_->compute_rdm12();
      trans_natorb();
    }
    
    auto mref = asd_dmrg_->multisite()->sref();
    const int nclosed = mref->nclosed();
    const int nact = mref->nact();
    const int nocc = nclosed + nact;

    shared_ptr<const Matrix> cfockao = nclosed ? make_shared<Fock<1>>(mref->geom(), mref->hcore(), nullptr, coeff_->slice(0, nclosed), true/*store*/, true/*rhf*/)
                                                : make_shared<Matrix>(*mref->hcore());
    shared_ptr<const Matrix> afockao = compute_active_fock(coeff_->slice(nclosed, nocc), asd_dmrg_->rdm1_av());
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const Matrix> qxr = compute_qvec(coeff_->slice(nclosed, nocc), asd_dmrg_->rdm2_av());

    shared_ptr<const ASD_DMRG_RotFile> grad = compute_gradient(cfock, afock, qxr);

    // check gradient and break if converged
    const double gradient = grad->rms();
    print_iteration(iter, asd_dmrg_->energies(), gradient);
    if (gradient < thresh_) {
      cout << endl << "    * Second-Order Optimization Converged. *" << endl << endl;
      break;
    }

    // half-transformed integrals (with JJ)
    shared_ptr<const DFHalfDist> half_1j = nclosed ? dynamic_pointer_cast<const Fock<1>>(cfockao)->half() : nullptr;
    shared_ptr<const DFHalfDist> half = nclosed ? half_1j->apply_J() : nullptr;
    shared_ptr<const DFHalfDist> halfa = mref->geom()->df()->compute_half_transform(coeff_->slice(nclosed, nocc));
    shared_ptr<const DFHalfDist> halfa_JJ = halfa->apply_JJ();

    // compute_denominator
    shared_ptr<const ASD_DMRG_RotFile> denom = compute_denom(half, half_1j, halfa, halfa_JJ, cfock, afock);

    AugHess<ASD_DMRG_RotFile> solver(max_micro_iter_, grad);

    // initial trial vector
    shared_ptr<ASD_DMRG_RotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    trot->normalize();

#ifdef DEBUG_Hess
  // build Hessian matrix
  const int nvirt = mref->nvirt();
  const int norb = coeff_->mdim();
  const int rotsize = denom->size();
  Matrix rdm1(nact, nact);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1.size(), rdm1.data());
  auto rdm2 = *asd_dmrg_->rdm2_av();
  Matrix mo2e(norb, norb);
  {
    const MatView coeff(*coeff_);
    shared_ptr<const DFHalfDist> halfx = mref->geom()->df()->compute_half_transform(coeff);
    shared_ptr<const DFFullDist> fullx_1j = halfx->compute_second_transform(coeff)->apply_J();
    mo2e = *fullx_1j->form_4index(fullx_1j, 1.0);
  const int vc_offset = va_offset + nvirt * nact;
#ifdef AAROT
  const int aa_offset = vc_offset + nvirt * nclosed;
#endif
  const Matrix fock = *cfock + *afock;
  shared_ptr<const Matrix> cfk_xt = cfock->get_submatrix(0, nclosed, norb, nact);
  const Matrix cfkd = *cfk_xt * rdm1;
  
  VectorB grad_check(rotsize);
  {
    // closed-active
    for (int t = 0; t != nact; ++t) {
      for (int i = 0; i != nclosed; ++i) {
        double value = 4.0 * fock(i, nclosed+t) - 2.0 * cfkd(i, t);
        for (int u = 0; u != nact; ++u) {
          for (int v = 0; v != nact; ++v) {
            for (int w = 0; w != nact; ++w) {
              value -= 2.0 * rdm2(t, u, v, w) * mo2e(nclosed+u+i*norb, nclosed+v+(nclosed+w)*norb);
            }
          }
        }
        grad_check(i+t*nclosed) = value;
      }
    }
    // virtual-active
    for (int t = 0; t != nact; ++t) {
      for (int a = 0; a != nvirt; ++a) {
        double value = 2.0 * cfkd(nocc+a, t);
        for (int u = 0; u != nact; ++u) {
          for (int v = 0; v != nact; ++v) {
            for (int w = 0; w != nact; ++w) {
              value += 2.0 * rdm2(t, u, v, w) * mo2e(nocc+a+(nclosed+u)*norb, nclosed+v+(nclosed+w)*norb);
            }
          }
        }
        grad_check(va_offset+a+t*nvirt) = value;
      }
    }
    // virtual-closed
    for (int i = 0; i != nclosed; ++i) {
      for (int a = 0; a != nvirt; ++a) {
        grad_check(vc_offset+a+i*nvirt) = 4.0 * fock(nocc+a, i);
      }
    }
#ifdef AAROT
    // active-active
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norbi;
      const int jnorb = block.norbj;
      const int offset = block.offset;

      for (int t = 0; t != inorb; ++t) {
        for (int u = 0; u != jnorb; ++u) {
          double value = 2.0 * (cfkd(nclosed+istart+t, jstart+u) - cfkd(nclosed+jstart+u, istart+t));
          for (int v = 0; v != nact; ++v) {
            for (int w = 0; w != nact; ++w) {
              for (int x = 0; x != nact; ++x) {
                value += 2.0 * rdm2(jstart+u, v, w, x) * mo2e(nclosed+istart+t+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb)
                       - 2.0 * rdm2(istart+t, v, w, x) * mo2e(nclosed+jstart+u+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb);
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
      for (int u = 0; u != nact; ++u) {
        for (int b = 0; b != nvirt; ++b) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 2.0 * rdm1(t, u) * cfock->element(nocc+a, nocc+b);
              if (a == b) value -= cfkd(nclosed+u, t) + cfkd(nclosed+t, u);
              hessian->element(va_offset+a+t*nvirt, va_offset+b+u*nvirt) = value;
            }
          }
        }
      }
      // (at, ui)
      for (int u = 0; u != nact; ++u) {
        for (int i = 0; i != nclosed; ++i) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = -2.0 * rdm1(t, u) * cfock->element(nocc+a, i);
              if (t == u) value += 2.0 * fock(nocc+a, i);
              hessian->element(va_offset+a+t*nvirt, i+u*nclosed) = value;
              hessian->element(i+u*nclosed, va_offset+a+t*nvirt) = value;
            }
          }
        }
      }
      // (bi, at)
      for (int t = 0; t != nact; ++t) {
        for (int a = 0; a != nvirt; ++a) {
          for (int i = 0; i != nclosed; ++i) {
            for (int b = 0; b != nvirt; ++b) {
              if (a == b) {
                double value = -2.0 * fock(i, nclosed+t) - cfkd(i, t);
                hessian->element(vc_offset+b+i*nvirt, va_offset+a+t*nvirt) = value;
                hessian->element(va_offset+a+t*nvirt, vc_offset+b+i*nvirt) = value;
              }
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact; ++u) {
        for (int j = 0; j != nclosed; ++j) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              double value = 2.0 * rdm1(t, u) * cfock->element(i, j);
              if (i == j) {
                value += 4.0 * fock(nclosed+t, nclosed+u) - cfkd(nclosed+t, u) - cfkd(nclosed+u, t);
              }
              if (t == u) {
                value -= 4.0 * fock(i, j);
              }
              hessian->element(i+t*nclosed, j+u*nclosed) = value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed; ++j) {
        for (int a = 0; a != nvirt; ++a) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              if (i == j) {
                double value = 4.0 * fock(nocc+a, nclosed+t) - cfkd(nocc+a, t);
                hessian->element(i+t*nclosed, vc_offset+a+j*nvirt) = value;
                hessian->element(vc_offset+a+j*nvirt, i+t*nclosed) = value;
              }
            }
          }
        }
      }
      // (ai, bj)
      for (int j = 0; j != nclosed; ++j) {
        for (int b = 0; b != nvirt; ++b) {
          for (int i = 0; i != nclosed; ++i) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 0;
              if (a == b) {
                value -= 4.0 * fock(i, j);
              }
              if (i == j) {
                value += 4.0 * fock(nocc+a, nocc+b);
              }
              hessian->element(vc_offset+a+i*nvirt, vc_offset+b+j*nvirt) = value;
            }
          }
        }
      }
#ifdef AAROT
      for (auto& block : act_rotblocks_) {
        const int istart = block.iorbstart;
        const int jstart = block.jorbstart;
        const int inorb = block.norbi;
        const int jnorb = block.norbj;
        const int offset = block.offset;

        // (ti, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact; ++t) {
              for (int i = 0; i != nclosed; ++i) {
                double value = 0.0;
                if (t == jstart+v) {
                  value += 2.0 * cfock->element(i, nclosed+istart+u) + cfkd(i, istart+u);
                }
                if (t == istart+u) {
                  value -= 2.0 * cfock->element(i, nclosed+jstart+v) + cfkd(i, jstart+v);
                }
                value += 2.0 * rdm1(t, istart+u) * cfock->element(i, nclosed+jstart+v) - 2.0 * rdm1(t, jstart+v) * cfock->element(i, nclosed+istart+u);
                hessian->element(i+t*nclosed, aa_offset+offset+u+v*inorb) = value;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed) = value;
              }
            }
          }
        }
        // (at, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact; ++t) {
              for (int a = 0; a != nvirt; ++a) {
                double value = 2.0 * rdm1(t, jstart+v) * cfock->element(nocc+a, nclosed+istart+u)
                             - 2.0 * rdm1(t, istart+u) * cfock->element(nocc+a, nclosed+jstart+v);
                if (t == istart+u) {
                  value += cfkd(nocc+a, jstart+v);
                }
                if (t == jstart+v) {
                  value -= cfkd(nocc+a, istart+u);
                }
                hessian->element(va_offset+a+t*nvirt, aa_offset+offset+u+v*inorb) = value;
                hessian->element(aa_offset+offset+u+v*inorb, va_offset+a+t*nvirt) = value;
              }
            }
          }
        }
        // (tu, vw)
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norbi;
          const int jnorb2 = block2.norbj;
          const int offset2 = block2.offset;

          for (int w = 0; w != jnorb2; ++w) {
            for (int v = 0; v != inorb2; ++v) {
              for (int u = 0; u != jnorb; ++u) {
                for (int t = 0; t != inorb; ++t) {
                  double value = 2.0 * rdm1(jstart2+w, jstart+u) * cfock->element(nclosed+istart+t, nclosed+istart2+v)
                               - 2.0 * rdm1(istart2+v, jstart+u) * cfock->element(nclosed+istart+t, nclosed+jstart2+w)
                               - 2.0 * rdm1(jstart2+w, istart+t) * cfock->element(nclosed+jstart+u, nclosed+istart2+v)
                               + 2.0 * rdm1(istart2+v, istart+t) * cfock->element(nclosed+jstart+u, nclosed+jstart2+w);
                  if (jstart+u == istart2+v) {
                    value += cfkd(nclosed+istart+t, jstart2+w) + cfkd(nclosed+jstart2+w, istart+t);
                  }
                  if (jstart+u == jstart2+w) {
                    value -= cfkd(nclosed+istart+t, istart2+v) + cfkd(nclosed+istart2+v, istart+t);
                  }
                  if (istart+t == istart2+v) {
                    value -= cfkd(nclosed+jstart+u, jstart2+w) + cfkd(nclosed+jstart2+w, jstart+u);
                  }
                  if (istart+t == jstart2+w) {
                    value += cfkd(nclosed+jstart+u, istart2+v) + cfkd(nclosed+istart2+v, jstart+u);
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
      for (int u = 0; u != nact; ++u) {
        for (int i = 0; i != nclosed; ++i) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 0;
              for (int v = 0; v != nact; ++v) {
                assert(mo2e(nocc+a+(nclosed+v)*norb, nclosed+u+i*norb) - mo2e(nclosed+v+(nocc+a)*norb, nclosed+u+i*norb) < 1.0e-14);
                value += 2.0 * rdm1(t, v) * (4.0 * mo2e(nocc+a+(nclosed+v)*norb, nclosed+u+i*norb) - mo2e(nocc+a+(nclosed+u)*norb, nclosed+v+i*norb)
                                                                                                        - mo2e(nocc+a+i*norb, nclosed+u+(nclosed+v)*norb));
              }
              hessian->element(va_offset+a+t*nvirt, i+u*nclosed) += value;
              hessian->element(i+u*nclosed, va_offset+a+t*nvirt) += value;
            }
          }
        }
      }
      // (at, bi)
      for (int i = 0; i != nclosed; ++i) {
        for (int b = 0; b != nvirt; ++b) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 0;
              for (int u = 0; u != nact; ++u) {
                value += 2.0 * rdm1(t, u) * (4.0 * mo2e(nocc+a+(nclosed+u)*norb, nocc+b+i*norb) - mo2e(nocc+a+(nocc+b)*norb, nclosed+u+i*norb)
                                                                                                     - mo2e(nocc+a+i*norb, nclosed+u+(nocc+b)*norb));
              }
              hessian->element(va_offset+a+t*nvirt, vc_offset+b+i*nvirt) += value;
              hessian->element(vc_offset+b+i*nvirt, va_offset+a+t*nvirt) += value;
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact; ++u) {
        for (int j = 0; j != nclosed; ++j) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              double value = 16.0 * mo2e(nclosed+t+i*norb, nclosed+u+j*norb) - 4.0 * mo2e(nclosed+t+j*norb, i+(nclosed+u)*norb)
                                                                                 - 4.0 * mo2e(nclosed+t+(nclosed+u)*norb, i+j*norb);
              for (int v = 0; v != nact; ++v) {
                value -= 2.0 * rdm1(u, v) * (4.0 * mo2e(nclosed+t+i*norb, nclosed+v+j*norb) - mo2e(nclosed+t+j*norb, i+(nclosed+v)*norb)
                                                                                                - mo2e(nclosed+t+(nclosed+v)*norb, i+j*norb));
                value -= 2.0 * rdm1(t, v) * (4.0 * mo2e(nclosed+v+i*norb, nclosed+u+j*norb) - mo2e(nclosed+v+j*norb, nclosed+u+i*norb)
                                                                                                - mo2e(nclosed+v+(nclosed+u)*norb, i+j*norb));
              }
              hessian->element(i+t*nclosed, j+u*nclosed) += value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed; ++j) {
        for (int a = 0; a != nvirt; ++a) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              double value = 16.0 * mo2e(nclosed+t+i*norb, nocc+a+j*norb) - 4.0 * mo2e(nclosed+t+j*norb, nocc+a+i*norb)
                                                                              - 4.0 * mo2e(nclosed+t+(nocc+a)*norb, i+j*norb);
              for (int u = 0; u != nact; ++u) {
                value -= 2.0 * rdm1(t, u) * (4.0 * mo2e(nclosed+u+i*norb, nocc+a+j*norb) - mo2e(nclosed+u+j*norb, nocc+a+i*norb)
                                                                                             - mo2e(nclosed+u+(nocc+a)*norb, i+j*norb));
              }
              hessian->element(i+t*nclosed, vc_offset+a+j*nvirt) += value;
              hessian->element(vc_offset+a+j*nvirt, i+t*nclosed) += value;
            }
          }
        }
      }
      // (ai, bj)
      for (int j = 0; j != nclosed; ++j) {
        for (int b = 0; b != nvirt; ++b) {
          for (int i = 0; i != nclosed; ++i) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 16.0 * mo2e(nocc+a+i*norb, nocc+b+j*norb) - 4.0 * mo2e(nocc+a+(nocc+b)*norb, i+j*norb)
                                                                           - 4.0 * mo2e(nocc+a+j*norb, nocc+b+i*norb);
              hessian->element(vc_offset+a+i*nvirt, vc_offset+b+j*nvirt) += value;
            }
          }
        }
      }
#ifdef AAROT
      for (auto& block : act_rotblocks_) {
        const int istart = block.iorbstart;
        const int jstart = block.jorbstart;
        const int inorb = block.norbi;
        const int jnorb = block.norbj;
        const int offset = block.offset;

        // (ti, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact; ++t) {
              for (int i = 0; i != nclosed; ++i) {
                double value = 0.0;
                for (int w = 0; w != nact; ++w) {
                  value += rdm1(w, jstart+v) * (8.0 * mo2e(nclosed+istart+u+(nclosed+w)*norb, nclosed+t+i*norb)
                                              - 2.0 * mo2e(nclosed+istart+u+(nclosed+t)*norb, nclosed+w+i*norb)
                                              - 2.0 * mo2e(nclosed+istart+u+i*norb, nclosed+w+(nclosed+t)*norb));
                  value -= rdm1(w, istart+u) * (8.0 * mo2e(nclosed+jstart+v+(nclosed+w)*norb, nclosed+t+i*norb)
                                              - 2.0 * mo2e(nclosed+jstart+v+(nclosed+t)*norb, nclosed+w+i*norb)
                                              - 2.0 * mo2e(nclosed+jstart+v+i*norb, nclosed+w+(nclosed+t)*norb));
                }
                if (t == jstart+v) {
                  for (int x = 0; x != nact; ++x) {
                    for (int y = 0; y != nact; ++y) {
                      value += rdm1(x, y) * (2.0 * mo2e(nclosed+x+(nclosed+y)*norb, nclosed+istart+u+i*norb)
                                                 - mo2e(nclosed+istart+u+(nclosed+y)*norb, nclosed+x+i*norb));
                    }
                  }
                }
                if (t == istart+u) {
                  for (int x = 0; x != nact; ++x) {
                    for (int y = 0; y != nact; ++y) {
                      value -= rdm1(x, y) * (2.0 * mo2e(nclosed+x+(nclosed+y)*norb, nclosed+jstart+v+i*norb)
                                                 - mo2e(nclosed+jstart+v+(nclosed+y)*norb, nclosed+x+i*norb));
                    }
                  }
                }
                hessian->element(i+t*nclosed, aa_offset+offset+u+v*inorb) += value;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed) += value;
              }
            }
          }
        }
        // (ai, tu)
        for (int u = 0; u != jnorb; ++u) {
          for (int t = 0; t != inorb; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              for (int a = 0; a != nvirt; ++a) {
                double value = 0.0;
                for (int v = 0; v != nact; ++v) {
                  value += rdm1(jstart+u, v) * (8.0 * mo2e(nclosed+istart+t+(nclosed+v)*norb, nocc+a+i*norb)
                                              - 2.0 * mo2e(nclosed+istart+t+(nocc+a)*norb, nclosed+v+i*norb)
                                              - 2.0 * mo2e(nclosed+istart+t+i*norb, nclosed+v+(nocc+a)*norb));
                  value -= rdm1(istart+t, v) * (8.0 * mo2e(nclosed+jstart+u+(nclosed+v)*norb, nocc+a+i*norb)
                                              - 2.0 * mo2e(nclosed+jstart+u+(nocc+a)*norb, nclosed+v+i*norb)
                                              - 2.0 * mo2e(nclosed+jstart+u+i*norb, nclosed+v+(nocc+a)*norb));
                }
                hessian->element(vc_offset+a+i*nvirt, aa_offset+offset+t+u*inorb) += value;
                hessian->element(aa_offset+offset+t+u*inorb, vc_offset+a+i*nvirt) += value;
              }
            }
          }
        }
      }
#endif
    } // end of 2-electron integral part

    // Qvec part
    {
      // (at, bu)
      for (int u = 0; u != nact; ++u) {
        for (int b = 0; b != nvirt; ++b) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 0.0;
              for (int v = 0; v != nact; ++v) {
                for (int w = 0; w != nact; ++w) {
                  value += 2.0 * rdm2(t, u, v, w) * mo2e(nocc+a+(nocc+b)*norb, nclosed+v+(nclosed+w)*norb)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nocc+b+(nclosed+v)*norb, nocc+a+(nclosed+w)*norb);
                  if (a == b) {
                    for (int x = 0; x != nact; ++x) {
                      value -= rdm2(t, v, w, x) * mo2e(nclosed+u+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb)
                             + rdm2(u, v, w, x) * mo2e(nclosed+t+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb);
                    }
                  }
                }
              }
              hessian->element(va_offset+a+t*nvirt, va_offset+b+u*nvirt) += value;
            }
          }
        }
      }
      // (at, ui)
      for (int u = 0; u != nact; ++u) {
        for (int i = 0; i != nclosed; ++i) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              double value = 0.0;
              for (int v = 0; v != nact; ++v) {
                for (int w = 0; w != nact; ++w) {
                  value -= 2.0 * rdm2(t, u, v, w) * mo2e(nocc+a+i*norb, nclosed+v+(nclosed+w)*norb)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nclosed+v+i*norb, nclosed+w+(nocc+a)*norb);
                }
              }
              hessian->element(va_offset+a+t*nvirt, i+u*nclosed) += value;
              hessian->element(i+u*nclosed, va_offset+a+t*nvirt) += value;
            }
          }
        }
      }
      // (at, bi)
      for (int i = 0; i != nclosed; ++i) {
        for (int b = 0; b != nclosed; ++b) {
          for (int t = 0; t != nact; ++t) {
            for (int a = 0; a != nvirt; ++a) {
              if (a == b) {
                double value = 0.0;
                for (int u = 0; u != nact; ++u) {
                  for (int v = 0; v != nact; ++v) {
                    for (int w = 0; w != nact; ++w) {
                      value -= rdm2(t, u, v, w) * mo2e(nclosed+u+i*norb, nclosed+v+(nclosed+w)*norb);
                    }
                  }
                }
                hessian->element(va_offset+a+t*nvirt, vc_offset+b+i*nvirt) += value;
                hessian->element(vc_offset+b+i*nvirt, va_offset+a+t*nvirt) += value;
              }
            }
          }
        }
      }
      // (ti, uj)
      for (int u = 0; u != nact; ++u) {
        for (int j = 0; j != nclosed; ++j) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              double value = 0.0;
              for (int v = 0; v != nact; ++v) {
                for (int w = 0; w != nact; ++w) {
                  value += 2.0 * rdm2(t, u, v, w) * mo2e(i+j*norb, nclosed+v+(nclosed+w)*norb)
                         + 2.0 * (rdm2(u, v, t, w) + rdm2(v, u, t, w)) * mo2e(nclosed+v+j*norb, nclosed+w+i*norb);
                  if (i == j) {
                    for (int x = 0; x != nact; ++x) {
                      value -= rdm2(u, v, w, x) * mo2e(nclosed+t+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb)
                             + rdm2(t, v, w, x) * mo2e(nclosed+u+(nclosed+v)*norb, nclosed+w+(nclosed+x)*norb);
                    }
                  }
                }
              }
              hessian->element(i+t*nclosed, j+u*nclosed) += value;
            }
          }
        }
      }
      // (ti, aj)
      for (int j = 0; j != nclosed; ++j) {
        for (int a = 0; a != nvirt; ++a) {
          for (int t = 0; t != nact; ++t) {
            for (int i = 0; i != nclosed; ++i) {
              double value = 0.0;
              if (i == j) {
                for (int u = 0; u != nact; ++u) {
                  for (int v = 0; v != nact; ++v) {
                    for (int w = 0; w != nact; ++w) {
                      value -= rdm2(t, u, v, w) * mo2e(nocc+a+(nclosed+u)*norb, nclosed+v+(nclosed+w)*norb);
                    }
                  }
                }
              }
              hessian->element(i+t*nclosed, vc_offset+a+j*nvirt) += value;
              hessian->element(vc_offset+a+j*nvirt, i+t*nclosed) += value;
            }
          }
        }
      }
#ifdef AAROT
      for (auto& block : act_rotblocks_) {
        const int istart = block.iorbstart;
        const int jstart = block.jorbstart;
        const int inorb = block.norbi;
        const int jnorb = block.norbj;
        const int offset = block.offset;

        // (ti, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact; ++t) {
              for (int i = 0; i != nclosed; ++i) {
                double value = 0.0;
                if (t == jstart+v) {
                  for (int w = 0; w != nact; ++w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        value += rdm2(istart+u, w, x, y) * mo2e(nclosed+w+i*norb, nclosed+x+(nclosed+y)*norb);
                      }
                    }
                  }
                }
                if (t == istart+u) {
                  for (int w = 0; w != nact; ++w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        value -= rdm2(jstart+v, w, x, y) * mo2e(nclosed+w+i*norb, nclosed+x+(nclosed+y)*norb);
                      }
                    }
                  }
                }
                hessian->element(i+t*nclosed, aa_offset+offset+u+v*inorb) += value;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed) += value;

                // Q' and Q''
                double value2 = 0.0;
                for (int x = 0; x != nact; ++x) {
                  for (int y = 0; y != nact; ++y) {
                    value2 += 2.0 * (rdm2(t, istart+u, x, y) * mo2e(nclosed+jstart+v+i*norb, nclosed+x+(nclosed+y)*norb)
                                  + (rdm2(istart+u, x, t, y) + rdm2(istart+u, x, y, t)) * mo2e(nclosed+jstart+v+(nclosed+x)*norb, nclosed+y+i*norb));
                    value2 -= 2.0 * (rdm2(t, jstart+v, x, y) * mo2e(nclosed+istart+u+i*norb, nclosed+x+(nclosed+y)*norb)
                                  + (rdm2(jstart+v, x, t, y) + rdm2(jstart+v, x, y, t)) * mo2e(nclosed+istart+u+(nclosed+x)*norb, nclosed+y+i*norb));
                  }
                }
                hessian->element(i+t*nclosed, aa_offset+offset+u+v*inorb) += value2;
                hessian->element(aa_offset+offset+u+v*inorb, i+t*nclosed) += value2;
              }
            }
          }
        }
        // (at, uv)
        for (int v = 0; v != jnorb; ++v) {
          for (int u = 0; u != inorb; ++u) {
            for (int t = 0; t != nact; ++t) {
              for (int a = 0; a != nvirt; ++a) {
                double value = 0.0;
                if (t == istart+u) {
                  for (int w = 0; w != nact; ++w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        value += rdm2(jstart+v, w, x, y) * mo2e(nclosed+w+(nocc+a)*norb, nclosed+x+(nclosed+y)*norb);
                      }
                    }
                  }
                }
                if (t == jstart+v) {
                  for (int w = 0; w != nact; ++w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        value -= rdm2(istart+u, w, x, y) * mo2e(nclosed+w+(nocc+a)*norb, nclosed+x+(nclosed+y)*norb);
                      }
                    }
                  }
                }
                hessian->element(va_offset+a+t*nvirt, aa_offset+offset+u+v*inorb) += value;
                hessian->element(aa_offset+offset+u+v*inorb, va_offset+a+t*nvirt) += value;

                // Q' and Q''
                double value2 = 0.0;
                for (int x = 0; x != nact; ++x) {
                  for (int y = 0; y != nact; ++y) {
                    value2 += 2.0 * (rdm2(t, jstart+v, x, y) * mo2e(nclosed+istart+u+(nocc+a)*norb, nclosed+x+(nclosed+y)*norb)
                                  + (rdm2(jstart+v, x, t, y) + rdm2(jstart+v, x, y, t)) * mo2e(nclosed+istart+u+(nclosed+x)*norb, nclosed+y+(nocc+a)*norb));
                    value2 -= 2.0 * (rdm2(t, istart+u, x, y) * mo2e(nclosed+jstart+v+(nocc+a)*norb, nclosed+x+(nclosed+y)*norb)
                                  + (rdm2(istart+u, x, t, y) + rdm2(istart+u, x, y, t)) * mo2e(nclosed+jstart+v+(nclosed+x)*norb, nclosed+y+(nocc+a)*norb));
                  }
                }
                hessian->element(va_offset+a+t*nvirt, aa_offset+offset+u+v*inorb) += value2;
                hessian->element(aa_offset+offset+u+v*inorb, va_offset+a+t*nvirt) += value2;
              }
            }
          }
        }
        // (tu, vw)
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norbi;
          const int jnorb2 = block2.norbj;
          const int offset2 = block2.offset;
          
          for (int w = 0; w != jnorb2; ++w) {
            for (int v = 0; v != inorb2; ++v) {
              for (int u = 0; u != jnorb; ++u) {
                for (int t = 0; t != inorb; ++t) {
                  double value = 0.0;
                  if (jstart+u == istart2+v) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        for (int z = 0; z != nact; ++z) {
                          value += rdm2(jstart2+w, x, y, z) * mo2e(nclosed+istart+t+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb)
                                 + rdm2(istart+t, x, y, z) * mo2e(nclosed+jstart2+w+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb);
                        }
                      }
                    }
                  }
                  if (istart+t == jstart2+w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        for (int z = 0; z != nact; ++z) {
                          value += rdm2(jstart+u, x, y, z) * mo2e(nclosed+istart2+v+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb)
                                 + rdm2(istart2+v, x, y, z) * mo2e(nclosed+jstart+u+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb);
                        }
                      }
                    }
                  }
                  if (jstart+u == jstart2+w) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        for (int z = 0; z != nact; ++z) {
                          value -= rdm2(istart2+v, x, y, z) * mo2e(nclosed+istart+t+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb)
                                 + rdm2(istart+t, x, y, z) * mo2e(nclosed+istart2+v+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb);
                        }
                      }
                    }
                  }
                  if (istart+t == istart2+v) {
                    for (int x = 0; x != nact; ++x) {
                      for (int y = 0; y != nact; ++y) {
                        for (int z = 0; z != nact; ++z) {
                          value -= rdm2(jstart+u, x, y, z) * mo2e(nclosed+jstart2+w+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb)
                                 + rdm2(jstart2+w, x, y, z) * mo2e(nclosed+jstart+u+(nclosed+x)*norb, nclosed+y+(nclosed+z)*norb);
                        }
                      }
                    }
                  }
                  hessian->element(aa_offset+offset+t+u*inorb, aa_offset+offset2+v+w*inorb2) += value;

                  // Q' and Q''
                  double value2 = 0.0;
                  for (int x = 0; x != nact; ++x) {
                    for (int y = 0; y != nact; ++y) {
                      value2 += rdm2(jstart2+w, jstart+u, x, y) * mo2e(nclosed+istart2+v+(nclosed+istart+t)*norb, nclosed+x+(nclosed+y)*norb)
                             + (rdm2(jstart2+w, x, jstart+u, y) + rdm2(x, jstart2+w, jstart+u, y)) * mo2e(nclosed+istart2+v+(nclosed+x)*norb, 
                                                                                                          nclosed+istart+t+(nclosed+y)*norb);
                      value2 -= rdm2(istart2+v, jstart+u, x, y) * mo2e(nclosed+jstart2+w+(nclosed+istart+t)*norb, nclosed+x+(nclosed+y)*norb)
                             + (rdm2(istart2+v, x, jstart+u, y) + rdm2(x, istart2+v, jstart+u, y)) * mo2e(nclosed+jstart2+w+(nclosed+x)*norb,
                                                                                                          nclosed+istart+t+(nclosed+y)*norb);
                      value2 -= rdm2(jstart2+w, istart+t, x, y) * mo2e(nclosed+istart2+v+(nclosed+jstart+u)*norb, nclosed+x+(nclosed+y)*norb)
                             + (rdm2(jstart2+w, x, istart+t, y) + rdm2(x, jstart2+w, istart+t, y)) * mo2e(nclosed+istart2+v+(nclosed+x)*norb,
                                                                                                          nclosed+jstart+u+(nclosed+y)*norb);
                      value2 += rdm2(istart2+v, istart+t, x, y) * mo2e(nclosed+jstart2+w+(nclosed+jstart+u)*norb, nclosed+x+(nclosed+y)*norb)
                             + (rdm2(istart2+v, x, istart+t, y) + rdm2(x, istart2+v, istart+t, y)) * mo2e(nclosed+jstart2+w+(nclosed+x)*norb,
                                                                                                          nclosed+jstart+u+(nclosed+y)*norb);
                    }
                  }
                  hessian->element(aa_offset+offset+t+u*inorb, aa_offset+offset2+v+w*inorb2) += 2.0 * value2;
                }
              }
            }
          }
        }
      }
#endif
    } // end of Qvec part

    // check for compute_gradient
    {
      cout << " * checking compute_grad" << endl;
      VectorB grad_vec(rotsize);
      copy_n(grad->data(), rotsize, grad_vec.data());
      VectorB grad_diff(rotsize);
      VectorB tmp(grad_vec - grad_check);
      for (int i = 0; i != rotsize; ++i) {
        if (fabs(grad_vec(i)) >= 1.0e-15) {
          grad_diff(i) = fabs(tmp(i) / grad_vec(i));
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
        if (fabs(denom_vec(i)) >= 1.0e-15) {
          denom_diff(i) = fabs(tmp(i) / denom_vec(i));
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
      auto hess_trial_rotfile = compute_hess_trial(rot, half, halfa, halfa_JJ, cfock, afock, qxr);
      VectorB hess_trial_vec(rotsize);
      copy_n(hess_trial_rotfile->data(), rotsize, hess_trial_vec.data());
      VectorB rot_vec(rotsize);
      assert(rot->size() == rotsize);
      copy_n(rot->data(), rotsize, rot_vec.data());
      VectorB hess_debug_vec(rotsize);
      dgemv_("N", hessian->ndim(), hessian->mdim(), 1.0, hessian->data(), rotsize, rot_vec.data(), 1, 0.0, hess_debug_vec.data(), 1);
      VectorB hess_t_diff(rotsize);
      VectorB tmp(hess_trial_vec - hess_debug_vec);
      for (int i = 0; i != rotsize; ++i) {
        if (fabs(hess_trial_vec(i)) >= 1.0e-15) {
          hess_t_diff(i) = fabs(tmp(i) / hess_trial_vec(i));
          if (hess_t_diff(i) > 1.0e-8) cout << i << " : " << hess_t_diff(i) << ", hess_t_vec = " << hess_trial_vec(i) << endl;
        }
      }
      cout << "diff hess_trial rms : " << hess_t_diff.rms() << endl;
    }
  }
#endif // end of DEBUG_HESS

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      
      shared_ptr<const ASD_DMRG_RotFile> sigma = compute_hess_trial(trot, half, halfa, halfa_JJ, cfock, afock, qxr);
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
      cout << endl << "    * Max iteration reached during the second-order optimization.  Convergence not reached! *   " << endl << endl;
    }
  
  } // end of macro iter
  // block diagonalize coeff_ in nclosed and nvirt
//  coeff_ = semi_canonical_orb();
  
  // TODO maybe one more ASD iteration
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_gradient(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto mref = asd_dmrg_->multisite()->sref();
  const int nclosed = mref->nclosed();
  const int nact = mref->nact();
  const int nocc = nclosed + nact;
  const int nvirt = mref->nvirt();

  auto grad = make_shared<ASD_DMRG_RotFile>(nclosed, nact, nvirt, naa_);
  shared_ptr<const RDM<1>> rdm1 = asd_dmrg_->rdm1_av();
  
  // closed-active section, closed runs first
  if (nclosed) {
    double* target = grad->ptr_ca();
    for (int t = 0; t != nact; ++t, target += nclosed) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(0, nclosed+t), nclosed, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(0, nclosed+t), nclosed, target);
      blas::ax_plus_y_n(-2.0, qxr->element_ptr(0, t), nclosed, target);
      for (int u = 0; u != nact; ++u)
        blas::ax_plus_y_n(-2.0*rdm1->element(u, t), cfock->element_ptr(0, nclosed+u), nclosed, target);
    }
  }
  // virtual-active section, virtual runs first
  {
    double* target = grad->ptr_va();
    for (int t = 0; t != nact; ++t, target += nvirt) {
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nocc, t), nvirt, target);
      for (int u = 0; u != nact; ++u)
        blas::ax_plus_y_n(2.0*rdm1->element(u, t), cfock->element_ptr(nocc, nclosed+u), nvirt, target);
    }
  }
  // virtual-closed asection, virtual runs firsgt
  if (nclosed){
    double* target = grad->ptr_vc();
    for (int i = 0; i != nclosed; ++i, target += nvirt) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(nocc, i), nvirt, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(nocc, i), nvirt, target);
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
      for (int v = 0; v != nact; ++v) {
        blas::ax_plus_y_n(2.0*rdm1->element(v, jstart+j), cfock->element_ptr(nclosed+istart, nclosed+v), inorb, target);
        blas::ax_plus_y_n(-2.0*cfock->element(nclosed+v, nclosed+jstart+j), rdm1->element_ptr(istart, v), inorb, target);
      }
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nclosed+istart, jstart+j), inorb, target);
      blas::ax_plus_y_n(-2.0, qxr->transpose()->element_ptr(istart, nclosed+jstart+j), inorb, target);
    }
  }
#endif

  return grad;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_Second::compute_denom(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> half_1j, shared_ptr<const DFHalfDist> halfa,
    shared_ptr<const DFHalfDist> halfa_JJ, shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock) const {
  auto mref = asd_dmrg_->multisite()->sref();
  const int nclosed = mref->nclosed();
  const int nact = mref->nact();
  const int nocc = nclosed + nact;
  const int nvirt = mref->nvirt();

  auto denom = make_shared<ASD_DMRG_RotFile>(nclosed, nact, nvirt, naa_);
  const MatView ccoeff = coeff_->slice(0, nclosed);
  const MatView acoeff = coeff_->slice(nclosed, nocc);
  const MatView vcoeff = coeff_->slice(nocc, nocc+nvirt);

  Matrix rdm1(nact, nact);
  copy_n(asd_dmrg_->rdm1_av()->data(), rdm1.size(), rdm1.data());

  // Fock related part
  const Matrix fcd = *cfock->get_submatrix(nclosed, nclosed, nact, nact) * rdm1;
  const Matrix fock = *cfock + *afock;
  {
    // closed-active
    for (int i = 0; i != nact; ++i) 
      for (int j = 0; j != nclosed; ++j)
        denom->ele_ca(j, i) = 4.0 * fock(i+nclosed, i+nclosed) - 4.0 * fock(j, j) - 2.0 * fcd(i, i) + 2.0 * (*cfock)(j, j) * rdm1(i, i);
  
    // virtual-active
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nvirt; ++j)
        denom->ele_va(j, i) = 2.0 * (*cfock)(j+nocc, j+nocc) * rdm1(i, i) - 2.0 * fcd(i, i);
  
    // virtual-closed
    for (int i = 0; i != nclosed; ++i)
      for (int j = 0; j != nvirt; ++j)
        denom->ele_vc(j, i) = 4.0 * fock(j+nocc, j+nocc) - 4.0 * fock(i, i);
  }

  const int nao = coeff_->ndim();

  // rdm-integral part  
  // [tt|pq] = \Gamma_{vw,tt}(vw|pq)
  shared_ptr<const DFFullDist> vaa = halfa_JJ->compute_second_transform(acoeff);
  const int nri = vaa->block(0)->asize();
  shared_ptr<const DFFullDist> vgaa = vaa->apply_2rdm(*asd_dmrg_->rdm2_av());
  {
    Matrix tmp_ao(nao, nao);
    for (int i = 0; i != nact; ++i) {
      dgemv_("T", nri, nao*nao, 1.0, mref->geom()->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact*i), 1, 0.0, tmp_ao.data(), 1);
      // tmp_ao.allreduce();
      Matrix tmp_virt = vcoeff % tmp_ao * vcoeff;
      blas::ax_plus_y_n(2.0, tmp_virt.diag().get(), nvirt, denom->ptr_va()+nvirt*i);
      if (nclosed) {
        Matrix tmp_clo = ccoeff % tmp_ao * ccoeff;
        blas::ax_plus_y_n(2.0, tmp_clo.diag().get(), nclosed, denom->ptr_ca()+nclosed*i);
      }
    }
  }
  
  // [t,t] = \Gamma_{vw,xt}(vw|xt)
  shared_ptr<const DFFullDist> vaa_exc = halfa->compute_second_transform(acoeff);
  shared_ptr<const Matrix> mo2e = vaa->form_4index(vaa_exc, 1.0);
  {
    for (int i = 0; i != nact; ++i) {
      const double e2 = -2.0 * blas::dot_product(mo2e->element_ptr(0, i*nact), nact*nact*nact, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,i));
      for (int j = 0; j != nvirt; ++j)
        denom->ele_va(j, i) += e2;
      for (int k = 0; k != nclosed; ++k)
        denom->ele_ca(k, i) += e2;
    }
  }

  // mixed rdm2
  Matrix rdmk(nact*nact, nact); // stores \Gamma_{k,i,j,i} + \Gamma_{k,i,i,j}
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        rdmk(k+nact*j, i) = asd_dmrg_->rdm2_av()->element(k, i, j, i) + asd_dmrg_->rdm2_av()->element(k, i, i, j);
  {
    shared_ptr<const DFFullDist> vav = halfa->compute_second_transform(vcoeff)->apply_J();
    denom->ax_plus_y_va(2.0, *(rdmk % *vav->form_4index_diagonal_part()).transpose());
    if (nclosed) {
      shared_ptr<const DFFullDist> vac = halfa->compute_second_transform(ccoeff)->apply_J();
      shared_ptr<const Matrix> mcaa = vac->form_4index_diagonal_part()->transpose();
      denom->ax_plus_y_ca(2.0, *mcaa * rdmk);
      shared_ptr<Matrix> mcaad = mcaa->copy();
      dgemm_("N", "N", nclosed*nact, nact, nact, -1.0, mcaa->data(), nclosed*nact, rdm1.data(), nact, 1.0, mcaad->data(), nclosed*nact);
      for (int i = 0; i != nact; ++i)
        blas::ax_plus_y_n(12.0, mcaad->element_ptr(0, i+nact*i), nclosed, denom->ptr_ca()+i*nclosed);
      
      Matrix tmp(nao, nao);
      shared_ptr<DFFullDist> vgaa = vaa->copy();
      vgaa = vgaa->transform_occ1(make_shared<Matrix>(rdm1));
      vgaa->ax_plus_y(-1.0, vaa);
      for (int i = 0; i != nact; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, mref->geom()->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact*i), 1, 0.0, tmp.data(), 1);
        // tmp.allreduce();
        Matrix tmp0 = ccoeff % tmp * ccoeff;
        blas::ax_plus_y_n(4.0, tmp0.diag().get(), nclosed, denom->ptr_ca()+nclosed*i);
      }
    }
  }

  // 4-index integral part
  {
    // virtual-closed
    if (nclosed) {
      auto vvc = half_1j->compute_second_transform(vcoeff)->form_4index_diagonal()->transpose();
      denom->ax_plus_y_vc(12.0, *vvc);

      shared_ptr<const DFFullDist> vgcc = half->compute_second_transform(ccoeff);
      const int nri = vgcc->block(0)->asize();
      Matrix tmp_ao(nao, nao);
      for (int i = 0; i != nclosed; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, mref->geom()->df()->block(0)->data(), nri, vgcc->block(0)->data()+nri*(i+nclosed*i), 1, 0.0, tmp_ao.data(), 1);
        // tmp_ao.allreduce();
        Matrix tmp_virt = vcoeff % tmp_ao * vcoeff;
        blas::ax_plus_y_n(-4.0, tmp_virt.diag().get(), nvirt, denom->ptr_vc()+nvirt*i);
      }
    }
  }

#ifdef AAROT
  // active-active part
  shared_ptr<const Matrix> maa = vaa_exc->apply_J()->form_4index_diagonal_part();
  Matrix mgaa(nact, nact);
  {
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        mgaa.element(j, i) = blas::dot_product(vaa_exc->block(0)->data()+nri*(j+nact*i), nri, vgaa->block(0)->data()+nri*(j+nact*i));
  }
  for (auto& block : act_rotblocks_) {
    const int istart = block.iorbstart;
    const int jstart = block.jorbstart;
    const int inorb = block.norb_i;
    const int jnorb = block.norb_j;
    const int offset = block.offset;
    
    // prepare for mixed rdm2
    btas::CRange<3> range(nact, nact, nact*nact);
    auto tmptensor = make_shared<btas::Tensor3<double>>(range, asd_dmrg_->rdm2_av()->storage());
    vector<double> buf(nact*nact);
    for (int i = 0; i != tmptensor->extent(2); ++i) {
      copy_n(&(*tmptensor)(0,0,i), buf.size(), buf.data());
      blas::transpose(buf.data(), nact, nact, &(*tmptensor)(0,0,i));
    }
    RDM<2> rdmmix = *asd_dmrg_->rdm2_av()->copy();
    blas::ax_plus_y_n(1.0, tmptensor->data(), tmptensor->size(), rdmmix.data());
    auto mat1 = make_shared<Matrix>(nact*nact, inorb*jnorb);
    auto mat2 = mat1->clone();
    
    for (int j = 0; j != jnorb; ++j) {
      // [t,t] = \Gamma_{vw,xt}(vw|xt)
      const double e2j = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact*(jstart+j)), nact*nact*nact, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,jstart+j));

      // Fock related part
      for (int i = 0; i != inorb; ++i) {
        const double e2i = -2.0 * blas::dot_product(mo2e->element_ptr(0, nact*(istart+i)), nact*nact*nact, asd_dmrg_->rdm2_av()->element_ptr(0,0,0,istart+i));

        denom->ele_aa_offset(i, inorb, j, offset) = 2.0 * (*cfock)(nclosed+istart+i, nclosed+istart+i) * rdm1(jstart+j, jstart+j)
                                                    + 2.0 * (*cfock)(nclosed+jstart+j, nclosed+jstart+j) * rdm1(istart+i, istart+i)
                                                    - 4.0 * (*cfock)(nclosed+istart+i, nclosed+jstart+j) * rdm1(istart+i, jstart+j)
                                                    - 2.0 * fcd(istart+i, istart+i) - 2.0 * fcd(jstart+j, jstart+j)
                                                    + e2j + e2i;
        
        // mixed rdm2
        for (int y = 0; y != nact; ++y) {
          copy_n(&rdmmix(0, istart+i, y, jstart+j), nact, mat1->element_ptr(y*nact, i+j*inorb));
          copy_n(mo2e->element_ptr((jstart+j)*nact, y+(istart+i)*nact), nact, mat2->element_ptr(y*nact, i+j*inorb));
        }
      }

      // [tt|pq] = \Gamma_{vw,tt}(vw|pq)
      Matrix tmp1_act(nact, nact);
      dgemv_("T", nri, nact*nact, 1.0, vaa_exc->block(0)->data(), nri, vgaa->block(0)->data()+nri*(jstart+j+nact*(jstart+j)), 1, 0.0, tmp1_act.data(), 1);
      dgemv_("T", nri, nact*nact, 1.0, vgaa->block(0)->data(), nri, vaa_exc->block(0)->data()+nri*(jstart+j+nact*(jstart+j)), 1, 1.0, tmp1_act.data(), 1);
      blas::ax_plus_y_n(2.0, tmp1_act.diag().get()+istart, inorb, denom->ptr_aa_offset(offset)+j*inorb);
      blas::ax_plus_y_n(-4.0, mgaa.data()+istart+nact*(jstart+j), inorb, denom->ptr_aa_offset(offset)+j*inorb);

      blas::ax_plus_y_n(2.0, ((rdmk % *maa) + (*maa % rdmk)).data()+istart+nact*(jstart+j), inorb, denom->ptr_aa_offset(offset)+j*inorb);

    } // end of looping over second active index

    // mixed rdm2
    const Matrix out = *mat1 % *mat2;
    blas::ax_plus_y_n(-4.0, out.diag().get(), inorb*jnorb, denom->ptr_aa_offset(offset));

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
  auto mref = asd_dmrg_->multisite()->sref();
  const int nclosed = mref->nclosed();
  const int nact = mref->nact();
  const int nocc = nclosed + nact;
  const int nvirt = mref->nvirt();

  shared_ptr<ASD_DMRG_RotFile> sigma = trot->clone();

  shared_ptr<const Matrix> va = trot->va_mat();
  shared_ptr<const Matrix> ca = nclosed ? trot->ca_mat() : nullptr;
  shared_ptr<const Matrix> vc = nclosed ? trot->vc_mat() : nullptr;

  const MatView ccoeff = coeff_->slice(0, nclosed);
  const MatView acoeff = coeff_->slice(nclosed, nocc);
  const MatView vcoeff = coeff_->slice(nocc, nocc+nvirt);

  Matrix rdm1(nact, nact);
  copy_n(asd_dmrg_->rdm1_av()->data(), nact*nact, rdm1.data());

  // lambda for computing g(D)
  auto compute_gd = [&, this] (shared_ptr<const DFHalfDist> halft, shared_ptr<const DFHalfDist> halfjj, const MatView pcoeff) {
    shared_ptr<const Matrix> pcoefft = make_shared<Matrix>(pcoeff)->transpose();
    shared_ptr<Matrix> gd = mref->geom()->df()->compute_Jop(halft, pcoefft);
    shared_ptr<Matrix> ex0 = halfjj->form_2index(halft, 1.0);
    ex0->symmetrize();
    gd->ax_plus_y(-0.5, ex0);
    return gd;
  };

  // g(t_vc) operator and g(t_ca) operator
  if (nclosed) {
    const Matrix tcoeff2c = vcoeff * *vc + acoeff * *ca->transpose();
    auto halft = mref->geom()->df()->compute_half_transform(tcoeff2c);
    const Matrix gt = *compute_gd(halft, half, ccoeff);
    sigma->ax_plus_y_va(16.0, vcoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_ca(32.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_ca(-16.0, ccoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_vc(32.0, vcoeff % gt * ccoeff);
  }

  const Matrix tcoeff2a = nclosed ? (vcoeff * *va - ccoeff * *ca) : (vcoeff *  *va);
  shared_ptr<const DFHalfDist> halfta = mref->geom()->df()->compute_half_transform(tcoeff2a);
  
  // g(t_va - t_ca)
  if (nclosed) {
    shared_ptr<DFHalfDist> halftad = halfta->copy();
    halftad = halftad->transform_occ(make_shared<Matrix>(rdm1));
    const Matrix gt = *compute_gd(halftad, halfa_JJ, acoeff);
    sigma->ax_plus_y_ca(16.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_vc(16.0, vcoeff % gt * ccoeff);
  }

#ifdef AAROT
  // active-active 2-electron integral part
  if (nclosed){
    for (auto& block : act_rotblocks_) {
      const int istart = block.iorbstart;
      const int jstart = block.jorbstart;
      const int inorb = block.norb_i;
      const int jnorb = block.norb_j;
      const int offset = block.offset;
      const int bsize = block.size;

      auto rotblock_aa = make_shared<Matrix>(inorb, jnorb);
      copy_n(trot->ptr_aa_offset(offset), bsize, rotblock_aa->data());

      const MatView acoeffi = coeff_->slice(nclosed+istart, nclosed+istart+inorb);
      const MatView acoeffj = coeff_->slice(nclosed+jstart, nclosed+jstart+jnorb);

      const MatView rdmxi = rdm1.slice(istart, istart+inorb);
      const MatView rdmxj = rdm1.slice(jstart, jstart+jnorb);

      const MatView cai = ca->slice(istart, istart+inorb);
      const MatView caj = ca->slice(jstart, jstart+jnorb);

      { // ai->tu
        const Matrix tcoeffv2c = vcoeff * *vc;
        auto halftv2c = mref->geom()->df()->compute_half_transform(tcoeffv2c);
        const Matrix gtv2c = *compute_gd(halftv2c, half, ccoeff);
        sigma->ax_plus_y_aa_offset(16.0, acoeffi % gtv2c * acoeff * rdmxj , offset);
        sigma->ax_plus_y_aa_offset(-16.0, rdmxi % (acoeff % gtv2c * acoeffj), offset);
      }

      shared_ptr<DFHalfDist> halfxy = halfa->copy();
      halfxy = halfxy->transform_occ(make_shared<Matrix>(rdm1));
      const Matrix gtaa = *compute_gd(halfxy, halfa_JJ, acoeff);

      { // ti->uv
        const Matrix tcoeffc2a = ccoeff * *ca;
        auto halftc2a = mref->geom()->df()->compute_half_transform(tcoeffc2a);
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
        auto halftix = mref->geom()->df()->compute_half_transform(tcoeffi2x);
        const Matrix gt1 = *compute_gd(halftix, halfa_JJ, acoeff);
        sigma->ax_plus_y_vc(16.0, vcoeff % gt1 * ccoeff);
        sigma->ax_plus_y_ca(16.0, ccoeff % gt1 * acoeff);
        // \gamma_{tu} and \gamma_{uw} part
        const Matrix tcoeffj2x = acoeffj ^ (rdmxi * *rotblock_aa);
        auto halftjx = mref->geom()->df()->compute_half_transform(tcoeffj2x);
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
    shared_ptr<const Matrix> qaa = qxr->cut(nclosed, nocc);
    shared_ptr<const Matrix> qva = qxr->cut(nocc, nocc+nvirt);
    shared_ptr<const Matrix> qca = qxr->cut(0, nclosed);
    sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
    sigma->ax_plus_y_va(-2.0, *va * *qaa);
    if (nclosed) {
      sigma->ax_plus_y_va(-2.0, *vc * *qca);
      sigma->ax_plus_y_ca(-2.0, *ca ^ *qaa);
      sigma->ax_plus_y_ca(-2.0, *ca * *qaa);
      sigma->ax_plus_y_ca(-2.0, *vc % *qva);
      sigma->ax_plus_y_vc(-2.0, *va ^ *qca);
      sigma->ax_plus_y_vc(-2.0, *qva ^ *ca);
    }

#ifdef AAROT
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
      const MatView qvai = qva->slice(istart, istart+inorb);
      const MatView vai = va->slice(istart, istart+inorb);
      const MatView vaj = va->slice(jstart, jstart+jnorb);
      sigma->ax_plus_y_va_offset(2.0, qvaj ^ *rotblock_aa, istart);
      sigma->ax_plus_y_va_offset(-2.0, qvai * *rotblock_aa, jstart);
      sigma->ax_plus_y_aa_offset(2.0, vai % qvaj, offset);
      sigma->ax_plus_y_aa_offset(-2.0, qvai % vaj, offset);

      if (nclosed) {
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
          auto tmpmat = make_shared<const Matrix>((*qaaiI * *rotblock2_aa) + (*qaaiI % *rotblock2_aa));
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
    if (nclosed)
      sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp + *qpp));

#ifdef AAROT
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

      const MatView acoeffi = coeff_->slice(nclosed+istart, nclosed+istart+inorb);
      const MatView acoeffj = coeff_->slice(nclosed+jstart, nclosed+jstart+jnorb);

      auto filti = make_shared<Matrix>(nact, nact);
      for (int i = istart; i != istart+inorb; ++i)
        filti->element(i, i) = 1.0;
      auto filtj = make_shared<Matrix>(nact, nact);
      for (int j = jstart; j != jstart+jnorb; ++j)
        filtj->element(j, j) = 1.0;
      Matrix rotmat_aa(nact, nact);
      rotmat_aa.copy_block(istart, jstart, inorb, jnorb, rotblock_aa->data());
      
      shared_ptr<const DFHalfDist> half1;
      shared_ptr<const DFFullDist> full1;
      shared_ptr<const DFFullDist> full2D;
      // uv->at and uv->it
      {
        half1 = halfa->transform_occ(make_shared<Matrix>((*filti) * rotmat_aa));
        full2D = fullaaD->transform_occ1(filtj);
        shared_ptr<const Matrix> qpp1 = half1->form_2index(full2D, 1.0);

        half1 = halfa->transform_occ(make_shared<Matrix>((*filtj) ^ rotmat_aa));
        full2D = fullaaD->transform_occ1(filti);
        shared_ptr<const Matrix> qpp2 = half1->form_2index(full2D, 1.0);

        half1 = halfa;
        auto tmpfull2 = fullaa->transform_occ1(make_shared<Matrix>(*filti * rotmat_aa));
        auto tmpfull2s = tmpfull2->swap();
        tmpfull2->ax_plus_y(1.0, tmpfull2s);
        full2D = tmpfull2->apply_2rdm(*asd_dmrg_->rdm2_av())->swap();
        shared_ptr<const Matrix> qp1 = half1->form_2index(full2D, 1.0);

        half1 = halfa;
        tmpfull2 = fullaa->transform_occ1(make_shared<Matrix>((*filtj) ^ rotmat_aa));
        tmpfull2s = tmpfull2->swap();
        tmpfull2->ax_plus_y(1.0, tmpfull2s);
        full2D = tmpfull2->apply_2rdm(*asd_dmrg_->rdm2_av())->swap();
        shared_ptr<const Matrix> qp2 = half1->form_2index(full2D, 1.0);

        sigma->ax_plus_y_va(4.0, vcoeff % (*qp1 - *qp2 + *qpp1 - *qpp2));
        if (nclosed)
          sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp1 - *qp2 + *qpp1 - *qpp2));
      }

      // at->uv and it->uv
      {
        Matrix Qp(inorb, jnorb);

        full1 = halfta->compute_second_transform(acoeff);
        full2D = fullaaD;
        shared_ptr<const Matrix> qpp = full1->form_2index(full2D, 1.0);
        Qp += *qpp->get_submatrix(istart, jstart, inorb, jnorb) - *(qpp->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

        full1 = fullaa;
        full2D = fullta->apply_2rdm(*asd_dmrg_->rdm2_av())->swap();
        shared_ptr<const Matrix> qp = full1->form_2index(full2D, 1.0);
        Qp += *qp->get_submatrix(istart, jstart, inorb, jnorb) - *(qp->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

        sigma->ax_plus_y_aa_offset(4.0, Qp, offset);
      }

      // (tu, vw) part
      {
        for (auto& block2 : act_rotblocks_) {
          const int istart2 = block2.iorbstart;
          const int jstart2 = block2.jorbstart;
          const int inorb2 = block2.norb_i;
          const int jnorb2 = block2.norb_j;
          const int offset2 = block2.offset;
          const int bsize2 = block2.size;
        
          Matrix out(inorb, jnorb);
          
          auto rotblock2_aa = make_shared<Matrix>(inorb2, jnorb2);
          copy_n(trot->ptr_aa_offset(offset2), bsize2, rotblock2_aa->data());

          auto filti2 = filti->clone();
          for (int i = istart2; i != istart2+inorb2; ++i)
            filti2->element(i, i) = 1.0;
          auto filtj2 = filtj->clone();
          for (int j = jstart2; j != jstart2+jnorb2; ++j)
            filtj2->element(j, j) = 1.0;
          Matrix rotmat2_aa(nact, nact);
          rotmat2_aa.copy_block(istart2, jstart2, inorb2, jnorb2, rotblock2_aa->data());

          shared_ptr<const DFFullDist> fullaa_0J = halfa->compute_second_transform(acoeff);
          shared_ptr<const DFFullDist> fulltaj2_0J = fullaa_0J->transform_occ1(make_shared<Matrix>((*filti2) * rotmat2_aa));
          shared_ptr<const DFFullDist> fulltai2_0J = fullaa_0J->transform_occ1(make_shared<Matrix>((*filtj2) ^ rotmat2_aa));

          full1 = fulltaj2_0J;
          full2D = fullaaD->transform_occ1(filtj2);
          shared_ptr<const Matrix> Qpp1 = full1->form_2index(full2D, 1.0);
          out += *Qpp1->get_submatrix(istart, jstart, inorb, jnorb) - *(Qpp1->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

          full1 = fulltai2_0J;
          full2D = fullaaD->transform_occ1(filti2);
          shared_ptr<const Matrix> Qpp2 = full2D->form_2index(full1, 1.0);
          out += *Qpp2->get_submatrix(istart, jstart, inorb, jnorb) - *(Qpp2->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

          full1 = fullaa;
          shared_ptr<DFFullDist> tmpfull2 = fulltaj2_0J->copy();
          auto tmpfull2s = tmpfull2->swap();
          tmpfull2->ax_plus_y(1.0, tmpfull2s);
          full2D = tmpfull2->apply_2rdm(*asd_dmrg_->rdm2_av())->swap();
          shared_ptr<const Matrix> Qp1 = full1->form_2index(full2D, 1.0);
          out += *Qp1->get_submatrix(istart, jstart, inorb, jnorb) - *(Qp1->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

          full1 = fullaa;
          tmpfull2 = fulltai2_0J->copy();
          tmpfull2s = tmpfull2->swap();
          tmpfull2->ax_plus_y(1.0, tmpfull2s);
          full2D = tmpfull2->apply_2rdm(*asd_dmrg_->rdm2_av())->swap();
          shared_ptr<const Matrix> Qp2 = full2D->form_2index(full1, 1.0);
          out += *Qp2->get_submatrix(istart, jstart, inorb, jnorb) - *(Qp2->get_submatrix(jstart, istart, jnorb, inorb)->transpose());

          sigma->ax_plus_y_aa_offset(4.0, out, offset);
        }
      }
    } // end of looping over block
#endif
  } // end of Q' and Q'' part

  // Fock related terms
  {
    // construct submatrices
    shared_ptr<const Matrix> fcaa = cfock->get_submatrix(nclosed, nclosed, nact, nact);
    shared_ptr<const Matrix> faaa = afock->get_submatrix(nclosed, nclosed, nact, nact);
    shared_ptr<const Matrix> fcva = cfock->get_submatrix(nocc, nclosed, nvirt, nact);
    shared_ptr<const Matrix> fava = afock->get_submatrix(nocc, nclosed, nvirt, nact);
    shared_ptr<const Matrix> fcvv = cfock->get_submatrix(nocc, nocc, nvirt, nvirt);
    shared_ptr<const Matrix> favv = afock->get_submatrix(nocc, nocc, nvirt, nvirt);
    shared_ptr<const Matrix> fccc = nclosed ? cfock->get_submatrix(0, 0, nclosed, nclosed) : nullptr;
    shared_ptr<const Matrix> facc = nclosed ? afock->get_submatrix(0, 0, nclosed, nclosed) : nullptr;
    shared_ptr<const Matrix> fcca = nclosed ? cfock->get_submatrix(0, nclosed, nclosed, nact) : nullptr;
    shared_ptr<const Matrix> faca = nclosed ? afock->get_submatrix(0, nclosed, nclosed, nact) : nullptr;
    shared_ptr<const Matrix> fcvc = nclosed ? cfock->get_submatrix(nocc, 0, nvirt, nclosed) : nullptr;
    shared_ptr<const Matrix> favc = nclosed ? afock->get_submatrix(nocc, 0, nvirt, nclosed) : nullptr;

    sigma->ax_plus_y_va( 4.0, *fcvv * *va * rdm1);
    sigma->ax_plus_y_va(-2.0, *va * (rdm1 * *fcaa + *fcaa * rdm1));
    if (nclosed) {
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
  
      if (nclosed) { // (ti, uv)
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
  auto mref = asd_dmrg_->multisite()->sref();
  const int nclosed = mref->nclosed();
  const int nact = mref->nact();

  auto trans = make_shared<Matrix>(nact, nact);
  trans->add_diag(2.0);
  blas::ax_plus_y_n(-1.0, asd_dmrg_->rdm1_av()->data(), nact*nact, trans->data());

  VectorB occup(nact);
  trans->diagonalize(occup);

  asd_dmrg_->rotate_rdms(trans);

  auto cnew = make_shared<Coeff>(*coeff_);
  cnew->copy_block(0, nclosed, cnew->ndim(), nact, coeff_->slice(nclosed, nclosed+nact) * *trans);
  coeff_ = cnew;
}


