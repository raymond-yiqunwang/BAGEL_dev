//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Date: July 28, 2014
//


#ifndef __ECP_TESTCODE_SOINTEGRAL_SOH_H
#define __ECP_TESTCODE_SOINTEGRAL_SO_H

#include <src/math/bessel.h>
#include "src/integral/ecp/sphharmonics.h"
#include "cartgauss.h"

using namespace bagel;
using namespace std;

class Angular {
  protected:
    array<double, 3> CA_;
    array<int, 3> angA_;

  public:

    Angular(const array<double, 3> CA, const array<int, 3> angA)
      : CA_(CA), angA_(angA) {}

    ~Angular() {}

    double integrate3USP(const int i, const int j, const int k) const {
      return (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) ? (4.0 * pi__ * dfact(i-1) * dfact(j-1) * dfact(k-1) / dfact(i+j+k+1)) : 0.0;
    }

    double integrate2SH1USP(const array<int, 2> lm1, const array<int, 2> lm2, const array<int, 3> ijk) const {
      vector<pair<double, int>> usp;

      auto zlm1 = make_shared<SphHarmonics>(lm1);
      const int l1 = lm1[0];
      int cnt = 0;
      int n1 = 0;
      for (int lz = 0; lz <= l1; ++lz) {
        for (int ly = 0; ly <= l1-lz; ++ly) {
          const int lx = l1-lz-ly;
          const double coeff = zlm1->sph_to_USP(lx, ly);
          if (coeff != 0.0) {
            ++n1;
            pair<double, int> c_usp(coeff, cnt);
            usp.push_back(c_usp);
          }
          ++cnt;
        }
      }

      auto zlm2 = make_shared<SphHarmonics>(lm2);
      const int l2 = lm2[0];
      cnt = 0;
      int n2 = 0;
      for (int lz = 0; lz <= l2; ++lz) {
        for (int ly = 0; ly <= l2-lz; ++ly) {
          const int lx = l2-lz-ly;
          const double coeff = zlm2->sph_to_USP(lx, ly);
          if (coeff != 0.0) {
            ++n2;
            pair<double, int> c_usp(coeff, cnt);
            usp.push_back(c_usp);
          }
          ++cnt;
        }
      }

      assert (n1+n2-usp.size()==0);
      double ans = 0.0;
      for (int i = 0; i != n1; ++i) {
        for (int j = n1; j != n1 + n2; ++j) {
          const double coeff = usp[i].first * usp[j].first;
          array<int, 3> ki;
          int id = usp[i].second;
          int kz = 0;
          for (int lp1 = l1+1; lp1 != 0; --lp1) {
            if (id - lp1 < 0) {
              ki[2] = kz;
              ki[1] = id;
              ki[0] = l1-ki[2]-ki[1];
              break;
            } else {
              kz++;
              id -= lp1;
            }
          }
          array<int, 3> kj;
          id = usp[j].second;
          kz = 0;
          for (int lp1 = l2+1; lp1 != 0; --lp1) {
            if (id - lp1 < 0) {
              kj[2] = kz;
              kj[1] = id;
              kj[0] = l2-kj[2]-kj[1];
              break;
            } else {
              kz++;
              id -= lp1;
            }
          }
          const int x = ki[0] + kj[0] + ijk[0];
          const int y = ki[1] + kj[1] + ijk[1];
          const int z = ki[2] + kj[2] + ijk[2];
          ans += coeff*integrate3USP(x, y, z);
        }
      }

      return ans;

    }

    double compute_omega(const int a, const int b, const int c, const int ld, const int l, const int m) const {
      double out = 0.0;
      const array<int, 3> ijk = {{a, b, c}};
      const array<int, 2> lm = {{l, m}};
      for (int mu = 0; mu != 2*ld; ++mu) {
        auto shAB = make_shared<SphHarmonics>(ld, mu-ld, CA_);
        const array<int, 2> ldmu = {{ld, mu}};
        out += shAB->zlm() * integrate2SH1USP(ldmu, lm, ijk);
      }

      return out;
    }

    double compute(const int h, const int ld, const int l, const int m) const {
      const static Comb comb;

      const int nA = angA_[0];
      const int lA = angA_[1];
      const int mA = angA_[2];

      double out = 0.0;
      const int amin = max(0, h-lA-mA);
      const int amax = min(nA, h);
      for (int a = amin; a != amax; ++a) {
        const int bmin = max(0, h-a-mA);
        const int bmax = min(lA, h-a);
        for (int b = bmin; b != bmax; ++b) {
          const double omega = compute_omega(a, b, h-a-b, ld, l, m);
          out += comb(nA, a) * comb(lA, b) * comb(mA, h-a-b) * pow(CA_[0], nA-a) * pow(CA_[1], lA-b) * pow(CA_[2], mA+a+b-h) * omega;
        }
      }

      return out;
    }

};

class SOIntegral { /* ACB */
  protected:
    shared_ptr<CartesianGauss> cartA_, cartB_;
    shared_ptr<SphHarmonics> sphC_;

    array<double, 3> AB_, CA_, CB_;
    double dAB_, dCA_, dCB_;

  public:
    SOIntegral(shared_ptr<CartesianGauss> cA, shared_ptr<CartesianGauss> cB, shared_ptr<SphHarmonics> sC)
      : cartA_(cA), cartB_(cB), sphC_(sC) { init(); }

    ~SOIntegral() {}

    array<double, 3> compute_mlm(const int m1, const int m2, const int l) {

       array<double, 3> out = {{0.0, 0.0, 0.0}};
       const double tau1 = (m1 < 0) ? 0.0 : 1.0;
       const double tau2 = (m2 < 0) ? 0.0 : 1.0;
       const complex<double> th1 = (m1==0) ? complex<double>(0.5, 0.0) : complex<double>(tau1/sqrt(2.0), (1.0-tau1)/sqrt(2.0));
       const complex<double> th2 = (m2==0) ? complex<double>(0.5, 0.0) : complex<double>(tau2/sqrt(2.0), (1.0-tau2)/sqrt(2.0));

       if (m1==-m2) out[0] = 0.5*m2;
       complex<double> mu = conj(th1)*th2;
       if (m1*m2==0) mu *= 2.0;
       const double dp = (abs(m1)==abs(m2)+1) ? 1.0 : 0.0;
       const double dm = (abs(m1)==abs(m2)-1) ? 1.0 : 0.0;
       const complex<double> ab = 0.5*sqrt((l+m1*m1-abs(m1*m2))*(l+m2*m2-abs(m1*m2)))*mu*(dp - pow(-1.0, tau1+tau2)*dm);
       out[1] = real(ab);
       out[2] = real(ab);

       return out;
    }

    double compute_radial(const int N, const int l1, const int l2, const double k1, const double k2, const double a, const double r) const {
    /* Q_{l1,l2}^N(k1,k2,a) (MD E25) */

      const static MSphBesselI msbessel;
      const double b1 = msbessel.compute(l1, k1*r);
      const double b2 = msbessel.compute(l2, k2*r);

      return pow(r, N) * exp(-a*r*r) * b1 * b2;
    }

    array<double, 3> compute(const double r) {
      array<double, 3> out = {{0.0, 0.0, 0.0}};

      const double expA = cartA_->exponent();
      const double expB = cartB_->exponent();

      const int jA = cartA_->total_ang();
      const int jB = cartB_->total_ang();
      const int l = sphC_->angular_momentum(0);

      const double coeff = exp(-expA*expB*dAB_*dAB_/(expA+expB));
      auto pA = make_shared<const Angular>(CA_, cartA_->angular_momentum());
      auto pB = make_shared<const Angular>(CB_, cartB_->angular_momentum());

      for (int ld = max(0, l-jA); ld != l+jA; ++ld) {
        for (int mu = max(0, l-jB); mu != l+jB; ++mu) {
          const int gmin = abs(ld-l)+abs(mu-l);
          const int gmax = jA-(jA-abs(ld-l))%2 + jB-(jB-abs(mu-l))%2;
          for (int g = gmin; g != gmax; g += 2) {
            const double Qldmu = compute_radial(g, ld, mu, 2.0*expA*dCA_*r, 2.0*expB*dCB_*r, expA+expB, r);
            const int hmin = max(abs(ld-l), g-(jA-abs(ld-l))%2);
            const int hmax = min(jA-(jA-abs(ld-l))%2, g-abs(mu-l));
            for (int h = hmin; h != hmax; h +=2) {
              for (int m1 = 0; m1 != 2*l; ++m1) {
                for (int m2 = 0; m2 != m1-1; ++m2) {
                  const array<double, 3> fmm = compute_mlm(m1-l, m2-l, l);
                  for (int i = 0; i != 3; ++i)
                    out[i] += coeff * Qldmu * (pA->compute(h, ld, l, m1-l)*pB->compute(g-h, mu, l, m2-l) - pA->compute(h, ld, l, m2-l)*pB->compute(g-h, mu, l, m1-l)) * fmm[i];
                }
              }
            }
          }
        }
      }

      return out;
    }

    void init() {
      for (int i = 0; i != 3; ++i) {
        AB_[i] = cartA_->centre(i) - cartB_->centre(i);
        CA_[i] = sphC_->centre(i) - cartA_->centre(i);
        CB_[i] = sphC_->centre(i) - cartB_->centre(i);
      }
      dAB_ = sqrt(pow(AB_[0], 2) + pow(AB_[1], 2) + pow(AB_[2], 2));
      dCA_ = sqrt(pow(CA_[0], 2) + pow(CA_[1], 2) + pow(CA_[2], 2));
      dCB_ = sqrt(pow(CB_[0], 2) + pow(CB_[1], 2) + pow(CB_[2], 2));
    }

};

#endif
