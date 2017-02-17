//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/dimer/dimer_metal.cc
// Copyright (C) 2015 Toru Shiozaki
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

#include<src/asd/dimer/dimer.h>

using namespace std;
using namespace bagel;

void Dimer::set_active_metal(shared_ptr<const PTree> idata) {
  auto actinfo = idata->get_child_optional("dimer_active");
  set<int> actset;
  if (actinfo) for (auto& i : *actinfo) { actset.insert(lexical_cast<int>(i->data())-1); }
  if (actset.empty())
    throw runtime_error("Active space of the dimer MUST be specified in dimer_active!");

  set<int> Alist, Blist, Llist;
 
  // Now dealing with localized orbitals, assign them to A and B fragments, briding active orbitals are assigned to Llist.
  {
    vector<pair<int, int>> bounds;
    vector<int> sizes = idata->get_vector<int>("region_sizes"); // [A, B, link]
    int nbasis = 0;
    int natoms = 0;
    // set bounds : first & last basis function indices for the given region
    for (int& region : sizes) {
      const int atomstart = natoms;
      const int basisstart = nbasis;
      for (int atom = atomstart; atom != atomstart + region; ++atom)
        nbasis += sgeom_->atoms()[atom]->nbasis();
      natoms += region;
      if (basisstart != nbasis)
        bounds.emplace_back(basisstart, nbasis);
    }
    if (natoms != count_if(sgeom_->atoms().begin(), sgeom_->atoms().end(), [](const shared_ptr<const Atom> a){ return !a->dummy();}))
      throw logic_error("All atoms must be assigned to regions");
    
    cout << "  o Assigning localized dimer active orbitals to fragment A and B" << endl;
    shared_ptr<Matrix> coeff = isolated_refs_.first->coeff()->copy();
    for (auto& amo : actset) {
      const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, amo), bounds[0].second - bounds[0].first, coeff->element_ptr(bounds[0].first, amo));
      const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, amo), bounds[1].second - bounds[1].first, coeff->element_ptr(bounds[1].first, amo));
      if (bounds.size() == 2) { // [A, B, 0]only two fragments, no linking atoms
        if (sum_A > sum_B && abs(sum_A - sum_B) > region_thresh_) {
          cout << "    - active orbital(" << amo + 1 << ") is assigned to monomer A." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          Alist.insert(amo);
        }
        if (sum_A < sum_B && abs(sum_A - sum_B) > region_thresh_) {
          cout << "    - active orbital(" << amo + 1 << ") is assigned to monomer B." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          Blist.insert(amo);
        } else {
          cout << "    - active orbital(" << amo + 1 << ") is assigned to linked active orbitals." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          Llist.insert(amo);
        }
      } else { // have [A, B, link]
        // TODO deal with such system 
      }
    }
    cout << "    - orbitals are assigned as : " << Alist.size() << "(A), " << Blist.size() << "(B) and " << Llist.size() << " bridging active orbitals." << endl;

    // Assigning bridging active orbitals to fragments A and B by projection
    if (!Llist.empty()) {
      const int dimerbasis = sgeom_->nbasis();
      const int nLact = Llist.size();
      auto Lcoeff = make_shared<Matrix>(dimerbasis, nLact);
      int iter = 0;
      for (auto& lmo : Llist) {
        Lcoeff->copy_block(0, iter, dimerbasis, 1, isolated_refs_.first->coeff()->slice(lmo, lmo+1));
        ++iter;
      }
    
      {
        Overlap S(sgeom_);
        const int abasis = bounds[0].second - bounds[0].first;
        const int bbasis = bounds[1].second - bounds[1].first;
        auto SA = make_shared<Matrix>(abasis, abasis);
        auto SB = make_shared<Matrix>(bbasis, bbasis);
        auto SAmix = make_shared<Matrix>(abasis, dimerbasis);
        auto SBmix = make_shared<Matrix>(bbasis, dimerbasis);
        SA->copy_block(0, 0, abasis, abasis, S.get_submatrix(0, 0, abasis, abasis));
        SB->copy_block(0, 0, bbasis, bbasis, S.get_submatrix(abasis, abasis, bbasis, bbasis));
        SAmix->copy_block(0, 0, abasis, dimerbasis, S.get_submatrix(0, 0, abasis, dimerbasis));
        SBmix->copy_block(0, 0, bbasis, dimerbasis, S.get_submatrix(abasis, 0, bbasis, dimerbasis));

        auto SA_inv = make_shared<Matrix>(*SA);
        SA_inv->inverse_symmetric();
        auto SB_inv = make_shared<Matrix>(*SB);
        SA_inv->inverse_symmetric();

        // Forming projected coeff to A and B
        auto projected_A = make_shared<const Matrix>(*SA_inv * *SAmix * *Lcoeff);
        auto projected_B = make_shared<const Matrix>(*SB_inv * *SBmix * *Lcoeff);

        // Get rid of redundancy in projected coeffs
        shared_ptr<Matrix> reduced_MO_A;
        shared_ptr<Matrix> reduced_MO_B;
        {
          auto CSC = make_shared<Matrix>(*projected_A % *SA * *projected_A);
          VectorB eig(projected_A->mdim());
          CSC->diagonalize(eig);
          auto P = CSC->get_submatrix(0, CSC->mdim()/2, CSC->ndim(), CSC->mdim()/2);
          auto reduced_A = make_shared<Matrix>(*projected_A * *P);
          reduced_MO_A = make_shared<Matrix>(dimerbasis, reduced_A->mdim());
          reduced_MO_A->copy_block(0, 0, abasis, reduced_A->mdim(), reduced_A->data());
        }
        {
          auto CSC = make_shared<Matrix>(*projected_B % *SB * *projected_B);
          VectorB eig(projected_B->mdim());
          CSC->diagonalize(eig);
          auto P = CSC->get_submatrix(0, CSC->mdim()/2, CSC->ndim(), CSC->mdim()/2);
          auto reduced_B = make_shared<Matrix>(*projected_B * *P);
          reduced_MO_B = make_shared<Matrix>(dimerbasis, reduced_B->mdim());
          reduced_MO_B->copy_block(abasis, 0, bbasis, reduced_B->mdim(), reduced_B->data());
        }
        auto reduced_MO_AB = reduced_MO_A->merge(reduced_MO_B);
        int iter = 0;
        for (auto& imo : Llist) {
          coeff->copy_block(0, imo, dimerbasis, 1, reduced_MO_AB->slice(iter, 1));
          ++iter;
        }
      
        // lowdin orthogonalization
        auto tildex = make_shared<Matrix>(*coeff % S * *coeff);
        tildex->inverse_half();
        coeff = make_shared<Matrix>(*coeff * *tildex);
      }
#if 0 // Compare new orthogonal coeff with original coeff
    auto compare = make_shared<Matrix>(*isolated_refs_.first->coeff() % S * *coeff);
    cout << "Overlap between original and new coeff :" << endl;
    for (int i = 0; i != compare->ndim(); ++i)
      cout << "(" << i << ") = " << *compare->element_ptr(i,i) << endl;
#endif
      
      for (auto& amo : Llist) {
        const double sum_A = blas::dot_product(coeff->element_ptr(bounds[0].first, amo), bounds[0].second - bounds[0].first,coeff->element_ptr(bounds[0].first, amo));
        const double sum_B = blas::dot_product(coeff->element_ptr(bounds[1].first, amo), bounds[1].second - bounds[1].first,coeff->element_ptr(bounds[1].first, amo));
        cout << "amo : " << amo + 1 << endl;
        cout << "sumA : " << sum_A << endl;
        cout << "sumB : " << sum_B << endl;
        if (sum_A > sum_B && abs(sum_A - sum_B) > region_thresh_) {
          cout << "    - projected active orbital(" << amo + 1 << ") is assigned to monomer A." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          Alist.insert(amo);
        } else if (sum_A < sum_B && abs(sum_A - sum_B) > region_thresh_) {
          cout << "    - projected active orbital(" << amo + 1 << ") is assigned to monomer B." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          Blist.insert(amo);
        } else 
          throw runtime_error("Wrong choice of active orbitals. The projected orbital(" + to_string(amo+1) + ") does not belong to any monomer.");
      }
    } // end of !Llist.empty() 
  } // Now we have new orthonormal coeff as well as active orbital lists.

  // Make new Reference; active orbitals are placed after closed orbitals
  auto tmpref = isolated_refs_.first->set_active_metal(Alist, Blist);
  active_refs_ = {tmpref, tmpref};

  // Update dimer info
  const int nclosedA = active_refs_.first->nclosed();
  const int nclosedB = active_refs_.second->nclosed();
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nact = nactA + nactB;
  const int nactvirtA = isolated_refs_.first->nvirt() - active_refs_.first->nvirt();
  const int nactvirtB = isolated_refs_.second->nvirt() - active_refs_.second->nvirt();
  const int dimerbasis = sgeom_->nbasis();
  assert(dimerbasis == geoms_.first->nbasis());
  assert(dimerbasis == geoms_.second->nbasis());
  const int nclosed_HF = sref_->nclosed();
  const int nvirt_HF = sref_->nvirt();
  assert(dimerbasis == nclosed_HF + nvirt_HF);
  assert(sref_->nact() == 0);
  const int nclosed = nclosed_HF - (nclosed_HF - nclosedA) - (nclosed_HF - nclosedB);
cout << "nclosedA = " << nclosedA << endl;
cout << "nclosedB = " << nclosedB << endl;
cout << "nactA = " << nactA << endl;
cout << "nactB = " << nactB << endl;
cout << "nact = " << nact << endl;
cout << "nactvirtA = " << nactvirtA << endl;
cout << "nactvirtB = " << nactvirtB << endl;
cout << "dimerbasis = " << dimerbasis << endl;
cout << "nclosed_HF = " << nclosed_HF << endl;
cout << "nvirt_HF = " << nvirt_HF << endl;
cout << "nclosed = " << nclosed << endl;


}
