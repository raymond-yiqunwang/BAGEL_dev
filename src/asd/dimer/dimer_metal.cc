//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dimer_metal.cc
// Copyright (C) 2017 Toru Shiozaki
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
#include<src/scf/hf/fock.h>

using namespace std;
using namespace bagel;

void Dimer::set_active_metal(shared_ptr<const PTree> idata) {
  auto actinfo = idata->get_child_optional("dimer_active");
  set<int> actset;
  if (actinfo) for (auto& i : *actinfo) { actset.insert(lexical_cast<int>(i->data())-1); }
  if (actset.empty())
    throw runtime_error("Active space of the dimer MUST be specified in dimer_active!");

  set<int> Alist, Blist, Llist;
 
  vector<pair<int, int>> bounds;
  vector<int> sizes = idata->get_vector<int>("region_sizes"); // [A, B, bridge]
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
    if (bounds.size() == 2) { // [A, B, 0]only two fragments, no bridging atoms
      if (sum_A > sum_B && fabs(sum_A - sum_B) > region_thresh_) {
        cout << "    - active orbital(" << amo + 1 << ") is assigned to monomer A." << endl;
        cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
        Alist.insert(amo);
      } else if (sum_A < sum_B && fabs(sum_A - sum_B) > region_thresh_) {
        cout << "    - active orbital(" << amo + 1 << ") is assigned to monomer B." << endl;
        cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
        Blist.insert(amo);
      } else {
        cout << "    - active orbital(" << amo + 1 << ") is assigned to linked active orbitals." << endl;
        cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
        Llist.insert(amo);
      }
    } else { // have [A, B, bridge]
      // TODO deal with such system 
    }
  }
  cout << "    - orbitals are assigned as : " << Alist.size() << "(A), " << Blist.size() << "(B) and " << Llist.size() << " (bridging) active orbitals." << endl;

  active_refs_ = {isolated_refs_.first->set_active_metal(Alist, Llist), isolated_refs_.second->set_active_metal(Blist, Llist)};

  // Update dimer info
  int nactA = active_refs_.first->nact();  
  const int nactcloA = active_refs_.first->nactclo();
  const int nactvirtA = active_refs_.first->nactvirt();
  int nactB = active_refs_.second->nact();
  const int nactcloB = active_refs_.second->nactclo();
  const int nactvirtB = active_refs_.second->nactvirt();
  const int nclosed = active_refs_.first->nclosed() - nactB;
  assert(nclosed == active_refs_.second->nclosed() - nactA);
  const int nlink = active_refs_.first->nlink();
  const int nact = nactA + nactB + nlink;
  const int nvirt = sref_->coeff()->mdim() - nclosed - nact;
  const int dimerbasis = sgeom_->nbasis();

  // active MO matrix
  auto activeMO = make_shared<Matrix>(dimerbasis, nactA + nactB + nlink);
  if (nactA) activeMO->copy_block(0, 0, dimerbasis, nactA, active_refs_.first->coeff()->get_submatrix(0, nclosed, dimerbasis, nactA));
  if (nactB) activeMO->copy_block(0, nactA, dimerbasis, nactB, active_refs_.second->coeff()->get_submatrix(0, nclosed, dimerbasis, nactB));
  if (nlink) activeMO->copy_block(0, nactA + nactB, dimerbasis, nlink, active_refs_.first->coeff()->get_submatrix(0, nclosed + nactA, dimerbasis, nlink));

  // pick active orbitals in localized MOs and reorder to closed - actcloA - actvirtA - actcloB - actvirtB - Link - virtual
  cout << endl << "  o Picking up active orbitals in localized MOs" << endl;
  shared_ptr<Matrix> new_coeff = pick_active(activeMO, sref_->coeff());
      
  auto out_coeff = new_coeff->clone();
  // build out_coeff except for actLA and actLB part
  out_coeff->copy_block(0, 0, dimerbasis, nclosed, new_coeff->get_submatrix(0, 0, dimerbasis, nclosed)); // closed
  out_coeff->copy_block(0, nclosed, dimerbasis, nactcloA, new_coeff->get_submatrix(0, nclosed, dimerbasis, nactcloA)); // actcloA
  out_coeff->copy_block(0, nclosed + nactcloA, dimerbasis, nactvirtA, new_coeff->get_submatrix(0, nclosed + nactcloA, dimerbasis, nactvirtA)); // actvirtA
  out_coeff->copy_block(0, nclosed + nactA + nlink/2, dimerbasis, nactcloB, new_coeff->get_submatrix(0, nclosed + nactA + nlink/2, dimerbasis, nactcloB)); // actcloB
  out_coeff->copy_block(0, nclosed + nactA + nlink + nactcloB, dimerbasis, nactvirtB, 
                         new_coeff->get_submatrix(0, nclosed + nactA + nlink/2 + nactcloB, dimerbasis, nactvirtB)); // actvirtB
  out_coeff->copy_block(0, nclosed + nactA + nactB + nlink, dimerbasis, nvirt, new_coeff->get_submatrix(0, nclosed + nactA + nactB + nlink, dimerbasis, nvirt)); // virtual

  Overlap S(sgeom_);

  // project Link(active) orbitals to fragments and construct new MO coeff
  if (nlink) {
    auto Lcoeff = new_coeff->get_submatrix(0, nclosed + nactA + nactB, dimerbasis, nlink);
  
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
    SB_inv->inverse_symmetric();

    // Forming projected coeff to A and B
    auto projected_A = make_shared<const Matrix>(*SA_inv * *SAmix * *Lcoeff);
    auto projected_B = make_shared<const Matrix>(*SB_inv * *SBmix * *Lcoeff);

    // Get rid of redundancy in projected coeffs
    shared_ptr<Matrix> reduced_MO_A;
    shared_ptr<Matrix> reduced_MO_B;
    {
      auto CC = make_shared<Matrix>(*projected_A % *projected_A);
      VectorB eig(projected_A->mdim());
      CC->diagonalize(eig);
      auto P = CC->get_submatrix(0, CC->mdim()/2, CC->ndim(), CC->mdim()/2);
      auto reduced_A = make_shared<Matrix>(*projected_A * *P);
      reduced_MO_A = make_shared<Matrix>(dimerbasis, reduced_A->mdim());
      reduced_MO_A->copy_block(0, 0, abasis, reduced_A->mdim(), reduced_A->data());
    }
    {
      auto CC = make_shared<Matrix>(*projected_B % *projected_B);
      VectorB eig(projected_B->mdim());
      CC->diagonalize(eig);
      auto P = CC->get_submatrix(0, CC->mdim()/2, CC->ndim(), CC->mdim()/2);
      auto reduced_B = make_shared<Matrix>(*projected_B * *P);
      reduced_MO_B = make_shared<Matrix>(dimerbasis, reduced_B->mdim());
      reduced_MO_B->copy_block(abasis, 0, bbasis, reduced_B->mdim(), reduced_B->data());
    }
    auto reduced_MO_AB = reduced_MO_A->merge(reduced_MO_B);
    assert(nlink == reduced_MO_AB->mdim());
      
    // normalization
    {
      auto csc = make_shared<Matrix>(*reduced_MO_AB % S * *reduced_MO_AB);
      for (int j = 0; j != nlink; ++j) 
        for_each(reduced_MO_AB->element_ptr(0, j), reduced_MO_AB->element_ptr(dimerbasis, j), [&j, &csc](double& p) { p /= sqrt(*csc->element_ptr(j, j)); });
    }

    int iactLA = nclosed + nactcloA;
    int iactLB = nclosed + nactA + nlink/2 + nactcloB;
    {
      int countA = 0; int countB = 0;
      for (int i = 0; i != nlink; ++i) {
        const double sum_A = blas::dot_product(reduced_MO_AB->element_ptr(bounds[0].first, i), bounds[0].second - bounds[0].first,
                                               reduced_MO_AB->element_ptr(bounds[0].first, i));
        const double sum_B = blas::dot_product(reduced_MO_AB->element_ptr(bounds[1].first, i), bounds[1].second - bounds[1].first,
                                               reduced_MO_AB->element_ptr(bounds[1].first, i));
        cout << "sumA : " << sum_A << ", sumB : " << sum_B << endl;
        if (sum_A > sum_B && fabs(sum_A - sum_B) > region_thresh_) {
          cout << "    - projected active orbital(" << i + 1 << ") is assigned to monomer A." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          copy_n(reduced_MO_AB->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, iactLA++));
          ++nactA; ++countA; 
        } 
        else if (sum_A < sum_B && fabs(sum_A - sum_B) > region_thresh_) {
          cout << "    - projected active orbital(" << i + 1 << ") is assigned to monomer B." << endl;
          cout << "      A(" << setw(6) << setprecision(3) << sum_A << "), B(" << setw(6) << setprecision(3) << sum_B << ")" << endl;
          copy_n(reduced_MO_AB->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, iactLB++));
          ++nactB; ++countB;
        } 
        else
          throw runtime_error("projected active orbital still cannot be assigned to either fragment");
      }
      assert(countA == countB); // to make sure projection to both sides
    }
  } 

  // lowdin orthogonalization
  auto tildex = make_shared<Matrix>(*out_coeff % S * *out_coeff);
  tildex->inverse_half();
  out_coeff = make_shared<Matrix>(*out_coeff * *tildex);

  sref_ = make_shared<Reference>(sgeom_, make_shared<Coeff>(*out_coeff), nclosed, nact, nvirt);
}


shared_ptr<Matrix> Dimer::pick_active(shared_ptr<const Matrix> control, shared_ptr<const Matrix> treatment) const {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();

  const int nactcloA = active_refs_.first->nactclo();
  const int nactcloB = active_refs_.second->nactclo();
  const int nactvirtA = active_refs_.first->nactvirt();
  const int nactvirtB = active_refs_.second->nactvirt();
  const int nlink = active_refs_.first->nlink();

  const int nclosed = active_refs_.first->nclosed() - nactB;
  const int dimerbasis = sgeom_->nbasis();
  const int nclosed_HF = isolated_refs_.first->nclosed();

# if 0 // debugging purpose
  cout << "nclosed = " << nclosed << endl;
  cout << "nactA = " << nactA << endl;
  cout << "nactB = " << nactB << endl;
  cout << "nactcloA = " << nactcloA << endl;
  cout << "nactcloB = " << nactcloB << endl;
  cout << "nactvirtA = " << nactvirtA << endl;
  cout << "nactvirtB = " << nactvirtB << endl;
  cout << "nlink = " << nlink << endl;
  cout << "dimerbasis = " << dimerbasis << endl;
#endif
  
  vector<tuple<shared_ptr<const Matrix>, pair<int, int>, int, string, bool>> ovl_info;

  if (nactA) {
    auto activeA = make_shared<Matrix>(dimerbasis, nactA);
    activeA->copy_block(0, 0, dimerbasis, nactA, control->get_submatrix(0, 0, dimerbasis, nactA));
    ovl_info.emplace_back(activeA, make_pair(0, nclosed_HF), nactcloA, "A", true);
    ovl_info.emplace_back(activeA, make_pair(nclosed_HF, dimerbasis), nactvirtA, "A", false);
  }

  if (nactB) {
    auto activeB = make_shared<Matrix>(dimerbasis, nactB);
    activeB->copy_block(0, 0, dimerbasis, nactB, control->get_submatrix(0, nactA, dimerbasis, nactB));
    ovl_info.emplace_back(activeB, make_pair(0, nclosed_HF), nactcloB, "B", true);
    ovl_info.emplace_back(activeB, make_pair(nclosed_HF, dimerbasis), nactvirtB, "B", false);
  }

  if (nlink) {
    auto Link = make_shared<Matrix>(dimerbasis, nlink);
    Link->copy_block(0, 0, dimerbasis, nlink, control->get_submatrix(0, 0, dimerbasis, nlink));
    ovl_info.emplace_back(Link, make_pair(0, dimerbasis), nlink, "L", true);
  }

  const Overlap S(sgeom_);

  shared_ptr<Matrix> out_coeff = treatment->clone();
  size_t active_position = nclosed;

  set<int> mask;
  for (int i = 0; i != out_coeff->mdim(); ++i) mask.insert(i);

  for (auto& subset : ovl_info) {
    const Matrix& active = *get<0>(subset);
    const pair<int , int> bounds = get<1>(subset);
    const int norb = get<2>(subset);
    const string set_name = get<3>(subset);
    const bool closed = get<4>(subset);

    shared_ptr<const Matrix> subcoeff = treatment->slice_copy(bounds.first, bounds.second);

    const Matrix overlap(active % S * *subcoeff);

    multimap<double, int> norms;

    for (int i = 0; i != overlap.mdim(); ++i) {
      const double norm = blas::dot_product(overlap.element_ptr(0, i), overlap.ndim(), overlap.element_ptr(0, i));
      norms.emplace(norm, i);
    }

    cout << endl << "  o Forming dimer's active orbitals arising from " << (closed ? "closed " : "virtual") << set_name
                 << " orbitals. Threshold for includsion in candidate space: " << setw(6) << setprecision(3) << active_thresh_ << endl;

    vector<int> active_list;
    double max_overlap, min_overlap;
    {
      auto end = norms.rbegin(); advance(end, norb);
      end = find_if(end, norms.rend(), [this] (const pair<const double, int>& p) { return p.first < active_thresh_; });
      for_each(norms.rbegin(), end, [&active_list] (const pair<const double, int>& p) { active_list.emplace_back(p.second); });
      auto mnmx = minmax_element(norms.rbegin(), end);
      tie(min_overlap, max_overlap) = make_tuple(mnmx.first->first, mnmx.second->first);
    }

    const int active_size = active_list.size();
    cout << "    - size of candidate space: " << active_size << endl;
    cout << "    - largest overlap with monomer space: " << max_overlap << ", smallest: " << min_overlap << endl;

    if (active_size != norb)
      throw runtime_error("active size != norb, SVD required or check for other reasons...");
    else {
      const set<int> active_set(active_list.begin(), active_list.end());
      for (size_t i = 0; i != subcoeff->mdim(); ++i) {
        if (active_set.count(i)) {
          const int imo = bounds.first + i;
          assert(mask.count(imo));
          mask.erase(imo);
          copy_n(subcoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, active_position++)); // Finally should be [cA, vA, cB, vB, L]
        }
      }
    }
  }

  // fill common closed and virtual subspaces
  size_t closed_position = 0;
  for (int i = 0; i != nclosed_HF; ++i)
    if (mask.count(i))
      copy_n(treatment->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, closed_position++));

  size_t virt_position = nclosed + nactA + nactB + nlink;
  for (int i = nclosed_HF; i != dimerbasis; ++i)
    if (mask.count(i))
      copy_n(treatment->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, virt_position++));

  return out_coeff;
}


shared_ptr<Matrix> Dimer::form_semi_canonical_metal(shared_ptr<const PTree> idata) const {
  const int nactA = active_refs_.first->nact();
  const int nactB = active_refs_.second->nact();
  const int nlink = active_refs_.first->nlink();
  const int nact = nactA + nactB + nlink;

  const int nclosed = active_refs_.first->nclosed() - nactB;
  const int nactcloA = active_refs_.first->nactclo();
  const int nactcloB = active_refs_.second->nactclo();
  const int nvirt = sref_->coeff()->mdim() - nclosed - nact;
//  const int nactvirtA = active_refs_.first->nactvirt();
//  const int nactvirtB = active_refs_.second->nactvirt();

  const int dimerbasis = sgeom_->nbasis();
  
  auto tmp_coeff = sref_->coeff()->copy();
  
  // form AO Fock
  shared_ptr<const Matrix> ofockao;
  {
    const int nocc = nclosed + nactcloA + nactcloB + nlink;
    auto occoeff = make_shared<Matrix>(dimerbasis, nocc);
    occoeff->copy_block(0, 0                        , dimerbasis, nclosed           , tmp_coeff->get_submatrix(0, 0, dimerbasis, nclosed)); // closed
    occoeff->copy_block(0, nclosed                  , dimerbasis, nactcloA + nlink/2, tmp_coeff->get_submatrix(0, nclosed, dimerbasis, nactcloA + nlink/2)); // occ_activeA
    occoeff->copy_block(0, nclosed + nactA + nlink/2, dimerbasis, nactcloB + nlink/2, tmp_coeff->get_submatrix(0, nclosed + nactA + nlink/2, dimerbasis,
                                                                                                                nactcloB + nlink/2)); // occ_activeB
    auto population = make_shared<Matrix>(nocc, nocc);
    population->add_diag(2.0);
    {
      const int LA_start = nclosed + nactcloA;
      const int LA_end = LA_start + nlink/2;
      const int LB_start = nclosed + nactA + nlink/2 + nactcloB;
      const int LB_end = LB_start + nlink/2;
      for (int i = LA_start; i != LA_end; ++i)
        *population->element_ptr(i, i) = 1.0;
      for (int i = LB_start; i != LB_end; ++i)
        *population->element_ptr(i,i) = 1.0;
    }
    auto dimerdensity = make_shared<const Matrix>(*occoeff * *population ^ *occoeff);
    ofockao = make_shared<Fock<1>>(sgeom_, sref_->hcore(), dimerdensity, occoeff, /*store*/false, /*rhf*/true);
  }

  auto out_coeff = sref_->coeff()->clone();
  // form MO Fock
  size_t position = 0;
  { // common closed
    VectorB eigs(nclosed);
    auto mocoeff = tmp_coeff->slice_copy(0, nclosed);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock;
    for (int i = 0; i != nclosed; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, position++));
  }
  { // activeA
    VectorB eigs(nactA + nlink/2);
    auto mocoeff = tmp_coeff->get_submatrix(0, nclosed, dimerbasis, nactA + nlink/2);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock;
    for (int i = 0; i != nactA + nlink/2; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, position++));
  }
  { // activeB
    VectorB eigs(nactB + nlink/2);
    auto mocoeff = tmp_coeff->get_submatrix(0, nclosed + nactA + nlink/2, dimerbasis, nactB + nlink/2);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock;
    for (int i = 0; i != nactB + nlink/2; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, position++));
  }
cout << "5" << endl;
  if (nvirt != 0) { // virtual
    VectorB eigs(nvirt);
    auto mocoeff = tmp_coeff->get_submatrix(0, nclosed + nact + nlink, dimerbasis, nvirt);
    auto fock = make_shared<Matrix>(*mocoeff % *ofockao * *mocoeff);
    fock->diagonalize(eigs);
    *mocoeff *= *fock;
    for (int i = 0; i != nvirt; ++i)
      copy_n(mocoeff->element_ptr(0, i), dimerbasis, out_coeff->element_ptr(0, position++));
  }

  return out_coeff;
}


