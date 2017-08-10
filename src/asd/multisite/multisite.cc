//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/multisite/multisite.h>
#include <src/mat1e/overlap.h>
#include <src/util/muffle.h>
#include <src/scf/hf/rhf.h>
#include <src/wfn/get_energy.h>
#include <src/util/io/moldenout.h>
#include <src/mat1e/mixedbasis.h>
#include <src/integral/os/overlapbatch.h>
#include <src/scf/hf/fock.h>
#include <src/ci/fci/harrison.h>

using namespace std;
using namespace bagel;

extern "C" {
  void dgels_(const char* trans, const int* m, const int* n, const int* nrhs, double* a, const int* lda,
              double* b, const int* ldb, double* work, const int* lwork, int* info);
}

namespace {
  void dgels_(const char* trans, const int m, const int n, const int nrhs, double* a, const int lda, 
              double* b, const int ldb, double* work, const int lwork, int& info) {
    ::dgels_(trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info); }
}

MultiSite::MultiSite(shared_ptr<const PTree> input, vector<shared_ptr<const Reference>> refs) : input_(input), nsites_(refs.size()) {
  // build the super geometry
  const shared_ptr<const PTree> mdata = input_->get_child_optional("molecule");
  for (auto& r : refs)
    geoms_.push_back(mdata ? make_shared<Geometry>(*r->geom(), mdata) : r->geom());

  auto envdata = input_->get_child_optional("environment");
  vector<shared_ptr<const Geometry>> geovec = geoms_;
  if (envdata) {
    Muffle hide_cout;
    geovec.push_back(make_shared<Geometry>(envdata));
  }

  auto sgeom = make_shared<Geometry>(geovec);

  // build combined references
  for (auto& r : refs)
    isolated_refs_.push_back(r->project_coeff(sgeom, false));
  const size_t norb = accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0ul, [] (size_t x, shared_ptr<const Reference> r) { return x + r->coeff()->mdim(); });
  Matrix coeff(sgeom->nbasis(), norb);

  // pull out all of the orbitals and put them in closed, active, virtual order
  int current = 0;
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,0),                      r->nclosed()*sgeom->nbasis(), coeff.element_ptr(0, current));
    closed_bounds_.push_back({current, current+r->nclosed()});
    current += r->nclosed();
  }
  for (auto& r : isolated_refs_) {
    copy_n(r->coeff()->element_ptr(0,r->nclosed()+r->nact()), r->nvirt()*sgeom->nbasis(),   coeff.element_ptr(0, current));
    virt_bounds_.push_back({current, current+r->nvirt()});
    current += r->nvirt();
  }

  Overlap S(sgeom);
  Matrix transform(coeff % S * coeff);
  transform.inverse_half();
  auto scoeff = make_shared<Coeff>(coeff * transform);

  sref_ = make_shared<Reference>(sgeom, scoeff,
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nclosed(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); }),
                                 accumulate(isolated_refs_.begin(), isolated_refs_.end(), 0, [] (int x, shared_ptr<const Reference> r) { return x + r->nvirt(); }));
}


MultiSite::MultiSite(shared_ptr<const PTree> itree, shared_ptr<const Reference> ref) : input_(itree), geom_(ref->geom()), hf_ref_(ref) {
  
  cout << " ===== Constructing MultiSite Molecular Orbitals ===== " << endl;
  
  // reorder MO coeff to closed - active - virtual
  set_active_metal();
  
  // project active orbitals to fragments
  project_active();
  
  // canonicalize active orbitals in sub spaces
  canonicalize();
}


// pick active orbitals from reference actorbs
void MultiSite::set_active_metal() {
  
  auto isp = input_->get_child("multisite_active");
  set<int> ActList;
  for (auto& s : *isp) { ActList.insert(lexical_cast<int>(s->data())-1); };
  
  active_electrons_ = input_->get_vector<int>("active_electrons");
  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  const int charge = input_->get<int>("charge", 0);

  const int multisitebasis = geom_->nbasis();
  const int nactive = ActList.size();
  const int nclosed = (geom_->nele() - charge - nactele) / 2;
  assert ((geom_->nele() - charge -nactele) % 2 == 0);
  const int nvirt = multisitebasis - nactive - nclosed;

  auto hf_coeff = hf_ref_->coeff();
  auto out_coeff = hf_coeff->clone();
  
  int iclosed = 0;
  int closed_position = 0;
  int active_position = nclosed;
  int virt_position = nclosed + nactive;
  for (int i = 0; i != multisitebasis; ++i) {
    if (ActList.find(i) != ActList.end()) copy_n(hf_coeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, active_position++));
    else copy_n(hf_coeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, ((iclosed < nclosed) ? (++iclosed, closed_position++) : virt_position++)));
  } 
  assert (virt_position == multisitebasis + 1);

  active_ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*out_coeff)), nclosed, nactive, nvirt);
}


void MultiSite::project_active() {

  // define fragments (regions) of the multisite
  vector<int> region_sizes = input_->get_vector<int>("region_sizes");
  nsites_ = region_sizes.size();
  vector<pair<int, int>> basis_bounds;
  {
    int nbasis = 0;
    int natoms = 0;
    for (int& rsize : region_sizes) {
      const int atombegin = natoms;
      const int basisbegin = nbasis;
      for (int iatom = atombegin; iatom != atombegin + rsize; ++iatom)
        nbasis += geom_->atoms()[iatom]->nbasis();

      natoms += rsize;
      if (basisbegin != nbasis)
        basis_bounds.emplace_back(basisbegin, nbasis);
    }
    if (natoms != count_if(geom_->atoms().begin(), geom_->atoms().end(), [](const shared_ptr<const Atom> a){return !a->dummy();}))
      throw logic_error("All atoms must be assigned to regions");
  }

  const int nclosed = active_ref_->nclosed();
  const int nact = active_ref_->nact();
  const int nvirt = active_ref_->nvirt();
  const int multimerbasis = active_ref_->geom()->nbasis();

  Overlap S(geom_);

  // project active orbitals to each fragement and do SVD
  active_sizes_ = input_->get_vector<int>("active_sizes");
  int basisoffset = 0;
  int orboffset = 0;
  auto actcoeff = active_ref_->coeff()->get_submatrix(0, nclosed, multimerbasis, nact);
  auto tmp = actcoeff->clone();
  for (int i = 0; i != nsites_; ++i) {
    const int nbasis = basis_bounds[i].second - basis_bounds[i].first;
    auto Sfrag = make_shared<Matrix>(nbasis, nbasis);
    auto Smix = make_shared<Matrix>(nbasis, multimerbasis);
    Sfrag->copy_block(0, 0, nbasis, nbasis, S.get_submatrix(basisoffset, basisoffset, nbasis, nbasis));
    Smix->copy_block(0, 0, nbasis, multimerbasis, S.get_submatrix(basisoffset, 0, nbasis, multimerbasis));

    auto S_inv = make_shared<Matrix>(*Sfrag);
    S_inv->inverse_symmetric();
 
    // forming projected orbitals of current fragment
    auto projected = make_shared<const Matrix>(*S_inv * *Smix * *actcoeff);

    // get rid of redundancy
    shared_ptr<Matrix> reduced;
    {
      auto CC = make_shared<Matrix>(*projected % *Sfrag *  *projected);
      VectorB eig(projected->mdim());
      CC->diagonalize(eig);
      auto P = CC->get_submatrix(0, (CC->mdim() - active_sizes_[i]), CC->ndim(), active_sizes_[i]);
      auto temp = make_shared<const Matrix>(*projected * *P);
      reduced = make_shared<Matrix>(multimerbasis, temp->mdim());
      reduced->copy_block(basisoffset, 0, nbasis, temp->mdim(), temp->data());
    }
    tmp->copy_block(0, orboffset, multimerbasis, active_sizes_[i], reduced->data());

    basisoffset += nbasis;
    orboffset += active_sizes_[i];
  }
  assert(orboffset == nact);
  
  // Now we want to solve min|| Act_coeff * x - tmp||, 
  // DGELS is used to obtain the transformation matrix, projected orbitals cannot have linear dependancy
  shared_ptr<Matrix> solution;
  {
    int N = actcoeff->ndim();
    int M = actcoeff->mdim();
    assert(N > M);
    assert(M == tmp->mdim());
    int lwork = 10 * N;
    unique_ptr<double[]> work(new double[lwork]);
    auto actcopy = actcoeff->copy();
    double* adata = actcopy->data();
    double* bdata = tmp->data();
    int info = 0;
    dgels_("N", N, M, M, adata, N, bdata, N, work.get(), lwork, info);
    if (info != 0) throw runtime_error("dgels failed in projecting coeff");
  
    solution = tmp->get_submatrix(0, 0, M, M);
  }

  // Lowdin orthogonalization for transformation matrix
  {
    auto tildeX = make_shared<Matrix>(*solution % *solution);
    tildeX->inverse_half();
    cout << "    *** If linear dependency is detected, you shall try to find better initial active orbital guess ***   " << endl;
    solution = make_shared<Matrix>(*solution * *tildeX);
  }

  actcoeff = make_shared<Matrix>(*actcoeff * *solution);

  auto new_coeff = active_ref_->coeff()->copy();
  new_coeff->copy_block(0, nclosed, multimerbasis, nact, actcoeff->data());

  // update multimer information
  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*new_coeff)), nclosed, nact, nvirt);
}


void MultiSite::canonicalize() {

  auto out_coeff = ref_->coeff()->copy();
  
  const int multimerbasis = ref_->geom()->nbasis();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nactele = ref_->geom()->nele() - 2*nclosed;

  // projected active orbitals are partially occupied, modify density matrix accordingly
  auto density_coeff = ref_->coeff()->slice_copy(0, nclosed + nact);
  blas::scale_n(sqrt(0.5*nactele/nact), density_coeff->element_ptr(0, nclosed), multimerbasis * nact);
  auto density = density_coeff->form_density_rhf(nclosed + nact);
  auto fock = make_shared<const Fock<1>>(geom_, ref_->hcore(), density, density_coeff);

  // canonicalize active orbitals in each fragment
  {
    int offset = nclosed;
    for (int i = 0; i != active_sizes_.size(); ++i) {
      const int norb = active_sizes_[0];
      auto subspace = ref_->coeff()->slice_copy(offset, offset + norb);
      auto subfock = make_shared<Matrix>(*subspace % *fock * *subspace);
      VectorB eigs(norb);
      subfock->diagonalize(eigs);
      subspace = make_shared<Matrix>(*subspace * *subfock);
      copy_n(subspace->data(), multimerbasis * norb, out_coeff->element_ptr(0, offset));
      offset += norb;
    }
  }

  ref_ = make_shared<Reference>(geom_, make_shared<Coeff>(move(*out_coeff)), nclosed, nact, ref_->nvirt());
}


shared_ptr<Reference> MultiSite::build_reference(const int site, const vector<bool> meanfield, bool metal) const {

  assert(meanfield.size()==nsites_ && site<nsites_ && site>=0);
  if (!metal) {
    vector<shared_ptr<const MatView>> closed_orbitals = {make_shared<MatView>(sref_->coeff()->slice(0, sref_->nclosed()))};
    const int act_start = accumulate(active_refs_.begin(), active_refs_.begin() + site, sref_->nclosed(),
                                      [] (int x, shared_ptr<const Reference> r) { return x + r->nact(); });
    const int act_fence = act_start + active_refs_[site]->nact();
    const MatView active_orbitals = sref_->coeff()->slice(act_start, act_fence);

    int current = sref_->nclosed();
    for (int i = 0; i < nsites_; ++i) {
      if (meanfield[i] && i!=site)
        closed_orbitals.push_back(make_shared<const MatView>(sref_->coeff()->slice(current, current+isolated_refs_[i]->nclosed()-active_refs_[i]->nclosed())));
      current += active_refs_[i]->nact();
    }

    const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [] (int x, shared_ptr<const MatView> m) { return x + m->mdim(); });
    const int nact = active_orbitals.mdim();

    auto out = make_shared<Matrix>(sref_->geom()->nbasis(), nclosed+nact);

    current = 0;
    closed_orbitals.push_back(make_shared<MatView>(active_orbitals));
    for (auto& orbitals : closed_orbitals) {
      copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
      current += orbitals->mdim();
    }

    return make_shared<Reference>(sref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);

  } else {// Raymond version

    const int act_start = accumulate(active_sizes_.begin(), active_sizes_.begin()+site, ref_->nclosed());
    const int nact = active_sizes_[site];
    const Matrix active_orbitals = ref_->coeff()->slice(act_start, act_start + nact);

    vector<shared_ptr<const Matrix>> closed_orbitals = {make_shared<Matrix>(ref_->coeff()->slice(0, ref_->nclosed()))};
    int current = ref_->nclosed();
    for (int i = 0; i != nsites_; ++i) {
      if (meanfield[i] && i != site) {
        const int nele = active_electrons_[i];
        const int norb = active_sizes_[i];
        auto scale_coeff = ref_->coeff()->slice_copy(current, current + active_sizes_[i]);
        blas::scale_n(sqrt(0.5*nele/norb), scale_coeff->data(), scale_coeff->size());
        closed_orbitals.push_back(scale_coeff);
      }
      current += active_sizes_[i];
    }

    const int nclosed = accumulate(closed_orbitals.begin(), closed_orbitals.end(), 0, [] (int x, shared_ptr<const Matrix> m) { return x + m->mdim(); });

    auto out = make_shared<Matrix>(ref_->geom()->nbasis(), nclosed + nact);
    current = 0;
    closed_orbitals.push_back(make_shared<Matrix>(active_orbitals));
    for (auto& orbitals : closed_orbitals) {
      copy_n(orbitals->data(), orbitals->size(), out->element_ptr(0, current));
      current += orbitals->mdim();
    }

    return make_shared<Reference>(ref_->geom(), make_shared<Coeff>(move(*out)), nclosed, nact, 0);
  }
}
