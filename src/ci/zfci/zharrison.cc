//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/ci/zfci/zharrison.h>
#include <src/ci/zfci/relspace.h>
#include <src/util/math/comb.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/giao/relhcore_london.h>
#include <src/mat1e/rel/reloverlap.h>
#include <src/mat1e/giao/reloverlap_london.h>
/**/
#include <src/mat1e/rel/general_small1e.h>
#include <src/integral/os/overlapbatch.h>
/**/

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::ZHarrison)

using namespace std;
using namespace bagel;

ZHarrison::ZHarrison(std::shared_ptr<const PTree> idat, shared_ptr<const Geometry> g, shared_ptr<const Reference> r, const int ncore, const int norb, const int nstate, std::shared_ptr<const RelCoeff_Block> coeff_zcas, const bool restricted)
 : Method(idat, g, r), ncore_(ncore), norb_(norb), nstate_(nstate), restarted_(false) {
  if (!ref_) throw runtime_error("ZFCI requires a reference object");

  auto rr = dynamic_pointer_cast<const RelReference>(ref_);
  if (!rr) throw runtime_error("ZFCI currently requires a relativistic reference object");

  const bool frozen = idata_->get<bool>("frozen", false);
  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_fci", max_iter_);
  davidson_subspace_ = idata_->get<int>("davidson_subspace", 20);
  thresh_ = idata_->get<double>("thresh", 1.0e-10);
  thresh_ = idata_->get<double>("thresh_fci", thresh_);
  print_thresh_ = idata_->get<double>("print_thresh", 0.05);
  restart_ = idata_->get<bool>("restart", false);

  states_ = idata_->get_vector<int>("state", 0);
  nstate_ = 0;
  for (int i = 0; i != states_.size(); ++i)
    nstate_ += states_[i] * (i+1); // 2S+1 for 0, 1/2, 1, ...

  // if nstate was specified on construction, it should match
  assert(nstate == nstate_ || nstate == -1);

  gaunt_ = idata_->get<bool>("gaunt", rr->gaunt());
  breit_ = idata_->get<bool>("breit", rr->breit());
  if (gaunt_ != rr->gaunt())
    geom_ = geom_->relativistic(gaunt_);

  // Invoke Kramer's symmetry for any case without magnetic field
  tsymm_ = !geom_->magnetism();

  if (ncore_ < 0)
    ncore_ = idata_->get<int>("ncore", (frozen ? geom_->num_count_ncore_only()/2 : 0));
  if (norb_  < 0)
    norb_ = idata_->get<int>("norb", (rr->relcoeff()->mdim()/2-ncore_));

  // additional charge
  charge_ = idata_->get<int>("charge", 0);

  nele_ = geom_->nele() - charge_ - ncore_*2;

  if (norb_ < 0 || norb_ + ncore_ > geom_->nbasis())
    throw runtime_error("Invalid number of active orbitals");
  if (nele_ < 0)
    throw runtime_error("Number of active electrons is less than zero.");

  energy_.resize(nstate_);

  print_header();
  cout << "    * nstate   : " << setw(6) << nstate_ << endl;
  cout << "    * nclosed  : " << setw(6) << ncore_ << endl;
  cout << "    * nact     : " << setw(6) << norb_ << endl;
  cout << "    * nvirt    : " << setw(6) << ((coeff_zcas ? coeff_zcas->mdim() : rr->relcoeff_full()->mdim())/2-ncore_-norb_) << endl;

  if (!geom_->dfs())
    geom_ = geom_->relativistic(gaunt_);

  space_ = make_shared<RelSpace>(norb_, nele_);
  int_space_ = make_shared<RelSpace>(norb_, nele_-2, /*mute*/true, /*link up*/true);
  assert((restricted && coeff_zcas) || (!restricted && !coeff_zcas));

  // obtain the coefficient matrix in striped format
  shared_ptr<const RelCoeff_Block> coeff;
  if (coeff_zcas) {
    coeff = coeff_zcas;
  } else {
    auto scoeff = make_shared<const RelCoeff_Striped>(*rr->relcoeff_full(), ncore_, norb_, rr->relcoeff_full()->mdim()/4-ncore_-norb_, rr->relcoeff_full()->mdim()/2);

    shared_ptr<const ZMatrix> hcore, overlap;
    if (!geom_->magnetism()) {
      overlap = make_shared<RelOverlap>(geom_);
      hcore = make_shared<RelHcore>(geom_);
    } else {
      overlap = make_shared<RelOverlap_London>(geom_);
      hcore = make_shared<RelHcore_London>(geom_);
    }

    // then compute Kramers adapated coefficient matrices
    scoeff = scoeff->init_kramers_coeff(geom_, overlap, hcore, geom_->nele()-charge_, tsymm_, gaunt_, breit_);

    // generate modified virtual orbitals, if requested
    const bool mvo = idata_->get<bool>("generate_mvo", false);
    if (mvo) {
      const bool hcore_mvo = idata_->get<bool>("hcore_mvo", false);
      const int ncore_mvo = idata_->get<int>("ncore_mvo", geom_->num_count_ncore_only());
      if (ncore_mvo == 2*rr->relcoeff_full()->nocc())
        cout << "    +++ Modified virtuals are Dirac-Fock orbitals with this choice of the core +++ "<< endl;
      else
        scoeff = scoeff->generate_mvo(geom_, overlap, hcore, ncore_mvo, rr->relcoeff_full()->nocc(), hcore_mvo, tsymm_, gaunt_, breit_);
    }

    // Reorder as specified in the input so frontier orbitals contain the desired active space
    const shared_ptr<const PTree> iactive = idata_->get_child_optional("active");
    if (iactive) {
      set<int> active_indices;
      // Subtracting one so that orbitals are input in 1-based format but are stored in C format (0-based)
      for (auto& i : *iactive)
        active_indices.insert(lexical_cast<int>(i->data()) - 1);
      scoeff = scoeff->set_active(active_indices, geom_->nele()-charge_, tsymm_);
    }
    coeff = scoeff->block_format();
  }

  update(coeff, restricted);

}


void ZHarrison::print_header() const {
  cout << "  ----------------------------" << endl;
  cout << "  Relativistic FCI calculation" << endl;
  cout << "  ----------------------------" << endl << endl;
  cout << "    * Correlation of " << nele_ << " active electrons in " << norb_ << " orbitals."  << endl;
  cout << "    * Time-reversal symmetry " << (tsymm_ ? "will be assumed." : "violation will be permitted.") << endl;
  cout << "    * gaunt    : " << (gaunt_ ? "true" : "false") << endl;
  cout << "    * breit    : " << (breit_ ? "true" : "false") << endl;
}


// generate initial vectors
void ZHarrison::generate_guess(const int nelea, const int neleb, const int nstate, std::shared_ptr<RelZDvec> out, const int offset) {
  shared_ptr<const Determinants> cdet = space_->finddet(nelea, neleb);
  int ndet = nstate*10;
  int oindex = offset;
  while (oindex < offset+nstate) {
    vector<pair<bitset<nbit__>, bitset<nbit__>>> bits = detseeds(ndet, nelea, neleb);

    // Spin adapt detseeds
    oindex = offset;
    vector<bitset<nbit__>> done;
    for (auto& it : bits) {
      bitset<nbit__> alpha = it.second;
      bitset<nbit__> beta = it.first;
      bitset<nbit__> open_bit = (alpha^beta);

      // This can happen if all possible determinants are checked without finding nstate acceptable ones.
      if (alpha.count() + beta.count() != nele_)
        throw logic_error("ZFCI::generate_guess produced an invalid determinant.  Check the number of states being requested.");

      // make sure that we have enough unpaired alpha
      const int unpairalpha = (alpha ^ (alpha & beta)).count();
      const int unpairbeta  = (beta ^ (alpha & beta)).count();
      if (unpairalpha-unpairbeta < nelea-neleb) continue;

      if (find(done.begin(), done.end(), open_bit) != done.end()) continue;

      done.push_back(open_bit);
      pair<vector<tuple<int, int, int>>, double> adapt;
      adapt = space_->finddet(nelea, neleb)->spin_adapt(nelea-neleb, alpha, beta);

      const double fac = adapt.second;
      for (auto& iter : adapt.first) {
        out->find(nelea, neleb)->data(oindex)->element(get<0>(iter), get<1>(iter)) = get<2>(iter)*fac;
      }
      cout << "     guess " << setw(3) << oindex << ":   closed " <<
            setw(20) << left << print_bit(alpha&beta, norb_) << " open " << setw(20) << print_bit(open_bit, norb_) << right << endl;

      ++oindex;
      if (oindex == offset+nstate) break;
    }

    if (oindex < offset+nstate) {
      for (int i = offset; i != offset+oindex; ++i) {
        out->find(nelea, neleb)->data(i)->zero();
      }
      ndet *= 4;
    }
  }
  assert(oindex == offset+nstate);
  cout << endl;
}


// returns seed determinants for initial guess
vector<pair<bitset<nbit__> , bitset<nbit__>>> ZHarrison::detseeds(const int ndet, const int nelea, const int neleb) const {
  shared_ptr<const Determinants> cdet = space_->finddet(nelea, neleb);

  multimap<double, pair<bitset<nbit__>,bitset<nbit__>>> tmp;
  for (int i = 0; i != ndet; ++i) tmp.emplace(-1.0e10*(1+i), make_pair(bitset<nbit__>(0),bitset<nbit__>(0)));

  double* diter = denom_->find(cdet->nelea(), cdet->neleb())->data();
  for (auto& aiter : cdet->string_bits_a()) {
    for (auto& biter : cdet->string_bits_b()) {
      const double din = -(*diter);
      if (tmp.begin()->first < din) {
        tmp.emplace(din, make_pair(biter, aiter));
        tmp.erase(tmp.begin());
      }
      ++diter;
    }
  }
  assert(tmp.size() == ndet || ndet > cdet->string_bits_a().size()*cdet->string_bits_b().size());
  vector<pair<bitset<nbit__> , bitset<nbit__>>> out;
  for (auto iter = tmp.rbegin(); iter != tmp.rend(); ++iter)
    out.push_back(iter->second);
  return out;
}


void ZHarrison::compute() {
  Timer pdebug(2);

  if (geom_->nirrep() > 1) throw runtime_error("ZFCI: C1 only at the moment.");

  if (!restarted_) {
    // Creating an initial CI vector
    cc_ = make_shared<RelZDvec>(space_, nstate_); // B runs first

    // TODO really we should check the number of states for each S value, rather than total number
    const static Comb combination;
    const size_t max_states = combination(2*norb_, nele_);
    if (nstate_ > max_states) {
      const string space = "(" + to_string(nele_) + "," + to_string(norb_) + ")";
      throw runtime_error("Wrong states specified - a " + space + " active space can only produce " + to_string(max_states) + " eigenstates.");
    }

    // find determinants that have small diagonal energies
    int offset = 0;
    for (int ispin = 0; ispin != states_.size(); ++ispin) {
      int nstate = 0;
      for (int i = ispin; i != states_.size(); ++i)
        nstate += states_[i];

      if (nstate == 0)
        continue;

      if ((geom_->nele()+ispin-charge_) % 2 == 1) {
        if (states_[ispin] == 0) {
          continue;
        } else {
          if ((geom_->nele()-charge_) % 2 == 0) throw runtime_error("Wrong states specified - only integer spins are allowed for even electron counts.");
          else throw runtime_error("Wrong states specified - only half-integer spins are allowed for odd electron counts.");
        }
      }

      const int nelea = (geom_->nele()+ispin-charge_)/2 - ncore_;
      const int neleb = (geom_->nele()-ispin-charge_)/2 - ncore_;
      if (neleb < 0) throw runtime_error("Wrong states specified - there are not enough active electrons for the requested spin state.");
      if (nelea > norb_) throw runtime_error("Wrong states specified - there are not enough active orbitals for the requested spin state.");

      generate_guess(nelea, neleb, nstate, cc_, offset);
      offset += nstate;
      if (nelea != neleb) {
        generate_guess(neleb, nelea, nstate, cc_, offset);
        offset += nstate;
      }
    }
    pdebug.tick_print("guess generation");

    // Davidson utility
    davidson_ = make_shared<DavidsonDiag<RelZDvec, ZMatrix>>(nstate_, davidson_subspace_);
  }

  // nuclear energy retrieved from geometry
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();

  // main iteration starts here
  cout << "  === Relativistic FCI iteration ===" << endl << endl;
  // 0 means not converged
  vector<int> conv(nstate_,0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    Timer fcitime;

#ifndef DISABLE_SERIALIZATION
    if (restart_) {
      stringstream ss; ss << "zfci_" << iter;
      OArchive ar(ss.str());
      ar << static_cast<Method*>(this);
    }
#endif

    // form a sigma vector given cc
    shared_ptr<RelZDvec> sigma = form_sigma(cc_, jop_, conv);
    pdebug.tick_print("sigma vector");

    const vector<double> energies = davidson_->compute(cc_->dvec(conv), sigma->dvec(conv));
    // get residual and new vectors
    vector<shared_ptr<RelZDvec>> errvec = davidson_->residual();
    for (auto& i : errvec)
      i->synchronize();
    pdebug.tick_print("davidson");

    // compute errors
    vector<double> errors;
    for (int i = 0; i != nstate_; ++i) {
      errors.push_back(errvec[i]->rms());
      conv[i] = static_cast<int>(errors[i] < thresh_);
    }
    pdebug.tick_print("error");

    if (!*min_element(conv.begin(), conv.end())) {
      // denominator scaling

      auto ctmp = errvec.front()->clone();

      for (int ist = 0; ist != nstate_; ++ist) {
        if (conv[ist]) continue;
        for (auto& ib : space_->detmap()) {
          const int na = ib.second->nelea();
          const int nb = ib.second->neleb();
          const size_t size = cc_->find(na, nb)->data(ist)->size();
          complex<double>* target_array = ctmp->find(na, nb)->data();
          complex<double>* source_array = errvec[ist]->find(na, nb)->data();
          double* denom_array = denom_->find(na, nb)->data();
          const double en = energies[ist];
          for (int i = 0; i != size; ++i) {
            target_array[i] = source_array[i] / min(en - denom_array[i], -0.1);
          }
        }
        ctmp->normalize();
        cc_->set_data(ist, ctmp);
      }
    }
    pdebug.tick_print("denominator");

    // printing out
    if (nstate_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstate_; ++i) {
      cout << setw(7) << iter << setw(4) << i << " " << setw(2) << (conv[i] ? "*" : " ")
                              << setw(17) << fixed << setprecision(8) << energies[i]+nuc_core << "   "
                              << setw(10) << scientific << setprecision(2) << errors[i] << fixed << setw(10) << setprecision(2)
                              << fcitime.tick() << endl;
      energy_[i] = energies[i]+nuc_core;
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }
  // main iteration ends here

  cc_ = make_shared<RelZDvec>(davidson_->civec());
  cc_->print(print_thresh_);

#if 0
  for (auto& iprop : properties_) {
    iprop->compute(cc_);
    iprop->print();
  }
#endif


  // TODO When the Property class is implemented, this should be one
  if (idata_->get<bool>("aniso", false)) {
    compute_pseudospin_hamiltonian();
  }
}


shared_ptr<const ZMatrix> ZHarrison::swap_pos_neg(shared_ptr<const ZMatrix> coeffin) const {
  auto out = coeffin->clone();
  const int n = coeffin->ndim();
  const int m = coeffin->mdim()/2;
  assert(n % 4 == 0 && m % 2 == 0 && m * 2 == coeffin->mdim());
  out->copy_block(0, 0, n, m, coeffin->get_submatrix(0, m, n, m));
  out->copy_block(0, m, n, m, coeffin->get_submatrix(0, 0, n, m));
  return out;
}


shared_ptr<const RelCIWfn> ZHarrison::conv_to_ciwfn() const {
  using PairType = pair<shared_ptr<const RelSpace>,shared_ptr<const RelSpace>>;
  return make_shared<RelCIWfn>(geom_, ncore_, norb_, nstate_, energy_, cc_, make_shared<PairType>(make_pair(space_, int_space_)));
}


void ZHarrison::compute_pseudospin_hamiltonian() const {

  /**  Part 1: Compute numerical pseudospin Hamiltonian by diagonalizing S_z matrix  **/

  // First, we create spin matrices in the atomic orbital basis
  // TODO Create a class for this
  const int n = geom_->nbasis();
  auto overlap = make_shared<Overlap>(geom_);
  auto spinx = make_shared<ZMatrix>(4*n, 4*n);
  auto spiny = make_shared<ZMatrix>(4*n, 4*n);
  auto spinz = make_shared<ZMatrix>(4*n, 4*n);
  const complex<double> imag(0.0, 1.0);

  // Large component
  spinx->add_real_block( 0.5,   n,   0, n, n, *overlap);
  spinx->add_real_block( 0.5,   0,   n, n, n, *overlap);

  spiny->add_real_block( 0.5*imag,   n,   0, n, n, *overlap);
  spiny->add_real_block(-0.5*imag,   0,   n, n, n, *overlap);

  spinz->add_real_block( 0.5,   0,   0, n, n, *overlap);
  spinz->add_real_block(-0.5,   n,   n, n, n, *overlap);

  // Small component
  // TODO Simplify this code
  // Commented out lines cancel at zero-field; can be replaced with field*overlap for GIAO-RMB
  const double w = 1.0/(8.0*c__*c__);
  auto smallints = make_shared<General_Small1e<OverlapBatch>>(geom_);

  // x^x contributions
  spinx->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[0]);
  spinx->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[0]);
  spiny->add_real_block( imag*w, 2*n, 3*n, n, n, (*smallints)[0]);
  spiny->add_real_block(-imag*w, 3*n, 2*n, n, n, (*smallints)[0]);
  spinz->add_real_block(     -w, 2*n, 2*n, n, n, (*smallints)[0]);
  spinz->add_real_block(      w, 3*n, 3*n, n, n, (*smallints)[0]);

  // y^y contributions
  spinx->add_real_block(     -w, 2*n, 3*n, n, n, (*smallints)[1]);
  spinx->add_real_block(     -w, 3*n, 2*n, n, n, (*smallints)[1]);
  spiny->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[1]);
  spiny->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[1]);
  spinz->add_real_block(     -w, 2*n, 2*n, n, n, (*smallints)[1]);
  spinz->add_real_block(      w, 3*n, 3*n, n, n, (*smallints)[1]);

  // z^z contributions
  spinx->add_real_block(     -w, 2*n, 3*n, n, n, (*smallints)[2]);
  spinx->add_real_block(     -w, 3*n, 2*n, n, n, (*smallints)[2]);
  spiny->add_real_block( imag*w, 2*n, 3*n, n, n, (*smallints)[2]);
  spiny->add_real_block(-imag*w, 3*n, 2*n, n, n, (*smallints)[2]);
  spinz->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[2]);
  spinz->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[2]);

  // x^y contributions
  spinx->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[3]);
  spinx->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[3]);
  spiny->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[3]);
  spiny->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[3]);
  //spinz->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[3]);
  //spinz->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[3]);

  // y^x contributions
  spinx->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[6]);
  spinx->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[6]);
  spiny->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[6]);
  spiny->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[6]);
  //spinz->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[6]);
  //spinz->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[6]);

  // y^z contributions
  //spinx->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[4]);
  //spinx->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[4]);
  spiny->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[4]);
  spiny->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[4]);
  spinz->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[4]);
  spinz->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[4]);

  // z^y contributions
  //spinx->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[7]);
  //spinx->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[7]);
  spiny->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[7]);
  spiny->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[7]);
  spinz->add_real_block(-imag*w, 2*n, 3*n, n, n, (*smallints)[7]);
  spinz->add_real_block( imag*w, 3*n, 2*n, n, n, (*smallints)[7]);

  // z^x contributions
  spinx->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[5]);
  spinx->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[5]);
  //spiny->add_real_block(-imag*w, 2*n, 2*n, n, n, (*smallints)[5]);
  //spiny->add_real_block(-imag*w, 3*n, 3*n, n, n, (*smallints)[5]);
  spinz->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[5]);
  spinz->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[5]);

  // x^z contributions
  spinx->add_real_block(      w, 2*n, 2*n, n, n, (*smallints)[8]);
  spinx->add_real_block(     -w, 3*n, 3*n, n, n, (*smallints)[8]);
  //spiny->add_real_block( imag*w, 2*n, 2*n, n, n, (*smallints)[8]);
  //spiny->add_real_block( imag*w, 3*n, 3*n, n, n, (*smallints)[8]);
  spinz->add_real_block(      w, 2*n, 3*n, n, n, (*smallints)[8]);
  spinz->add_real_block(      w, 3*n, 2*n, n, n, (*smallints)[8]);


  const int ncol = 2*(ncore_+norb_);
  // TODO Probably we can ignore the closed part, but for now include it
  shared_ptr<const ZMatrix> closed_aodensity = jop_->coeff_input()->form_density_rhf(2*ncore_, 0, 1.0);
  const ZMatView active_coeff = jop_->coeff_input()->slice(2*ncore_, ncol);

  // S value of spin manifold to be mapped
  cout << endl << endl;
  const int nspin = idata_->get<int>("aniso_spin", states_.size()-1);
  const int nspin1 = nspin + 1;
  cout << "    Modeling Pseudospin Hamiltonian for S = " << nspin/2 << (nspin % 2 == 0 ? "" : " 1/2") << endl;

  // By default, just use the ground states
  vector<int> aniso_state;
  aniso_state.resize(nspin1);
  for (int i=0; i!=nspin1; ++i)
    aniso_state[i] = i;

  // aniso_state can be used to request mapping excited states instead
  const shared_ptr<const PTree> exstates = idata_->get_child_optional("aniso_state");
  if (exstates) {
    aniso_state = {};
    for (auto& i : *exstates)
      aniso_state.push_back(lexical_cast<int>(i->data()) - 1);
    if (aniso_state.size() != nspin1)
      throw runtime_error("Aniso:  Wrong number of states requested for this S value (should be " + to_string(nspin1) + ")");
    for (int i=0; i!=nspin1; ++i)
      if (aniso_state[i] < 0 || aniso_state[i] >= nstate_)
        throw runtime_error("Aniso:  Invalid state requested (should be between 1 and " + to_string(nstate_) + ")");
    cout << "    For the following states:  ";
    for (int i=0; i!=nspin1; ++i)
      cout << aniso_state[i] << "  ";
    cout << endl;
  } else {
    cout << "    For the ground spin-manifold" << endl;
  }
  cout << endl;

  // Compute spin matrices in the basis of ZFCI Hamiltonian eigenstates
  array<shared_ptr<ZMatrix>,3> spinmat;
  for (int i=0; i!=3; ++i) {
    spinmat[i] = make_shared<ZMatrix>(nspin1, nspin1);
  }
  for (int i=0; i!=nspin1; ++i) {
    for (int j=0; j!=nspin1; ++j) {
      shared_ptr<Kramers<2,ZRDM<1>>> temprdm = rdm1(aniso_state[i], aniso_state[j]);
      if (!temprdm->exist({1,0})) {
        cout << " * Need to generate an off-diagonal rdm of zeroes." << endl;
        temprdm->add({1,0}, temprdm->at({0,0})->clone());
      }
      shared_ptr<const ZRDM<1>> tmp = expand_kramers<1,complex<double>>(temprdm, norb_);
      auto rdmmat = make_shared<ZMatrix>(norb_*2, norb_*2);
      copy_n(tmp->data(), tmp->size(), rdmmat->data());

      ZMatrix modensity (2*norb_, 2*norb_);
      modensity.copy_block(0, 0, 2*norb_, 2*norb_, rdmmat);

      ZMatrix aodensity = (active_coeff * modensity ^ active_coeff) + *closed_aodensity;

      spinmat[0]->element(i,j) = aodensity.dot_product(*spinx);
      spinmat[1]->element(i,j) = aodensity.dot_product(*spiny);
      spinmat[2]->element(i,j) = aodensity.dot_product(*spinz);
    }
  }

  // Diagonalize S_z to get pseudospin eigenstates as combinations of ZFCI Hamiltonian eigenstates
  auto transform = spinmat[2]->copy();
  VectorB zeig(nspin1);
  transform->diagonalize(zeig);

  { // Reorder eigenvectors so positive M_s come first
    shared_ptr<ZMatrix> tempm = transform->clone();
    VectorB tempv = *zeig.clone();
    for (int i=0; i!=nspin1; ++i) {
      tempv[i] = zeig[nspin-i];
      tempm->copy_block(0, i, nspin1, 1, transform->slice(nspin-i, nspin-i+1));
    }
    transform = tempm;
    zeig = tempv;
  }

  for (int i=0; i!=nspin1; ++i)
    cout << "    Spin-z eigenvalue " << i+1 << " = " << zeig[i] << endl;

  // We will subtract out average energy so the pseudospin Hamiltonian is traceless
  complex<double> energy_avg = 0;
  for (int i=0; i!=nspin1; ++i)
    energy_avg += energy_[aniso_state[i]];
  energy_avg /= nspin1;

  // Now build up the numerical pseudospin Hamiltonian!
  auto energies = make_shared<ZMatrix>(nspin1,nspin1);
  for (int i=0; i!=nspin1; ++i) {
    energies->element(i,i) = energy_[aniso_state[i]] - energy_avg;
  }
  auto spinham = make_shared<ZMatrix>(*transform % *energies * *transform);

  auto diagspinx = make_shared<ZMatrix>(*transform % *spinmat[0] * *transform);
  auto diagspiny = make_shared<ZMatrix>(*transform % *spinmat[1] * *transform);
  auto diagspinz = make_shared<ZMatrix>(*transform % *spinmat[2] * *transform);

  cout << endl;
  spinham->print("Pseudospin Hamiltonian!");
/*
  diagspinx->print("Spin matrix - x-component");
  diagspiny->print("Spin matrix - y-component");
  diagspinz->print("Spin matrix - z-component");

  cout << endl;
  energies->print("Hamiltonian over the spin manifold, in ZFCI eigenstate basis");
  spinmat[0]->print("Spin matrix (in ZFCI eigenstate basis) - x-component");
  spinmat[1]->print("Spin matrix (in ZFCI eigenstate basis) - y-component");
  spinmat[2]->print("Spin matrix (in ZFCI eigenstate basis) - z-component");
*/
  cout << endl;

  /**  Part 2: Build up symbolic pseudospin Hamiltonian  **/

  // S_x, S_y, and S_z operators in pseudospin basis
  array<shared_ptr<ZMatrix>,3> pspinmat;
  for (int i=0; i!=3; ++i)
    pspinmat[i] = make_shared<ZMatrix>(nspin1, nspin1);

  auto spin_plus = make_shared<ZMatrix>(nspin1, nspin1);
  auto spin_minus = make_shared<ZMatrix>(nspin1, nspin1);
  const double sval = nspin/2.0;
  const double ssp1 = sval*(sval+1.0);

  for (int i=0; i!=nspin1; ++i) {
    const double ml1 = sval - i;
    pspinmat[2]->element(i,i) = ml1;
    if (i < nspin) {
      spin_plus->element(i,i+1) = std::sqrt(ssp1 - ml1*(ml1-1.0));
    }
    if (i > 0) {
      spin_minus->element(i,i-1) = std::sqrt(ssp1 - ml1*(ml1+1.0));
    }
  }

  pspinmat[0]->add_block( 0.5, 0, 0, nspin1, nspin1, spin_plus);
  pspinmat[0]->add_block( 0.5, 0, 0, nspin1, nspin1, spin_minus);
  pspinmat[1]->add_block( complex<double>( 0.0, -0.5), 0, 0, nspin1, nspin1, spin_plus);
  pspinmat[1]->add_block( complex<double>( 0.0,  0.5), 0, 0, nspin1, nspin1, spin_minus);

  // Transformation matrix to build pseudospin Hamiltonian from D-tensor
  // Rows correspond to pairs of pseudospins (SS, S-1S, S-2S...)
  // Columns correspond to elements of the D-tensor (Dxx, Dyx, Dzx, Dxy...)
  // Note that we look over the first indices before the second, so we can copy data between vectors and matrices and have it come out in the right order
  auto d2h = make_shared<ZMatrix>(nspin1*nspin1,9);
  for (int i=0; i!=3; ++i) {
    for (int j=0; j!=3; ++j) {
      ZMatrix temp = *pspinmat[i] * *pspinmat[j];
      d2h->copy_block(0, 3*j+i, nspin1*nspin1, 1, temp.data());
    }
  }

  /**  Part 3: Extract D-tensor from the numerical pseudospin Hamiltonian **/

  // h2d is the left-inverse of d2h
  // It converts from the pseudospin Hamiltonian to the D-tensor.
  ZMatrix d2h_sqinv = *d2h % *d2h;
  d2h_sqinv.inverse();
  ZMatrix h2d = d2h_sqinv ^ *d2h;

  auto spinham_vec = make_shared<ZMatrix>(nspin1*nspin1,1);
  spinham_vec->copy_block(0, 0, nspin1*nspin1, 1, spinham->element_ptr(0,0));
  ZMatrix Dtensor_vec = h2d * *spinham_vec;
  auto Dtensor = make_shared<ZMatrix>(3,3);
  Dtensor->copy_block(0, 0, 3, 3, Dtensor_vec.element_ptr(0,0));

#ifndef NDEBUG
  // Watch for numerical instability - might happen if higher spins have linear dependency, for example
  {
    ZMatrix test = h2d * *d2h;
    ZMatrix unit = *test.clone();
    unit.unit();
    const double error = (test-unit).rms();
    if (error > 1.0e-8)
      cout << "  **  RMS Error in left-inverse computation for extraction of D-tensor: " << scientific << setprecision(4) << error << endl;
  }
#endif

  Dtensor->print("D tensor");
  VectorB Ddiag(3);
  Dtensor->diagonalize(Ddiag);
  for (int i=0; i!=3; ++i)
    cout << "Diagonalized D-tensor value " << i << " = " << Ddiag[i] << endl;

  // Compute Davg so that it works even if D is not traceless (which shouldn't happen on accident)
  const double Davg = 1.0 / 3.0 * (Ddiag[0] + Ddiag[1] + Ddiag[2]);

  int jmax = 0;
  const array<int,3> fwd = {{ 1, 2, 0 }};
  const array<int,3> bck = {{ 2, 0, 1 }};
  if (std::abs(Ddiag[1]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 1;
  if (std::abs(Ddiag[2]-Davg) > std::abs(Ddiag[jmax]-Davg)) jmax = 2;
  const double Dval = Ddiag[jmax] - 0.5*(Ddiag[fwd[jmax]] + Ddiag[bck[jmax]]);
  const double Eval = 0.5*(Ddiag[fwd[jmax]] - Ddiag[bck[jmax]]);
  cout << " ** D = " << Dval << " E_h, or " << Dval * au2wavenumber__ << " cm-1" << endl;
  cout << " ** E = " << std::abs(Eval) << " E_h, or " << std::abs(Eval * au2wavenumber__) << " cm-1" << endl;

  ZMatrix check_spinham_vec = *d2h * Dtensor_vec;
  auto checkham = make_shared<ZMatrix>(nspin1,nspin1);
  checkham->copy_block(0, 0, nspin1, nspin1, check_spinham_vec.element_ptr(0,0));
  checkham->print("Pseudospin Hamiltonian, recomputed", 30);
  cout << "  Error in recomputation of spin Hamiltonian from D = " << (*checkham - *spinham).rms() << endl << endl;

  VectorB shenergies(nspin1);
  checkham->diagonalize(shenergies);

  if (nspin == 2) {
    cout << "  ** Relative energies expected from diagonalized D parameters: " << endl;
    if (Dval > 0.0) {
      cout << "     2  " << Dval + std::abs(Eval) << " E_h  =  " << (Dval + std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     1  " << Dval - std::abs(Eval) << " E_h  =  " << (Dval - std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
    } else {
      cout << "     2  " << -Dval + 0.5*std::abs(Eval) << " E_h  =  " << (-Dval + 0.5*std::abs(Eval))*au2wavenumber__ << " cm-1" << endl;
      cout << "     1  " << std::abs(Eval) << " E_h  =  " << std::abs(Eval)*au2wavenumber__ << " cm-1" << endl;
      cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
    }
  } else if (nspin == 3) {
    cout << "  ** Relative energies expected from diagonalized D parameters: " << endl;
    const double energy32 = 2.0*std::sqrt(Dval*Dval + 3.0*Eval*Eval);
    cout << "     3  " << energy32 << " E_h  =  " << energy32*au2wavenumber__ << " cm-1" << endl;
    cout << "     2  " << energy32 << " E_h  =  " << energy32*au2wavenumber__ << " cm-1" << endl;
    cout << "     1  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl;
    cout << "     0  " << 0.0 << " E_h  =  " << 0.0 << " cm-1" << endl << endl;
  }

  cout << "  ** Relative energies expected from the recomputed Pseudospin Hamiltonian: " << endl;
  for (int i=nspin; i>=0; --i)
    cout << "     " << i << "  " << shenergies[i] - shenergies[0] << " E_h  =  " << (shenergies[i] - shenergies[0])*au2wavenumber__ << " cm-1" << endl;
  cout << endl;

  cout << "  ** Relative energies observed by relativistic configuration interaction: " << endl;
  for (int i=nspin; i>=0; --i)
    cout << "     " << i << "  " << energy_[aniso_state[i]] - energy_[aniso_state[0]] << " E_h  =  " << (energy_[aniso_state[i]] - energy_[aniso_state[0]])*au2wavenumber__ << " cm-1" << endl;
  cout << endl;

}
