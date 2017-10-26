//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_rotfile.cc
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

#include <src/asd/dmrg/orbopt/asd_dmrg_rotfile.h>

using namespace std;
using namespace bagel;

// constructors
ASD_DMRG_RotFile::ASD_DMRG_RotFile(const int iclo, const int iact, const int ivirt)
  : nclosed_(iclo), nact_(iact), nvirt_(ivirt), size_(iclo*iact+ivirt*iact+ivirt*iclo/*TODO+iact*iact*/), data_(new double[size_]) {
  zero();  
}


ASD_DMRG_RotFile::ASD_DMRG_RotFile(const ASD_DMRG_RotFile& o)
  : nclosed_(o.nclosed_), nact_(o.nact_), nvirt_(o.nvirt_), size_(o.size_), data_(new double[o.size_]) {
  *this = o;
}


ASD_DMRG_RotFile::ASD_DMRG_RotFile(shared_ptr<const ASD_DMRG_RotFile> o)
  : nclosed_(o->nclosed_), nact_(o->nact_), nvirt_(o->nvirt_), size_(o->size_), data_(new double[o->size_]) {
  *this = *o;
}


// overloaded operators
ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator=(const ASD_DMRG_RotFile& o) {
  copy_n(o.data(), size(), data());
  return *this;
}


ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator+=(const ASD_DMRG_RotFile& o) {
  ax_plus_y(1.0, o);
  return *this;
}


ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator-=(const ASD_DMRG_RotFile& o) {
  ax_plus_y(-1.0, o);
  return *this;
}


ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator*=(const ASD_DMRG_RotFile& o) {
  for (int i = 0; i != size(); ++i)
    data(i) *= o.data(i);
  return *this;
}


ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator/=(const ASD_DMRG_RotFile& o) {
  for (int i = 0; i != size(); ++i)
    data(i) /= o.data(i);
  return *this;
}


ASD_DMRG_RotFile ASD_DMRG_RotFile::operator+(const ASD_DMRG_RotFile& o) const {
  ASD_DMRG_RotFile out(*this);
  out.ax_plus_y(1.0, o);
  return out;
}


ASD_DMRG_RotFile ASD_DMRG_RotFile::operator-(const ASD_DMRG_RotFile& o) const {
  ASD_DMRG_RotFile out(*this);
  out.ax_plus_y(-1.0, o);
  return out;
}


ASD_DMRG_RotFile ASD_DMRG_RotFile::operator*(const ASD_DMRG_RotFile& o) const {
  ASD_DMRG_RotFile out(*this);
  return out *= o;
}


ASD_DMRG_RotFile ASD_DMRG_RotFile::operator/(const ASD_DMRG_RotFile& o) const {
  ASD_DMRG_RotFile out(*this);
  return out /= o;
}


ASD_DMRG_RotFile& ASD_DMRG_RotFile::operator*=(const double a) {
  scale(a);
  return *this;
}


double ASD_DMRG_RotFile::normalize() {
  const double scal = 1.0 / this->norm();
  scale(scal);
  return 1.0/scal;
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_RotFile::clone() const {
  return make_shared<ASD_DMRG_RotFile>(nclosed_, nact_, nvirt_);
}


shared_ptr<ASD_DMRG_RotFile> ASD_DMRG_RotFile::copy() const {
  return make_shared<ASD_DMRG_RotFile>(*this);
}


void ASD_DMRG_RotFile::ax_plus_y_ca(const double a, const MatView mat) {
  assert(mat.ndim() == nclosed_ && mat.mdim() == nact_);
  blas::ax_plus_y_n(a, mat.data(), nclosed_*nact_, ptr_ca());
}


void ASD_DMRG_RotFile::ax_plus_y_va(const double a, const MatView mat) {
  assert(mat.ndim() == nvirt_ && mat.mdim() == nact_);
  blas::ax_plus_y_n(a, mat.data(), nvirt_*nact_, ptr_va());
}


void ASD_DMRG_RotFile::ax_plus_y_vc(const double a, const MatView mat) {
  assert(mat.ndim() == nvirt_ && mat.mdim() == nclosed_);
  blas::ax_plus_y_n(a, mat.data(), nvirt_*nclosed_, ptr_vc());
}


/*
void ASD_DMRG_RotFile::ax_plus_y_aa(const double a, const MatView mat) {
  assert(mat.ndim() == nact_ && mat.mdim() == nact_);
  blas::ax_plus_y_n(a, mat.data(), nact_*nact_, ptr_aa());
}
*/


shared_ptr<Matrix> ASD_DMRG_RotFile::ca_mat() const {
  auto out = make_shared<Matrix>(nclosed_, nact_);
  copy_n(ptr_ca(), nclosed_*nact_, out->data());
  return out;
}


shared_ptr<Matrix> ASD_DMRG_RotFile::va_mat() const {
  auto out = make_shared<Matrix>(nvirt_, nact_);
  copy_n(ptr_va(), nvirt_*nact_, out->data());
  return out;
}


shared_ptr<Matrix> ASD_DMRG_RotFile::vc_mat() const {
  auto out = make_shared<Matrix>(nvirt_, nclosed_);
  copy_n(ptr_vc(), nvirt_*nclosed_, out->data());
  return out;
}


void ASD_DMRG_RotFile::print(const string input) const {
  cout << " +++++ " + input + " +++++" <<endl;
  if (nclosed_ && nact_) {
    cout << " printing closed-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nclosed_; ++j)
        cout << setw(10) << setprecision(6) << ele_ca(j, i);
      cout << endl;
    }
  }
  if (nvirt_ && nact_) {
    cout << " printing virtual-active block" << endl;
    for (int i = 0; i != nact_; ++i) {
      for (int j = 0; j != nvirt_; ++j)
        cout << setw(10) << setprecision(6) << ele_va(j, i);
      cout << endl;
    }
  }
  if (nvirt_ && nclosed_) {
    cout << " printing virtual-closed block" << endl;
    for (int i = 0; i != nclosed_; ++i) {
      for (int j = 0; j != nvirt_; ++j)
        cout << setw(10) << setprecision(6) << ele_vc(j, i);
      cout << endl;
    }
  }
/*
  if (nact_) {
    cout << " printing active-active block" << endl;
    cout << " not implemented yet.." << endl;
  }
*/
}


