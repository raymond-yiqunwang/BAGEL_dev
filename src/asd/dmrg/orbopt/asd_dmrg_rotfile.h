//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg_rotfile.h
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

#ifndef __ASD_DMRG_ROTFILE_H
#define __ASD_DMRG_ROTFILE_H

#include <src/util/math/matrix.h>

namespace bagel {

class ASD_DMRG_RotFile {
  protected:
    const int nclosed_;
    const int nact_;
    const int nvirt_;
    const int size_;
    std::unique_ptr<double[]> data_;

  public:
    // constructors
    ASD_DMRG_RotFile(const int iclo, const int iact, const int ivirt);
    ASD_DMRG_RotFile(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile(std::shared_ptr<const ASD_DMRG_RotFile> o);

    // overloaded operators
    ASD_DMRG_RotFile& operator=(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile& operator+=(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile& operator-=(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile& operator*=(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile& operator/=(const ASD_DMRG_RotFile& o);
    ASD_DMRG_RotFile operator+(const ASD_DMRG_RotFile& o) const;
    ASD_DMRG_RotFile operator-(const ASD_DMRG_RotFile& o) const;
    ASD_DMRG_RotFile operator*(const ASD_DMRG_RotFile& o) const;
    ASD_DMRG_RotFile operator/(const ASD_DMRG_RotFile& o) const;
    ASD_DMRG_RotFile& operator*=(const double a);


    // util functions
    std::shared_ptr<ASD_DMRG_RotFile> clone() const;
    std::shared_ptr<ASD_DMRG_RotFile> copy() const;

    int size() const { return size_; }
    void zero() { fill(0.0); }
    void fill(const double& a) { std::fill_n(data(), size_, a); }
    double dot_product(const ASD_DMRG_RotFile& o) const { return blas::dot_product(data(), size_, o.data()); }
    double dot_product(const std::shared_ptr<const ASD_DMRG_RotFile> o) const { return dot_product(*o); }
    void scale(const double& a) { std::for_each(data(), data()+size_, [&a](double& p) { p *= a; }); } // scale function
    double norm() const { return std::sqrt(detail::real(dot_product(*this))); } // returns norm of the vector
    double rms() const { return norm() / std::sqrt(static_cast<double>(size())); }
    
    void ax_plus_y(const double& a, const ASD_DMRG_RotFile& o) { blas::ax_plus_y_n(a, o.data(), size_, data()); }
    void ax_plys_y(const double& a, const std::shared_ptr<const ASD_DMRG_RotFile> o) { ax_plus_y(a, *o); }

    // return data
    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }
    double& data(const size_t i) { return data_[i]; }
    const double& data(const size_t i) const { return data_[i]; }
    double* begin() { return data(); }
    double* end() { return data() + size_; }

    // closed-active block, closed runs first
    double* ptr_ca() { return data(); }
    double& ele_ca(const int ic, const int ia) { return data_[ic + ia*nclosed_]; }
    // virtual-active block, virtual runs first
    double* ptr_va() { return data() + nclosed_*nact_; }
    double& ele_va(const int iv, const int ia) { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    //virtual-closed block, virtual runs first
    double* ptr_vc() { return data() + (nclosed_+nvirt_)*nact_; }
    double& ele_vc(const int iv, const int ic) { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    // active-active block, the first active runs first
    double* ptr_aa() { return data() + (nclosed_+nvirt_)*nact_ + nvirt_*nclosed_; }
    double& ele_aa(const int ia1, const int ia2) { return data_[(nclosed_+nvirt_)*nact_ + nvirt_*nclosed_ + ia1 + ia2*nact_]; }
    // const references and pointers
    const double* ptr_ca() const { return data(); }
    const double* ptr_va() const { return data() + nclosed_*nact_; }
    const double* ptr_vc() const { return data() + (nclosed_+nvirt_)*nact_; }
    const double* ptr_aa() const { return data() + (nclosed_+nvirt_)*nact_ + nvirt_*nclosed_; }
    
    const double& ele_ca(const int ic, const int ia) const { return data_[ic + ia*nclosed_]; }
    const double& ele_va(const int iv, const int ia) const { return data_[nclosed_*nact_ + iv + ia*nvirt_]; }
    const double& ele_vc(const int iv, const int ic) const { return data_[(nclosed_+nvirt_)*nact_ + iv + ic*nvirt_]; }
    const double& ele_aa(const int ia1, const int ia2) const { return data_[(nclosed_+nvirt_)*nact_ + nvirt_*nclosed_ + ia1 + ia2*nact_]; }

    void  ax_plus_y_ca(const double a, const MatView mat);
    void  ax_plus_y_va(const double a, const MatView mat);
    void  ax_plus_y_vc(const double a, const MatView mat);
    void  ax_plus_y_aa(const double a, const MatView mat);

    // unpack to Matrix
    // to be implemented..

    void print(const std::string input = "") const;
    

};

}

#endif
