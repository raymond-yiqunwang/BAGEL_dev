//
// BAGEL - Parallel electron correlation program.
// Filename: lattice.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#ifndef __BAGEL_SRC_PERIODIC_LATTICE_H
#define __BAGEL_SRC_PERIODIC_LATTICE_H

#include <src/periodic/pdfdist.h>
#include <src/periodic/tree.h>

namespace bagel {

class Lattice {
  protected:
    std::shared_ptr<const Geometry> primitive_cell_;
    int ndim_;
    int k_parameter_;
    int extent_;
    int num_lattice_vectors_;
    int num_lattice_pts_;
    // real lattice vectors g
    std::vector<std::array<double, 3>> lattice_vectors_;
    std::map<int, std::array<int, 3>> lattice_map_;

    double nuclear_repulsion_;
    double compute_nuclear_repulsion() const;
    double thresh_;

    // ``volume'' of a unit cell
    double volume_;
    // primitive reciprocal lattice vectors
    std::vector<std::array<double, 3>> primitive_kvectors_;
    // recriprocal lattice vectors k
    std::vector<std::array<double, 3>> lattice_kvectors_;
    int num_lattice_kvectors_;
    int gamma_point_;

    double dot(const std::array<double, 3>& b, const std::array<double, 3>& c) const;
    std::array<double, 3> cross(const std::array<double, 3>& b, const std::array<double, 3>& c, double s = 1.0) const;

    int nele_;

    void init();
    std::shared_ptr<const Tree> fmmtree_;
    std::shared_ptr<const Geometry> supergeom_;
    std::vector<double> schwarz_;
    double schwarz_thresh_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & ndim_ & extent_ & num_lattice_vectors_ & num_lattice_pts_ & primitive_cell_ & lattice_vectors_
         & nuclear_repulsion_ & volume_ & primitive_kvectors_ & lattice_kvectors_ & k_parameter_ & num_lattice_kvectors_ & nele_;
    }

  public:
    Lattice() { }
    Lattice(const std::shared_ptr<const Geometry> g, const int k_parameter, const int extent = 0, const bool dofmm = false,
            const std::tuple<int, int, bool, bool, double>& fmmp = std::tuple<int, int, bool, bool, double>());
    virtual ~Lattice() { }

    int ndim() const { return ndim_; }
    int extent() const {return extent_; }
    int num_lattice_pts() const { return num_lattice_pts_; }
    int num_lattice_vectors() const { return num_lattice_vectors_; }
    int num_lattice_kvectors() const { return num_lattice_kvectors_; }

    std::shared_ptr<const Geometry> primitive_cell() const { return primitive_cell_; }
    const std::vector<std::array<double, 3>>& lattice_vectors() const { return lattice_vectors_; }
    const std::array<double, 3>& lattice_vectors(const int i) const { return lattice_vectors_[i]; }

    double nuclear_repulsion() const { return nuclear_repulsion_; };
    double volume() const { return volume_; }
    int nele() const { return nele_; }

    const std::vector<std::array<double, 3>>& primitive_kvectors() const { return primitive_kvectors_; }
    const std::array<double, 3>& primitive_kvectors(const int i) const { return primitive_kvectors_[i]; }
    const std::vector<std::array<double, 3>>& lattice_kvectors() const { return lattice_kvectors_; }
    const std::array<double, 3>& lattice_kvectors(const int i) const { return lattice_kvectors_[i]; }
    void generate_kpoints();
    int gamma_point() const { return gamma_point_; }
    int central_cell() const;

    int find_lattice_vector(const int i, const int j) const;
    void print_primitive_vectors() const;
    void print_primitive_kvectors() const;
    void print_lattice_vectors() const;
    void print_lattice_kvectors() const;
    void print_lattice_coordinates() const; // write .XYZ file
    void print_atoms() const;
    std::array<double, 3> centre() const { return move(primitive_cell_->charge_center()); }
    double centre(const int i) const { return primitive_cell_->charge_center()[i]; }
    std::array<double, 3> cell_centre(const int icell) const;
    std::shared_ptr<const Geometry> supergeom() const { return supergeom_; }
    double thresh() const { return thresh_; }

    // density fitting
    std::shared_ptr<const PDFDist> form_df() const;
    // PFMM
    std::shared_ptr<const Tree> fmmtree() const { return fmmtree_; }
    const std::vector<double>& schwarz() const { return schwarz_; }
    double schwarz_thresh() const { return schwarz_thresh_; }
    void build_tree(const std::tuple<int, int, bool, bool, double>& fmm_param);
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Lattice)

#endif

