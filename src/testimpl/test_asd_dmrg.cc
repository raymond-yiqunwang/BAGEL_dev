//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: test_asd_dmrg.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/asd/dmrg/rasd.h>


double asd_dmrg_energy(std::string inp) {
  Muffle hide_cout(inp + ".testout");

  // a bit ugly to hardwire an input file, but anyway...
  std::stringstream ss; ss << location__ << inp << ".json";
  auto idata = std::make_shared<const PTree>(ss.str());
  auto keys = idata->get_child("bagel");
  std::shared_ptr<Geometry> geom;

  std::shared_ptr<const Reference> ref;
  std::shared_ptr<MultiSite> multisite;

  std::map<std::string, std::shared_ptr<const Reference>> saved;

  for (auto& itree : *keys) {
    const std::string method = to_lower(itree->get<std::string>("title", ""));

    if (method == "molecule") {
      geom = std::make_shared<Geometry>(itree);

    } else if (method == "hf") {
      auto scf = std::make_shared<RHF>(itree, geom, ref);
      scf->compute();
      ref = scf->conv_to_ref();
    } else if (method == "multisite") {
      const int nsites = itree->get<int>("nsites");
      multisite = std::make_shared<MultiSite>(itree, ref, nsites);
      multisite->compute();
      ref = multisite->conv_to_ref();
      *geom = *ref->geom();
    } else if (method == "asd_dmrg") {
        if (!multisite)
          throw std::runtime_error("multisite must be called before asd_dmrg");
        auto asd = std::make_shared<RASD>(itree, multisite);
        asd->compute();
        return asd->energies(0);
    }
  }
  assert(false);
  return 0.0;
}

BOOST_AUTO_TEST_SUITE(TEST_ASD_DMRG)

BOOST_AUTO_TEST_CASE(RASD_DMRG) {
    BOOST_CHECK(compare(asd_dmrg_energy("he3_svp_asd-dmrg"), -8.56338715, 1.0e-8));
}

BOOST_AUTO_TEST_SUITE_END()
