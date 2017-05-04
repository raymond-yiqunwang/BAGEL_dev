//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: process.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <cassert>
#include <src/util/parallel/process.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

Process::Process() : print_level_(3), muted_(false) {
  if (mpi__->world_rank() != 0) {
    cout_orig = cout.rdbuf();
    cout.rdbuf(ss_.rdbuf());
    muted_ = true;
  }
}


Process::~Process() {
  if (muted_)
    cout.rdbuf(cout_orig);
}


void Process::cout_on() {
  if (mpi__->world_rank() != 0) {
    assert(muted_);
    cout.rdbuf(cout_orig);
    muted_ = false;
  }
}


void Process::cout_off() {
  if (mpi__->world_rank() != 0) {
    assert(!muted_);
    cout.rdbuf(ss_.rdbuf());
    muted_ = true;
  }
}
