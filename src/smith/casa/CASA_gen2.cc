//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASA_gen2.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/casa/CASA_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASA;

Task50::Task50(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x2, x1, x0}}, in, t[0], range));
}

Task51::Task51(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task52::Task52(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task53::Task53(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task54::Task54(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task55::Task55(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a3, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a3, c1, x3}}, in, t[0], range));
}

Task56::Task56(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task57::Task57(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, c1, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, c1, x4, x3}}, in, t[0], range));
}

Task58::Task58(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task59::Task59(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range, e));
}

Task60::Task60(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task61::Task61(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task62::Task62(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task63::Task63(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x0, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x0, c1, a2}}, in, t[0], range));
}

Task64::Task64(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task65::Task65(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, c3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, c3, x0}}, in, t[0], range));
}

Task66::Task66(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, c3, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, c3, x1}}, in, t[0], range));
}

Task67::Task67(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task68::Task68(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, c1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, c1, x0}}, in, t[0], range));
}

Task69::Task69(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, c1, x1}}, in, t[0], range));
}

Task70::Task70(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task71::Task71(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x1, x0}}, in, t[0], range));
}

Task72::Task72(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x1, x0}}, in, t[0], range));
}

Task73::Task73(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task74::Task74(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c1, x0, x1}}, in, t[0], range));
}

Task75::Task75(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a2, c1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a2, c1, x2}}, in, t[0], range));
}

Task76::Task76(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task77::Task77(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x0}}, in, t[0], range));
}

Task78::Task78(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task79::Task79(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x0}}, in, t[0], range));
}

Task80::Task80(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task81::Task81(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, x0}}, in, t[0], range));
}

Task82::Task82(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task83::Task83(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, x0}}, in, t[0], range));
}

Task84::Task84(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task85::Task85(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x1, x0}}, in, t[0], range));
}

Task86::Task86(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task87::Task87(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, c1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, x0}}, in, t[0], range));
}

Task88::Task88(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, x1}}, in, t[0], range));
}

Task89::Task89(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task90::Task90(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range, e));
}

Task91::Task91(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task92::Task92(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x3}}, in, t[0], range));
}

Task93::Task93(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x0, c1}}, in, t[0], range));
}

Task94::Task94(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x5}}, in, t[0], range));
}

Task95::Task95(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x1, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, x0, a1}}, in, t[0], range));
}

Task96::Task96(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, c2}}, in, t[0], range));
}

Task97::Task97(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, c3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, c3, x0, x1}}, in, t[0], range));
}

Task98::Task98(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, c3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, c3, x0, x1}}, in, t[0], range));
}

Task99::Task99(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, c2}}, in, t[0], range));
}

#endif
