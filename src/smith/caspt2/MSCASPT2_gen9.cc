//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen9.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task400::Task400(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c1}}, in, t[0], range));
}

Task401::Task401(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task402::Task402(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task403::Task403(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task404::Task404(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c1}}, in, t[0], range));
}

Task405::Task405(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task406::Task406(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task407::Task407(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task408::Task408(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x1}}, in, t[0], range));
}

Task409::Task409(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x1}}, in, t[0], range));
}

Task410::Task410(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task411::Task411(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c4}}, in, t[0], range));
}

Task412::Task412(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c4}}, in, t[0], range));
}

Task413::Task413(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task414::Task414(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a3, a1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, a1, c2}}, in, t[0], range));
}

Task415::Task415(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a3, c4, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, c4, a1}}, in, t[0], range));
}

Task416::Task416(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a3, a1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, a1, c2}}, in, t[0], range));
}

Task417::Task417(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a4, c2, a3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a4, c2, a3}}, in, t[0], range));
}

Task418::Task418(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a3, a1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, a1, c2}}, in, t[0], range));
}

Task419::Task419(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a4, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a4, c2, a1}}, in, t[0], range));
}

Task420::Task420(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[5]->is_local(x5, x2, x4, x3, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x2, x4, x3, x1, x0}}, in, out, range));
}

Task421::Task421(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              if (t[0]->is_local(x1, x0, x2, x5, x4, x3))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, x2, x5, x4, x3}}, in, t[0], range));
}

Task422::Task422(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, x2}}, in, t[0], range));
}

Task423::Task423(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[5]->is_local(x5, x0, x3, x4, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x3, x4, x2, x1}}, in, out, range));
}

Task424::Task424(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x1, x0}}, in, t[0], range));
}

Task425::Task425(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, a1, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, a1, x4, x3}}, in, t[0], range));
}

Task426::Task426(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[5]->is_local(x5, x4, x3, x0, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x0, x2, x1}}, in, out, range));
}

Task427::Task427(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x1, x0}}, in, t[0], range));
}

Task428::Task428(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x5, x4, x3}}, in, t[0], range));
}

Task429::Task429(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x7 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x6 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& x4 : *range[1])
                for (auto& x3 : *range[1])
                  if (t[5]->is_local(x7, x0, x6, x5, x2, x1))
                    subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x5, x2, x1, x4, x3}}, in, out, range));
}

Task430::Task430(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x2, x1, x0, x7, x6, x5))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x7, x6, x5}}, in, t[0], range));
}

Task431::Task431(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[5]->is_local(x5, x0, x4, x3, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, out, range));
}

Task432::Task432(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x1, x0}}, in, t[0], range));
}

Task433::Task433(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, x4, x3, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x4, x3, a1}}, in, t[0], range));
}

Task434::Task434(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, x4, x3, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x4, x3, a1}}, in, t[0], range));
}

Task435::Task435(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x1, x0}}, in, t[0], range));
}

Task436::Task436(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a1, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a1, x0, x2}}, in, t[0], range));
}

Task437::Task437(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, out, range));
}

Task438::Task438(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task439::Task439(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x3 : *range[1])
      if (t[0]->is_local(x3, a1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, a1}}, in, t[0], range));
}

Task440::Task440(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a3, c2, a1}}, in, t[0], range));
}

Task441::Task441(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task442::Task442(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a3 : *range[2])
      if (t[0]->is_local(a3, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, x0}}, in, t[0], range));
}

Task443::Task443(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task444::Task444(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, x0}}, in, t[0], range));
}

Task445::Task445(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task446::Task446(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, a1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a1, x0, x1}}, in, t[0], range));
}

Task447::Task447(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task448::Task448(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& a1 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a1, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a1, a2, x2}}, in, t[0], range));
}

Task449::Task449(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& a1 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a1, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a1, a2, x2}}, in, t[0], range));
}

#endif
