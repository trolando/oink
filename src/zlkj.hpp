/*
 * Copyright 2017-2020 Tom van Dijk, Johannes Kepler University Linz
 * Copyright 2019-2020 Ruben Lapauw, KU Leuven
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef ZLKJ_HPP
#define ZLKJ_HPP

#include <queue>
#include <vector>
#include <unordered_set>
#include "solver.hpp"
#include "lace.h"
#include "uintqueue.hpp"

namespace pg {

class ZLKJSolver : public Solver
{
public:
    ZLKJSolver(Oink *oink, Game *game);
    virtual ~ZLKJSolver();

    virtual void run();

    int iterations;

protected:
    int *inverse;
    int max_prio;

    int *region;
    int *winning;
    int *strategy;
    std::vector<std::unordered_set<int> > just;


    bool to_inversion = true;
    bool only_recompute_when_attracted = true;

    uintqueue Q;
    uintqueue D;

    int attractExt(int i, int r, std::vector<int> *R);
    std::pair<int, int> attractLosing(int i, int r, std::vector<int> *S, int next_r);
};
}
#endif
