/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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

#ifndef ZLK_HPP
#define ZLK_HPP

#include <queue>
#include <lace.h>

#include "oink/solver.hpp"

namespace pg {

class ZLKSolver : public Solver
{
public:
    ZLKSolver(Oink& oink, Game& game);
    virtual ~ZLKSolver();

    virtual void run();

    int iterations;

    int *inverse;
    int max_prio;

    int *region;
    int *winning;
    int *strategy;

    bool to_inversion = true;
    bool only_recompute_when_attracted = true;

    uintqueue Q;

    int attractExt(int i, int r, std::vector<int> *R);
    int attractLosing(int i, int r, std::vector<int> *S, std::vector<int> *R);
    void attractParT(WorkerP*, Task*, int pl, int cur, int r);
    int attractPar(WorkerP*, Task*, int i, int r, std::vector<int>* R);
    // friend void updateOutcount_WORK(WorkerP*, Task*, int, int, ZLKSolver*);
};

class UnoptimizedZLKSolver : public ZLKSolver
{
public:
    UnoptimizedZLKSolver(Oink& oink, Game& game) : ZLKSolver(oink, game) { to_inversion = false; only_recompute_when_attracted = false; }
    virtual ~UnoptimizedZLKSolver() { }
};

}

#endif 
