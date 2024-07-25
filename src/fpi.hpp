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

#ifndef FPI_HPP
#define FPI_HPP

#include <queue>

#include "oink/solver.hpp"
#include "lace.h"

namespace pg {

class FPISolver : public Solver
{
public:
    FPISolver(Oink& oink, Game& game);
    virtual ~FPISolver();

    virtual void run();

    int updateBlock(int i, int n);
    void freezeThawReset(int i, int n, int p);
    void runPar(void);
    void runSeq(void);

    int update_block_rec(WorkerP*, Task*, int i, int n);
    void run_par(WorkerP*, Task*);

private:
    unsigned long long iterations = 0;
    int *frozen;
    int *strategy;
    bitset parity;
    bitset distraction;
};

}

#endif 
