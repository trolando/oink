/*
 * Copyright 2020-2024 Tom van Dijk
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

#ifndef FPJM_HPP
#define FPJM_HPP

#include "oink/solver.hpp"

namespace pg {

/**
 * Multi-strategy variant of the fixpoint-iteration-with-justifications solver
 * FPJ (sequential, non-greedy). Identical algorithm, but instead of recording a
 * single winning move per vertex it records *all* moves that are winning at the
 * point the vertex is justified, in the solution's multi-strategy. A justified
 * vertex is invalidated only when its *last* recorded strategy edge is pruned.
 * Won vertices report the -1 sentinel single strategy; their moves are read from
 * the multi-strategy.
 */
class FPJMSolver : public Solver
{
public:
    FPJMSolver(Oink& oink, Game& game);
    virtual ~FPJMSolver();

    unsigned long long iterations = 0;

    void runSeq(void);

    virtual void run();
};

}

#endif
