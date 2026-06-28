/*
 * Copyright 2017-2024 Tom van Dijk
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

#ifndef FPIM_HPP
#define FPIM_HPP

#include "oink/solver.hpp"

namespace pg {

/**
 * Multi-strategy variant of the (sequential) distraction fixpoint iteration
 * solver FPI. Identical algorithm, but instead of selecting a single winning
 * move per vertex it records *all* moves that are winning at the point the
 * vertex is decided, in the game's MultiStrategy. A representative single
 * strategy is still reported so the standard solution stays valid.
 */
class FPIMSolver : public Solver
{
public:
    FPIMSolver(Oink& oink, Game& game);
    virtual ~FPIMSolver();

    virtual void run();

    int updateBlock(int i, int n);
    void freezeThawReset(int i, int n, int p);
    void runSeq(void);

private:
    unsigned long long iterations = 0;
    int *frozen;
    bitset parity;
    bitset distraction;
};

}

#endif
