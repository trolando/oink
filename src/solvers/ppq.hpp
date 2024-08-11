/*
 * Copyright 2021 Tom van Dijk, University of Twente
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

#ifndef PPQ_HPP
#define PPQ_HPP

#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * Implementation of ZLKQ but accelerated with priority promotion.
 */

class PPQSolver : public Solver
{
public:
    PPQSolver(Oink& oink, Game& game);
    virtual ~PPQSolver();

    virtual void run();

protected:
    unsigned long long iterations = 0;

    uintqueue Q;
    int *str;

    bitset G;
    bitset R;
    bitset L;
    bitset U;
    int *r, *u;

    bitset W0;
    bitset W1;

    inline void attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y);
    bool maybePromote(int pr, int pl, bitset &R, bitset &L);
    void solve(int vtop, int pr, int prc, int pe, int po);
};

}

#endif
