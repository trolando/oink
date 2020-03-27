/*
 * Copyright 2019 Tom van Dijk, University of Twente
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

#ifndef ZLKQ_HPP
#define ZLKQ_HPP

#include <vector>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

/**
 * This is an implementation of Pawel Parys's Zielonka variant
 * as presented by himself on his website.
 * https://www.mimuw.edu.pl/~parys/publications/2018-parity-algorithm.pdf
 *
 * I modified the algorithm to produce a strategy as well.
 */

class ZLKQSolver : public Solver
{
public:
    ZLKQSolver(Oink *oink, Game *game);
    virtual ~ZLKQSolver();

    virtual void run();

protected:
    unsigned long long iterations = 0;

    uintqueue Q;
    uintqueue Hs;
    int *str;

    bitset W0, W1; // current approximation of winning areas

    inline void attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y);
    void solve(bitset &Subgame, int vtop, int pr, int pe, int po);

    // version of the Liverpools
    void solve2(bitset &Subgame, int vtop, int pr, int pl, int pe, int po);
};

}

#endif
