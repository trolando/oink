/*
 * Copyright 2019-2020 Tom van Dijk, University of Twente
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

#include "oink/solver.hpp"

namespace pg {

/**
 * This is an implementation of the universal recursive algorithm run on a quasi-polynomial tree
 * The original version was based on https://www.mimuw.edu.pl/~parys/publications/2018-parity-algorithm.pdf
 * I modified the algorithm to produce a strategy as well.
 * I changed it with the optimizations from Lehtinen et al (arXiv)
 * Furthermore I added some optimizations, i.e., shortcuts in the tree
 */

class ZLKQSolver : public Solver
{
public:
    ZLKQSolver(Oink& oink, Game& game);
    virtual ~ZLKQSolver();

    virtual void run();

protected:
    unsigned long long iterations = 0;

    uintqueue Q;
    int *str;

    bitset W0, W1; // current approximation of winning areas

    inline void attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y);
    void solve(bitset &Subgame, int vtop, int pe, int po);
};

}

#endif
