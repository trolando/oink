/*
 * Copyright 2024-2025 Tom van Dijk, University of Twente
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

#ifndef TLQ_HPP
#define TLQ_HPP

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

class TLQSolver : public Solver
{
public:
    TLQSolver(Oink& oink, Game& game);
    virtual ~TLQSolver();

    virtual void run();

protected:
    unsigned long long iterations = 0;
    unsigned int dominions = 0;
    unsigned int tangles = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle

    uintqueue pea_state; // v,i,...
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto;
    bitset escapes; // which escapes we considered

    uintqueue Q;
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    std::vector<bitset*> Rs;
    std::vector<bitset*> Hs;
    std::vector<int> Ps; // priorities
    std::vector<int> heads; // record heads for the tangle learning
    int *str;

    bitset W0, W1; // current approximation of winning areas

    bool attracts(const int pl, const int v, bitset &Z, bitset &R);
    void attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y);
    bool attractTangle(int t, int pl, bitset &R, bitset &Z, bitset &G);
    void attractTangles(int pl, int v, bitset &R, bitset &Z, bitset &G);
    bool extractTangles(int startvertex, bitset &R);
    void solve(bitset &Subgame, int vtop, int pe, int po);
};

}

#endif
