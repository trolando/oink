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

#ifndef PTL_HPP
#define PTL_HPP

#include "oink/solver.hpp"

namespace pg {

class PTLSolver : public Solver
{
public:
    PTLSolver(Oink& oink, Game& game);
    virtual ~PTLSolver();

    virtual void run();

protected:
    bool multiplayer = true;

    int iterations = 0;
    int dominions = 0;
    int tangles = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle

    uintqueue Q;
    uintqueue SolvedQ0;
    uintqueue SolvedQ1;

    uintqueue Zvec;
    int *str;

    uintqueue pea_vS; // vS
    uintqueue pea_iS; // iS
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index
    std::vector<int> tangle; // stores the new tangle

    uintqueue tangleto;
    bitset bs_exits;

    bitset H; // currently halted vertices
    bitset G; // the unsolved game
    bitset S0; // solved vertices (in queue Q)
    bitset S1; // solved vertices (in queue Q)

    std::string path; // for logging

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z, bitset &G);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G);
    inline int attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G);

    bool search(bitset &R, int top, int player);
    bool search_rec(bitset &R, int top, int player, bitset &XX);
    bool extractTangles(int i, bitset &R, int *str);

    void solve(void);
};

class SPPTLSolver : public PTLSolver
{
public:
    SPPTLSolver(Oink& oink, Game& game) : PTLSolver(oink, game) { multiplayer = false; }
    virtual ~SPPTLSolver() { }
};

}

#endif
