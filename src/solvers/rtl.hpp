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

#ifndef RTL_HPP
#define RTL_HPP

#include <stack>
#include <set>
#include <map>
#include <tuple>

#include "oink/solver.hpp"

namespace pg {

#define CACHEOUT 0

class RTLSolver : public Solver
{
public:
    RTLSolver(Oink& oink, Game& game);
    virtual ~RTLSolver();

    virtual void run();

protected:
    int max_prio;
    // int *inverse; // reverse lookup priority->vertex

    int iterations = 0;
    int dominions = 0;
    int tangles = 0;

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

    bitset R; // remaining region (in rtl)
    bitset Z; // current region (in rtl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset V; // heads 1
    bitset W; // heads 2

    uintqueue Q;
    int *str;

    bool onesided = false;

    void attractVertices(int pl, int v, bitset &R, bitset &Z, bitset &G, int maxpr);
    bool attractTangle(int t, int pl, bitset &R, bitset &Z, bitset &G, int maxpr);
    void attractTangles(int pl, int v, bitset &R, bitset &Z, bitset &G, int maxpr);
    bool extractTangles(int startvertex, bitset &R);
    bool rtl(bitset &R, int only_player, int depth);
};

class ORTLSolver : public RTLSolver
{
public:
    ORTLSolver(Oink& oink, Game& game) : RTLSolver(oink, game) { onesided = true; }
    virtual ~ORTLSolver() { }
};

}

#endif
