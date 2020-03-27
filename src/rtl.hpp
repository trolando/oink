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

#include "oink.hpp"
#include "solver.hpp"
#include "npp.hpp"

namespace pg {

#define CACHEOUT 0

class RTLSolver : public Solver
{
public:
    RTLSolver(Oink *oink, Game *game);
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

    uintqueue Q;
    uintqueue Qtar;
    uintqueue Zvec;
    uintqueue heads;
    int *str;
    int *tarj;
    int pre;
    uintqueue tarres;
    uintqueue tangleto;
    bitset bs_exits;

    int *region;
    bitset H; // currently halted vertices
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    uintqueue open_heads;

    bool onesided = false;

    int *val; // for cascader

    inline void attractVertices(const int pr, const int pl, int v, bitset &R, int *str, bitset &Z);
    bool attractTangle(const int t, const int pl, bitset &R, int *str, bitset &Z);
    inline int attractTangles(const int pr, const int pl, int v, bitset &R, int *str, bitset &Z);
    inline int attractTangles(const int pr, const int pl, bitset &R, int *str, bitset &Z);
    void search(bitset &R, int top, int only_player, int depth);
    void extractTangles(int i, bitset &R, int *str);
    void loop(int only_player);
};

class ORTLSolver : public RTLSolver
{
public:
    ORTLSolver(Oink *oink, Game *game) : RTLSolver(oink, game) { onesided = true; }
    virtual ~ORTLSolver() { }
};

}

#endif
