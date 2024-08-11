/*
 * Copyright 2017-2019 Tom van Dijk
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

#ifndef DTL_HPP
#define DTL_HPP

#include "oink/solver.hpp"

namespace pg {

class DTLSolver : public Solver
{
public:
    DTLSolver(Oink& oink, Game& game);
    virtual ~DTLSolver();

    virtual void run();

protected:
    bool interleaved = false;

    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    int steps = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle

    uintqueue Q; // main queue when attracting vertices
    uintqueue SQ; // auxiliary queue for solved vertices (in dominions)

    std::vector<int> dom_vector; // stores the solved vertices (in dominions)

    uintqueue Zvec; // stores current region when using trace level 2 or higher
    int *str; // stores currently assigned strategy of each vertex
    int *Candidates; // list of current distraction candidates

    int *dvalue; // for trace
    int dvalue_idx;

    int *reg; // for trace
    int regidx;
    bitset regflag;

    uintqueue pea_vS; // vS
    uintqueue pea_iS; // iS
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset Z; // current region (in sptl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset W; // already did tangle this iteration

    bitset Player; // vertices of player 0/1

    int cur_depth; // set before extractTangles

    inline bool attracts(const int pl, const int v, bitset &Z, bitset &R);
    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);

    bool search(const int player);
    void search_rec(bitset &R, const int player, const int depth);
    bool sptl(bitset &V, bitset &W, bitset &R, const int player);
    bool extractTangles(int i, bitset &R, int *str);
    bool computeDistanceValues(bitset &V, bitset &R, bitset &G, int* val, const int player);
    int best_vertices(bitset &V, bitset &G, const int player);
};

class IDTLSolver : public DTLSolver
{
public:
    IDTLSolver(Oink& oink, Game& game) : DTLSolver(oink, game) { interleaved = true; }
    virtual ~IDTLSolver() { }
};


}

#endif
