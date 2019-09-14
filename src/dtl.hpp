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

#include <boost/dynamic_bitset.hpp>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class DTLSolver : public Solver
{
public:
    DTLSolver(Oink *oink, Game *game);
    virtual ~DTLSolver();

    virtual void run();

    typedef boost::dynamic_bitset<unsigned long long> bitset;

protected:
    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    int steps = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle

    uintqueue Q;

    uintqueue Zvec;
    int *str;
    int *Candidates;
    bitset Distractions;

    uintqueue pea_vS; // vS
    uintqueue pea_iS; // iS
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index
    std::vector<int> tangle; // stores the new tangle

    uintqueue tangleto;
    bitset bs_exits;

    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset Z; // current region

    bitset Even, Odd, CurG, SubEven, SubOdd;

    inline bool attracts(const int pl, const int v, bitset &Z, bitset &R);
    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z, bitset &G);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G);
    inline int attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G);
    inline void attractVerticesM(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);
    bool attractTangleM(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio);
    inline int attractTanglesM(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);

    void partition(bitset &R, int top, bitset &Even, bitset &Odd, bool check_distractions);
    bool sptl(bitset &R, int top, int player, bitset &Even, bitset &Odd);
    bool extractTangles(int i, bitset &R, int *str);

    void go(const int player);
    bool sptl_loop(const int player);
    bool prune(const int player);
};

}

#endif
