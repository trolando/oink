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

#ifndef PMTL_HPP
#define PMTL_HPP

#include <queue>

#include "oink/solver.hpp"
#include "gpm.hpp"

namespace pg {

class PMTLSolver : public Solver
{
public:
    PMTLSolver(Oink& oink, Game& game, MeasureKind kind = MeasureKind::Ordered);
    virtual ~PMTLSolver();

    virtual void run();
    virtual bool parseOptions(std::string&);

private:
    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    int lifts = 0;
    MeasureKind measure_kind = MeasureKind::Ordered;
    bool use_tangles = true;
    bool go_up = false;
    bool use_shortcuts = true;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle
    std::vector<std::vector<int>> tangles_of_pr; // for each priority, its tangles (in order of creation)
    std::vector<int> Zlist; // the vertices added to <Z> (some may have been removed again), so that
                            // work per region is proportional to its size, not to the whole game

    bitset G; // remaining unsolved game
    uintqueue Q; // main queue when attracting vertices
    int *str; // stores currently assigned strategy of each vertex
    // for each player, its unsolved vertices by measure as sorted in its last update, and the
    // vertices whose measure changed since (so that only these need to be sorted again)
    std::vector<int> order[2];
    bitset Changed[2];
    bool sorted[2] = {false, false};
    std::vector<int> moved, merged;

    bitset Z;
    bitset R;
    bitset E;
    bitset Updated;
    bitset OpenHead; // heads of the current region that can escape it
    bitset WasHead; // vertices that were an extra head in this round (at most once)
    bitset TZ; // the part of a partially closed region that does not depend on open heads
    std::vector<int> heads; // the heads of the current region

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);

    template<bool GoUp, bool UseTangles> bool update(Measures &pm, Measures &target_pm, const int player);
    void solve(Measures &pm, const int player);
    void shortcuts(const int player, Measures &src, Measures &target1, Measures &target2);

    /** stuff for extractTangles **/
    uintqueue pea_vS; // vS
    uintqueue pea_iS; // iS
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index
    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered
    bitset S; // solved vertices (in queue Q)
    std::vector<int> dom_vector; // stores the solved vertices (in dominions)
    bool extractTangles(int i, bitset &R, int *str);

};

}

#endif
