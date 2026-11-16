/*
 * Copyright 2026 Tom van Dijk
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

#ifndef TLJ_HPP
#define TLJ_HPP

#include "oink/solver.hpp"

namespace pg {

/**
 * Tangle learning with justifications.
 *
 * Regions persist in region[] together with their implicit justification graph.
 * A vertex assigned to its owner is justified by its selected strategy edge;
 * a vertex assigned to the opponent is justified by every residual successor.
 * Changes that break these dependencies reset the backwards justification
 * closure. Original heads group consecutive enabled priorities of one parity.
 *
 * Visits caused only by promoted tangles reuse the surviving closure and process
 * newly attracted vertices. New targets and invalidated heads require complete
 * region replay. Current TL also retains unaffected regions; persistence alone
 * does not determine how much work the two solvers perform.
 *
 * Reference: Lapauw, Bruynooghe, Denecker, "Improving Parity Game Solvers with
 * Justifications", VMCAI 2020 (Algorithm 7).
 */
class TLJSolver : public Solver
{
public:
    TLJSolver(Oink& oink, Game& game);
    virtual ~TLJSolver();

    virtual void run();

protected:
    int iterations = 0;
    int dominions = 0;
    int tangles = 0;

    size_t promoted = 0; // the tangles before this one are linked into tin (see promote)
    std::vector<int> tvdata; // the vertex-strategy pairs of the tangles (each ending with -1)
    std::vector<int> toutdata; // the escapes of the tangles (each ending with -1)
    std::vector<size_t> tvpos; // for each tangle, its vertex-strategy pairs in <tvdata>
    std::vector<size_t> toutpos; // for each tangle, its escapes in <toutdata>
    struct Backlinks { std::vector<int> by_player[2]; };
    std::vector<Backlinks> tin; // allocated only for vertices with backlinks
    std::vector<int> tin_index;
    bitset has_tin[2]; // vertices with a backlink from each player
    bitset tangled; // the vertices of the opponent attracted into their region with a tangle
                    // (the mark can remain when it moves on: see solveDominion)
    std::vector<int> tpr; // priority of a tangle (-1 = deleted)
    std::vector<int> tcnt; // for attracting: the escapes of a tangle (-1: not counted)
    std::vector<unsigned> tcpos; // the size of <zlist> when <tcnt> was counted
    std::vector<int> tblock; // remembered escape, revalidated before counting
    std::vector<int> ttouched; // the tangles with a counter

    uintqueue Q; // work queue when attracting to a dominion
    uintqueue Dq; // justification-invalidation closure queue

    int *str; // currently assigned strategy of each vertex (implicit justification)
    int *region; // persistent region (priority) of each vertex, -1 = unassigned

    uintqueue pea_state; // v,i,...
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    std::vector<int> pea_seen; // the vertices with a nonzero rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset Z; // current region being grown
    bitset C; // the vertices in G with a region below <cpr> or unassigned (the candidates
              // for attraction when growing region <cpr>; see setCandidates)
    int cpr = -1; // see <C> (-1: to be computed again)
    bitset G; // the unsolved game
    bitset S; // dominion vertices found during extraction

    bitset dirty; // regions (by priority) that still need to be (re)examined
    struct RegionWork {
        int first = -1, last = -1; // promoted tangles targeting this region
        bool replay = true; // old members must be replayed
    };
    std::vector<RegionWork> region_work;
    std::vector<int> pending_next; // next promoted tangle at the same target region
    bitset extracted; // priorities whose (closed) region has been extracted at its current size

    std::vector<int> gdeg; // the number of successors in G of each vertex in G
    std::vector<int> rcount; // the number of vertices in each region
    bitset nonempty[2]; // the nonempty regions of each parity
    bitset anyregion; // the nonempty regions
    std::vector<int> cnt; // for attracting: the escapes of a vertex of the opponent (-1: not
                          // counted; with a region of the player above: 1 while it has one)
    std::vector<unsigned> cpos; // for attracting: the out-edge of a vertex of the opponent to
                                // an escape (with a region of the player above)
    uintqueue touched; // the vertices with a counter
    std::vector<int> rprio; // for each vertex, the (compressed) priority of the region it heads
    std::vector<int> mpool; // the members of the regions (in segments, see setMembers)
    std::vector<unsigned> mbeg, mlen; // for each priority p, the segment in <mpool> with the
                                      // vertices of the last region p (a superset of the
                                      // current one)
    std::vector<int> gfirst; // for each region priority, the first of its heads (vertices)
    std::vector<int> gend; // for each region priority, one past the last of its heads (0: none)
    uintqueue zlist; // the vertices of <Z>, in the order of attraction (its queue)
    int zlast = -1; // the region whose vertices (its members) are in <Z>, if any
    std::vector<int> heads; // the heads of the current region, in order
    int sweep_lo = 0; // the nonempty regions from <sweep_lo> up to the current one are
                      // marked dirty (see run)
    bool restart = false; // a dominion was solved: restart the sweep (see run)
    int pull_lo = 0; // the lowest priority of a region with a new <pl>-target (see pull)
    int pull_count = 0; // the number of vertices pulled into the current region
    uintqueue flipped; // the vertices pulled from a region of the other player

    inline void setRegion(const int v, const int r);
    inline int findEscape(const int from, int i, const int p, const int pl);
    inline void invalidate(const int m);
    inline void pull(const int v, const int s, const int p, const int pl);
    bool attractTangleInto(const int t, const unsigned vpos, const int p, const int pl);
    void invalidateOpponent(const int start, const int p, const int pl);
    void invalidateSolved();
    void solveDominion(bitset &Dom, const int pl, const bool attract);
    void setMembers(const int p);
    void setCandidates(const int p);
    void improveRegion(int p);
    void promote();
    bool extractTangles(int i, bitset &R, int *str);
};

}

#endif
