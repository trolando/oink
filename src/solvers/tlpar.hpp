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

#ifndef TLPAR_HPP
#define TLPAR_HPP

#include <vector>

#include "oink/solver.hpp"
#include "lace.h"

namespace pg {

/**
 * Parallel tangle learning: TLSolver with parallel attraction (using Lace).
 *
 * The attractor computations (regions, and extending dominions in the whole game) are
 * where TL spends its time. Here each attraction drains its queue in rounds: while few
 * vertices are left to process, one vertex at a time as in TL; once at least <PAR_MIN> are
 * left, all of them form a frontier that the Lace workers process in parallel, collecting
 * the newly attracted vertices as the next frontier. Vertices are claimed by atomically
 * setting their bit in the region; escape counters (of vertices and tangles) are decremented
 * with a compare-and-swap that also does their lazy initialization.
 *
 * The regions (as sets) are the same as TL's, but attractor strategies depend on the order
 * of attraction (and so on scheduling), so the tangles learned can differ. Leak checks and
 * extracting tangles are sequential.
 *
 * Parallel rounds only pay off on large games: Oink solves games below the threshold of
 * tlpar (Solvers::parallelMinVertices) with TL instead, and runs with fewer than 2 Lace
 * workers are solved exactly as by TL.
 */

class TLParSolver : public Solver
{
public:
    TLParSolver(Oink& oink, Game& game);
    virtual ~TLParSolver();

    virtual void run();
    virtual bool parseOptions(std::string& opts);

    int par_min = 2048; // smallest queue that is processed in parallel ("parmin=N")
    int par_chunk = 256; // vertices per task in a parallel round ("chunk=N")

    // the Lace task bodies
    void runPar(WorkerP* __lace_worker, Task* __lace_dq_head);
    void roundChunk(int start, int n);

protected:
    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    long par_rounds = 0; // number of parallel rounds
    long par_vertices = 0; // vertices processed in parallel rounds

    std::vector<int> tvdata; // the vertex-strategy pairs of the tangles (each ending with -1)
    std::vector<int> toutdata; // the escapes of the tangles (each ending with -1)
    std::vector<size_t> tvpos; // for each tangle, its vertex-strategy pairs in <tvdata>
    std::vector<size_t> toutpos; // for each tangle, its escapes in <toutdata>
    // for each vertex, the tangles that escape to it (if <has_tin>), in the order they were
    // found: <tin_size> of them, in a block of <tin_cap> at <tin_off> in <tin_pool>
    std::vector<int> tin_pool;
    size_t *tin_off;
    unsigned int *tin_size;
    unsigned int *tin_cap;
    bitset has_tin; // the vertices with tangles that escape to it (to skip reading <tin>)
    std::vector<int> tpr; // priority of a tangle
    std::vector<unsigned int> tescs; // number of escapes of a tangle

    // main queue when attracting vertices (during a pass: the vertices of its regions); a
    // parallel round writes the vertices it attracts after its end (see roundChunk)
    struct Queue : uintqueue { void grow(uint n) { pointer += n; } };
    Queue Q;

    int *str; // stores currently assigned strategy of each vertex
    unsigned int *escs; // stores remaining escapes of each vertex
    unsigned int *gdeg; // number of successors in <G> of each vertex in <G>, or 0 if not known yet
    bool escs_is_gdeg = false; // after a dominion: for vertices left in <G>, escs (if not 0) is the new gdeg

    // the regions of the current pass, from the top (see tl): their first vertex in <Q>, first
    // head in <heads>, priority, and (for the dominion check) whether there is no region of
    // each player above it
    struct Region { size_t begin, hbegin; int pr; bool highest0, highest1; };
    std::vector<Region> regions;
    std::vector<int> heads; // the heads of these regions, region by region
    int *regidx; // for each vertex in <G>, the index of its region in <regions>
    int keep = -1; // the next pass keeps the regions before this one and extends it (see tl)
    std::vector<int> seeds; // the new tangles that go into region <keep> (see tl)
    bitset Zx, Vx; // the vertices and heads of region <zx> of the current pass, if not -1 (see tl)
    int zx = -1;
    size_t pass_tangles = 0; // the tangles from this one on were found in the current pass
    bool count_in_R = false; // count the escapes in <R> instead of <G> (see attractVertices)

    // region tangles (see learnRegion): closed regions learned as a whole, and decomposed into
    // their tangles when broken (see analyseLazies)
    std::vector<char> tlazy; // for each tangle: 1 if a region tangle, 2 if a broken one
    std::vector<int> lreg; // for each region tangle: the region of the current pass it is in, or -1
    std::vector<int> lazies; // the region tangles
    int cur_ridx = -1; // the region being attracted (see attractTangle)
    int *astr; // the strategy of a region tangle being analysed
    bitset LB; // the vertices of a region tangle being analysed
    bitset H[2]; // the highest region of each player in the current pass

    uintqueue pea_state; // v,i,...
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset C; // during an attraction: the vertices in <R> but not in <Z>

    bitset R; // remaining region (in tl)
    bitset Z; // current region (in tl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset V; // heads 1
    bitset W; // heads 2 (only with PARTIALLY_CLOSED)

    // parallel attraction
    WorkerP* lw = nullptr; // the Lace worker running this solver, if any
    Task* lw_head = nullptr; // and its deque head (runPar has no tasks outstanding when it starts a round)
    // per Lace worker: the vertices attracted by the current chunk of a round (each in its own
    // cache line, as the workers update them all the time)
    struct alignas(64) WorkerBuf { std::vector<int> v; };
    std::vector<WorkerBuf> pbuf;
    // during a parallel round: the end of <Q> (see roundChunk). Every task adds to it, and every
    // task reads the members after it (and others before it) all the time, so it has a cache
    // line of its own; otherwise each add takes the line away from all the other workers.
    alignas(64) unsigned int q_tail = 0;
    alignas(64) int par_pl = 0; // the current parallel round attracts for this player
    int par_maxpr = 0; // ... tangles up to this priority
    int par_ridx = -1; // ... into this region (or -1, see drain)
    bitset *par_R = nullptr; // ... from this subgame
    bitset *par_Z = nullptr; // ... to this set

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, int maxpr);
    int tangleTarget(const int t);
    void addTangle(const int pr, const int *strat);
    bool learnRegion(const int ridx, const int pl, const int pr);
    bool analyseLazies(const bool broken);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, int maxpr);
    void setCandidates(const bitset &R, const bitset &Z);

    inline void attractVerticesPar(const int v, std::vector<int> &out);
    void attractTanglePar(const int t, std::vector<int> &out);
    void drain(const int pl, bitset &R, bitset &Z, const int maxpr, unsigned int &qh, const int ridx);

    bool tl(void);
    bool extractTangles(int i, bitset &R, int *str, int tangle_pr = -1);
};

}

#endif
