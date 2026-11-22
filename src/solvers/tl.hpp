/*
 * Copyright 2020 Tom van Dijk
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

#ifndef TL_HPP
#define TL_HPP

#include <memory>

#include "oink/solver.hpp"

namespace pg {

class TLSolver : public Solver
{
public:
    TLSolver(Oink& oink, Game& game);
    virtual ~TLSolver();

    virtual void run();

protected:
    uint64_t iterations = 0;
    uint64_t dominions = 0;
    uint64_t tangles = 0;
    uint64_t steps = 0;

    std::vector<int> tvdata; // the vertex-strategy pairs of the tangles (each ending with -1)
    std::vector<int> toutdata; // the escapes of the tangles (each ending with -1)
    std::vector<size_t> tvpos; // for each tangle, its vertex-strategy pairs in <tvdata>
    std::vector<size_t> toutpos; // for each tangle, its escapes in <toutdata>
    // The tangles that escape to a vertex, in discovery order in <tin_pool>.
    // Only pages with an escape destination have descriptors, and only vertices
    // marked in <has_tin> have an initialized descriptor.
    struct Backlinks { size_t offset; unsigned int size, capacity; };
    static constexpr unsigned int TIN_PAGE_SIZE = 256;
    std::vector<int> tin_pool;
    std::vector<std::unique_ptr<Backlinks[]>> tin_pages;
    bitset has_tin;
    std::vector<int> tpr; // priority of a tangle
    std::vector<unsigned int> tescs; // number of escapes of a tangle

    uintqueue Q; // main queue when attracting vertices (during a pass: the vertices of its regions)

    const int *inedges; // the in-edges of the game (each list ends with -1)
    const int *firstin; // offsets into <inedges>, owned by the game
    std::unique_ptr<int[]> str; // stores currently assigned strategy of each vertex
    std::unique_ptr<unsigned int[]> escs; // stores remaining escapes of each vertex
    std::unique_ptr<unsigned int[]> gdeg; // number of successors in <G> of each vertex in <G>, or 0 if not known yet
    bool escs_is_gdeg = false; // after a dominion: for vertices left in <G>, escs (if not 0) is the new gdeg

    // the regions of the current pass, from the top (see tl): their first vertex in <Q>, first
    // head in <heads>, priority, and (for the dominion check) whether there is no region of
    // each player above it
    struct Region { size_t begin, hbegin; int pr; bool highest0, highest1; };
    std::vector<Region> regions;
    std::vector<int> heads; // the heads of these regions, region by region
    std::unique_ptr<int[]> regidx; // for each vertex in <G>, the index of its region in <regions>
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
    std::unique_ptr<int[]> astr; // the strategy of a region tangle being analysed
    bitset LB; // the vertices of a region tangle being analysed
    bitset H[2]; // the highest region of each player in the current pass

    struct PearceFrame { unsigned int vertex, edge; };
    std::vector<PearceFrame> pea_state; // DFS stack
    std::vector<unsigned int> pea_S;
    std::unique_ptr<unsigned int[]> pea_vidx; // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    std::vector<int> tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset C; // during an attraction: the vertices in <R> but not in <Z>
    unsigned int ncand = 0; // the number of vertices in <C>
    bool prefetch = false; // prefetch the in-edges of the next vertices to attract from (see prefetchAhead)

    bitset R; // remaining region (in tl)
    bitset Z; // current region (in tl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset V; // heads 1
    bitset W; // heads 2

    bitset Even; // accumulate even regions
    bitset Odd; // accumulate odd regions

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z);
    inline void prefetchAhead(const unsigned int qh);
    bool attractTangle(const int t, const int v0, const int pl, bitset &R, bitset &Z, int maxpr);
    int tangleTarget(const int t);
    void addTangle(const int pr, const int *strat);
    bool learnRegion(const int ridx, const int pl, const int pr);
    bool analyseLazies(const bool broken);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, int maxpr);
    void setCandidates(const bitset &R, const bitset &Z);

    bool tl(void);
    void prepareExtraction();
    bool extractTangles(int i, bitset &R, int *str, int tangle_pr = -1);
};

}

#endif
