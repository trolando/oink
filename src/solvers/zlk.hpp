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

#ifndef ZLK_HPP
#define ZLK_HPP

#include <atomic>
#include <climits>
#include <cstdint>
#include <memory>
#include <queue>
#include <vector>
#include <lace.h>

#include "oink/solver.hpp"

namespace pg {

class ZLKSolver : public Solver
{
public:
    ZLKSolver(Oink& oink, Game& game);
    virtual ~ZLKSolver();

    virtual void run();

    uint64_t iterations = 0;

    int *inverse;
    int max_prio;

    // Shared synchronisation word per vertex (see zlk.cpp for the DIS/BOT/counter/r
    // encoding). Mutated concurrently by the parallel attractor, so it must be atomic.
    std::unique_ptr<std::atomic<int>[]> region;
    // Winner and strategy per vertex. Each vertex is claimed (in region[]) by exactly
    // one worker, which is then the sole writer here; reads happen only in the
    // sequential phases after the Lace tasks have fully joined, so these stay plain.
    // The winner is a byte (-1, 0 or 1) to keep the array small for the cache.
    std::unique_ptr<int8_t[]> winning;
    std::unique_ptr<int[]> strategy;

    // A superset of the open vertices: those in the current subgame that are not yet attracted
    // to a region, that is, with region BOT or an escape counter. The sequential attractor
    // removes the vertices it attracts, the parallel attractor does not; vertices become open
    // again when the region is reset. The head search skips the vertices not in <open>, and
    // drops the vertices that it finds closed.
    bitset open;

    bool to_inversion = true;
    bool only_recompute_when_attracted = true;

    // Per-worker scratch buffer collecting the vertices a worker attracts during one
    // parallel attractor pass. One entry per Lace worker; cache-line aligned so that
    // concurrent push_back from different workers does not false-share the headers.
    struct alignas(64) par_buffer {
        std::vector<int> items;
    };
    std::vector<par_buffer> pvec;

    // The region ids stay below this limit: they are renumbered before they reach it (see
    // renumber() in zlk.cpp). A test may set it low, to renumber often. (After the other members,
    // which keep their offsets in the hot code.)
    int region_limit = INT_MAX;
    uint64_t renumbered = 0; // how often the region ids were renumbered

    int attractExt(int i, int r, std::vector<int> *R, std::vector<int> *H);
    int attractLosing(int i, int r, std::vector<int> *H, std::vector<int> *R);
    void attractParT(lace_worker*, int pl, int cur, int r);
    int attractPar(lace_worker*, int i, int r, std::vector<int>* R, std::vector<int>* H);
    [[gnu::cold, gnu::noinline]] int renumber(std::vector<int> &history);
    // friend void updateOutcount_WORK(WorkerP*, Task*, int, int, ZLKSolver*);
};

class UnoptimizedZLKSolver : public ZLKSolver
{
public:
    UnoptimizedZLKSolver(Oink& oink, Game& game) : ZLKSolver(oink, game) { to_inversion = false; only_recompute_when_attracted = false; }
    virtual ~UnoptimizedZLKSolver() { }
};

}

#endif 
