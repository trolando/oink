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
 *
 * Zielonka's algorithm with justifications, after Algorithm 5 of Lapauw, Bruynooghe and
 * Denecker, "Improving Parity Game Solvers with Justifications" (VMCAI 2020), and Ruben
 * Lapauw's zlkjb in pull request trolando/oink#4.
 */

#ifndef ZLKJ_HPP
#define ZLKJ_HPP

#include <climits>
#include <cstdint>
#include <memory>
#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * Zielonka's algorithm with justifications (sequential).
 *
 * As ZLK, ZLKJ attracts the vertices of the highest priorities of a subgame (the heads) to a
 * region and solves the rest of the subgame recursively. When the opponent attracts vertices
 * of the region or of its part of the rest, ZLK solves the whole remaining subgame again.
 * ZLKJ keeps the vertices whose justification is still valid, resets only the vertices whose
 * justification depends on the vertices the opponent attracted, attracts these back to the
 * region where it can, and solves only the remaining ones again, in a loop.
 *
 * The justification is implicit: a vertex won by its owner is justified by its strategy edge,
 * a vertex won by the opponent of its owner by all its edges. Heads have no justification
 * (strategy -2) until their region is final. Unlike ZLK, the heads of a region are only the
 * highest run of priorities of one parity in the subgame (no skipping of attracted vertices
 * of the other parity), because vertices that are attracted back must not form a cycle with
 * a lower head through a vertex of the other parity with a higher priority.
 */
class ZLKJSolver : public Solver
{
public:
    ZLKJSolver(Oink& oink, Game& game);
    virtual ~ZLKJSolver();

    virtual void run();

protected:
    uint64_t iterations = 0, rounds = 0, resets = 0, reattracted = 0;

    // One allocation for the arrays per vertex: region, strategy and the winner (bytes).
    std::unique_ptr<int[]> workspace;
    // DIS (disabled), BOT (open: not yet in a region), a negative escape counter (open), or
    // the region (frame) the vertex was attracted to.
    int *region;
    // The justification: the successor for a vertex won by its owner, -1 for a vertex won by
    // the opponent of its owner (all its edges), -2 for a head that has no strategy yet.
    int *strategy;
    int8_t *winning;

    // A superset of the open vertices; the head search skips the other vertices with it.
    bitset open;
    // Marks for removing duplicates from the vertices won by a player (see finish).
    bitset seen;

    // The vertices won by each player. A frame owns the entries from its start on: its region,
    // the vertices it gave to the opponent and the results of its subgames. Entries of vertices
    // that were reset or lost are removed when the frame finishes.
    std::vector<int> W[2];
    // Scratch vectors of the opponent attractor: vertices whose justification edge was lost,
    // and the vertices that are reset.
    std::vector<int> lost, reset;

    // The region ids stay below this limit: they are renumbered before they reach it (see
    // renumber() in zlkj.cpp). A test may set it low, to renumber often. (After the other members,
    // which keep their offsets in the hot code.)
    int region_limit = INT_MAX;
    uint64_t renumbered = 0; // how often the region ids were renumbered

    int attractHeads(int i, int r, std::vector<int> &A, std::vector<int> &H);
    int attractLosing(int r, int pl, std::vector<int> &O, size_t start, std::vector<int> &H);
    void resetClosure(int r, int pl);
    int reattract(int r, int r2, int pl, std::vector<int> &A);
    void finish(int r, int pl, std::vector<int> &H);
    void compact(int p, size_t start);
    [[gnu::cold, gnu::noinline]] int renumber(std::vector<int> &live);
};

}

#endif
