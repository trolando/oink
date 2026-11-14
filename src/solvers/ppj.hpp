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
#ifndef PPJ_HPP
#define PPJ_HPP

#include <climits>
#include <memory>
#include "oink/solver.hpp"

namespace pg {

/**
 * Priority promotion with justifications, after Algorithm 4.4 of Ruben Lapauw's
 * thesis, Reconstructing and Improving Parity Game Solvers with Justifications.
 *
 * A region is a lower bound on justification levels. Its native heads remain
 * parameters until the region closes; all other members are justified by one
 * strategy edge (owner wins) or every residual successor (opponent wins).
 * Promotion keeps lower justifications. A winner change invalidates only the
 * backwards justification closure, with resets delayed until attraction ends.
 */
class PPJSolver : public Solver
{
public:
    PPJSolver(Oink& oink, Game& game) : Solver(oink, game) {}
    void run() override;

private:
    static constexpr int SOLVED = INT_MAX;
    bitset G, dirty, nonempty[2], dom;
    std::unique_ptr<int[]> workspace;
    int *region, *str, *native, *degree, *cnt, *witness;
    std::vector<int> count, headfirst, headend, members;
    std::vector<unsigned> begin, length;
    uintqueue work, invalid, flipped, touched;
    uint64_t iterations = 0, promotions = 0, dominions = 0, resets = 0;

    void move(int v, int p);
    void storeMembers(int p);
    void invalidate(int v);
    void invalidateClosure(int p, int pl);
    int findEscape(int v, int i, int p, int pl);
    void grow(int p);
    void solveDominion(int p, bool attract);
};
}
#endif
