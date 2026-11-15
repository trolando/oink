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

    // A stack of vertices in a buffer of the workspace (at most every vertex once).
    struct Stack {
        int *data = nullptr;
        unsigned top = 0;
        void push(int v) { data[top++] = v; }
        int pop() { return data[--top]; }
        int operator[](unsigned i) const { return data[i]; }
        unsigned size() const { return top; }
        bool empty() const { return top == 0; }
        bool nonempty() const { return top != 0; }
        void clear() { top = 0; }
    };

    bitset G, dirty, nonempty[2];
    // One allocation holds all arrays of a vertex or a priority and the stacks: on small games,
    // separate allocations cost about as much as the solving.
    std::unique_ptr<int[]> workspace;
    int *region, *str, *native, *degree, *cnt, *witness;
    int *count, *headfirst, *headend;
    unsigned *begin, *length, priorities = 0;
    std::vector<int> members;
    // <flipped> is also the queue of invalidated vertices: it is empty after every closure.
    Stack work, flipped, touched;
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
