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

#include <cassert>
#include <cstring>
#include <vector>

#include "fpj.hpp"

namespace pg {

FPJSolver::FPJSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

FPJSolver::~FPJSolver()
{
}

/**
 * Fixpoint iteration with justifications.
 *
 * Vertices are processed in blocks of consecutive (enabled) vertices with the same parity.
 * A vertex is justified when its current evaluation (winner and strategy) is known to be
 * consistent; justified vertices are skipped. When a block contains new distractions, all
 * justified vertices that depend on them are reset, and solving continues at the lowest
 * reset vertex.
 *
 * - greedy: every evaluated vertex is justified immediately (variation 1);
 * - otherwise: only new distractions are justified immediately, and the other evaluated
 *   vertices of a block are justified when the block ends without new distractions.
 */
void
FPJSolver::solveSeq(bool greedy)
{
    const int n = nodecount();

    /**
     * Allocate and initialize data structures
     */

    std::vector<int> strategy(n, -1); // the current strategy
    bitset justified(n);   // whether a vertex is justified
    bitset distraction(n); // whether a vertex is won by the opponent
    bitset parity(n);      // the parity of every vertex's priority
    for (int v=0; v<n; v++) parity[v] = priority(v)&1;
    const bitset enabled = ~disabled;

    // won[pl][v]: v is enabled and currently won by pl (one bit test per edge)
    bitset won[2];
    won[1] = parity;
    won[1] -= disabled;
    won[0] = ~parity;
    won[0] -= disabled;

    auto makeDistraction = [&](int v) {
        const int par = parity.test(v);
        distraction.set(v);
        won[par].reset(v);
        won[1-par].set(v);
    };
    auto clearDistraction = [&](int v) {
        const int par = parity.test(v);
        distraction.reset(v);
        won[1-par].reset(v);
        won[par].set(v);
    };

    // blockend[v]: for an enabled vertex v, the first enabled vertex after v with the other
    // parity (or n); disabled vertices do not end a block
    std::vector<int> blockend(n, n);
    for (int v=n-1, next=n; v>=0; v--) {
        if (disabled[v]) continue;
        if (next == n) blockend[v] = n;
        else if (parity[next] != parity[v]) blockend[v] = next;
        else blockend[v] = blockend[next];
        next = v;
    }

    uintqueue Q(n);

    const uint64_t *J = justified.data();
    const uint64_t *EN = enabled.data();

    // find the first enabled vertex
    int i = 0;
    while (i < n and disabled[i]) i++;

    bool blockchanged = false; // did the current block change (new distractions)
    int blockstart = i;        // first vertex of the current block
    int cur_parity = i < n ? parity.test(i) : 0;

    while (i < n) {
        const int end = blockend[i];

        /**
         * Evaluate the unjustified vertices in [i, end), in order, 64 at a time
         */
        const int last_word = (end - 1) >> 6;
        for (int w = i >> 6; w <= last_word; w++) {
            uint64_t bits = EN[w] & ~J[w];
            if (w == (i >> 6)) bits &= ~uint64_t(0) << (i & 63);
            if (w == last_word and (end & 63) != 0) bits &= (uint64_t(1) << (end & 63)) - 1;

            while (bits) {
                const int v = (w << 6) + intrinsics::countr_zero64(bits);
                bits &= bits - 1;

                if (greedy) justified.set(v);

                // compute one step winner of <v> and update the strategy
                const int pl = owner(v);
                const uint64_t *good = won[pl].data();
                int onestep_winner = 1 - pl, str = -1;
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    const int to = *curedge;
                    if ((good[to >> 6] >> (to & 63)) & 1) {
                        onestep_winner = pl;
                        str = to;
                        break;
                    }
                }
                strategy[v] = str;

                if (cur_parity != onestep_winner) {
                    Q.push(v); // add to the Queue, because all justified predecessors are now invalid
                    makeDistraction(v);
                    if (!greedy) justified.set(v);
                    blockchanged = true;
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[38;5;165;1mjustified*\033[m " << label_vertex(v);
                        if (strategy[v] != -1) logger << " => " << label_vertex(strategy[v]);
                        logger << std::endl;
                    }
#endif
                }
#ifndef NDEBUG
                else if (greedy and trace >= 2) {
                    logger << "\033[38;5;165;1mjustified\033[m " << label_vertex(v);
                    if (strategy[v] != -1) logger << " => " << label_vertex(strategy[v]);
                    logger << std::endl;
                }
#endif
            }
        }

        /**
         * The block ended
         */
        i = end;

        if (blockchanged) {
            /**
             * We have to reset all justified vertices that are now no longer justified...
             */
            while (!Q.empty()) {
                const int v = Q.pop();
                for (auto curedge = ins(v); *curedge != -1; curedge++) {
                    const int from = *curedge;
                    if (disabled[from]) continue;
                    if (justified[from]) {
                        if (strategy[from] == -1 or strategy[from] == v) {
#ifndef NDEBUG
                            if (trace >= 2) logger << "\033[31;1mresetting\033[m " << label_vertex(from) << std::endl;
#endif
                            justified.reset(from);
                            clearDistraction(from);
                            Q.push(from);
                            if (from < i) i = from; // afterwards, continue at lowest unjustified vertex
                        }
                    }
                }
            }
#ifndef NDEBUG
            if (trace) logger << "restarting after finding distractions of prio " << priority(end-1) << std::endl;
#endif
            iterations++;
            blockchanged = false;
            // non-greedy: the unjustified vertices of the block must be evaluated again
            if (!greedy and i > blockstart) i = blockstart;
        } else if (!greedy) {
            /**
             * We now know that the current strategy of all unjustified vertices of the block is
             * justified: justified |= enabled on [blockstart, end)
             */
            uint64_t *JW = justified.data();
            const int last = (end - 1) >> 6;
            for (int w = blockstart >> 6; w <= last; w++) {
                uint64_t mask = EN[w];
                if (w == (blockstart >> 6)) mask &= ~uint64_t(0) << (blockstart & 63);
                if (w == last and (end & 63) != 0) mask &= (uint64_t(1) << (end & 63)) - 1;
#ifndef NDEBUG
                if (trace >= 2) {
                    for (uint64_t m = mask & ~JW[w]; m; m &= m - 1) {
                        const int v = (w << 6) + intrinsics::countr_zero64(m);
                        logger << "\033[38;5;165;1mjustified\033[m " << label_vertex(v);
                        if (strategy[v] != -1) logger << " => " << label_vertex(strategy[v]);
                        logger << std::endl;
                    }
                }
#endif
                JW[w] |= mask;
            }
        }

        if (i == n) break;

        // continue with the block of vertex i
        cur_parity = parity.test(i);
        blockstart = i;
    }

    // done
    for (int v=0; v<n; v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    logger << "solved with " << iterations << " iterations." << std::endl;
}

/**
 * Variation 1: greedy justification
 */
void
FPJSolver::runSeqGreedy()
{
    solveSeq(true);
}

void
FPJSolver::runSeq()
{
    solveSeq(false);
}

void
FPJSolver::run()
{
    if (greedy) runSeqGreedy();
    else runSeq();
}

}
