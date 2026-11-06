/*
 * Copyright 2020-2024 Tom van Dijk
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

#include <vector>

#include "fpjm.hpp"
#include "oink/uintqueue.hpp"

namespace pg {

FPJMSolver::FPJMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

FPJMSolver::~FPJMSolver()
{
}

void
FPJMSolver::runSeq()
{
    const int n = nodecount();

    /**
     * Allocate and initialize data structures
     */

    std::vector<int> nstrat(n, 0); // number of recorded winning moves per vertex
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

    // allocate the multi-strategy (set of all winning moves), unless a previous
    // call (e.g. on a different subgame) already did so
    if (!hasMultiStrategy()) initMultiStrategy();

    // build the in-edge array plus the in->out edge map, so that pruning the
    // strategy edge from->v while walking ins(v) is O(1) (no scan of from's edges)
    game.build_in_to_out();
    const int* in_base = game.inedges();
    const int* in_to_out = game.inToOut();

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

                // compute one step winner of <v> and update the strategy. Unlike FPJ, we
                // do not stop at the first winning move: we record *every* outgoing edge
                // to a co-winning vertex (the maximal set of moves valid at this point).
                clearStrategyEdges(v); // drop any edges recorded on a previous pass
                const int o = owner(v);
                const uint64_t *good = won[o].data();
                int onestep_winner = 1 - o; // default: owner cannot reach a vertex good for itself
                int cnt = 0;
                int k = 0;
                for (auto curedge = outs(v); *curedge != -1; curedge++, k++) {
                    const int to = *curedge;
                    if ((good[to >> 6] >> (to & 63)) & 1) {
                        onestep_winner = o;
                        addStrategyEdge(v, k);
                        cnt++;
                    }
                }
                nstrat[v] = cnt;

                // evaluation stays the same
                if (cur_parity != onestep_winner) {
                    Q.push(v); // add to the Queue, because all justified predecessors are now invalid
                    makeDistraction(v);
                    justified.set(v);
                    blockchanged = true;
#ifndef NDEBUG
                    if (trace >= 2) logger << "\033[38;5;165;1mjustified*\033[m " << label_vertex(v) << std::endl;
#endif
                }
            }
        }

        /**
         * The block ended
         */
        i = end;

        if (blockchanged) {
            /**
             * Propagate the new distractions: a justified predecessor loses its
             * justification only when its *last* winning move is pruned. A vertex
             * with no winning move of its own (nstrat == 0, won by the opponent of
             * its owner) depends on all successors and is reset immediately.
             */
            while (!Q.empty()) {
                const int v = Q.pop();
                for (auto curedge = ins(v); *curedge != -1; curedge++) {
                    const int from = *curedge;
                    if (disabled[from]) continue;
                    if (!justified[from]) continue;
                    bool reset = false;
                    if (nstrat[from] == 0) {
                        // <from> is won by the opponent of its owner: reset
                        reset = true;
                    } else {
                        // prune the strategy edge from->v in O(1) via the in->out map
                        const int oidx = in_to_out[curedge - in_base];
                        if (isStrategyEdgeIndex(oidx)) {
                            removeStrategyEdgeIndex(oidx);
                            if (--nstrat[from] == 0) reset = true; // last move pruned
                        }
                        // otherwise from->v was not a winning move: still justified
                    }
                    if (reset) {
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
#ifndef NDEBUG
            if (trace) logger << "restarting after finding distractions of prio " << priority(end-1) << std::endl;
#endif
            iterations++;
            blockchanged = false;
            // the unjustified vertices of the block must be evaluated again
            if (i > blockstart) i = blockstart;
        } else {
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
                        logger << "\033[38;5;165;1mjustified\033[m " << label_vertex((w << 6) + intrinsics::countr_zero64(m)) << std::endl;
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

    // done: report the winner. Won vertices use the -1 sentinel; their winning
    // moves live in the multi-strategy (see strategyTargets/the verifier).
    for (int v=0; v<n; v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, -1);
    }

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPJMSolver::run()
{
    // multi-strategy variant is sequential only
    runSeq();
}

}
