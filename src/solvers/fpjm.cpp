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
    /**
     * Allocate and initialize data structures
     */

    int *nstrat = new int[nodecount()]; // number of recorded winning moves per vertex
    bitset justified(nodecount()); // whether a vertex is justified
    bitset distraction(nodecount()); // whether a vertex is won by the opponent

    bitset parity(nodecount()); // optimization: precompute the parity of every vertex's priority
    for (int v=0; v<nodecount(); v++) parity[v] = priority(v)&1;

    // allocate the multi-strategy (set of all winning moves), unless a previous
    // call (e.g. on a different subgame) already did so
    if (!hasMultiStrategy()) initMultiStrategy();

    // build the in-edge array plus the in->out edge map, so that pruning the
    // strategy edge from->v while walking ins(v) is O(1) (no scan of from's edges)
    game.build_in_to_out();
    const int* in_base = game.inedges();
    const int* in_to_out = game.inToOut();

    uintqueue Q(nodecount());

    /**
     * Initialize loop
     */

    bool blockchanged = false; // did the current block change (new distractions)
    int cur_parity = parity[0]; // parity of current block
    int i = 0; // the current vertex
    int blockstart = 0; // first vertex of the current block

    for (;;) {
        /**
         * First detect if we are at the end of a block (vertices of same parity)
         */
        bool blockended = false;
        if (i == nodecount()) {
            blockended = true;
        } else {
            if (disabled[i]) { i++; continue; }
            if (parity[i] != cur_parity) blockended = true;
        }

        if (blockended) {
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
                            justified[from] = false;
                            distraction[from] = false;
                            Q.push(from);
                            if (from < i) i = from; // afterwards, continue at lowest unjustified vertex
                        }
                    }
                }
#ifndef NDEBUG
                if (trace) logger << "restarting after finding distractions of prio " << priority(i-1) << std::endl;
#endif
                iterations++;
                blockchanged = false;
                if (i > blockstart) i = blockstart;
            } else {
                // We now know that the current strategy of all unjustified vertices of the block is justified
                for (int v=blockstart; v<i; v++) {
                    if (!disabled[v] and !justified[v]) {
                        justified[v] = true;
#ifndef NDEBUG
                        if (trace >= 2) logger << "\033[38;5;165;1mjustified\033[m " << label_vertex(v) << std::endl;
#endif
                    }
                }
            }

            if (i == nodecount()) break; // in case the last block didn't result in unjustifications

            // continue with the new block, update the current parity
            cur_parity = parity[i];
            blockstart = i;
        }

        // if the current vertex is justified, we don't need to check it
        if (justified[i]) {
            i++;
            continue;
        }

        // compute one step winner of <i> and update the strategy. Unlike FPJ, we
        // do not stop at the first winning move: we record *every* outgoing edge
        // to a co-winning vertex (the maximal set of moves valid at this point).
        clearStrategyEdges(i); // drop any edges recorded on a previous pass
        const int o = owner(i);
        int onestep_winner = 1 - o; // default: owner cannot reach a vertex good for itself
        int cnt = 0;
        int k = 0;
        for (auto curedge = outs(i); *curedge != -1; curedge++, k++) {
            const int to = *curedge;
            if (disabled[to]) continue;
            const int winner_to = parity[to] ^ distraction[to];
            if (winner_to == o) {
                onestep_winner = o;
                addStrategyEdge(i, k);
                cnt++;
            }
        }
        nstrat[i] = cnt;

        // evaluation stays the same
        if (cur_parity != onestep_winner) {
            Q.push(i); // add to the Queue, because all justified predecessors are now invalid
            distraction[i] = true;
            justified[i] = true;
            blockchanged = true;
#ifndef NDEBUG
            if (trace >= 2) logger << "\033[38;5;165;1mjustified*\033[m " << label_vertex(i) << std::endl;
#endif
        }

        i++;
    }

    // done: report the winner. Won vertices use the -1 sentinel; their winning
    // moves live in the multi-strategy (see strategyTargets/the verifier).
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, -1);
    }

    // free allocated data structures
    delete[] nstrat;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPJMSolver::run()
{
    // multi-strategy variant is sequential only
    runSeq();
}

}
