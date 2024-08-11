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
#include <queue>
#include <stack>

#include "fpj.hpp"

namespace pg {

FPJSolver::FPJSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

FPJSolver::~FPJSolver()
{
}

/**
 * Variation 1: greedy justification
 */
void
FPJSolver::runSeqGreedy()
{
    /**
     * Allocate and initialize data structures
     */

    int *strategy = new int[nodecount()]; // the current strategy
    bitset justified(nodecount()); // whether a vertex is justified
    bitset distraction(nodecount()); // whether a vertex is won by the opponent

    bitset parity(nodecount()); // optimization: precompute the parity of every vertex's priority
    for (int v=0; v<nodecount(); v++) parity[v] = priority(v)&1;

    uintqueue Q(nodecount());

    /**
     * Initialize loop
     */

    bool blockchanged = false; // did the current block change (new distractions)
    int cur_parity = parity[0]; // parity of current block
    int i = 0; // the current vertex

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
                 * We have to reset all justified vertices that are now no longer justified...
                 */
                while (!Q.empty()) {
                    const int v = Q.pop();
                    for (auto curedge = ins(v); *curedge != -1; curedge++) {
                        const int from = *curedge;
                        if (justified[from]) {
                            if (strategy[from] == -1 or strategy[from] == v) {
#ifndef NDEBUG
                                if (trace >= 2) logger << "\033[31;1mresetting\033[m " << label_vertex(from) << std::endl;
#endif
                                justified[from] = false;
                                distraction[from] = false; // reset it
                                Q.push(from);
                                if (from < i) i = from; // afterwards, continue at lowest unjustified vertex
                            }
                        }
                    }
                }
#ifndef NDEBUG
                if (trace) logger << "restarting after finding distractions of prio " << priority(i-1) << std::endl;
#endif
                iterations++;
                blockchanged = false;
            }

            if (i == nodecount()) break; // in case the last block didn't result in unjustifications

            // continue with the new block, update the current parity
            cur_parity = parity[i];
        }

        // if the current vertex is justified, we don't need to check it
        if (justified[i]) { i++; continue; }
        justified[i] = true;

        // compute one step winner of <i> and update the strategy
        int onestep_winner, str = -1;
        if (owner(i) == 0) {
            // see if player Even can go to a vertex currently good for Even
            onestep_winner = 1;
            for (auto curedge = outs(i); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue;
                const int winner_to = parity[to] ^ distraction[to];
                if (winner_to == 0) {
                    // good for player Even
                    onestep_winner = 0;
                    // and set the strategy
                    str = to;
                    break;
                }
            }
        } else {
            // see if player Odd can go to a vertex currently good for Odd
            onestep_winner = 0;
            for (auto curedge = outs(i); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue;
                const int winner_to = parity[to] ^ distraction[to];
                if (winner_to == 1) {
                    // good for player Odd
                    onestep_winner = 1;
                    // and set the strategy
                    str = to;
                    break;
                }
            }
        }

        strategy[i] = str;

        // evaluation stays the same
        if (cur_parity != onestep_winner) {
            Q.push(i); // add to the Queue, because all justified predecessors are now invalid
            distraction[i] = true;
            blockchanged = true;
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[38;5;165;1mjustified*\033[m " << label_vertex(i);
                if (strategy[i] != -1) logger << " => " << label_vertex(strategy[i]);
                logger << std::endl;
            }
#endif
        } else {
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[38;5;165;1mjustified\033[m " << label_vertex(i);
                if (strategy[i] != -1) logger << " => " << label_vertex(strategy[i]);
                logger << std::endl;
            }
#endif
        }

        i++;
    }

    // done
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    // free allocated data structures
    delete[] strategy;

    logger << "solved with " << iterations << " iterations." << std::endl;
}


void
FPJSolver::runSeq()
{
    /**
     * Allocate and initialize data structures
     */

    int *strategy = new int[nodecount()]; // the current strategy
    bitset justified(nodecount()); // whether a vertex is justified
    bitset distraction(nodecount()); // whether a vertex is won by the opponent

    bitset parity(nodecount()); // optimization: precompute the parity of every vertex's priority
    for (int v=0; v<nodecount(); v++) parity[v] = priority(v)&1;

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
                                justified[from] = false;
                                distraction[from] = false; // reset it
                                Q.push(from);
                                if (from < i) i = from; // afterwards, continue at lowest unjustified vertex
                            }
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
                        if (trace >= 2) {
                            logger << "\033[38;5;165;1mjustified\033[m " << label_vertex(v);
                            if (strategy[v] != -1) logger << " => " << label_vertex(strategy[v]);
                            logger << std::endl;
                        }
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

        // compute one step winner of <i> and update the strategy
        int onestep_winner, str = -1;
        if (owner(i) == 0) {
            // see if player Even can go to a vertex currently good for Even
            onestep_winner = 1;
            for (auto curedge = outs(i); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue;
                const int winner_to = parity[to] ^ distraction[to];
                if (winner_to == 0) {
                    // good for player Even
                    onestep_winner = 0;
                    // and set the strategy
                    str = to;
                    break;
                }
            }
        } else {
            // see if player Odd can go to a vertex currently good for Odd
            onestep_winner = 0;
            for (auto curedge = outs(i); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue;
                const int winner_to = parity[to] ^ distraction[to];
                if (winner_to == 1) {
                    // good for player Odd
                    onestep_winner = 1;
                    // and set the strategy
                    str = to;
                    break;
                }
            }
        }

        strategy[i] = str;

        // evaluation stays the same
        if (cur_parity != onestep_winner) {
            Q.push(i); // add to the Queue, because all justified predecessors are now invalid
            distraction[i] = true;
            justified[i] = true;
            blockchanged = true;
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[38;5;165;1mjustified*\033[m " << label_vertex(i);
                if (strategy[i] != -1) logger << " => " << label_vertex(strategy[i]);
                logger << std::endl;
            }
#endif
        }

        i++;
    }

    // done
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    // free allocated data structures
    delete[] strategy;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPJSolver::run()
{
    if (greedy) runSeqGreedy();
    else runSeq();
}

}
