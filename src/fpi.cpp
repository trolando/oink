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

#include <cassert>
#include <cstring>
#include <unistd.h>

#include "fpi.hpp"
#include "uintqueue.hpp"

#define ANTIFREEZE 0

namespace pg {

FPISolver::FPISolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

FPISolver::~FPISolver()
{
}

TASK_3(int, update_block_rec, FPISolver*, solver, int, i, int, n)
{
    if (n>128) {
        // because dynamic bitset is not thread safe, work in blocks of 64...
        int N = (n/128)*64;
        SPAWN(update_block_rec, solver, i+N, n-N);
        int a = CALL(update_block_rec, solver, i, N);
        int b = SYNC(update_block_rec);
        return a+b;
    } else {
        return solver->updateBlock(i, n);
    }
}

VOID_TASK_4(freeze_thaw_reset_rec, FPISolver*, solver, int, i, int, n, int, p)
{
    if (n>128) {
        // because dynamic bitset is not thread safe, work in blocks of 64...
        int N = (n/128)*64;
        SPAWN(freeze_thaw_reset_rec, solver, i+N, n-N, p);
        CALL(freeze_thaw_reset_rec, solver, i, N, p);
        SYNC(freeze_thaw_reset_rec);
    } else {
        solver->freezeThawReset(i, n, p);
    }
}

int
FPISolver::updateBlock(int i, int n)
{
    int res = 0;
    for (; n != 0; i++, n--) {
        if (disabled[i]) continue;
        if (ANTIFREEZE == 0 and frozen[i]) continue;
        if (distraction[i]) continue;

        // update whether current vertex <i> is a distraction by computing the one step winner
        int onestep_winner;
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
                    if (!frozen[i]) strategy[i] = to;
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
                    if (!frozen[i]) strategy[i] = to;
                    break;
                }
            }
        }
        if (parity[i] != onestep_winner) {
            distraction[i] = true;
            res++;
#ifndef NDEBUG
            if (trace >= 2) logger << "vertex " << label_vertex(i) << " is now a distraction (won by " << onestep_winner << ")" << std::endl;
#endif
        }
    }
    return res;
}

/**
 * Called after a vertex of priority p is flipped
 * Check for vertex i until (i+n) to update.
 */
void
FPISolver::freezeThawReset(int i, int n, int p)
{
    const int pl = p&1;
    for (; n != 0; i++, n--) {
        if (disabled[i]) continue; // not in the subgame
        if (ANTIFREEZE==0 and frozen[i] >= p) continue; // already frozen

        if (frozen[i]) {
            if ((frozen[i]&1) == pl) {
                frozen[i] = p;
            } else {
                frozen[i] = 0;
                distraction[i] = 0;
#ifndef NDEBUG
                if (trace >= 2) logger << "\033[38;5;202;1mthaw\033[m " << label_vertex(i) << std::endl;
#endif
            }
        } else if ((parity[i] ^ distraction[i]) != pl) {
            frozen[i] = p;
#ifndef NDEBUG
            if (trace >= 2) logger << "\033[38;5;51;1mfreeze\033[m " << label_vertex(i) << " at priority " << p << std::endl;
#endif
        } else if (distraction[i]) {
            distraction[i] = 0;
#ifndef NDEBUG
            if (trace >= 2) logger << "\033[31;1mresetting\033[m " << label_vertex(i) << std::endl;
#endif
        }
    }
}

void
FPISolver::runPar()
{
    /**
     * Allocate and initialize data structures
     */

    parity.resize(nodecount());
    distraction.resize(nodecount());
    parity.reset();
    distraction.reset();
    strategy = new int[nodecount()]; // the current strategy for winning the game
    frozen = new int[nodecount()]; // records for every vertex at which level it is frozen (or 0 if not frozen)

    memset(frozen, 0, sizeof(int[nodecount()])); // initially no vertex is frozen (we don't freeze at level 0)

    int d = priority(nodecount()-1);
    int *p_start = new int[d+1];
    int *p_len = new int[d+1];

    // initialize p_start, p_len, parity
    {
        int v=0;
        for (int p=0; p<=d; p++) {
            if (priority(v) == p) {
                p_start[p] = v;
                while (v < nodecount() and priority(v) == p) {
                    parity[v] = p&1;
                    v++;
                }
                p_len[p] = v - p_start[p];
            } else {
                p_start[p] = -1;
                p_len[p] = 0;
            }
        }
    }

    LACE_ME;

    iterations = 1;
    int p = 0;
    while (p <= d) {
        if (p_len[p] == 0) {
            p++;
        } else if (CALL(update_block_rec, this, p_start[p], p_len[p])) {
            // something changed, freeze and reset
            if (p != 0) {
                CALL(freeze_thaw_reset_rec, this, 0, p_start[p], p);
                p = 0;
            }
            iterations++;
#ifndef NDEBUG
            if (trace >= 2) logger << "restarting after finding distractions" << std::endl;
#endif
        } else {
            // nothing changed
            p++;
        }
    }

    // done
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        oink->solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    // free allocated data structures
    delete[] strategy;
    delete[] frozen;
    delete[] p_start;
    delete[] p_len;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPISolver::runSeq2()
{
    /**
     * Allocate and initialize data structures
     */
    distraction.resize(nodecount());
    strategy = new int[nodecount()]; // the current strategy for winning the game
    frozen = new int[nodecount()]; // records for every vertex at which level it is frozen (or 0 if not frozen)
    memset(frozen, 0, sizeof(int[nodecount()])); // initially no vertex is frozen (we don't freeze at level 0)

    int d = priority(nodecount()-1);
    int *p_start = new int[d+1];
    int *p_len = new int[d+1];
    parity.resize(nodecount());

    /**
     * Initialize p_start, p_len, parity
     */
    int v=0;
    for (int p=0; p<=d; p++) {
        if (priority(v) == p) {
            p_start[p] = v;
            while (v < nodecount() and priority(v) == p) {
                parity[v] = p&1;
                v++;
            }
            p_len[p] = v - p_start[p];
        } else {
            p_start[p] = -1;
            p_len[p] = 0;
        }
    }

    /**
     * The main loop
     */
    iterations = 1;
    int p = 0;
    while (p <= d) {
        if (p_len[p] == 0 or updateBlock(p_start[p], p_len[p]) == 0) {
            p++;
            continue;
        } 

        if (p != 0) {
            // actually we don't freeze at priority 0 :-)
            freezeThawReset(0, p_start[p] /* +p_len[p] */, p);
            p = 0;
        }

        iterations++;
#ifndef NDEBUG
        if (trace >= 2) logger << "restarting after finding distractions" << std::endl;
#endif
    }

    /**
     * Done, now tell Oink the solution
     */
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        oink->solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    /**
     * Free allocated data structures
     */
    delete[] strategy;
    delete[] frozen;
    delete[] p_start;
    delete[] p_len;

    logger << "solved with " << iterations << " iterations." << std::endl;

#ifndef NDEBUG
    if (trace) {
        logger << "Distractions:\n";
        for (int v=0; v<nodecount(); v++) {
            if (distraction[v]) logger << "\033[1;38:5:124m" << label_vertex(v) << "\033[m";
            else logger << "\033[1;38:5:34m" << label_vertex(v) << "\033[m";
            logger << std::endl;
        }
    }
#endif
}

void
FPISolver::runSeq()
{
    /**
     * Allocate and initialize data structures
     */

    bitset distraction(nodecount()); // the main data structure: is a vertex a distraction or not?

    int *strategy = new int[nodecount()]; // the current strategy for winning the game
    memset(strategy, -1, sizeof(int[nodecount()])); // initially set no strategy

    int *frozen = new int[nodecount()]; // records for every vertex at which level it is frozen (or 0 if not frozen)
    memset(frozen, 0, sizeof(int[nodecount()])); // initially no vertex is frozen (we don't freeze at level 0)

    bitset parity(nodecount()); // optimization: precompute the parity of every vertex's priority
    for (int v=0; v<nodecount(); v++) parity[v] = priority(v)&1;

    unsigned int iterations = 0; // record the number of iterations that we needed
    for (;;) {
        iterations++;

        /**
         * We update distractions from low to high.
         * We process per block of priorities.
         * When a block changes, this means we found more distractions.
         * We then freeze all (lower) vertices that the opponent (who wins the distraction) wins
         *   and reset all (lower) vertices that the current player <cur_pl> wins
         * When the block no longer changes, we thaw the frozen vertices of the current block
         *
         * The algorithm refines a partition into two winning regions.
         * FPI terminates when the winning regions are inductive, i.e.,
         * no player escapes in one step.
         */
        bool changed = false; // did the current block change (new distractions)
        int cur_pr = priority(0); // priority of the current block
        int cur_pl = cur_pr & 1;
        int i;

        for (i=0; i<nodecount();i++) {
            if (disabled[i]) continue;
            // if (priority(i) != cur_pr) { // no on-the-fly compression
            if (parity[i] != cur_pl) { // on-the-fly compression ((like Verver's Zielonka trick))
                // new block!
                if (changed) {
                    // break the loop to handle the end of a block
                    break;
                } else {
                    // the previous block did not change, continue with new block
                    cur_pr = priority(i);
                    cur_pl = cur_pr & 1;
                    // also thaw all vertices below the new <cur_pr>
                    if (cur_pr > 1) {
                        for (int v=0; v<i; v++) {
                            if (frozen[v] and frozen[v] < cur_pr) {
                                frozen[v] = 0;
#ifndef NDEBUG
                                if (trace >= 2) logger << "\033[38;5;202;1mthaw\033[m " << label_vertex(v) << std::endl;
#endif
                            }
                        }
                    }
                }
            }
            if (!frozen[i] and !distraction[i]) {
                // update whether current vertex <i> is a distraction by computing the one step winner
                int onestep_winner;
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
                            strategy[i] = to;
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
                            strategy[i] = to;
                            break;
                        }
                    }
                }
                if (cur_pl != onestep_winner) {
                    distraction[i] = true;
                    changed = true;
#ifndef NDEBUG
                    if (trace >= 2) logger << "vertex " << label_vertex(i) << " is now a distraction (won by " << onestep_winner << ")" << std::endl;
#endif
                }
            }
        }

        /**
         * If we are here, either we are at the end of a changed block
         * or we are at the end of the entire game (highest block, possibly unchanged)
         */

        if (!changed) break; // nothing changed, we're done

        // we found distractions, so we freeze or reset the lower game
        // (we can safely freeze also the vertices of the block)
        for (int v=0; v<i; v++) {
            if (disabled[v]) continue; // not in the subgame
            if (frozen[v]) continue; // already frozen
            int winner_v = parity[v] ^ distraction[v];
            if (winner_v != cur_pl) {
                // freeze (monotonically increases, btw)
                frozen[v] = cur_pr;
#ifndef NDEBUG
                if (trace >= 2) logger << "\033[38;5;51;1mfreeze\033[m " << label_vertex(v) << " at priority " << cur_pr << std::endl;
#endif
            } else {
                // reset
#ifndef NDEBUG
                if (trace >= 2 and distraction[v]) logger << "\033[31;1mresetting\033[m " << label_vertex(v) << std::endl;
#endif
                if (distraction[v]) distraction[v] = 0;
            }
        }

#ifndef NDEBUG
        if (trace) logger << "restarting after finding distractions of prio " << cur_pr << std::endl;
#endif
    }

    // done
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        oink->solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    // free allocated data structures
    delete[] strategy;
    delete[] frozen;

    logger << "solved with " << iterations << " iterations." << std::endl;
}


void
FPISolver::run()
{
    if (lace_workers() != 0) {
        runPar();
    } else {
        runSeq2();
    }
}

}
