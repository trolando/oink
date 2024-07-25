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

namespace pg {

FPISolver::FPISolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

FPISolver::~FPISolver()
{
}

TASK_3(int, update_block_rec, FPISolver*, solver, int, i, int, n)
{
    return solver->update_block_rec(__lace_worker, __lace_dq_head, i, n);
}

int
FPISolver::update_block_rec(WorkerP* __lace_worker, Task* __lace_dq_head, int i, int n)
{
    if (n>128) {
        // because dynamic bitset is not thread safe, work in blocks of 64...
        if (i&127) {
            // start not yet aligned
            int N = 128 - (i&127);
            SPAWN(update_block_rec, this, i+N, n-N);
            int a = this->updateBlock(i, N);
            int b = SYNC(update_block_rec);
            return a+b;
        } else {
            int N = (n/128)*64;
            SPAWN(update_block_rec, this, i+N, n-N);
            int a = CALL(update_block_rec, this, i, N);
            int b = SYNC(update_block_rec);
            return a+b;
        }
    } else {
        return this->updateBlock(i, n);
    }
}

VOID_TASK_4(freeze_thaw_reset_rec, FPISolver*, solver, int, i, int, n, int, p)
{
    if (n>128) {
        // because dynamic bitset is not thread safe, work in blocks of 64...
        if (i&127) {
            // start not yet aligned
            int N = 128 - (i&127);
            SPAWN(freeze_thaw_reset_rec, solver, i+N, n-N, p);
            solver->freezeThawReset(i, N, p);
            SYNC(freeze_thaw_reset_rec);
        } else {
            int N = (n/128)*64;
            SPAWN(freeze_thaw_reset_rec, solver, i+N, n-N, p);
            CALL(freeze_thaw_reset_rec, solver, i, N, p);
            SYNC(freeze_thaw_reset_rec);
        }
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
        if (frozen[i]) continue;
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
        if (frozen[i] >= p) continue; // already frozen

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
        } else if (distraction[i]) {
            if (parity[i] == pl) {
                frozen[i] = p;
#ifndef NDEBUG
                if (trace >= 2) logger << "\033[38;5;51;1mfreeze\033[m " << label_vertex(i) << " at priority " << p << std::endl;
#endif
            } else {
                distraction[i] = 0;
#ifndef NDEBUG
                if (trace >= 2) logger << "\033[31;1mresetting\033[m " << label_vertex(i) << std::endl;
#endif
            }
        } else if (parity[i] != pl) {
            frozen[i] = p;
#ifndef NDEBUG
            if (trace >= 2) logger << "\033[38;5;51;1mfreeze\033[m " << label_vertex(i) << " at priority " << p << std::endl;
#endif
        }
    }
}

VOID_TASK_1(fpi_run_par, FPISolver*, _this)
{
    _this->run_par(__lace_worker, __lace_dq_head);
}

void
FPISolver::run_par(WorkerP* __lace_worker, Task* __lace_dq_head)
{
    int d = this->priority(this->nodecount()-1);
    int *p_start = new int[d+1];
    int *p_len = new int[d+1];

    // initialize p_start, p_len, parity
    {
        int v=0;
        for (int p=0; p<=d; p++) {
            if (this->priority(v) == p) {
                p_start[p] = v;
                while (v < this->nodecount() and this->priority(v) == p) {
                    this->parity[v] = p&1;
                    v++;
                }
                p_len[p] = v - p_start[p];
            } else {
                p_start[p] = -1;
                p_len[p] = 0;
            }
        }
    }

    this->iterations = 1;
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
            this->iterations++;
#ifndef NDEBUG
            if (this->trace >= 2) this->logger << "restarting after finding distractions" << std::endl;
#endif
        } else {
            // nothing changed
            p++;
        }
    }

    delete[] p_start;
    delete[] p_len;
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

    RUN(fpi_run_par, this);

    // done
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, winner == owner(v) ? strategy[v] : -1);
    }

    // free allocated data structures
    delete[] strategy;
    delete[] frozen;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPISolver::runSeq()
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
        Solver::solve(v, winner, winner == owner(v) ? strategy[v] : -1);
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
FPISolver::run()
{
    if (lace_workers() != 0) {
        runPar();
    } else {
        runSeq();
    }
}

}
