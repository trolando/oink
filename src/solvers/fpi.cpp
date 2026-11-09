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
    if (n <= 0) return 0;

    int res = 0;
    const int end = i + n;
    const int last_word = (end - 1) >> 6;
    const uint64_t *cand = candidate.data();

    // visit only the candidates (enabled, not frozen, not a distraction) in [i, end),
    // 64 vertices at a time
    for (int w = i >> 6; w <= last_word; w++) {
        uint64_t bits = cand[w];
        if (w == (i >> 6)) bits &= ~uint64_t(0) << (i & 63);
        if (w == last_word and (end & 63) != 0) bits &= (uint64_t(1) << (end & 63)) - 1;

        while (bits) {
            const int v = (w << 6) + intrinsics::countr_zero64(bits);
            bits &= bits - 1;

            // update whether vertex <v> is a distraction by computing the one step winner:
            // can the owner go to a vertex currently won by the owner?
            const int pl = owner(v);
            const uint64_t *good = won[pl].data();
            int onestep_winner = 1 - pl;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                const int to = *curedge;
                if ((good[to >> 6] >> (to & 63)) & 1) {
                    onestep_winner = pl;
                    strategy[v] = to;
                    break;
                }
            }
            if (parity[v] != onestep_winner) {
                makeDistraction(v);
                res++;
#ifndef NDEBUG
                if (trace >= 2) logger << "vertex " << label_vertex(v) << " is now a distraction (won by " << onestep_winner << ")" << std::endl;
#endif
            }
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
    if (n <= 0) return;

    // Process 64 vertices (one bitset word) at a time. For every enabled vertex v in
    // [i, i+n) with frozen[v] < p, with pl = p&1:
    //  - frozen at a level of parity pl:       freeze at p
    //  - frozen at a level of the other parity: thaw (unfreeze, reset distraction)
    //  - not frozen, distraction, parity pl:    freeze at p
    //  - not frozen, distraction, other parity: reset distraction
    //  - not frozen, no distraction, other parity: freeze at p
    // The lane predicates on frozen[] are simple loops the compiler can vectorize; the
    // case analysis is done with 64-bit masks. frozen[] is padded to a multiple of 64.
    const uint64_t PL = (p & 1) ? ~uint64_t(0) : 0;
    const uint64_t *DIS = disabled.data();
    uint64_t *DIST = distraction.data();
    uint64_t *CAND = candidate.data();
    uint64_t *WON0 = won[0].data();
    uint64_t *WON1 = won[1].data();
    const uint64_t *PAR = parity.data();

    const int end = i + n;
    const int last_word = (end - 1) >> 6;
    for (int w = i >> 6; w <= last_word; w++) {
        uint64_t en = ~DIS[w];
        if (w == (i >> 6)) en &= ~uint64_t(0) << (i & 63);
        if (w == last_word and (end & 63) != 0) en &= (uint64_t(1) << (end & 63)) - 1;
        if (en == 0) continue;

        int *f = frozen + (w << 6);
        uint64_t below = 0, nz = 0, odd = 0;
        for (int j = 0; j < 64; j++) {
            const int fj = f[j];
            below |= uint64_t(fj < p) << j;
            nz |= uint64_t(fj != 0) << j;
            odd |= uint64_t(fj & 1) << j;
        }
        const uint64_t act = en & below;
        if (act == 0) continue;

        const uint64_t D = DIST[w];
        const uint64_t P = PAR[w];
        const uint64_t frozen_same = ~(odd ^ PL); // frozen at a level of parity pl
        const uint64_t par_same = ~(P ^ PL);      // vertex parity is pl

        const uint64_t refreeze = act & nz & frozen_same;
        const uint64_t thaw = act & nz & ~frozen_same;
        const uint64_t freeze_d = act & ~nz & D & par_same;
        const uint64_t reset_d = act & ~nz & D & ~par_same;
        const uint64_t freeze_c = act & ~nz & ~D & ~par_same;

        const uint64_t setp = refreeze | freeze_d | freeze_c; // frozen = p
        const uint64_t clear = thaw | reset_d;                // no longer a distraction, candidate again

        DIST[w] = D & ~clear;
        WON1[w] = (WON1[w] & ~clear) | (clear & P);
        WON0[w] = (WON0[w] & ~clear) | (clear & ~P);
        CAND[w] = (CAND[w] | clear) & ~freeze_c;

        if (setp | thaw) {
            for (int j = 0; j < 64; j++) {
                const int fj = f[j];
                f[j] = ((setp >> j) & 1) ? p : (((thaw >> j) & 1) ? 0 : fj);
            }
        }

#ifndef NDEBUG
        if (trace >= 2) {
            for (uint64_t m = thaw; m; m &= m - 1) logger << "\033[38;5;202;1mthaw\033[m " << label_vertex((w << 6) + intrinsics::countr_zero64(m)) << std::endl;
            for (uint64_t m = freeze_d | freeze_c; m; m &= m - 1) logger << "\033[38;5;51;1mfreeze\033[m " << label_vertex((w << 6) + intrinsics::countr_zero64(m)) << " at priority " << p << std::endl;
            for (uint64_t m = reset_d; m; m &= m - 1) logger << "\033[31;1mresetting\033[m " << label_vertex((w << 6) + intrinsics::countr_zero64(m)) << std::endl;
        }
#endif
    }
}

/**
 * Initialize the derived state from <parity> and <disabled>, with nothing frozen
 * and no distractions: every enabled vertex is a candidate, won by its parity.
 */
void
FPISolver::initState()
{
    won[1] = parity;
    won[1] -= disabled;
    won[0] = ~parity;
    won[0] -= disabled;
    candidate = ~disabled;
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
    this->initState();

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
    // records for every vertex at which level it is frozen (or 0 if not frozen), padded to a
    // multiple of 64 for freezeThawReset; initially no vertex is frozen (we don't freeze at level 0)
    const size_t frozen_size = (size_t(nodecount()) + 63) & ~size_t(63);
    frozen = new int[frozen_size]();

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
    // records for every vertex at which level it is frozen (or 0 if not frozen), padded to a
    // multiple of 64 for freezeThawReset; initially no vertex is frozen (we don't freeze at level 0)
    const size_t frozen_size = (size_t(nodecount()) + 63) & ~size_t(63);
    frozen = new int[frozen_size]();

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
    initState();

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
    if (parallel()) {
        runPar();
    } else {
        runSeq();
    }
}

}
