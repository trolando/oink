/*
 * Copyright 2017-2024 Tom van Dijk
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

#include "fpim.hpp"

namespace pg {

FPIMSolver::FPIMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

FPIMSolver::~FPIMSolver()
{
}

/**
 * Update the block [i, i+n). For each undecided vertex, compute its one-step
 * winner. Unlike FPI, we do not stop at the first winning move: we record
 * *every* outgoing edge to a co-winning vertex in the solution's multi-strategy.
 * With <TrackChanges>, only the affected candidates are evaluated: the winning
 * moves of the others are unchanged, as none of their successors changed winner.
 */
template<bool TrackChanges>
int
FPIMSolver::updateBlockImpl(int i, int n, uint64_t &evaluated)
{
    if (n <= 0) return 0;

    int res = 0;
    const int end = i + n;
    const int last_word = (end - 1) >> 6;
    const uint64_t *cand = candidate.data();

    const int first_word = i >> 6;
    auto prev_word = [&](int word) {
        if constexpr (TrackChanges) {
            if (word < first_word) return word;
            const auto prev = dirty_words[word] ? size_t(word) : dirty_words.find_prev(word);
            return prev == bitset::npos || int(prev) < first_word ? first_word - 1 : int(prev);
        } else {
            return word;
        }
    };

    // visit only the candidates (enabled, not frozen, not a distraction) in [i, end),
    // 64 vertices at a time, from the last to the first: in the file order of generated
    // games (states in the order of exploration) a vertex mostly comes before its
    // successors, so a change then reaches its predecessors in the same pass
    for (int w = prev_word(last_word); w >= first_word; w = prev_word(w - 1)) {
        uint64_t bits = cand[w];
        if constexpr (TrackChanges) bits &= dirty.data()[w];
        if (w == (i >> 6)) bits &= ~uint64_t(0) << (i & 63);
        if (w == last_word and (end & 63) != 0) bits &= (uint64_t(1) << (end & 63)) - 1;

        while (bits) {
            const int v = (w << 6) + 63 - intrinsics::countl_zero64(bits);
            bits &= ~(uint64_t(1) << (v & 63));
            if constexpr (TrackChanges) dirty.reset(v);
            evaluated++;

            // Recompute the winning-edge mask while evaluating <v>.
            const int o = owner(v);
            const uint64_t* good = won[o].data();
            const int* edges = outs(v);
            const int degree = game.outcount(v);
            bool any_winning = false;
            auto winningBits = [&](int offset, unsigned count) {
                uint64_t bits = 0;
                for (unsigned i = 0; i < count; ++i) {
                    const int to = edges[offset + i];
                    bits |= ((good[to >> 6] >> (to & 63)) & 1) << i;
                }
                return bits;
            };
            // Most vertices fit in one mask; avoid the outer chunk loop.
            if (degree <= 64) {
                const uint64_t bits = winningBits(0, degree);
                replaceStrategyBits(v, 0, degree, bits);
                any_winning = bits != 0;
            } else {
                for (int offset = 0; offset < degree;) {
                    const unsigned count = std::min(64, degree - offset);
                    const uint64_t bits = winningBits(offset, count);
                    replaceStrategyBits(v, offset, count, bits);
                    any_winning |= bits != 0;
                    offset += count;
                }
            }
            const int onestep_winner = any_winning ? o : 1 - o;

            if (parity[v] != onestep_winner) {
                makeDistraction(v);
                if constexpr (TrackChanges) {
                    // Only predecessors can get other winning moves. Visit newly
                    // affected lower lanes in this word in the usual order.
                    for (auto edge = ins(v); *edge != -1; ++edge) {
                        const int from = *edge;
                        markDirty(from);
                        if ((from >> 6) == w && from < v && from >= i && candidate[from])
                            bits |= uint64_t(1) << (from & 63);
                    }
                }
                res++;
#ifndef NDEBUG
                if (trace >= 2) logger << "vertex " << label_vertex(v) << " is now a distraction (won by " << onestep_winner << ")" << std::endl;
#endif
            }
        }
        if constexpr (TrackChanges) {
            // Ineligible vertices are dirtied again if thaw/reset reactivates
            // them. Retain candidates outside this priority block in the word.
            dirty.data()[w] &= cand[w];
            if (!dirty.data()[w]) dirty_words.reset(w);
        }
    }
    return res;
}

int
FPIMSolver::updateBlock(int i, int n)
{
    uint64_t evaluated = 0;
    return updateBlockImpl<false>(i, n, evaluated);
}

/**
 * Called after a vertex of priority p is flipped
 * Check for vertex i until (i+n) to update.
 */
void
FPIMSolver::freezeThawReset(int i, int n, int p)
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
        if (track_changes) {
            // Thawed and reset vertices must be evaluated again even when none
            // of their successors changed. Only actual winner flips change the
            // winning moves of their predecessors.
            dirty.data()[w] = (dirty.data()[w] | clear) & CAND[w];
            if (dirty.data()[w]) dirty_words.set(w);
            else dirty_words.reset(w);
            for (uint64_t changed = D & clear; changed; changed &= changed - 1) {
                const int v = (w << 6) + intrinsics::countr_zero64(changed);
                for (auto edge = ins(v); *edge != -1; ++edge) markDirty(*edge);
            }
        }

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
FPIMSolver::initState()
{
    won[1] = parity;
    won[1] -= disabled;
    won[0] = ~parity;
    won[0] -= disabled;
    candidate = ~disabled;
}

/**
 * Run the main loop from priority <p> until the fixpoint is reached. Without
 * <TrackChanges>, stop and return true at a restart once change tracking is
 * expected to pay off. (Checking this for every block costs too much when the
 * blocks are small.)
 */
template<bool TrackChanges>
bool
FPIMSolver::iterate(int &p, const int d, const int *p_start, const int *p_len)
{
    const uint64_t track_work = 16 * uint64_t(nodecount());
    uint64_t evaluated = 0, changes = 0;
    while (p <= d) {
        const int res = p_len[p] == 0 ? 0 : updateBlockImpl<TrackChanges>(p_start[p], p_len[p], evaluated);
        if (res == 0) {
            p++;
            continue;
        }
        changes += res;

        if (p != 0) {
            // actually we don't freeze at priority 0 :-)
            freezeThawReset(0, p_start[p] /* +p_len[p] */, p);
            p = 0;
        }

        iterations++;
#ifndef NDEBUG
        if (trace >= 2) logger << "restarting after finding distractions" << std::endl;
#endif
        if constexpr (!TrackChanges) {
            if (evaluated >= track_work && changes * 8 < evaluated) return true;
        }
    }
    return false;
}

void
FPIMSolver::runSeq()
{
    /**
     * Allocate and initialize data structures
     */
    distraction.resize(nodecount());
    // records for every vertex at which level it is frozen (or 0 if not frozen), padded to a
    // multiple of 64 for freezeThawReset; initially no vertex is frozen (we don't freeze at level 0)
    const size_t frozen_size = (size_t(nodecount()) + 63) & ~size_t(63);
    frozen = new int[frozen_size]();

    // allocate the multi-strategy (set of all winning moves), unless a previous
    // call (e.g. on a different subgame) already did so
    if (!hasMultiStrategy()) initMultiStrategy();

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
    track_changes = false;

    /**
     * The main loop. As in FPI: once FPIM has evaluated many vertices with few
     * changes, it evaluates only the affected vertices, which needs the in-edges.
     */
    iterations = 1;
    int p = 0;
    if (iterate<false>(p, d, p_start, p_len)) {
        game.build_in_array(false);
        dirty = candidate;
        dirty_words.resize(candidate.block_count());
        dirty_words.reset();
        for (size_t word = 0; word < candidate.block_count(); ++word)
            if (candidate.data()[word]) dirty_words.set(word);
        track_changes = true;
        iterate<true>(p, d, p_start, p_len);
    }

    /**
     * Done, now tell Oink the solution. Won vertices report the -1 sentinel: the
     * winning moves live in the multi-strategy (see strategyTargets/the verifier).
     */
    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        const int winner = parity[v] ^ distraction[v];
        Solver::solve(v, winner, -1);
    }

    /**
     * Free allocated data structures
     */
    delete[] frozen;
    delete[] p_start;
    delete[] p_len;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

void
FPIMSolver::run()
{
    // multi-strategy variant is sequential only
    runSeq();
}

}
