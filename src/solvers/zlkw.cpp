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

#include <algorithm>
#include <cassert>
#include <climits>
#include <cstring>

#include "zlkw.hpp"
#include "printf.hpp"

namespace pg {

/**
 * The Warsaw quasi-polynomial version of Zielonka's algorithm.
 *
 * This follows ZLK_WARSAW in zlkpp.cpp (the implementation by Pawel Parys that accompanies
 * Lehtinen, Parys, Schewe, Wojtczak 2021). A call solve(G, pc, po) for the player <pl> of the
 * top vertex of G, with precision <pc> for <pl> and <po> for the opponent, is Zielonka's loop of
 * steps. A step attracts to the heads, solves the rest H with a recursive call, lets the
 * opponent attract to the part of H that the opponent won, and removes that attractor:
 *
 *   repeat: step with opponent precision po/2, until the opponent attracts nothing new
 *   if |H| > po/2 and the opponent's precision ran out in the last step:
 *     step with opponent precision po; if the opponent attracted something new:
 *       repeat: step with opponent precision po/2, until the opponent attracts nothing new
 *
 * A step with opponent precision p < min_dominion does nothing: there is no opponent dominion
 * of size p, so the rest of G is presumed won by <pl>. The flags rb_cur and rb_opo record
 * whether that happened to the precision of <pl> or of the opponent anywhere below. If the
 * opponent's precision never ran out in the last step, the result is the exact Zielonka result
 * and the call stops (the reached-bottom exit); if |H| <= po/2, the halved precision already
 * covers every opponent dominion in H (the size test of Lehtinen et al.). The precision
 * bookkeeping, the flags and the loop control are exactly those of ZLKPP-WAW.
 *
 * The difference with ZLKPP-WAW is in how the steps are computed, with the techniques of ZLK:
 *
 * [I] Data structures of ZLK. Every vertex has a region id or is open (with an escape counter
 *     for the opponent's vertices during an attractor); region ids grow, so a call's subgame
 *     is the set of open vertices when it starts, and the vertices with an id >= its first id
 *     afterwards. Attraction is FIFO over in-edges with escape counters (each in-edge visited
 *     once), the head search skips closed vertices with the <open> bitset, and the results of
 *     a call are passed up in two vectors W[0], W[1] that are swapped rather than copied. The
 *     sets computed are the same as with Parys's linked list of the subgame.
 *
 * [A] On-the-fly compression of priorities. The heads of a step are the open vertices from the
 *     top down to the first open vertex of the other parity, attracted together, and the run
 *     continues below a vertex of the other parity that the heads attract (Maks Verver's
 *     optimisation, as in ZLK). One level of the recursion thus covers a run of priorities of
 *     one parity, where ZLKPP-WAW makes one level per priority (also for absent priorities),
 *     each level halving a precision. This changes the recursion tree, so it is an algorithmic
 *     change; it keeps the bound (never more levels) and the correctness: the lemmas behind
 *     the precision argument use only that (1) H is the complement of an attractor of the
 *     player of the call, so an x-dominion of G restricted to H is an x-dominion of H, and
 *     (2) if D is a nonempty dominion of the opponent in G, then H contains a nonempty
 *     opponent dominion inside D, because the player of the call wins every play that visits
 *     the region infinitely often. (2) holds for these heads: every vertex of the other parity
 *     above the lowest head is in the region, attracted by heads of a higher priority, so a
 *     play that keeps visiting the region but eventually avoids the attractor of the heads
 *     above some point sees only priorities of the call's parity at that level.
 *     The same compression is used by ZLKQ.
 *
 * [O] When the opponent attracts nothing new, the step leaves the subgame as it is and the
 *     loop ends; when the opponent won nothing in H, there is nothing to attract. Both are
 *     part of ZLKPP-WAW (do_step returns whether the attractor grew).
 *
 * [O] The opponent's attractor starts from the heads only (as in ZLK): every other vertex of
 *     the region was attracted to the heads, and the part of H that the call's player won is a
 *     trap for the opponent there, so only a head can be attracted directly. This gives the
 *     same attractor with less work.
 *
 * [G] min_dominion and the initial precision use the enabled vertices only (ZLKPP uses all
 *     vertices). A dominion of the subgame that Oink asks us to solve consists of enabled
 *     vertices, so this is as sound. It changes the precisions, and so possibly intermediate
 *     results, but only when Oink's preprocessing disabled vertices.
 *
 * [I] After a step, the vertices that stay in the subgame keep their region ids; they are
 *     opened again at the start of the next step that does work. If the precision runs out
 *     instead, they are the result as they are.
 *
 * Strategies are computed as in ZLK: attractor strategies, the strategies of the recursive
 * calls, and at the end of a call every vertex of the player of the call whose strategy is not
 * in the call's region gets a successor in the region (the heads, and vertices whose successor
 * the opponent attracted). Where a step presumed a region won (precision ran out), that region
 * has no computed strategy, and gets a successor in the region; the precision argument proves
 * the winning regions, not these strategies. As for ZLKQ, the solutions are verified.
 */

static const int DIS = 0x80000000; // permanently disabled vertex
static const int BOT = 0x80000001; // open vertex (no escape counter)

/**
 * The highest vertex in <b> that is at most <v>, or -1 if there is none.
 */
static inline int highest_at_most(const bitset &b, int v)
{
    if (v < 0) return -1;
    const uint64_t *p = b.data();
    long w = v/64;
    uint64_t m = p[w] & (~uint64_t(0) >> (63 - v%64));
    while (m == 0) {
        if (--w < 0) return -1;
        m = p[w];
    }
    return (int)(w*64 + bsr(m));
}

ZLKWSolver::ZLKWSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

ZLKWSolver::~ZLKWSolver()
{
}

/**
 * Attract for <pl> to the heads at and below <i>, with region id <r>, as attractExt in ZLK.
 * The heads are the open vertices from <i> down to the first open vertex of the other parity
 * that is not attracted; consecutive heads are added at once and attracted together, in FIFO
 * order with <A> as the queue. The heads are also added to <H>.
 * Returns the highest open vertex after the attraction (the top of the rest), or -1.
 */
int
ZLKWSolver::attract(int i, const int pl, const int r, std::vector<int> &A, std::vector<int> &H)
{
    while (true) {
        size_t next = A.size();

        while (true) {
            // With no open vertices left, the FIFO cannot attract anything else.
            if (i < 0) return -1;

            const int ri = region[i];
            if (ri == DIS or ri >= 0) {
                // not open: continue at the next open vertex
                open.reset(i);
                i = highest_at_most(open, i);
                continue;
            }

            if ((priority(i)&1) != pl) break; // until parity inversion

            region[i] = r;
            open.reset(i);
            winning[i] = pl;
            strategy[i] = -1; // heads get a strategy when the region is final
            H.push_back(i);
            A.push_back(i);

#ifndef NDEBUG
            if (trace >= 2) logger << "\033[36;1mhead\033[m " << label_vertex(i) << std::endl;
#endif
            i--;
        }

        // No new heads: <i> ends the run.
        if (next == A.size()) break;

        while (next < A.size()) {
            const int cur = A[next++];

            for (auto curedge = ins(cur); *curedge != -1; curedge++) {
                const int from = *curedge;
                if (from > i) continue; // above <i>: nothing is open there
                const int rf = region[from];
                if (rf == DIS or rf >= 0) continue; // not open

                if (owner(from) == pl) {
                    region[from] = r;
                    open.reset(from);
                    winning[from] = pl;
                    strategy[from] = cur;
                    A.push_back(from);
                } else {
                    // the escape counter: minus the number of successors that are not yet attracted
                    int count = rf;
                    if (count == BOT) {
                        count = 1;
                        for (auto e = outs(from); *e != -1; e++) {
                            const int rt = region[*e];
                            if (rt == DIS) continue;
                            if (rt >= 0 and rt < r) continue; // removed from the subgame
                            count--;
                        }
                    } else {
                        count++;
                    }
                    if (count == 0) {
                        region[from] = r;
                        open.reset(from);
                        winning[from] = pl;
                        strategy[from] = -1;
                        A.push_back(from);
                    } else {
                        region[from] = count;
                    }
                }
            }
        }

        // the heads may have attracted <i>; if so, the run of heads continues below it
    }

    return i;
}

/**
 * The opponent of <pl> attracts, in the subgame of the step with region id <r>, the vertices
 * won by <pl> to the vertices of <R> (won by the opponent), as attractLosing in ZLK.
 * <H> holds the heads of the region of the step. The attracted vertices are added to <R>.
 * Returns the number of attracted vertices.
 */
long
ZLKWSolver::attractLosing(const int r, const int pl, const std::vector<int> &H, std::vector<int> &R)
{
    // does <v> (owned by <pl>) have a successor in the subgame that is still won by <pl>?
    auto can_escape = [&](int v) {
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            const int to = *curedge;
            if (region[to] < r) continue; // not in the subgame, or disabled
            if (winning[to] != pl) continue; // not an escape
            return true;
        }
        return false;
    };

    // move <v> to the opponent with strategy <str> and add it to <R> (the FIFO queue)
    auto lose = [&](int v, int str) {
        region[v] = r;
        winning[v] = 1-pl;
        strategy[v] = str;
        R.push_back(v);
    };

    // attract the predecessors of <cur> (won by the opponent) that are won by <pl>
    auto attract_to = [&](int cur) {
        for (auto curedge = ins(cur); *curedge != -1; curedge++) {
            const int from = *curedge;
            if (region[from] < r) continue; // not in the subgame, or disabled
            if (winning[from] != pl) continue; // already lost
            if (owner(from) != pl) lose(from, cur);
            else if (!can_escape(from)) lose(from, -1);
        }
    };

    /**
     * Only a head can be attracted directly: every other vertex of the region was attracted to
     * the heads, so it has a strategy into the region (if <pl> owns it) or all its successors in
     * the subgame are in the region (if the opponent owns it); and the part of the rest won by
     * <pl> is a trap for the opponent there. So either check the predecessors of <R> or the
     * successors of the heads; take the smaller set.
     */
    const size_t initial = R.size();
    if (initial <= H.size()) {
        for (size_t k = 0; k < initial; k++) attract_to(R[k]);
    } else {
        for (int v : H) {
            if (owner(v) == pl) {
                if (!can_escape(v)) lose(v, -1);
            } else {
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    const int to = *curedge;
                    if (region[to] < r) continue; // not in the subgame, or disabled
                    if (winning[to] == pl) continue; // not attracting
                    lose(v, to);
                    break;
                }
            }
        }
    }

    for (size_t k = initial; k < R.size(); k++) attract_to(R[k]);

    return (long)(R.size() - initial);
}

/**
 * Solve the subgame of the open vertices (all at or below <top>; there are <n> > 0 of them) for
 * the player <pl>, with precision <pc> for <pl> and <po> for the opponent, as solve() with
 * ZLK_WARSAW in zlkpp.cpp. <rb_cur> and <rb_opo> are set when the precision of <pl> or of the
 * opponent runs out below. On return, W[0] and W[1] hold the subgame split by the winner.
 */
void
ZLKWSolver::solve(const int top, const int pl, long n, const long pc, const long po, bool &rb_cur, bool &rb_opo)
{
    calls++;

    const int d = depth++;
    struct DepthGuard { int &depth; ~DepthGuard() { depth--; } } guard{depth};
    if ((int)levels.size() <= d) levels.emplace_back(); // references to the other levels stay valid
    Level &L = levels[d];

    auto &Wm = W[pl];   // after a step: the part of the rest won by <pl>
    auto &Wo = W[1-pl]; // after a step: the part of the rest won by the opponent

#ifndef NDEBUG
    if (!Wm.empty() or !Wo.empty() or !L.G.empty() or !L.Y.empty()) LOGIC_ERROR;
#endif

    // Whether L.G holds the subgame. Initially the subgame is the set of open vertices; after a
    // step it is the list L.G, whose vertices keep their region ids until the next step needs
    // them open again (a step that runs out of precision does not).
    bool listed = false;
    bool changed = false; // whether the opponent attracted vertices in the last step
    int r = -1; // the region id of the last step

    /**
     * One step of the loop with precision <p> for the opponent, as do_step in zlkpp.cpp.
     * Sets <num_h> to the size of the rest, or LONG_MAX if the step did nothing.
     * Returns whether the opponent attracted vertices beyond what it won in the rest.
     */
    auto step = [&](const long p, long &num_h, bool &rb_step) -> bool {
        num_h = LONG_MAX;
        if (n == 0) return false;
        if (p < min_dominion) {
            rb_step = true;
            presumed++;
            return false;
        }
        iterations++;

        if (listed) {
            for (int v : L.G) {
                region[v] = BOT;
                open.set(v);
            }
            L.G.clear();
            listed = false;
        }

        r = new_region();
        L.A.clear();
        L.H.clear();
        const int j = attract(top, pl, r, L.A, L.H);
        const long nh = n - (long)L.A.size();

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "depth %d player %d step pc=%ld po=%ld: region of %zu vertices, rest %ld\n", d, pl, pc, p, L.A.size(), nh);
        if ((nh == 0) != (j == -1)) LOGIC_ERROR;
#endif

        // the opponent solves the rest with its precision <p> and ours <pc>
        if (nh > 0) solve(j, 1-pl, nh, p, pc, rb_step, rb_cur);
        num_h = nh;

        // the opponent attracts to what it won in the rest
        const long count = Wo.empty() ? -1 : attractLosing(r, pl, L.H, Wo);

        // remove the opponent's attractor from the subgame (merge the smaller vector into the larger)
        n -= (long)Wo.size();
        if (Wo.size() > L.Y.size()) L.Y.swap(Wo);
        L.Y.insert(L.Y.end(), Wo.begin(), Wo.end());
        Wo.clear();

        // the subgame is now the region and the part of the rest won by <pl>, minus what the
        // opponent attracted
        L.G.swap(Wm);
        L.G.insert(L.G.end(), L.A.begin(), L.A.end());
        if (count > 0) {
            L.G.erase(std::remove_if(L.G.begin(), L.G.end(), [&](int v) { return winning[v] != pl; }), L.G.end());
        }
        listed = true;
        changed = count > 0;

#ifndef NDEBUG
        if ((long)L.G.size() != n) LOGIC_ERROR;
#endif

        return changed;
    };

    /**
     * The loop of ZLK_WARSAW.
     */
    long num_h;
    bool atr_greater = true, rb_opo_last = false;
    while (atr_greater) {
        rb_opo_last = false;
        atr_greater = step(po/2, num_h, rb_opo_last);
        rb_opo |= rb_opo_last;
    }
    if (!(num_h <= po/2 or !rb_opo_last)) {
        if (step(po, num_h, rb_opo)) {
            while (step(po/2, num_h, rb_opo)) { }
        }
    }

    /**
     * The rest of the subgame is won by <pl>. Give every vertex of <pl> in it a strategy in it.
     * If the last step attracted nothing for the opponent and the loop ended there, the rest is
     * its region and the part of the lower subgame won by <pl>, and only the heads of the region
     * lack a strategy (as in ZLK); otherwise a step ran out of precision, and all vertices of
     * <pl> are checked (the opponent may have attracted their successors).
     */
    if (!listed) {
        // Every step ran out of precision: the subgame (the open vertices) is presumed won by <pl>.
        r = new_region();
        int v = top;
        for (long k = 0; k < n; k++) {
            v = highest_at_most(open, v);
#ifndef NDEBUG
            if (v < 0) LOGIC_ERROR;
#endif
            open.reset(v);
            region[v] = r;
            winning[v] = pl;
            strategy[v] = -1;
            L.G.push_back(v);
        }
        changed = true;
    }

    for (int v : changed ? L.G : L.A) {
        if (owner(v) != pl) continue;
        const int s = strategy[v];
        if (s != -1 and region[s] >= r and winning[s] == pl) continue; // in the region
        strategy[v] = -1;
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            const int to = *curedge;
            if (region[to] >= r and winning[to] == pl) {
                strategy[v] = to;
                break;
            }
        }
#ifndef NDEBUG
        if (strategy[v] == -1) LOGIC_ERROR;
#endif
    }

#ifndef NDEBUG
    if (trace) fmt::printf(logger, "depth %d player %d: won %zu, lost %zu\n", d, pl, L.G.size(), L.Y.size());
#endif

    // pass the result up
    Wm.swap(L.G);
    Wo.swap(L.Y);

    // keep the vectors of this depth, unless they are large
    for (auto *vec : {&L.A, &L.H, &L.G, &L.Y}) {
        if (vec->capacity() > 65536) std::vector<int>().swap(*vec);
        else vec->clear();
    }
}

void
ZLKWSolver::run()
{
    const long N = nodecount();

    region = std::make_unique<int[]>(N);
    winning = std::make_unique<int8_t[]>(N);
    strategy = std::make_unique<int[]>(N);

    memset(winning.get(), -1, N);

    open = disabled;
    open.flip();

    long n = 0;
    int top = -1;
    for (long v = 0; v < N; v++) {
        strategy[v] = -1;
        if (disabled[v]) {
            region[v] = DIS;
        } else {
            region[v] = BOT;
            top = (int)v;
            n++;
        }
    }
    if (n == 0) return;

    // the smallest possible dominion has one vertex if there is a self-loop, otherwise two
    min_dominion = 2;
    for (long v = 0; v < N and min_dominion == 2; v++) {
        if (disabled[v]) continue;
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            if (*curedge == v) {
                min_dominion = 1;
                break;
            }
        }
    }

    // the precision 2^k-1 >= n, as in ZLKPP
    long prec = 1;
    while (prec < n) prec = 2*prec+1;

    bool rb_cur = false, rb_opo = false;
    solve(top, priority(top)&1, n, prec, prec, rb_cur, rb_opo);

    for (long v = 0; v < N; v++) {
        if (disabled[v]) continue;
#ifndef NDEBUG
        if (winning[v] == -1) LOGIC_ERROR;
#endif
        Solver::solve((int)v, winning[v], strategy[v]);
    }

    W[0].clear();
    W[1].clear();
    levels.clear();
    region.reset();
    winning.reset();
    strategy.reset();

    logger << "solved with " << iterations << " iterations (" << calls << " calls, " << presumed << " presumed)." << std::endl;
}

}
