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
#include <climits>
#include <cstring>

#include "zlkwj.hpp"
#include "printf.hpp"

namespace pg {

/**
 * The Warsaw quasi-polynomial version of Zielonka's algorithm with justifications.
 *
 * ZLKWJ is ZLKW (see zlkw.cpp) with the implicit justifications of ZLKJ (see zlkj.cpp), after
 * Lapauw, Bruynooghe and Denecker, "Improving Parity Game Solvers with Justifications" (VMCAI
 * 2020). The loop of a call, the precisions, the reached-bottom flags and the size test are
 * those of ZLKW (and ZLKPP-WAW); what changes is what a step solves.
 *
 * The justification is the strategy: a vertex won by its owner is justified by its strategy
 * edge, a vertex won by the opponent of its owner by all its edges in the subgame. The heads of
 * a call have no justification (they stay in the region of the call). In ZLKW every step of a
 * call attracts the heads again and solves the whole rest. In ZLKWJ only the first step does;
 * a later step resets the pending vertices of the player <pl> of the call and the vertices won
 * by <pl> whose justification reaches them, attracts the reset vertices back to the vertices
 * still won by <pl> where it can, and solves only the rest of the reset vertices. Pending are
 * the vertices whose justification edge the opponent attracted (they follow another escape
 * meanwhile, as in ZLKJ), and the vertices whose winner is presumed: where a step ran out of
 * precision, ZLKW presumes the rest won by the player of that call, without computing a
 * strategy. Those vertices are marked <unjust>, and a call adds the unjust vertices of <pl> in
 * the result of its recursive call to its pending vertices. A step that runs out of precision
 * keeps the pending vertices won by <pl> and marks them unjust. So the extra step at precision
 * po, which ZLKW takes because the precision of the opponent ran out below, solves again
 * exactly the part where it ran out and what depends on it.
 *
 * Why the precision argument still holds: a step of ZLKWJ is a step of Zielonka's algorithm
 * that attracts to the heads and to the set J of vertices won by <pl> that are not reset. J is
 * won by <pl> in the subgame: its justification stays in J and avoids every presumed vertex,
 * so it is a winning strategy computed by exact steps. The two lemmas the argument uses (see
 * zlkw.cpp) hold for this target: the rest is the complement of an attractor of <pl>, and <pl>
 * wins every play that visits the attractor infinitely often (it visits the heads infinitely
 * often, or it reaches J and stays there). The size test uses the size of the rest that the
 * step solved, which contains every opponent dominion of the subgame.
 *
 * Two choices differ from ZLKW:
 * - The heads are the open vertices from the top down to the first open vertex of the other
 *   parity, as in ZLKJ: there is no Verver extension past an attracted vertex of the other
 *   parity. With it, a reset vertex attracted back to a lower head could close a cycle through
 *   such a vertex, which has a higher priority than the lower head and the parity of the
 *   opponent.
 * - A recursive call is for the player of the top vertex of its subgame. After the first step,
 *   or when the heads attracted the first vertex of the other parity, that can be <pl>; the
 *   call then gets the precision of <pl> for <pl> and that of the step for the opponent, which
 *   gives the guarantees the step needs.
 *
 * ZLKWJ reuses the results of calls on the same subgame with the same precisions as ZLKW does,
 * with the same filter, slots and budget (the computation of a call still depends only on its
 * subgame and its precisions); the stored result includes the unjust marks.
 *
 * Strategies are those of the justification; the heads get a successor in the region at the
 * end of their call. Presumed vertices get a successor in the subgame. As for ZLKW and ZLKQ,
 * the solutions are verified.
 */

static const int DIS = 0x80000000; // permanently disabled vertex
static const int BOT = 0x80000001; // open vertex (no escape counter)
static const int HEAD = -2;        // strategy of a head that has no strategy yet

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

static inline bool is_open(int rv)
{
    return rv < 0 and rv != DIS;
}

ZLKWJSolver::ZLKWJSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

ZLKWJSolver::~ZLKWJSolver()
{
}

/**
 * Attract for <pl> to the heads at and below <i>, with region id <r>, as attractHeads in ZLKJ:
 * the heads are the open vertices from <i> down to the first open vertex of the other parity.
 * Attracted vertices are added to <A> (also the FIFO queue), heads also to <H>.
 * Returns the highest open vertex after the attraction (the top of the rest), or -1.
 */
int
ZLKWJSolver::attractHeads(int i, const int pl, const int r, std::vector<int> &A, std::vector<int> &H)
{
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
        if ((priority(i)&1) != pl) break;

        region[i] = r;
        open.reset(i);
        winning[i] = pl;
        strategy[i] = HEAD;
        unjust.reset(i);
        H.push_back(i);
        A.push_back(i);
#ifndef NDEBUG
        if (trace >= 2) logger << "\033[36;1mhead\033[m " << label_vertex(i) << std::endl;
#endif
        i--;
    }

    // <i> is the first open vertex of the other parity; every vertex above it is now in a
    // region, outside the subgame or disabled
    for (; next < A.size(); next++) {
        const int cur = A[next];
        for (const int *e = ins(cur); *e != -1; e++) {
            const int from = *e;
            if (from > i) continue;
            const int rf = region[from];
            if (rf == DIS or rf >= 0) continue; // not open

            if (owner(from) == pl) {
                strategy[from] = cur;
            } else {
                // the escape counter: minus the number of successors that are not yet attracted
                int count = rf;
                if (count == BOT) {
                    count = 1;
                    for (const int *f = outs(from); *f != -1; f++) {
                        const int rt = region[*f];
                        if (rt == DIS) continue;
                        if (rt >= 0 and rt < r) continue; // removed from the subgame
                        count--;
                    }
                } else {
                    count++;
                }
                if (count != 0) {
                    region[from] = count;
                    continue;
                }
                strategy[from] = -1;
            }
            region[from] = r;
            open.reset(from);
            winning[from] = pl;
            unjust.reset(from);
            A.push_back(from);
        }
    }

    // the top of the rest
    while (i >= 0 and !is_open(region[i])) {
        open.reset(i);
        i = highest_at_most(open, i);
    }
    return i;
}

/**
 * The opponent of <pl> attracts, in the subgame of the call with region id <r>, the vertices
 * won by <pl> to its new vertices <O> (the vertices it won in the rest), as attractLosing in
 * ZLKJ. Attracted vertices are added to <O>. Vertices of <pl> whose justification edge goes to
 * an attracted vertex, but that have another escape, follow that escape and are added to <U>.
 * <H> holds the heads of the call. Returns the number of attracted vertices.
 */
long
ZLKWJSolver::attractLosing(const int r, const int pl, std::vector<int> &O, const std::vector<int> &H, std::vector<int> &U)
{
    // a successor of <v> in the subgame that is still won by <pl>, or -1
    auto escape = [&](int v) {
        for (const int *e = outs(v); *e != -1; e++) {
            const int to = *e;
            if (winning[to] == pl and region[to] >= r) return to;
        }
        return -1;
    };

    auto lose = [&](int v, int str) {
        region[v] = r;
        winning[v] = 1-pl;
        strategy[v] = str;
        unjust.reset(v);
        O.push_back(v);
    };

    // attract the vertices of <pl> that depend on <cur>, which is won by the opponent
    auto attract_to = [&](int cur) {
        for (const int *e = ins(cur); *e != -1; e++) {
            const int from = *e;
            if (winning[from] != pl) continue; // not won by <pl>
            if (region[from] < r) continue; // not in the subgame
            if (owner(from) != pl) {
                lose(from, cur); // justified by all its edges, so lost
                continue;
            }
            const int s = strategy[from];
            if (s >= 0 and s != cur) continue; // its justification edge is still won by <pl>
            const int to = escape(from);
            if (to == -1) {
                lose(from, -1);
            } else if (s != HEAD) {
                // lost its justification edge: pending, but follow another escape now
                strategy[from] = to;
                U.push_back(from);
            }
        }
    };

    /**
     * Only a head can be attracted directly: every other vertex won by <pl> is justified by
     * vertices won by <pl>, and the new vertices of the opponent were open before they were
     * solved. So either check the successors of the heads, or the predecessors of the new
     * vertices of the opponent; take the smaller set.
     */
    const size_t initial = O.size();
    if (initial <= H.size()) {
        for (size_t k = 0; k < initial; k++) attract_to(O[k]);
    } else {
        for (int v : H) {
            if (winning[v] != pl) continue; // lost in an earlier step
            if (owner(v) == pl) {
                if (escape(v) == -1) lose(v, -1);
            } else {
                for (const int *e = outs(v); *e != -1; e++) {
                    const int to = *e;
                    if (winning[to] == 1-pl and region[to] >= r) {
                        lose(v, to);
                        break;
                    }
                }
            }
        }
    }

    for (size_t k = initial; k < O.size(); k++) attract_to(O[k]);

    return (long)(O.size() - initial);
}

/**
 * Reset the pending vertices <U> that are still won by <pl> in the subgame of the call with
 * region id <r>, and the vertices won by <pl> whose justification reaches a reset vertex
 * (following the justification edges backwards, as resetClosure in ZLKJ). Heads have no
 * justification and stay. Fills <reset> and clears <U>.
 */
void
ZLKWJSolver::resetClosure(const int r, const int pl, std::vector<int> &U)
{
    auto do_reset = [&](int v) {
        winning[v] = -1;
        region[v] = BOT;
        strategy[v] = -1;
        unjust.reset(v);
        open.set(v);
        reset.push_back(v);
    };

    reset.clear();
    for (int v : U) {
        if (winning[v] == pl and region[v] >= r and strategy[v] != HEAD) do_reset(v);
    }
    U.clear();

    for (size_t k = 0; k < reset.size(); k++) {
        const int v = reset[k];
        for (const int *e = ins(v); *e != -1; e++) {
            const int from = *e;
            if (winning[from] != pl) continue;
            if (region[from] < r) continue; // not in the subgame
            const int s = strategy[from];
            if (s == HEAD) continue;
            if (owner(from) == pl and s != v) continue; // justified by another edge
            do_reset(from);
        }
    }
    resets += reset.size();
}

/**
 * Attract the reset vertices back to the vertices won by <pl> in the subgame of the call with
 * region id <r>, into the new region <r2>, as reattract in ZLKJ. Attracted vertices are added
 * to <A>. Returns the highest reset vertex that is still open, or -1.
 */
int
ZLKWJSolver::reattract(const int r, const int r2, const int pl, std::vector<int> &A)
{
    const size_t start = A.size();

    auto attract = [&](int v, int str) {
        region[v] = r2;
        open.reset(v);
        winning[v] = pl;
        strategy[v] = str;
        A.push_back(v);
    };

    // reset vertices with an edge to (if <pl> owns them) or only edges to (otherwise) vertices
    // won by <pl>; the others have a reset successor and are found by the FIFO
    for (int v : reset) {
        if (region[v] != BOT) continue;
        if (owner(v) == pl) {
            for (const int *e = outs(v); *e != -1; e++) {
                const int to = *e;
                if (winning[to] == pl and region[to] >= r) {
                    attract(v, to);
                    break;
                }
            }
        } else {
            bool all = true;
            for (const int *e = outs(v); *e != -1; e++) {
                if (is_open(region[*e])) {
                    all = false;
                    break;
                }
            }
            if (all) attract(v, -1);
        }
    }

    for (size_t next = start; next < A.size(); next++) {
        const int cur = A[next];
        for (const int *e = ins(cur); *e != -1; e++) {
            const int from = *e;
            const int rf = region[from];
            if (!is_open(rf)) continue;
            if (owner(from) == pl) {
                attract(from, cur);
            } else {
                // the escape counter is the negative number of successors that are open or in
                // region <r2> and have not yet attracted <from>
                int count = rf;
                if (count == BOT) {
                    count = 1;
                    for (const int *f = outs(from); *f != -1; f++) {
                        const int rt = region[*f];
                        if (is_open(rt) or rt == r2) count--;
                    }
                } else {
                    count++;
                }
                if (count == 0) attract(from, -1);
                else region[from] = count;
            }
        }
    }
    reattracted += A.size() - start;

    int top = -1;
    for (int v : reset) {
        if (is_open(region[v]) and v > top) top = v;
    }
    return top;
}

/**
 * Take the result of an earlier call on the same subgame, as reuse in ZLKW, with the unjust
 * marks of its vertices.
 */
bool
ZLKWJSolver::reuse(const long n, const long pc, const long po, const uint64_t key, bool &rb_cur, bool &rb_opo)
{
    auto range = index.equal_range(key);
    for (auto it = range.first; it != range.second; ++it) {
        Slot &M = slots[it->second];
        if ((long)M.n0 + (long)M.n1 != n) continue;

        const Prec *P = nullptr;
        for (const Prec &q : M.precs) {
            if (pc != q.pc and (pc < q.pc or q.rbc)) continue;
            if (po != q.po and (po < q.po or q.rbo)) continue;
            P = &q;
            break;
        }
        if (P == nullptr) continue;

        const int *V = M.V.data();
        long k = 0;
        while (k < n) {
            const int v = V[k] < 0 ? ~V[k] : V[k];
            if (!is_open(region[v])) break;
            k++;
        }
        if (k != n) continue; // another subgame with the same hash

        const int r = new_region();
        const int *S = V + n;
        W[0].resize(M.n0);
        W[1].resize(M.n1);
        for (k = 0; k < n; k++) {
            const int v = V[k] < 0 ? ~V[k] : V[k];
            const int w = k < (long)M.n0 ? 0 : 1;
            region[v] = r;
            open.reset(v);
            winning[v] = w;
            strategy[v] = S[k];
            unjust[v] = V[k] < 0;
            W[w][w ? k - M.n0 : k] = v;
        }
        rb_cur |= P->rbc;
        rb_opo |= P->rbo;
        res_hash[0] = M.h0;
        res_hash[1] = key - M.h0;
        touch(it->second);
        reused++;
        return true;
    }
    return false;
}

/**
 * Store the result of the call at depth <d> that just ended, as store in ZLKW.
 */
void
ZLKWJSolver::store(const int d, const long pc, const long po, const uint64_t key, const bool rbc, const bool rbo)
{
    const size_t n0 = W[0].size(), n1 = W[1].size(), n = n0 + n1;
    auto enc = [&](int v) { return unjust.test(v) ? ~v : v; };

    // the same result as an earlier call on this subgame: add the precisions to its slot
    auto range = index.equal_range(key);
    for (auto it = range.first; it != range.second; ++it) {
        Slot &M = slots[it->second];
        if (M.n0 != n0 or M.n1 != n1) continue;
        const int *V = M.V.data();
        size_t k = 0;
        while (k < n0 and V[k] == enc(W[0][k])) k++;
        while (k >= n0 and k < n and V[k] == enc(W[1][k-n0])) k++;
        if (k != n) continue;
        k = 0;
        while (k < n and V[n+k] == strategy[V[k] < 0 ? ~V[k] : V[k]]) k++;
        if (k != n) continue;
        M.precs.push_back(Prec{pc, po, rbc, rbo});
        touch(it->second);
        return;
    }

    if (depth_head.size() <= (size_t)d) {
        depth_head.resize(d+1, UINT32_MAX);
        depth_tail.resize(d+1, UINT32_MAX);
        depth_count.resize(d+1, 0);
    }

    // Over the budget, a depth with more than SLOTS slots gives up its least recently used
    // slots until the result fits, and a depth with SLOTS slots replaces its least recently
    // used one; otherwise it takes a new slot.
    uint32_t s = UINT32_MAX;
    while (depth_count[d] >= SLOTS and pool + 2*n > BUDGET) {
        const uint32_t t = depth_head[d];
        auto old = index.equal_range(slots[t].key);
        for (auto it = old.first; it != old.second; ++it) {
            if (it->second == t) {
                index.erase(it);
                break;
            }
        }
        const bool replace = depth_count[d] == SLOTS or pool - slots[t].V.size() + 2*n <= BUDGET;
        unlink(t);
        if (replace) {
            s = t;
            break;
        }
        pool -= slots[t].V.size();
        std::vector<int>().swap(slots[t].V);
        free_slots.push_back(t);
    }
    if (s == UINT32_MAX) {
        if (free_slots.empty()) {
            s = (uint32_t)slots.size();
            slots.emplace_back();
        } else {
            s = free_slots.back();
            free_slots.pop_back();
        }
    }
    Slot &M = slots[s];
    pool -= M.V.size();
    pool += 2*n;

    M.key = key;
    M.h0 = res_hash[0];
    M.n0 = (uint32_t)n0;
    M.n1 = (uint32_t)n1;
    M.V.resize(2*n);
    int *V = M.V.data();
    for (size_t k = 0; k < n0; k++) V[k] = enc(W[0][k]);
    for (size_t k = 0; k < n1; k++) V[n0+k] = enc(W[1][k]);
    for (size_t k = 0; k < n; k++) V[n+k] = strategy[k < n0 ? W[0][k] : W[1][k-n0]];
    M.precs.assign(1, Prec{pc, po, rbc, rbo});
    index.emplace(key, (uint32_t)s);
    link(s, d);
}

/**
 * The hash of the subgame of the <n> open vertices at or below <top>.
 */
uint64_t
ZLKWJSolver::hash_open(int v, const long n)
{
    uint64_t h = 0;
    long k = 0;
    while (k < n) {
        v = highest_at_most(open, v);
        if (is_open(region[v])) {
            h += vhash(v);
            k++;
        }
        v--;
    }
    return h;
}

/**
 * Solve the subgame of the open vertices (all at or below <top>; there are <n> > 0 of them) for
 * the player <pl> of <top>, with precision <pc> for <pl> and <po> for the opponent, as solve()
 * in ZLKW. <rb_cur> and <rb_opo> are set when the precision of <pl> or of the opponent runs
 * out below. If <hashed>, <hg> is the hash of the subgame. On return, W[0] and W[1] hold the
 * subgame split by the winner.
 */
void
ZLKWJSolver::solve(const int top, long n, const long pc, const long po, bool &rb_cur_out, bool &rb_opo_out, uint64_t hg, bool hashed)
{
    // as in ZLKW, subgames are hashed only after a precision ran out
    if (memo_on) {
        if (!hashed) {
            hg = hash_open(top, n);
            hashed = true;
        }
        if (reuse(n, pc, po, hg, rb_cur_out, rb_opo_out)) return;
    }

    calls++;
    const int pl = priority(top)&1;
    const uint64_t key = hg;
    uint64_t hy = 0; // the hash of what the opponent won (L.Y)
    bool rb_cur = false, rb_opo = false; // the flags of this call; passed up at the end

    const int d = depth++;
    struct DepthGuard { int &depth; ~DepthGuard() { depth--; } } guard{depth};
    if ((int)levels.size() <= d) levels.emplace_back(); // references to the other levels stay valid
    Level &L = levels[d];

    auto &Wm = W[pl];   // after a recursive call: the part of the rest won by <pl>
    auto &Wo = W[1-pl]; // after a recursive call: the part of the rest won by the opponent

#ifndef NDEBUG
    if (!Wm.empty() or !Wo.empty() or !L.G.empty() or !L.Y.empty() or !L.U.empty()) LOGIC_ERROR;
#endif

    const int r0 = new_region(); // the region of the heads; the subgame has region ids >= r0
    bool started = false; // whether a step attracted the heads

    // remove the stale entries and duplicates from L.G
    auto compact = [&]() {
        size_t k = 0;
        for (size_t t = 0; t < L.G.size(); t++) {
            const int v = L.G[t];
            if (winning[v] != pl or region[v] < r0 or seen.test(v)) continue;
            seen.set(v);
            L.G[k++] = v;
        }
        L.G.resize(k);
        for (size_t t = 0; t < k; t++) seen.reset(L.G[t]);
    };

    /**
     * One step of the loop with precision <p> for the opponent, as in ZLKW.
     * Sets <num_h> to the size of the rest that was solved, or LONG_MAX if the step ran out of
     * precision. Returns whether the opponent attracted vertices beyond what it won in the rest.
     */
    auto step = [&](const long p, long &num_h, bool &rb_step) -> bool {
        num_h = LONG_MAX;
        if (n == 0) return false;
        if (p < min_dominion) {
            rb_step = true;
            presumed++;
            memo_on = true;
            // the pending vertices stay won by <pl>, presumed
            for (int v : L.U) {
                if (winning[v] == pl and region[v] >= r0) unjust.set(v);
            }
            return false;
        }

        int j;
        long nh = 0;
        uint64_t hr = 0;
        bool hrv = false; // whether hr is the hash of the rest
        if (!started) {
            iterations++;
            started = true;
            j = attractHeads(top, pl, r0, L.G, L.H);
            nh = n - (long)L.G.size();
            if (hashed) {
                hr = hg;
                for (int v : L.G) hr -= vhash(v);
                hrv = true;
            }
        } else {
            resetClosure(r0, pl, L.U);
            if (reset.empty()) {
                num_h = 0;
                return false;
            }
            iterations++;
            j = reattract(r0, new_region(), pl, L.G);
            for (int v : reset) {
                if (is_open(region[v])) {
                    nh++;
                    if (memo_on) hr += vhash(v);
                }
            }
            hrv = memo_on;
        }

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "depth %d player %d step pc=%ld po=%ld: rest %ld\n", d, pl, pc, p, nh);
        if ((nh == 0) != (j == -1)) LOGIC_ERROR;
#endif

        num_h = nh;
        if (nh == 0) return false;

        // the player of the top of the rest solves it, with its precision for itself
        bool sc = false, so = false;
        const int x = priority(j)&1;
        if (x != pl) solve(j, nh, p, pc, sc, so, hr, hrv);
        else solve(j, nh, pc, p, sc, so, hr, hrv);
        const bool rbp = x == pl ? sc : so; // the precision of <pl> ran out below
        const bool rbq = x == pl ? so : sc; // the precision of the opponent ran out below
        rb_cur |= rbp;
        rb_step |= rbq;

        // the presumed vertices won by <pl> are pending
        if (rbp or rbq) {
            for (int v : Wm) {
                if (unjust.test(v)) L.U.push_back(v);
            }
        }
        L.G.insert(L.G.end(), Wm.begin(), Wm.end());
        Wm.clear();

        // the opponent attracts to what it won in the rest
        const size_t w0 = Wo.size();
        const long count = Wo.empty() ? -1 : attractLosing(r0, pl, Wo, L.H, L.U);

        // remove the opponent's attractor from the subgame
        n -= (long)Wo.size();
        if (hashed) {
            // the hash of the opponent's part of the rest comes with the result of the call
            uint64_t hw = res_hash[1-pl]; // (a call with a hash gives its recursive calls one)
            for (size_t k = w0; k < Wo.size(); k++) hw += vhash(Wo[k]);
            hg -= hw;
            hy += hw;
        }
        if (Wo.size() > L.Y.size()) L.Y.swap(Wo);
        L.Y.insert(L.Y.end(), Wo.begin(), Wo.end());
        Wo.clear();
        if (L.G.size() > 2*(size_t)n + 1024) compact();

        return count > 0;
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

    if (!started) {
        // Every step ran out of precision: the subgame (the open vertices) is presumed won by pl.
        int v = top;
        long k = 0;
        while (k < n) {
            v = highest_at_most(open, v);
#ifndef NDEBUG
            if (v < 0) LOGIC_ERROR;
#endif
            open.reset(v);
            if (!is_open(region[v])) continue;
            region[v] = r0;
            winning[v] = pl;
            strategy[v] = -1;
            unjust.set(v);
            L.G.push_back(v);
            k++;
        }
        for (int u : L.G) {
            if (owner(u) != pl) continue;
            for (const int *e = outs(u); *e != -1; e++) {
                const int to = *e;
                if (region[to] >= r0 and winning[to] == pl) {
                    strategy[u] = to;
                    break;
                }
            }
#ifndef NDEBUG
            if (strategy[u] == -1) LOGIC_ERROR;
#endif
        }
    } else {
        // pending vertices that no step solved again are presumed
        for (int v : L.U) {
            if (winning[v] == pl and region[v] >= r0) unjust.set(v);
        }

        // the heads get a strategy
        for (int v : L.H) {
            if (winning[v] != pl or strategy[v] != HEAD) continue;
            strategy[v] = -1;
            if (owner(v) != pl) continue;
            for (const int *e = outs(v); *e != -1; e++) {
                const int to = *e;
                if (winning[to] == pl and region[to] >= r0) {
                    strategy[v] = to;
                    break;
                }
            }
            if (strategy[v] == -1) LOGIC_ERROR;
        }

        // the vertices won by pl, once
        compact();
    }

#ifndef NDEBUG
    if (trace) fmt::printf(logger, "depth %d player %d: won %zu, lost %zu\n", d, pl, L.G.size(), L.Y.size());
    if ((long)L.G.size() != n) LOGIC_ERROR;
#endif

    // pass the result up
    Wm.swap(L.G);
    Wo.swap(L.Y);
    if (hashed) {
        res_hash[1-pl] = hy;
        res_hash[pl] = key - hy;
    }
    if (hashed and (rb_cur or rb_opo) and seen_before(key, pc, po)) store(d, pc, po, key, rb_cur, rb_opo);
    rb_cur_out |= rb_cur;
    rb_opo_out |= rb_opo;

    // keep the vectors of this depth, unless they are large
    for (auto *vec : {&L.H, &L.G, &L.Y, &L.U}) {
        if (vec->capacity() > 65536) std::vector<int>().swap(*vec);
        else vec->clear();
    }
}

void
ZLKWJSolver::run()
{
    const long N = nodecount();

    region = std::make_unique<int[]>(N);
    winning = std::make_unique<int8_t[]>(N);
    strategy = std::make_unique<int[]>(N);

    memset(winning.get(), -1, N);

    open = disabled;
    open.flip();
    unjust.resize(N);
    seen.resize(N);

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

    filter.assign(FILTER, 0);
    memo_on = false;

    bool rb_cur = false, rb_opo = false;
    solve(top, n, prec, prec, rb_cur, rb_opo, 0, false);

    for (long v = 0; v < N; v++) {
        if (disabled[v]) continue;
#ifndef NDEBUG
        if (winning[v] == -1) LOGIC_ERROR;
        if (owner(v) == winning[v] and strategy[v] < 0) LOGIC_ERROR;
#endif
        Solver::solve((int)v, winning[v], owner(v) == winning[v] ? strategy[v] : -1);
    }

    W[0].clear();
    W[1].clear();
    levels.clear();
    slots.clear();
    depth_head.clear();
    depth_tail.clear();
    depth_count.clear();
    free_slots.clear();
    index.clear();
    std::vector<uint64_t>().swap(filter);
    region.reset();
    winning.reset();
    strategy.reset();

    logger << "solved with " << iterations << " iterations (" << calls << " calls, " << presumed << " presumed, " << reused << " reused, " << resets << " resets, " << reattracted << " attracted back)." << std::endl;
}

}
