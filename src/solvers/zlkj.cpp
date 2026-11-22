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
 *
 * Zielonka's algorithm with justifications, after Algorithm 5 of Lapauw, Bruynooghe and
 * Denecker, "Improving Parity Game Solvers with Justifications" (VMCAI 2020), and Ruben
 * Lapauw's zlkjb in pull request trolando/oink#4.
 */

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "zlkj.hpp"
#include "printf.hpp"

namespace pg {

static const int DIS = 0x80000000; // permanently disabled vertex
static const int BOT = 0x80000001; // open vertex without an escape counter
static const int HEAD = -2;        // strategy of a head that has no strategy yet

#define KC "\033[36;1m"

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

ZLKJSolver::ZLKJSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

ZLKJSolver::~ZLKJSolver()
{
}

/**
 * Attract the heads of a new frame to region <r>: the open vertices from <i> down to the first
 * open vertex of the other parity. Attracted vertices are added to <A> (also the FIFO queue),
 * heads also to <H>. Returns the first open vertex of the other parity, or -1.
 */
int
ZLKJSolver::attractHeads(int i, int r, std::vector<int> &A, std::vector<int> &H)
{
    const int pl = priority(i)&1;
    size_t next = A.size();

    while (true) {
        // With no open vertices left, every predecessor of a head fails the test below
        // (from > i): the heads cannot attract anything else.
        if (i < 0) return -1;

        const int ri = region[i];
        if (ri == DIS or ri >= 0) {
            // not open: drop it from <open>, continue at the next open vertex
            open.reset(i);
            i = highest_at_most(open, i);
            continue;
        }
        if ((priority(i)&1) != pl) break;
        region[i] = r;
        open.reset(i);
        winning[i] = pl;
        strategy[i] = HEAD;
        H.push_back(i);
        A.push_back(i);
#ifndef NDEBUG
        if (trace >= 2) logger << KC"head\033[m " << label_vertex(i) << std::endl;
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
            if (rf == DIS or rf >= 0) continue; // cannot be attracted

            if (owner(from) == pl) {
                strategy[from] = cur;
            } else {
                // the escape counter is the negative number of successors that are open or in
                // this region and have not yet attracted <from>
                int count = rf;
                if (count == BOT) {
                    count = 1;
                    for (const int *f = outs(from); *f != -1; f++) {
                        const int rt = region[*f];
                        if (rt == DIS) continue;
                        if (rt >= 0 and rt < r) continue;
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
            A.push_back(from);
#ifndef NDEBUG
            if (trace >= 2) logger << KC"attracted\033[m " << label_vertex(from) << std::endl;
#endif
        }
    }

    return i;
}

/**
 * The opponent of <pl> attracts from the subgame of region <r>, starting with its new vertices
 * O[start..]. Attracted vertices are added to <O>. Vertices of <pl> whose justification edge
 * goes to an attracted vertex are collected in <lost> (and get another escape as strategy, so
 * that a later loss of that escape is noticed); they are reset afterwards.
 * Returns the number of attracted vertices.
 */
int
ZLKJSolver::attractLosing(int r, int pl, std::vector<int> &O, size_t start, std::vector<int> &H)
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
        O.push_back(v);
    };

    // attract the vertices of <pl> that depend on <cur>, which is won by the opponent
    auto attract_to = [&](int cur) {
        for (const int *e = ins(cur); *e != -1; e++) {
            const int from = *e;
            if (winning[from] != pl) continue; // not won by <pl>, or open
            if (region[from] < r) continue; // not in the subgame
            if (owner(from) != pl) {
                // justified by all its edges, so lost
#ifndef NDEBUG
                if (trace >= 2) logger << KC"attracted\033[m " << label_vertex(from) << std::endl;
#endif
                lose(from, cur);
                continue;
            }
            const int s = strategy[from];
            if (s >= 0 and s != cur) continue; // its justification edge is still won by <pl>
            const int to = escape(from);
            if (to == -1) {
#ifndef NDEBUG
                if (trace >= 2) logger << KC"forced\033[m " << label_vertex(from) << std::endl;
#endif
                lose(from, -1);
            } else if (s != HEAD) {
                // lost its justification edge: reset it later, but follow another escape now
                strategy[from] = to;
                lost.push_back(from);
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
    if (initial - start <= H.size()) {
        for (size_t k = start; k < initial; k++) {
            // (skip stale entries of vertices that were reset or lost in the subgame)
            if (winning[O[k]] == 1-pl) attract_to(O[k]);
        }
    } else {
        for (int v : H) {
            if (winning[v] != pl) continue; // lost in an earlier round
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

    return O.size() - initial;
}

/**
 * Reset the vertices won by <pl> in the subgame of region <r> whose justification reaches a
 * vertex that the opponent attracted: start with <lost>, then follow the justification edges
 * backwards. Heads have no justification and stay. Fills <reset>.
 */
void
ZLKJSolver::resetClosure(int r, int pl)
{
    auto do_reset = [&](int v) {
        winning[v] = -1;
        region[v] = BOT;
        strategy[v] = -1;
        open.set(v);
        reset.push_back(v);
    };

    reset.clear();
    for (int v : lost) {
        if (winning[v] == pl) do_reset(v); // (not attracted by the opponent after all)
    }
    lost.clear();

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
 * Attract the reset vertices back to the vertices won by <pl> in the subgame of region <r>,
 * into the new region <r2>. Attracted vertices are added to <A>. Returns the highest reset
 * vertex that is still open, or -1.
 */
int
ZLKJSolver::reattract(int r, int r2, int pl, std::vector<int> &A)
{
    const size_t start = A.size();

    auto attract = [&](int v, int str) {
        region[v] = r2;
        open.reset(v);
        winning[v] = pl;
        strategy[v] = str;
        A.push_back(v);
    };

    // an opponent vertex was reset because a successor was reset, so it is not attracted yet
    for (int v : reset) {
        if (owner(v) != pl or region[v] != BOT) continue;
        for (const int *e = outs(v); *e != -1; e++) {
            const int to = *e;
            if (winning[to] == pl and region[to] >= r) {
                attract(v, to);
                break;
            }
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
    if (A.size() == start) {
        for (int v : reset) top = std::max(top, v);
    } else {
        for (int v : reset) {
            if (is_open(region[v]) and v > top) top = v;
        }
    }
    return top;
}

/**
 * The frame of region <r> is done: give the heads a strategy.
 */
void
ZLKJSolver::finish(int r, int pl, std::vector<int> &H)
{
    for (int v : H) {
        if (winning[v] != pl or strategy[v] != HEAD) continue;
        if (owner(v) != pl) {
            strategy[v] = -1;
            continue;
        }
        strategy[v] = -1;
        for (const int *e = outs(v); *e != -1; e++) {
            const int to = *e;
            if (winning[to] == pl and region[to] >= r) {
                strategy[v] = to;
                break;
            }
        }
        if (strategy[v] == -1) LOGIC_ERROR;
    }
}

/**
 * Remove the entries of W[p] from <start> on of vertices not won by <p>, and duplicates.
 */
void
ZLKJSolver::compact(int p, size_t start)
{
    std::vector<int> &L = W[p];
    size_t k = start;
    for (size_t j = start; j < L.size(); j++) {
        const int v = L[j];
        if (winning[v] != p or seen[v]) continue;
        seen[v] = true;
        L[k++] = v;
    }
    L.resize(k);
    for (size_t j = start; j < k; j++) seen[L[j]] = false;
}

/**
 * Renumber the region ids, before they reach <region_limit>. <live> holds the ids of the frames
 * on the stack, from the bottom up, and gets their new ids; returns the next id.
 *
 * Every frame takes a new region id when it is pushed, and so does every round that attracts
 * reset vertices back (<r2> in reattract), so the ids count the steps of the whole run, not the
 * depth. After 2^31 steps they would overflow into the negative values of region[] (DIS, BOT
 * and the escape counters), and the new regions would look open or disabled.
 *
 * Only the ids of the frames on the stack (F.r in <frames>) are live: every comparison of a
 * region id is with 0 (ids versus the negative values), with DIS or BOT, or with the id r of
 * the top frame: rt >= 0 and rt < r in attractHeads, region < r and region >= r in
 * attractLosing, resetClosure, reattract and finish; the one exception is rt == r2 in
 * reattract, which takes its id just before and gives it to no frame. A frame lower on the
 * stack compares with its id again when it is the top frame. The other ids, of the frames that
 * are done and of the vertices attracted back, are dead: they stay in region[], but only to be
 * compared with the live ids. The ids on the stack increase from the bottom up, and a new frame
 * has the highest id. The justifications (strategy) and the vectors W hold vertices, not ids.
 *
 * So the k-th live id from the bottom (k = 1..L) becomes 2k, and every other id becomes the odd
 * number between the live ids around it: 2k+1 for an id between the k-th and the (k+1)-th, 1
 * below the first, 2L+1 above the last. Every comparison of an id with a live id (<, >=, ==)
 * then has the same answer as before; the negative values stay; and the next id is 2L+2, above
 * all ids. So the solver takes the same steps, with the same regions and strategies.
 * This runs at the top of the main loop of run(), before it reads the frame.
 */
int
ZLKJSolver::renumber(std::vector<int> &live)
{
#ifndef NDEBUG
    for (size_t k = 1; k < live.size(); k++) if (live[k-1] >= live[k]) LOGIC_ERROR;
#endif
    if (live.size() >= (size_t)INT_MAX/2) throw std::runtime_error("zlkj: out of region ids");

    auto map = [&](int x) {
        const long k = std::upper_bound(live.begin(), live.end(), x) - live.begin();
        return (int)((k > 0 and live[k-1] == x) ? 2*k : 2*k+1);
    };
    for (long v = 0; v < nodecount(); v++) {
        if (region[v] >= 0) region[v] = map(region[v]);
    }
    for (size_t k = 0; k < live.size(); k++) live[k] = (int)(2*k+2);

    renumbered++;
    return (int)(2*live.size()+2);
}

void
ZLKJSolver::run()
{
    const int n = nodecount();

    workspace.reset(new int[size_t(n)*2 + size_t(n+3)/4]);
    region = workspace.get();
    strategy = region + n;
    winning = reinterpret_cast<int8_t*>(strategy + n);

    memset(winning, -1, n);
    memset(strategy, -1, sizeof(int)*size_t(n));
    int top = -1;
    for (int v = 0; v < n; v++) {
        region[v] = disabled[v] ? DIS : BOT;
        if (!disabled[v]) top = v;
    }
    if (top == -1) LOGIC_ERROR; // unexpected empty game

    // all vertices that are not disabled may be open
    open = disabled;
    open.flip();

    seen.resize(n);

    /**
     * The stack of frames: the top vertex and region of each frame, its phase (0: attract the
     * heads, 1: a subgame was solved), where its entries and those of its current subgame start
     * in W, and how many of its entries are stale (reset or lost vertices, duplicates). The
     * heads of each frame are kept in a vector per depth, reused by the next frame at that
     * depth (only large ones are freed).
     */
    struct Frame {
        int i, r, phase;
        size_t start[2], sub[2], stale[2];
    };
    std::vector<Frame> frames;
    std::vector<std::vector<int>> heads;
    auto push = [&](int i, int r) {
        if (frames.size() == heads.size()) heads.emplace_back();
        frames.push_back(Frame{i, r, 0, {W[0].size(), W[1].size()}, {0, 0}, {0, 0}});
    };
    // finish the top frame: compact its entries when at least half are stale, otherwise
    // leave the stale entries to the parent frame
    auto pop = [&]() {
        Frame &F = frames.back();
        finish(F.r, priority(F.i)&1, heads[frames.size()-1]);
        size_t stale[2];
        for (int p = 0; p < 2; p++) {
            stale[p] = F.stale[p];
            if (stale[p] != 0 and 2*stale[p] >= W[p].size() - F.start[p]) {
                compact(p, F.start[p]);
                stale[p] = 0;
            }
        }
        std::vector<int> &H = heads[frames.size()-1];
        if (H.capacity() > 65536) std::vector<int>().swap(H);
        else H.clear();
        frames.pop_back();
        if (!frames.empty()) {
            frames.back().stale[0] += stale[0];
            frames.back().stale[1] += stale[1];
        }
    };

    int next_r = 0;
    const int limit = region_limit; // see renumber
    push(top, next_r++);

    while (!frames.empty()) {
        // An iteration takes at most two region ids. Renumber them first when they run out: the
        // ids of the frames go to renumber in a vector of their own, so that <frames> stays
        // local to run().
        if (next_r >= limit - 1) {
            std::vector<int> live;
            for (const Frame &G : frames) live.push_back(G.r);
            next_r = renumber(live);
            for (size_t k = 0; k < frames.size(); k++) frames[k].r = live[k];
        }

        Frame &F = frames.back();
        const int r = F.r;
        const int pl = priority(F.i)&1;
        std::vector<int> &H = heads[frames.size()-1];

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "\n\033[1;37mDepth %zu phase %d\033[m: node %d priority %d\n", frames.size(), F.phase, F.i, priority(F.i));
#endif

        if (F.phase == 0) {
            const int j = attractHeads(F.i, r, W[pl], H);
            ++iterations;
            if (j == -1) {
                // the region is the whole subgame
                pop();
            } else {
                F.phase = 1;
                F.sub[0] = W[0].size();
                F.sub[1] = W[1].size();
                push(j, next_r++);
            }
            continue;
        }

        /**
         * A subgame was solved: its vertices won by each player are in W from F.sub on.
         */
        ++rounds;
        const int count = W[1-pl].size() == F.sub[1-pl] ? 0 : attractLosing(r, pl, W[1-pl], F.sub[1-pl], H);
        if (count == 0) {
            // the opponent attracts nothing from <pl>
            lost.clear();
            pop();
            continue;
        }

        resetClosure(r, pl);
        // the entries of the lost and reset vertices are stale (or become duplicates)
        F.stale[pl] += count + reset.size();
        const int new_i = reset.empty() ? -1 : reattract(r, next_r++, pl, W[pl]);

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "%zu vertices reset, new top %d\n", reset.size(), new_i);
#endif

        if (new_i == -1) {
            // every reset vertex was attracted back
            pop();
        } else {
            // solve the reset vertices that are still open (and stay in phase 1)
            F.sub[0] = W[0].size();
            F.sub[1] = W[1].size();
            push(new_i, next_r++);
        }
    }

    for (int v = 0; v < n; v++) {
        if (disabled[v]) continue;
#ifndef NDEBUG
        if (winning[v] == -1) LOGIC_ERROR;
        if (owner(v) == winning[v] and strategy[v] < 0) LOGIC_ERROR;
#endif
        Solver::solve(v, winning[v], owner(v) == winning[v] ? strategy[v] : -1);
    }

    W[0].clear();
    W[1].clear();
    workspace.reset();

    if (renumbered) logger << "renumbered the region ids " << renumbered << " times." << std::endl;
    logger << "solved with " << iterations << " iterations, " << rounds << " rounds, " << resets << " resets and " << reattracted << " attracted back." << std::endl;
}

}
