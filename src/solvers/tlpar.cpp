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

#include <climits>

#include "tlpar.hpp"


#define PARTIALLY_CLOSED 0 // find tangles in partially locally closed regions, but costs time on practical games
                           // (its search uses <Q> as a stack, which needs a change for the kept regions, see tl)

namespace pg {

TLParSolver::TLParSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}


TLParSolver::~TLParSolver()
{
}


/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
TLParSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        const int from = *curedge;
        // (only vertices in <R> but not in <Z>; the heads of <pl> get their strategy after the
        // attraction, see tl)
        if (C[from]) {
            if (owner(from) != pl) {
                // check if opponent can escape
                unsigned int e = escs[from];
                if (e == 0) {
                    if (count_in_R) {
                        // the higher regions are kept from the last pass, so their vertices
                        // were not attracted in this one: count the successors that are not
                        // attracted yet or in the current region (in <R>, see drain), as
                        // those processed before <v> would have counted already
                        e = 0;
                        for (auto curedge = outs(from); *curedge != -1; curedge++) {
                            if (R[*curedge]) e++;
                        }
                    } else {
                        // first time in this attraction: start from the number of successors in G
                        e = gdeg[from];
                        if (e == 0) {
                            for (auto curedge = outs(from); *curedge != -1; curedge++) {
                                if (G[*curedge]) e++;
                            }
                            gdeg[from] = e;
                        }
                    }
                }
                escs[from] = --e;
                if (e > 0) continue; // escapes
            }
            // attract
            Z[from] = true;
            C[from] = false;
            str[from] = owner(from) == pl ? v : -1;
            Q.push(from);
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << "\033[m by \033[1;36m" << pl << "\033[m";
                if (owner(from) == pl) logger << " (via " << label_vertex(v) << ")" << std::endl;
                else logger << " (forced)" << std::endl;
            }
#endif
        }
    }
}


/**
 * Try to attract tangle <t> for player <pl> to attractor set <Z>.
 * All vertices in tangle <t> must be in <Z+R> and the opponent may not escape to subgame <R\Z>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
bool
TLParSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, int maxpr)
{
    /**
     * Check if tangle is won by player <pl> and not deleted.
     */
    {
        const int tangle_pr = tpr[t];
        if (tangle_pr == -1) return false; // deleted tangle
        if (tlazy[t] == 2) return false; // a broken region tangle (analysed later)
        if (tangle_pr > maxpr) return false; // too high...
        if (pl != (tangle_pr&1)) return false; // not of desired parity
    }

    /**
     * Check if the tangle can escape to R\Z.
     */
    {
        int e = tescs[t];
        if (e == 0) {
            int v, *ptr = &toutdata[toutpos[t]];
            if (count_in_R) {
                // The higher regions are kept from the last pass (see tl), so their vertices
                // were not attracted in this pass: count the escapes not attracted yet or in
                // the current region (in <R>, see drain). As in a full pass, where the
                // counter starts at the first escape attracted into a region the tangle can
                // be attracted into, an escape in another region (of the other player or
                // lower than the tangle) blocks it.
                const int tpr_t = tpr[t];
                while ((v=*ptr++) != -1) {
                    if (!G[v]) continue;
                    if (R[v]) {
                        e++;
                    } else {
                        const Region &r = regions[regidx[v]];
                        if ((r.pr&1) != pl or r.pr < tpr_t) { e = INT_MAX/2; break; }
                    }
                }
            } else {
                while ((v=*ptr++) != -1) {
                    if (G[v]) e++;
                }
            }
        }
        tescs[t] = --e;
        if (e > 0) return false;
    }

    /**
     * On-the-fly detect out-of-game tangles
     */
    {
        int *ptr = &tvdata[tvpos[t]];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            ptr++; // skip strategy
            if (!G[v]) {
                // on-the-fly detect out-of-game tangles
                if (tlazy[t]) { tlazy[t] = 2; return false; } // (a region tangle: analysed later)
                tpr[t] = -1; // delete the tangle
                return false; // is now a deleted tangle
            }
        }
    }
    if (tlazy[t]) lreg[t] = cur_ridx;

    if (maxpr == INT_MAX) {
        // attracted to a dominion, so delete the tangle
        tpr[t] = -1;
    }

    /**
     * Attract!
     */
    {
        int *ptr = &tvdata[tvpos[t]];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            const int s = *ptr++;
            if (!C[v]) continue; // already in <Z>, or not in <R> (might be higher region)
            Z[v] = true;
            C[v] = false;
            str[v] = s;
            Q.push(v);

#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[36m" << label_vertex(v) << "\033[m by \033[1;36m" << pl << "\033[m";
                logger << " (via tangle " << t << ")" << std::endl;
            }
#endif
        }
    }

    return true;
}


/**
 * Attract vertices that are in tangles.
 * Current attracting set is <Z>.
 * Current subgame is <R>.
 * Only tangles that are contained in <R>+<Z>
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
__attribute__((always_inline)) void
TLParSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, int maxpr)
{
    if (!has_tin[v]) return;
    const int *ts = &tin_pool[tin_off[v]];
    for (unsigned int i = 0, n = tin_size[v]; i < n; i++) attractTangle(ts[i], pl, R, Z, maxpr);
}


/**
 * Before attracting to <Z> from <R>: the candidates for attraction are the vertices in <R>
 * but not in <Z> (in <C>). Every vertex added to <Z> is then removed from <C>. (Testing
 * one bitset instead of <Z> and then <R> makes the attraction about twice as fast on games
 * with many edges.)
 */
void
TLParSolver::setCandidates(const bitset &R, const bitset &Z)
{
    C = R;
    C -= Z;
}


/**
 * Atomic access to bitsets that the Lace workers share in a parallel round.
 */
static inline bool
atomic_test(const bitset &B, const int v)
{
    return (__atomic_load_n(B.data() + (v>>6), __ATOMIC_RELAXED) >> (v&63)) & 1;
}

/**
 * Atomically add <v> to <B>; true if this call added it.
 */
static inline bool
atomic_claim(bitset &B, const int v)
{
    const uint64_t mask = uint64_t(1) << (v&63);
    return (__atomic_fetch_or(B.data() + (v>>6), mask, __ATOMIC_ACQ_REL) & mask) == 0;
}

/**
 * Atomically remove <v> from <B> (as other workers may change other bits of its block).
 */
static inline void
atomic_remove(bitset &B, const int v)
{
    __atomic_fetch_and(B.data() + (v>>6), ~(uint64_t(1) << (v&63)), __ATOMIC_RELAXED);
}


/**
 * Parallel version of attractVertices, for the round in <par_pl>, <par_R>, <par_Z>:
 * attract via <v>, adding attracted vertices to <out>. <par_R> is not changed in a round.
 */
inline void
TLParSolver::attractVerticesPar(const int v, std::vector<int> &out)
{
    const int pl = par_pl;
    bitset &R = *par_R;
    bitset &Z = *par_Z;
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        const int from = *curedge;
        // (see attractVertices; a vertex that another worker attracts in this round may still
        // look like a candidate here, but then the claim below fails)
        if (atomic_test(C, from)) {
            if (owner(from) != pl) {
                // count down the escapes; the first visit in this attraction (0) initializes
                unsigned int e = __atomic_load_n(&escs[from], __ATOMIC_RELAXED);
                unsigned int next;
                do {
                    if (e == 0) {
                        unsigned int g = 0;
                        if (count_in_R) {
                            // (see attractVertices)
                            for (auto e2 = outs(from); *e2 != -1; e2++) {
                                if (R[*e2]) g++;
                            }
                        } else {
                            g = __atomic_load_n(&gdeg[from], __ATOMIC_RELAXED);
                            if (g == 0) {
                                for (auto e2 = outs(from); *e2 != -1; e2++) {
                                    if (G[*e2]) g++;
                                }
                                __atomic_store_n(&gdeg[from], g, __ATOMIC_RELAXED);
                            }
                        }
                        next = g-1;
                    } else {
                        next = e-1;
                    }
                } while (!__atomic_compare_exchange_n(&escs[from], &e, next, true, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED));
                // (a vertex already attracted by another worker may count below zero; then
                // <next> wraps around, and the claim below would fail anyway)
                if (next != 0) continue; // escapes
            }
            if (!atomic_claim(Z, from)) continue; // another worker attracted it
            atomic_remove(C, from);
            str[from] = owner(from) == pl ? v : -1;
            out.push_back(from);
        }
    }
}


/**
 * Parallel version of attractTangle (see there), adding attracted vertices to <out>.
 */
void
TLParSolver::attractTanglePar(const int t, std::vector<int> &out)
{
    const int pl = par_pl;
    const int maxpr = par_maxpr;
    bitset &R = *par_R;
    bitset &Z = *par_Z;

    const int tangle_pr = __atomic_load_n(&tpr[t], __ATOMIC_RELAXED);
    if (tangle_pr == -1) return; // deleted tangle
    const char lazy = __atomic_load_n(&tlazy[t], __ATOMIC_RELAXED);
    if (lazy == 2) return; // a broken region tangle (analysed later)
    if (tangle_pr > maxpr) return; // too high...
    if (pl != (tangle_pr&1)) return; // not of desired parity

    {
        unsigned int *ptr = &tescs[t];
        unsigned int e = __atomic_load_n(ptr, __ATOMIC_RELAXED);
        unsigned int next;
        do {
            if (e == 0) {
                unsigned int c = 0;
                int w, *p = &toutdata[toutpos[t]];
                if (count_in_R) {
                    // (see attractTangle)
                    while ((w=*p++) != -1) {
                        if (!G[w]) continue;
                        if (R[w]) {
                            c++;
                        } else {
                            const Region &r = regions[regidx[w]];
                            if ((r.pr&1) != pl or r.pr < tangle_pr) { c = INT_MAX/2; break; }
                        }
                    }
                } else {
                    while ((w=*p++) != -1) {
                        if (G[w]) c++;
                    }
                }
                next = c-1;
            } else {
                next = e-1;
            }
        } while (!__atomic_compare_exchange_n(ptr, &e, next, true, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED));
        if (next != 0) return;
    }

    for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
        if (!G[*ptr]) {
            // on-the-fly delete out-of-game tangle (a region tangle is analysed later)
            if (lazy) __atomic_store_n(&tlazy[t], (char)2, __ATOMIC_RELAXED);
            else __atomic_store_n(&tpr[t], -1, __ATOMIC_RELAXED);
            return;
        }
    }
    if (lazy) __atomic_store_n(&lreg[t], cur_ridx, __ATOMIC_RELAXED);

    if (maxpr == INT_MAX) __atomic_store_n(&tpr[t], -1, __ATOMIC_RELAXED); // attracted to a dominion

    for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
        const int v = *ptr;
        if (!atomic_test(C, v)) continue; // in <Z>, or not in <R>
        if (!atomic_claim(Z, v)) continue;
        atomic_remove(C, v);
        str[v] = *(ptr+1);
        out.push_back(v);
    }
}


/**
 * Process vertices [start, start+n) of the current parallel round (in <Q>), and append the
 * vertices attracted by them to <Q>, at the end of the round so far (<q_tail>).
 */
void
TLParSolver::roundChunk(int start, int n)
{
    auto &out = pbuf[lace_worker_id()].v;
    out.clear();
    for (int i = start; i < start+n; i++) {
        const int v = Q[i];
        if (par_ridx >= 0) regidx[v] = par_ridx;
        attractVerticesPar(v, out);
        if (has_tin[v]) {
            const int *ts = &tin_pool[tin_off[v]];
            for (unsigned int j = 0, m = tin_size[v]; j < m; j++) attractTanglePar(ts[j], out);
        }
    }
    if (!out.empty()) {
        const unsigned int off = __atomic_fetch_add(&q_tail, (unsigned int)out.size(), __ATOMIC_RELAXED);
        std::copy(out.begin(), out.end(), &Q[off]);
    }
}

VOID_TASK_3(tlpar_round, TLParSolver*, s, int, start, int, n)
void tlpar_round_CALL(lace_worker* lw, TLParSolver* s, int start, int n)
{
    if (n > s->par_chunk) {
        const int half = n/2;
        tlpar_round_SPAWN(lw, s, start+half, n-half);
        tlpar_round_CALL(lw, s, start, half);
        tlpar_round_SYNC(lw);
    } else {
        s->roundChunk(start, n);
    }
}


/**
 * Attract as player <pl> to <Z> from subgame <R>, processing the vertices of <Q> from <qh>
 * on until all are processed: one vertex at a time (in FIFO order) while few are left,
 * otherwise in parallel rounds over all that are left (the frontier), appending the newly
 * attracted vertices to <Q>. If <ridx> >= 0, the region of the processed vertices is <ridx>
 * (see the end of tl).
 *
 * Unlike in TL, the vertices of <Z> are not removed from <R> one at a time (when processed),
 * but by the caller afterwards, so that a parallel round does not change <R>, and counting
 * the escapes in <R> (see attractVertices) is the same in a round as one vertex at a time.
 */
void
TLParSolver::drain(const int pl, bitset &R, bitset &Z, const int maxpr, unsigned int &qh, const int ridx)
{
    while (qh < Q.size()) {
        if (lw == nullptr or (int)(Q.size()-qh) < par_min) {
            const int v = Q[qh++];
            if (ridx >= 0) regidx[v] = ridx;
            attractVertices(pl, v, R, Z);
            attractTangles(pl, v, R, Z, maxpr);
        } else {
            const unsigned int end = Q.size();
            par_pl = pl;
            par_maxpr = maxpr;
            par_ridx = ridx;
            par_R = &R;
            par_Z = &Z;
            q_tail = end;
            tlpar_round_CALL(lw, this, qh, end-qh);
            Q.grow(q_tail-end); // (the vertices attracted in the round, see roundChunk)
            par_rounds++;
            par_vertices += end-qh;
            qh = end;
        }
    }
}


/**
 * Store the tangle in <tangle> (with the strategy in <strat>) of priority <pr>, with the
 * escapes in <tangleto>, and link it into <tin>.
 */
void
TLParSolver::addTangle(const int pr, const int *strat)
{
    // new tangle idx
    const int tidx = tpr.size();

    // add back links to all normal vertices in our [out]
    for (unsigned int x = 0; x < tangleto.size(); x++) {
        const int e = tangleto[x];
        if (!has_tin[e]) {
            has_tin[e] = true;
            tin_size[e] = 0;
            tin_cap[e] = 0;
        }
        if (tin_size[e] == tin_cap[e]) {
            // move the list of <e> to a new block of twice the size at the end of the pool
            const unsigned int cap = tin_cap[e] ? 2*tin_cap[e] : 2;
            const size_t off = tin_pool.size();
            tin_pool.resize(off+cap);
            if (tin_size[e]) std::copy(&tin_pool[tin_off[e]], &tin_pool[tin_off[e]]+tin_size[e], &tin_pool[off]);
            tin_off[e] = off;
            tin_cap[e] = cap;
        }
        tin_pool[tin_off[e]+tin_size[e]++] = tidx;
    }

    // store the tangle: its vertex-strategy pairs and its escapes (each ending with -1)
    tvpos.push_back(tvdata.size());
    for (const int v : tangle) {
        tvdata.push_back(v);
        tvdata.push_back(strat[v]);
    }
    tvdata.push_back(-1);
    toutpos.push_back(toutdata.size());
    toutdata.insert(toutdata.end(), &tangleto[0], &tangleto[0]+tangleto.size());
    toutdata.push_back(-1);

    // and set p to pr
    tpr.push_back(pr);
    tescs.push_back(0);
    tlazy.push_back(0);
    lreg.push_back(-1);
}


/**
 * Compute SCCs in subgraph induced by <R> and <str>.
 * Start the SCC computation at vertex <startvertex>.
 * Every SCC is then processed as a tangle.
 * If the tangle is closed, it is a dominion and added to <S> and <Q>.
 */
bool
TLParSolver::extractTangles(int startvertex, bitset &R, int *str, int tangle_pr)
{
    bool new_tangles = false;
    const int pr = tangle_pr >= 0 ? tangle_pr : priority(startvertex);
    const int pl = pr&1;

    /**
     * The following is the nonrecursive implementation of David Pearce,
     * "A space-efficient algorithm for finding strongly connected components" (IPL, 2016),
     * modified to on-the-fly restrict the graph by <R> and <str>.
     */

    // beginVisiting
    pea_state.push(startvertex);
    pea_state.push(0);
    pea_root[startvertex] = true;
    pea_vidx[startvertex] = pea_curidx++;
    while (pea_state.nonempty()) {
pearce_again:
        // visitLoop
        const unsigned int n = pea_state.back2();
        unsigned int i = pea_state.back();

        if (owner(n) != pl) {
            auto edges = outs(n);
            if (i>0) {
                // finishEdge
                const int w = edges[i-1];
                if (pea_vidx[w] < pea_vidx[n]) {
                    pea_vidx[n] = pea_vidx[w];
                    pea_root[n] = false;
                }
            }
            for (;;) {
                const int to = edges[i];
                if (to == -1) break; // done
                // beginEdge
                if (R[to]) {
                    if (pea_vidx[to] == 0) {
                        pea_state.back() = i+1;
                        // beginVisiting
                        pea_state.push(to);
                        pea_state.push(0);
                        pea_root[to] = true;
                        pea_vidx[to] = pea_curidx++;
                        goto pearce_again; // break; continue;
                    } else {
                        // finishEdge
                        if (pea_vidx[to] < pea_vidx[n]) {
                            pea_vidx[n] = pea_vidx[to];
                            pea_root[n] = false;
                        }
                    }
                }
                i++;
            }
        } else {
            const int s = str[n];
            if (i == 0) {
                // beginEdge
                if (pea_vidx[s] == 0) {
                    pea_state.back() = 1;
                    // beginVisiting
                    pea_state.push(s);
                    pea_state.push(0);
                    pea_root[s] = true;
                    pea_vidx[s] = pea_curidx++;
                    goto pearce_again; // break; continue;
                }
            }
            // finishEdge
            if (pea_vidx[s] < pea_vidx[n]) {
                pea_vidx[n] = pea_vidx[s];
                pea_root[n] = false;
            }
        }
        // finishVisiting
        pea_state.pop2();
        if (pea_root[n]) {
            pea_curidx -= 1;
            tangle.push_back(n);
            while (pea_S.nonempty()) {
                const int t = pea_S.back();
                if (pea_vidx[n]>pea_vidx[t]) break;
                pea_S.pop();
                pea_curidx -= 1;
                pea_vidx[t] = (unsigned int)-1;
                tangle.push_back(t);
            }
            pea_vidx[n] = (unsigned int)-1;
        } else {
            pea_S.push(n);
            continue;
        }

        /**
         * End of Pearce's algorithm.
         * At this point, we have an SCC in <tangle>.
         *
         * Now check if the SCC is nontrivial, i.e., contains a cycle, i.e., either 2+ vertices
         * or a self-loop.
         */

        const bool is_not_tangle = (tangle.size() == 1) and
            ((unsigned int)str[n] != n) and (str[n] != -1 or !game.has_edge(n, n));
        if (is_not_tangle) {
            // std::cout << "not a tangle: " << pea_vidx[n] << std::endl;
            tangle.clear();
            continue;
        }

        /**
         * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
         */

        for (const int v : tangle) escapes[v] = true;

        for (const int v : tangle) {
            if (owner(v) != pl) {
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (G[to] and !S[to] && !escapes[to]) {
                        escapes[to] = true;
                        tangleto.push(to);
                    }
                }
            }
        }

        if ((tangle.size()+tangleto.size()) < (escapes.size()/64)) {
            for (unsigned int x = 0; x < tangleto.size(); x++) escapes[tangleto[x]] = false;
            for (const int v : tangle) escapes[v] = false;
        } else {
            escapes.reset(); // this is very expensive on large practical games with many small tangles
        }

        /**
         * If there are no outgoing edges, then we have found a dominion.
         */

        if (tangleto.empty()) {
            // dominion
#ifndef NDEBUG
            if (trace) {
                logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m";
                if (trace >= 2) {
                    for (const int v : tangle) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
                logger << std::endl;
            }
#endif
            for (const int v : tangle) S[v] = true;
            dominions++;
            tangle.clear();
            new_tangles = true;
            continue;
        }

        /**
         * We're not a dominion, we're a tangle.
         */

#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;38;5;198mnew tangle " << pr << "\033[m (" << tpr.size() << ")";
            if (trace >= 2) {
                for (const int v : tangle) {
                    logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                }
            }
            logger << " with " << tangleto.size() << " escape vertices.";
            logger << std::endl;
        }
#endif

        addTangle(pr, str);
        tangles++;
        tangle.clear();
        tangleto.clear();
        new_tangles = true;
    }

    pea_S.clear();
    return new_tangles;
}


/**
 * Learn the closed region <ridx> (of player <pl> and priority <pr>, in <Z>) as one region
 * tangle, if all its escapes are in the highest region of <pl> (in <H>): then none of its
 * tangles can be attracted into a higher region than the region as a whole, so it is only
 * decomposed into its tangles once it is broken (see analyseLazies). This avoids learning
 * many small tangles that all go into the highest region. Returns whether it was learned.
 */
bool
TLParSolver::learnRegion(const int ridx, const int pl, const int pr)
{
    bool single = true;
    for (size_t i = regions[ridx].begin; single and i < Q.size(); i++) {
        const int u = Q[i];
        if (owner(u) == pl) continue;
        for (auto curedge = outs(u); *curedge != -1; curedge++) {
            const int to = *curedge;
            if (!G[to] or Z[to] or escapes[to]) continue;
            if (!H[pl][to]) { single = false; break; }
            escapes[to] = true;
            tangleto.push(to);
        }
    }
    for (unsigned int x = 0; x < tangleto.size(); x++) escapes[tangleto[x]] = false;
    if (!single or tangleto.empty()) { // (without escapes, extracting finds the dominion)
        tangleto.clear();
        return false;
    }
    for (size_t i = regions[ridx].begin; i < Q.size(); i++) tangle.push_back(Q[i]);
    addTangle(pr, str);
    tlazy.back() = 1;
    lreg.back() = ridx;
    lazies.push_back((int)tpr.size()-1);
    tangles++;
    tangle.clear();
    tangleto.clear();
    return true;
}


/**
 * Decompose region tangles into their tangles (see learnRegion): all of them (after a
 * dominion), or (if <broken>) those of earlier passes that this pass did not attract again
 * as a whole, into the highest region of their player, or that lost vertices. Their tangles
 * are learned from their vertices that are still in <G>, with their strategy, and dominions
 * found are solved. Returns whether dominions were found.
 */
bool
TLParSolver::analyseLazies(const bool broken)
{
    bool found = false;
    size_t keepn = 0;
    for (size_t li = 0; li < lazies.size(); li++) {
        const int t = lazies[li];
        if (tpr[t] == -1 or tlazy[t] == 0) continue;
        if (broken and (size_t)t >= pass_tangles) { lazies[keepn++] = t; continue; } // (new)
        const int pr = tpr[t], pl = pr&1;
        // (still delayed if it is (as a whole) in the highest region of its player)
        if (broken and tlazy[t] == 1 and lreg[t] >= 0 and (size_t)lreg[t] < regions.size() and
            (pl == 0 ? regions[lreg[t]].highest0 : regions[lreg[t]].highest1)) { lazies[keepn++] = t; continue; }
        // its vertices in <G> and their strategy, without the vertices of <pl> whose strategy
        // leaves them
        for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
            if (G[*ptr]) { LB[*ptr] = true; astr[*ptr] = *(ptr+1); pea_vidx[*ptr] = 0; }
        }
        for (bool changed = true; changed; ) {
            changed = false;
            for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
                const int v = *ptr;
                if (LB[v] and owner(v) == pl and (astr[v] < 0 or !LB[astr[v]])) { LB[v] = false; changed = true; }
            }
        }
        const auto D = dominions;
        pea_curidx = 1;
        for (auto v = LB.find_last(); v != bitset::npos; v = LB.find_prev(v)) {
            if (pea_vidx[v] == 0) extractTangles((int)v, LB, astr, pr);
        }
        for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) LB[*ptr] = false;
        tpr[t] = -1;
        tlazy[t] = 0;
        if (D != dominions) {
            // solve the dominions (in <S>) with their attractor in <G>
            for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) str[v] = astr[v];
            count_in_R = false;
            std::fill(escs, escs+nodecount(), '\0');
            std::fill(tescs.begin(), tescs.end(), '\0');
            Q.clear();
            for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) Q.push(v);
            setCandidates(G, S);
            unsigned int qh = 0;
            drain(pl, G, S, INT_MAX, qh, -1);
            Q.clear();
            for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) Solver::solve(v, pl, str[v]);
            G -= S;
            S.reset();
            std::fill(gdeg, gdeg+nodecount(), 0); // (counted again when needed)
            escs_is_gdeg = false;
            found = true;
        }
    }
    lazies.resize(keepn);
    return found;
}


/**
 * The region of the current pass that tangle <t> would be attracted into in the next pass:
 * the lowest region of its escapes, as its last escape is attracted there, if all its
 * escapes are in regions of its player that are not lower than the tangle. Returns -1 if
 * there is none.
 */
int
TLParSolver::tangleTarget(const int t)
{
    if (tpr[t] == -1) return -1;
    const int tpl = tpr[t]&1;
    int k = -1;
    for (int *ptr = &toutdata[toutpos[t]]; *ptr != -1; ptr++) {
        if (!G[*ptr]) continue;
        const int r = regidx[*ptr];
        if ((regions[r].pr&1) != tpl or regions[r].pr < tpr[t]) return -1;
        if (r > k) k = r;
    }
    return k;
}


bool
TLParSolver::tl(void)
{
    bool new_tangles = false;

    if (escs_is_gdeg) {
        // The last dominion attraction left in <escs>, for each vertex still in <G> that has an
        // edge into the dominion, its number of successors that are not in the dominion.
        for (int v=0; v<nodecount(); v++) {
            if (escs[v] != 0 and G[v]) gdeg[v] = escs[v];
            escs[v] = 0;
        }
        escs_is_gdeg = false;
    } else {
        std::fill(escs, escs+nodecount(), '\0');
    }
    std::fill(tescs.begin(), tescs.end(), '\0');

    bool highest0 = true;
    bool highest1 = true;
    bool first_region = true; // the first region is attracted in the full unsolved game <G>

    pass_tangles = tpr.size();
    int extend = -1; // the region of the last pass that this pass extends (see below)
    const bool kept = keep >= 0;
    // (<H> keeps the highest region of each player if that region is kept, see below)
    for (int p = 0; p < 2; p++) {
        bool k = false;
        for (int i = 0; kept and i < keep; i++) if ((regions[i].pr&1) == p) { k = true; break; }
        if (!k) H[p].reset();
    }
    const int xreg = keep; // (the region extended by this pass, or -1)
    if (kept) {
        // The regions above region <keep> of the last pass did not change: none of its new
        // tangles can be attracted into them, so this pass would compute them again exactly
        // (and find no new tangles in them). Keep them. Region <keep> itself is computed in
        // the same subgame again, with more tangles, so it only grows: start from its
        // vertices (and heads), attract the new tangles that go into it, and attract from
        // there. Then continue with the lower regions, in the subgame of their vertices.
        const Region &k = regions[keep];
        highest0 = k.highest0;
        highest1 = k.highest1;
        first_region = keep == 0;
        const bool last = keep+1 == (int)regions.size();
        const size_t kend = last ? Q.size() : regions[keep+1].begin;
        const size_t hend = last ? heads.size() : regions[keep+1].hbegin;
        // (the vertices and heads of the region are in <Zx> and <Vx> if the last pass ended it
        // there, see the end of the region below)
        if (zx != keep) {
            Zx.reset();
            Vx.reset();
            for (size_t i = k.begin; i < kend; i++) Zx[Q[i]] = true;
            for (size_t i = k.hbegin; i < hend; i++) Vx[heads[i]] = true;
        }
        Z.swap(Zx);
        V.swap(Vx);
        zx = -1;
        if (keep == 0) {
            R = G;
            R -= Z;
        } else {
            R.reset();
            for (size_t i = kend; i < Q.size(); i++) R[Q[i]] = true;
        }
        while (Q.size() > kend) Q.pop();
        heads.resize(hend);
        regions.resize(keep+1);
        count_in_R = true;
        extend = keep;
        // (region tangles in the recomputed regions are attracted again, or analysed at the end)
        for (const int t : lazies) if (tlazy[t] == 1 and lreg[t] > keep) lreg[t] = -1;
    } else {
        // After a dominion: analyse the region tangles (see the end of tl).
        if (analyseLazies(false)) { keep = -1; return true; }
        R = G;
        Q.clear();
        heads.clear();
        regions.clear();
        count_in_R = false;
        zx = -1;
    }
    setCandidates(R, Z);
    // During a pass, <Q> holds the vertices of its regions, region by region, in the order of
    // attraction; the first <qn> have been processed (their predecessors attracted).
    unsigned int qn = Q.size();
    auto top = R.find_last();
    while (top != bitset::npos or extend >= 0) {
        int pr, pl, ridx;
        unsigned int vcnt = 0;

        if (extend >= 0) {
            ridx = extend;
            cur_ridx = ridx;
            extend = -1;
            pr = regions[ridx].pr;
            pl = pr&1;
            // attract the new tangles that go into this region (see the end of tl)
            for (const int t : seeds) {
                bool ingame = true;
                for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
                    if (!G[*ptr]) { ingame = false; break; }
                }
                if (!ingame) {
                    if (tlazy[t]) tlazy[t] = 2; // (a region tangle: analysed later)
                    else tpr[t] = -1;
                    continue;
                }
                if (tlazy[t]) lreg[t] = ridx;
                for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
                    const int u = *ptr;
                    if (!C[u]) continue; // (in <Z>, or not in <R>)
                    Z[u] = true;
                    C[u] = false;
                    str[u] = *(ptr+1);
                    Q.push(u);
                }
            }
            drain(pl, R, Z, pr, qn, kept ? ridx : -1); // (see the end of tl)
            R -= Z; // (see drain)
            // (more heads, of lower priorities, if no vertex of the other parity is left)
            top = R.find_last();
        } else {
            pr = priority(top);
            pl = priority(top)&1;
            ridx = (int)regions.size();
            cur_ridx = ridx;
            regions.push_back({Q.size(), heads.size(), pr, highest0, highest1});
        }

        // attract from all heads with priority <pr> that are in <R>
        while (top != bitset::npos and (priority(top)&1) == (pr&1)) { // otf compress
            // first add all heads of the same parity that are in <R> to <Z>, then attract
            // once: the region is the same, but the vertices of these priorities that the
            // earlier heads would attract are now heads, so the opponent's need no escape
            // counter (attracting from each head separately did several times more counter
            // updates, even on games where every vertex is a head)
            for (; top != bitset::npos; top = R.find_prev(top)) {
                if ((priority(top)&1) != (pr&1)) break;
                // if (priority(top) != pr) break; // disable otf compress

                vcnt++;
                V[top] = true; // heads
                heads.push_back(top);
                Z[top] = true; // add to <Z>
                C[top] = false;
                str[top] = -1;
                Q.push(top);

#ifndef NDEBUG
                // maybe report event
                if (trace >= 3) {
                    logger << "\033[1;37mattracting towards \033[36m" << label_vertex(top) << "\033[m for \033[1;36m" << pl << "\033[m" << std::endl;
                }
#endif
            }

            drain(pl, R, Z, pr, qn, kept ? ridx : -1); // (see the end of tl)
            R -= Z; // (see drain)

            // (more heads, of lower priorities, if no vertex of the other parity is left; the
            // vertices above <top> are no longer in <R>, so search down from <top>)
            if (top != bitset::npos and !R[top]) top = R.find_prev(top);
        }

        // the heads of <pl> without a strategy play to <Z> if they can (doing this during the
        // attraction costs a test of each predecessor in <Z> of each attracted vertex)
        for (size_t i = regions[ridx].hbegin; i < heads.size(); i++) {
            const int h = heads[i];
            if (owner(h) != pl or str[h] != -1) continue;
            for (auto curedge = outs(h); *curedge != -1; curedge++) {
                if (Z[*curedge]) { str[h] = *curedge; break; }
            }
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m ";
            logger << "\033[1;36m" << pr << "\033[m";
            for (auto v = Z.find_last(); v != bitset::npos; v = Z.find_prev(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif

#if PARTIALLY_CLOSED
        bool leaks = false;
        W.reset();

        if (top != bitset::npos) { // if not lowest region
            // figure out what are good heads.
            for (auto v = V.find_last(); v != bitset::npos; v = V.find_prev(v)) {
                // check if open
                if (owner(v) == pl) {
                    if (str[v] == -1) {
                        vcnt--;
                        W[v] = true;
                        Z[v] = false;
                        leaks = true;
                    }
                } else {
                    for (auto curedge = outs(v); *curedge != -1; curedge++) {
                        if (R[*curedge]) {
                            vcnt--;
                            W[v] = true;
                            Z[v] = false;
                            leaks = true;
                            break;
                        }
                    }
                }
            }
        }
#else
        bool leaks = false;

        if (top != bitset::npos) { // if not lowest region
            // figure out what are good heads.
            for (auto v = V.find_last(); v != bitset::npos; v = V.find_prev(v)) {
                // check if open
                if (owner(v) == pl) {
                    if (str[v] == -1) {
                        leaks = true;
                        break;
                    }
                } else {
                    for (auto curedge = outs(v); *curedge != -1; curedge++) {
                        if (R[*curedge]) {
                            leaks = true;
                            break;
                        }
                    }
                    if (leaks) break;
                }
            }
        }
#endif

#if PARTIALLY_CLOSED
        // if any heads in V left, we need to refine
        if (vcnt > 0) {
            if (leaks) {
                // then do a backwards search removing all vertices that cannot stay
                for (auto v = W.find_first(); v != bitset::npos; v = W.find_next(v)) {
                    W[v] = false; // reset W
                    V[v] = false; // remove from V
                    Q.push(v);
                }
                while (Q.nonempty()) {
                    int v = Q.pop();
#ifndef NDEBUG
                    if (trace >= 2) logger << "leaking: \033[1;36m" << label_vertex(v) << "\033[m" << std::endl;
#endif
                    for (auto curedge = ins(v); *curedge != -1; curedge++) {
                        auto from = *curedge;
                        if (Z[from] && (str[from] == -1 or str[from] == v)) {
                            Z[from] = false;
                            V[from] = false;
                            Q.push(from);
                        }
                    }
                }
            }
            // V &= Z;
#else
        // if any leaks, just skip
        if (!leaks) {
#endif
            // now Z is the vertices remaining

            if (Z.any()) {
                /**
                 * Extract tangles by computing bottom SCCs starting at each top vertex.
                 * Note: each bottom SCC must contain a top vertex.
                 */

                bool highest = pl == 0 ? highest0 : highest1;
                if (highest) {
                    // dominion
#ifndef NDEBUG
                    if (trace) {
                        logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m";
                        if (trace >= 2) {
                            for (const int v : tangle) {
                                logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                            }
                        }
                        logger << std::endl;
                    }
#endif
                    dominions++;

                    // Extend the dominion to its attractor in the full unsolved game <G>. Not
                    // needed for the first region: it was attracted in <G> with every tangle
                    // of <pl> (none has a higher priority), so it is already closed under
                    // attraction. Its tangles are deleted lazily once found out of the game.
                    // (If this pass extended the first region of the last one, its vertices of
                    // the last pass were not attracted again, so <escs> does not have the number
                    // of successors left in <G> of all vertices next to it: count those again.)
                    const bool recount = first_region and count_in_R;
                    if (!first_region) {
                        // start anew with escs
                        count_in_R = false;
                        std::fill(escs, escs+nodecount(), '\0');
                        std::fill(tescs.begin(), tescs.end(), '\0');

                        Q.clear();
                        for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) Q.push(v);

                        setCandidates(G, Z);
                        unsigned int qh = 0;
                        drain(pl, G, Z, INT_MAX, qh, -1);
                        Q.clear();
                    }

                    for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) {
                        Solver::solve(v, pl, str[v]);
                    }

                    G -= Z; // remove from G
                    if (recount) {
                        std::fill(gdeg, gdeg+nodecount(), 0); // (counted again when needed)
                        escs_is_gdeg = false;
                    } else {
                        escs_is_gdeg = true;
                    }
                    V.reset();
                    Z.reset();
                    keep = -1;
                    return true; // end tl so that we can restart...
                } else {
                    const auto D = dominions;
                    if (learnRegion(ridx, pl, pr)) {
                        new_tangles = true;
                    } else {
                        std::fill(pea_vidx, pea_vidx+nodecount(), '\0');
                        pea_curidx = 1;
                        for (auto v = V.find_last(); v != bitset::npos; v = V.find_prev(v)) {
                            if (extractTangles(v, Z, str)) new_tangles = true;
                        }
                    }

                    // Extend any dominions that were found.
                    // (Any solved vertices are now in <S>.)
                    if (D != dominions) {
                        // start anew with escs
                        count_in_R = false;
                        std::fill(escs, escs+nodecount(), '\0');
                        std::fill(tescs.begin(), tescs.end(), '\0');

                        Q.clear();
                        for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) Q.push(v);

                        setCandidates(G, S);
                        unsigned int qh = 0;
                        drain(pl, G, S, INT_MAX, qh, -1);
                        Q.clear();

                        for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) {
                            Solver::solve(v, pl, str[v]);
                        }

                        G -= S; // remove from G
                        escs_is_gdeg = true;
                        S.reset();
                        V.reset();
                        Z.reset();
                        keep = -1;
                        return true; // end tl so that we can restart...
                    }
                }
            }
        }

        if (pl == 0 ? regions[ridx].highest0 : regions[ridx].highest1) H[pl] = Z;
        if (ridx == 0 or ridx == xreg) {
            // (keep the vertices and heads of the first region and of an extended region, for
            // the next pass if it extends the region, see above)
            Z.swap(Zx);
            V.swap(Vx);
            zx = ridx;
        }
        V.reset();
        Z.reset();

        if (pr&1) highest1 = false;
        else highest0 = false;
        first_region = false;
    }

    // A new tangle is attracted (in the next pass) into the lowest region of its escapes, if
    // that region is of its player and not lower than the tangle, as its last escape is
    // attracted there. The next pass keeps the regions above the highest such region.
    // A region tangle that was not attracted again (as a whole) in a recomputed region, or
    // that lost vertices, is analysed now: its tangles are learned.
    if (analyseLazies(true)) { keep = -1; return true; }
    if (tpr.size() > pass_tangles) new_tangles = true;

    keep = -1;
    if (new_tangles) {
        // (a pass that keeps regions needs <regidx> while it runs, see attractTangle)
        if (!kept) {
            for (size_t k = 0; k < regions.size(); k++) {
                const size_t end = k+1 < regions.size() ? regions[k+1].begin : Q.size();
                for (size_t i = regions[k].begin; i < end; i++) regidx[Q[i]] = (int)k;
            }
        }
        int q = INT_MAX;
        seeds.clear();
        for (size_t t = pass_tangles; t < tpr.size(); t++) {
            const int k = tangleTarget((int)t);
            if (k < 0 or k > q) continue;
            if (k < q) { q = k; seeds.clear(); }
            seeds.push_back((int)t);
        }
        if (q != INT_MAX) keep = q;
    }

    return new_tangles;
}


/**
 * Options (comma separated): parmin=N (smallest queue processed in parallel),
 * chunk=N (vertices per task in a parallel round), pargame=N (smaller games are solved
 * sequentially, since the parallel rounds cost more than they gain on them).
 */
bool
TLParSolver::parseOptions(std::string& opts)
{
    std::string remaining = opts;
    while (!remaining.empty()) {
        const size_t comma = remaining.find(',');
        const std::string opt = remaining.substr(0, comma);
        remaining = comma == std::string::npos ? "" : remaining.substr(comma+1);
        const size_t eq = opt.find('=');
        const std::string key = opt.substr(0, eq);
        int value = 0;
        try {
            value = eq == std::string::npos ? 0 : std::stoi(opt.substr(eq+1));
        } catch (...) {
            value = 0;
        }
        if (key == "parmin" and value >= 1) {
            par_min = value;
        } else if (key == "chunk" and value >= 1) {
            par_chunk = value;
        } else if (key == "pargame" and value >= 0) {
            par_game = value;
        } else {
            logger << "unrecognised option for tlpar: " << opt << " (options: parmin=N, chunk=N, pargame=N)" << std::endl;
            return false;
        }
    }
    return true;
}

VOID_TASK_1(tlpar_run, TLParSolver*, s)
void tlpar_run_CALL(lace_worker* lw, TLParSolver* s)
{
    s->runPar(lw);
}

void
TLParSolver::run()
{
    // parallel attraction pays off only on large games; small ones are solved as by TL
    if (lace_worker_count() > 1 and nodecount() >= par_game) {
        pbuf.resize(lace_worker_count());
        tlpar_run(this); // runs runPar on a Lace worker
    } else {
        runPar(nullptr);
    }
}

void
TLParSolver::runPar(lace_worker* lw)
{
    this->lw = lw;

    tin_off = new size_t[nodecount()];
    tin_size = new unsigned int[nodecount()];
    tin_cap = new unsigned int[nodecount()];
    has_tin.resize(nodecount());
    str = new int[nodecount()];
    escs = new unsigned int[nodecount()];
    gdeg = new unsigned int[nodecount()];
    regidx = new int[nodecount()];
    astr = new int[nodecount()];
    LB.resize(nodecount());
    C.resize(nodecount());
    H[0].resize(nodecount());
    H[1].resize(nodecount());
    // if no vertex is disabled, every successor is in <G>; otherwise count lazily
    if (disabled.none()) {
        for (int v=0; v<nodecount(); v++) gdeg[v] = game.outcount(v);
    } else {
        std::fill(gdeg, gdeg+nodecount(), 0);
    }

    V.resize(nodecount());
    W.resize(nodecount());
    R.resize(nodecount());
    Z.resize(nodecount());
    Zx.resize(nodecount());
    Vx.resize(nodecount());
    S.resize(nodecount());
    G = disabled;
    G.flip();

    Q.resize(nodecount());

    tangleto.resize(nodecount());
    escapes.resize(nodecount());

    pea_state.resize(nodecount()*2);
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()];
    pea_root.resize(nodecount());

    while (G.any()) {
#ifndef NDEBUG
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "\033[m" << std::endl;
#endif
        iterations++;

        if (!tl()) break;
    }

    logger << "found " << dominions << " dominions and "<< tangles << " tangles." << std::endl;
    logger << "solved in " << iterations << " iterations." << std::endl;
    if (lw != nullptr) logger << "attracted " << par_vertices << " vertices in " << par_rounds << " parallel rounds." << std::endl;

#ifndef NDEBUG
    // Check if the whole game is now solved
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { THROW_ERROR("search was incomplete!"); }
    }
#endif

    // Free all explicitly allocated memory
    delete[] tin_off;
    delete[] tin_size;
    delete[] tin_cap;
    delete[] regidx;
    delete[] astr;
    delete[] str;
    delete[] pea_vidx;
    delete[] escs;
    delete[] gdeg;
}

}
