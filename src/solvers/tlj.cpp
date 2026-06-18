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
 * Tangle learning with justifications, after Algorithm 7 of Lapauw, Bruynooghe and
 * Denecker, "Improving Parity Game Solvers with Justifications" (VMCAI 2020).
 */

#include <algorithm>
#include <climits>

#include "tlj.hpp"

namespace pg {

static const int SOLVED = INT_MAX; // the region of a solved vertex (above every region)

TLJSolver::TLJSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}


TLJSolver::~TLJSolver()
{
}


/**
 * Compute SCCs in subgraph induced by <R> and <str>. Identical to the base
 * tangle-learning extraction; produces tangles and dominions (in <S>).
 */
bool
TLJSolver::extractTangles(int startvertex, bitset &R, int *str)
{
    bool new_tangles = false;
    const int pr = priority(startvertex);
    const int pl = pr&1;

    pea_state.push(startvertex);
    pea_state.push(0);
    pea_root[startvertex] = true;
    pea_vidx[startvertex] = pea_curidx++;
    pea_seen.push_back(startvertex);
    while (pea_state.nonempty()) {
pearce_again:
        const unsigned int n = pea_state.back2();
        unsigned int i = pea_state.back();

        if (owner(n) != pl) {
            auto edges = outs(n);
            if (i>0) {
                const int w = edges[i-1];
                if (pea_vidx[w] < pea_vidx[n]) {
                    pea_vidx[n] = pea_vidx[w];
                    pea_root[n] = false;
                }
            }
            for (;;) {
                const int to = edges[i];
                if (to == -1) break;
                if (R[to]) {
                    if (pea_vidx[to] == 0) {
                        pea_state.back() = i+1;
                        pea_state.push(to);
                        pea_state.push(0);
                        pea_root[to] = true;
                        pea_vidx[to] = pea_curidx++;
                        pea_seen.push_back(to);
                        goto pearce_again;
                    } else {
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
                if (pea_vidx[s] == 0) {
                    pea_state.back() = 1;
                    pea_state.push(s);
                    pea_state.push(0);
                    pea_root[s] = true;
                    pea_vidx[s] = pea_curidx++;
                    pea_seen.push_back(s);
                    goto pearce_again;
                }
            }
            if (pea_vidx[s] < pea_vidx[n]) {
                pea_vidx[n] = pea_vidx[s];
                pea_root[n] = false;
            }
        }
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

        const bool is_not_tangle = (tangle.size() == 1) and
            ((unsigned int)str[n] != n) and (str[n] != -1 or !game.has_edge(n, n));
        if (is_not_tangle) {
            tangle.clear();
            continue;
        }

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
            escapes.reset();
        }

        if (tangleto.empty()) {
            // dominion
            for (const int v : tangle) S[v] = true;
            dominions++;
            tangle.clear();
            new_tangles = true;
            continue;
        }

        // NOTE: the tangle is NOT linked into tin here; it is promoted (made
        // available for attraction) only after the current round settles (see run()).

        // store the tangle: its vertex-strategy pairs and its escapes (each ending with -1)
        tvpos.push_back(tvdata.size());
        for (const int v : tangle) {
            tvdata.push_back(v);
            tvdata.push_back(str[v]);
        }
        tvdata.push_back(-1);
        toutpos.push_back(toutdata.size());
        toutdata.insert(toutdata.end(), &tangleto[0], &tangleto[0]+tangleto.size());
        toutdata.push_back(-1);

        tpr.push_back(pr);
        tcnt.push_back(-1);
        tcpos.push_back(0);

        tangles++;
        tangle.clear();
        tangleto.clear();
        new_tangles = true;
    }

    pea_S.clear();
    return new_tangles;
}


/**
 * Find an escape of vertex <from> of the opponent of <pl> from region <p>, starting at its
 * out-edge <i>: a successor in G that is not already won by <pl> at a priority >= p, i.e.,
 * that is unassigned or in a lower region or in a region of the other player. Returns the
 * index of the out-edge, or -1 if there is none.
 */
inline int
TLJSolver::findEscape(const int from, int i, const int p, const int pl)
{
    const int *edges = outs(from);
    for (int to; (to = edges[i]) != -1; i++) {
        const int r = region[to];
        if (r < p) return i;                          // lower region or unassigned (-1 < p)
        if ((r&1) != pl && r != SOLVED) return i;     // higher region of the other player
    }
    return -1;
}


/**
 * Set the region of <v> to <r> (-1: unassigned), and keep track of the nonempty regions.
 */
inline void
TLJSolver::setRegion(const int v, const int r)
{
    const int old = region[v];
    if (old >= 0 && --rcount[old] == 0) { nonempty[old&1][old] = false; anyregion[old] = false; }
    if (r >= 0 && rcount[r]++ == 0) { nonempty[r&1][r] = true; anyregion[r] = true; }
    region[v] = r;
}


/**
 * The justification of <m> (in a region) no longer holds: reset it to unassigned, mark its
 * priority dirty and add it to <Dq> (to invalidate the justifications that reach it), unless
 * it is a head of its region, which is in it regardless: then it only loses its strategy.
 */
inline void
TLJSolver::invalidate(const int m)
{
    const int r = region[m];
    dirty[rprio[m]] = true;
    if (rprio[m] == r) {
        if (owner(m) == (r&1)) str[m] = -1;
        extracted[r] = false; // (its tangles change: each contains a head)
        return;
    }
    setRegion(m, -1);
    str[m] = -1;
    tangled[m] = false;
    Dq.push(m);
}


/**
 * Add <v> (not in region <p>, from a lower region or unassigned) to region <p> of <pl>
 * with strategy <s>, and record what this changes for the other regions.
 */
__attribute__((always_inline)) inline void
TLJSolver::pull(const int v, const int s, const int p, const int pl)
{
    const int oldr = region[v];
    // the lowest priority of a region of which <v> is now a new <pl>-target
    const int lo = (oldr >= 0 && (oldr&1) == pl) ? oldr+1 : 0;
    if (lo < pull_lo) pull_lo = lo;
    if (oldr != -1 && (oldr&1) != pl) flipped.push(v);
    pull_count++;
    setRegion(v, p);
    str[v] = s;
    Z[v] = true;
    zlist.push(v);
}


/**
 * Try to attract tangle <t> into region <p> (player <pl>) after its escape at <vpos> in
 * <zlist> has been popped: possible if the tangle is of the right parity, not deleted,
 * still entirely in G, and all of its escapes are already won by <pl> at priority >= p.
 * Adds attracted vertices to region <p> (see pull).
 */
bool
TLJSolver::attractTangleInto(const int t, const unsigned vpos, const int p, const int pl)
{
    if (tpr[t] == -1) return false;       // deleted
    if (tpr[t] > p) return false;         // cannot be attracted this high
    if ((tpr[t]&1) != pl) return false;   // wrong parity

    // all escapes must be won by <pl> at priority >= p: count the others when first seen;
    // each that joins region p later is subtracted when it is popped (-2: done for now)
    if (tcnt[t] == -2) return false;
    if (tcnt[t] == -1) {
        int c = 0;
        for (int *ptr = &toutdata[toutpos[t]]; *ptr != -1; ptr++) {
            const int e = *ptr;
            if (!G[e]) continue;
            if (region[e] < p || (region[e]&1) != pl) c++;
        }
        tcnt[t] = c;
        tcpos[t] = (unsigned)zlist.size();
        ttouched.push_back(t);
    } else if (vpos >= tcpos[t]) {
        tcnt[t]--;
    }
    if (tcnt[t] != 0) return false;
    tcnt[t] = -2;

    // verify the tangle is still entirely in G (on-the-fly delete otherwise)
    for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
        if (!G[*ptr]) { tpr[t] = -1; return false; }
    }

    bool attracted = false;
    for (int *ptr = &tvdata[tvpos[t]]; *ptr != -1; ptr += 2) {
        const int u = *ptr;
        const int s = *(ptr+1);
        if (region[u] == p) continue;     // already in region p
        if (region[u] > p) continue;      // in a higher region (should not happen)
        pull(u, s, p, pl);
        if (owner(u) != pl) tangled[u] = true;
        attracted = true;
    }
    return attracted;
}


/**
 * Invalidate the justification closure of <start>: every vertex m won by the
 * opponent of <pl> (region parity != pl) whose justification reaches <start>
 * (str[m] == the popped vertex, or m is a "loser" with str == -1) is reset to
 * unassigned and its priority marked dirty. Used after pulling a wrong-parity
 * vertex up into region <p>.
 */
void
TLJSolver::invalidateOpponent(const int start, const int p, const int pl)
{
    Dq.push(start);
    while (Dq.nonempty()) {
        const int v = Dq.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            const int m = *curedge;
            if (!G[m]) continue;
            if (region[m] == -1) continue;     // already unassigned
            if (region[m] >= p) continue;      // not in a lower region
            if ((region[m]&1) == pl) continue; // won by pl: handled by attraction, keep
            if (str[m] == v || str[m] == -1) invalidate(m); // justification reaches v
        }
    }
}


/**
 * Invalidate the justification closure of the vertices in <Dq> (reset to unassigned after a
 * dominion was solved): every vertex whose justification reaches them is reset to unassigned
 * and marked dirty, regardless of parity.
 */
void
TLJSolver::invalidateSolved()
{
    while (Dq.nonempty()) {
        const int v = Dq.pop();
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            const int m = *curedge;
            if (!G[m]) continue;
            if (region[m] == -1) continue;
            if (str[m] == v || str[m] == -1) invalidate(m);
        }
    }
}


/**
 * Solve dominion <Dom> of player <pl>: attract maximally to it over G (unless <attract> is
 * false: then it is already closed under attraction in G, and <touched> holds the vertices
 * in G with an edge into it, with in <cnt> their successors in G outside it), solve it and
 * remove it from G.
 */
void
TLJSolver::solveDominion(bitset &Dom, const int pl, const bool attract)
{
    if (attract) {
        // attract maximally over G (the opponent's vertices with a counter of their
        // successors in G that are not yet in Dom)
        for (auto v = Dom.find_first(); v != bitset::npos; v = Dom.find_next(v)) Q.push(v);
        while (Q.nonempty()) {
            const int v = Q.pop();
            for (auto curedge = ins(v); *curedge != -1; curedge++) {
                const int from = *curedge;
                if (!G[from] || Dom[from]) continue;
                if (owner(from) == pl) {
                    Dom[from] = true; str[from] = v; Q.push(from);
                } else {
                    // opponent: attracted iff all edges (in G) lead into Dom
                    if (cnt[from] < 0) {
                        cnt[from] = gdeg[from];
                        touched.push(from);
                    }
                    if (--cnt[from] == 0) { Dom[from] = true; str[from] = -1; Q.push(from); }
                }
            }
        }
    }

    // solve and remove from G
    G -= Dom;
    for (auto v = Dom.find_first(); v != bitset::npos; v = Dom.find_next(v)) {
        Solver::solve(v, pl, str[v]);
        setRegion(v, -1);
        region[v] = SOLVED;
    }

    // the vertices left in G with an edge into Dom (all of the opponent, as the others
    // are attracted) are in <touched>: update their number of successors in G, and
    // invalidate the strategies that play to Dom. A vertex justified by all its successors
    // in G (str == -1) stays justified when some of them leave G, as they were all in its
    // region (or higher) before it, unless a tangle attracted it (with vertices that may
    // have left G; a needless reset if it was attracted again since); its region may have
    // new tangles though.
    for (unsigned i = 0; i < touched.size(); i++) {
        const int m = touched[i];
        const int c = cnt[m];
        cnt[m] = -1;
        if (!G[m]) continue;
        gdeg[m] = c;
        if (region[m] == -1) continue;
        if (str[m] == -1 && !tangled[m]) extracted[region[m]] = false;
        else if (str[m] == -1 || Dom[str[m]]) invalidate(m);
    }
    touched.clear();
    invalidateSolved();
}


/**
 * Store the vertices of region <p> (in <zlist>) as its members: append them to <mpool>,
 * after removing the vertices that left their regions from it when it is full.
 */
void
TLJSolver::setMembers(const int p)
{
    const unsigned zn = zlist.size();
    if (mpool.size() + zn > mpool.capacity()) {
        // compact: keep for each region only its current vertices (at most all vertices)
        std::vector<int> fresh;
        fresh.reserve(std::max(mpool.capacity(), 2*(size_t)nodecount()+zn));
        for (int q = 0; q < (int)mlen.size(); q++) {
            const unsigned b = (unsigned)fresh.size();
            for (unsigned i = mbeg[q], e = i+mlen[q]; i < e; i++) {
                if (region[mpool[i]] == q) fresh.push_back(mpool[i]);
            }
            mbeg[q] = b;
            mlen[q] = (unsigned)fresh.size()-b;
        }
        mpool.swap(fresh);
    }
    mbeg[p] = (unsigned)mpool.size();
    mlen[p] = zn;
    mpool.insert(mpool.end(), &zlist[0], &zlist[0]+zn);
}


/**
 * Grow region <p> (player pl=p&1): re-attract from the current region-p vertices
 * and any unassigned vertices of the priorities of region p, pulling in vertices from
 * lower regions and attracting applicable tangles. The existing region p is never reset
 * (growth is monotone). After growth, reset the justification closure of any
 * wrong-parity vertices that were pulled up, then solve the region if it is a dominion,
 * or else extract its tangles (and dominions) if it is closed.
 */
void
TLJSolver::improveRegion(int p)
{
    const int pl = p&1;

    if (zlast >= 0) {
        for (unsigned i = mbeg[zlast], e = i+mlen[zlast]; i < e; i++) Z[mpool[i]] = false;
        zlast = -1;
    }
    zlist.clear();
    heads.clear();
    pull_lo = p;
    pull_count = 0;
    flipped.clear();

    // seed: the heads (the vertices of the priorities of region p that are in region p or
    // unassigned, in order), then the other vertices of region p. Region p only gains
    // vertices here, so these are among the vertices of the last call for p (its members).
    bool new_head = false;
    for (int v = gfirst[p]; v < gend[p]; v++) {
        if (region[v] == -1) {
            setRegion(v, p);
            str[v] = -1;
            new_head = true;
        } else if (region[v] != p) {
            continue; // in a higher region, or solved
        }
        Z[v] = true;
        zlist.push(v);
        heads.push_back(v);
    }
    for (unsigned i = mbeg[p], e = i+mlen[p]; i < e; i++) {
        const int v = mpool[i];
        if (region[v] != p || rprio[v] == p) continue; // (not in region p, or a head)
        Z[v] = true;
        zlist.push(v);
    }

    if (zlist.empty()) {
        mlen[p] = 0;
        return; // nothing at this priority
    }

    // without a region of <pl> above p, the escapes of a vertex of the opponent are all its
    // successors in G that are not in region p (and every vertex in region p is popped once)
    const bool pl_above = nonempty[pl].find_next(p) != bitset::npos;

    // grow region p (<zlist> is the queue)
    for (unsigned qi = 0; qi < zlist.size(); qi++) {
        const int v = zlist[qi];
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            const int from = *curedge;
            if (region[from] >= p) continue; // in region p or a higher region, or solved
            // candidate (region < p or unassigned)
            if (owner(from) == pl) {
                pull(from, v, p, pl);
            } else {
                // the opponent: attracted when it has no escape (a successor in G that is
                // not in region p or in a higher region of <pl>). Without a region of <pl>
                // above, count them (all successors not in region p; each is subtracted
                // when it is popped); otherwise keep one escape as a witness, and look for
                // the next one when it joins region p (during the attraction, a successor
                // that is not an escape remains one)
                if (cnt[from] < 0) {
                    touched.push(from);
                    if (pl_above) {
                        const int e = findEscape(from, 0, p, pl);
                        cnt[from] = e < 0 ? 0 : 1;
                        cpos[from] = (unsigned)e;
                    } else {
                        cnt[from] = gdeg[from]-1; // (<v> is one)
                    }
                } else if (!pl_above) {
                    cnt[from]--;
                } else if (region[outs(from)[cpos[from]]] == p) {
                    const int e = findEscape(from, cpos[from]+1, p, pl);
                    if (e < 0) cnt[from] = 0;
                    else cpos[from] = (unsigned)e;
                }
                if (cnt[from] == 0) pull(from, -1, p, pl);
            }
        }
        // tangle attraction (tangles escaping to <v> may now be attractable)
        if (has_tin[v]) for (int t : tin[v]) attractTangleInto(t, qi, p, pl);
    }
    for (const int t : ttouched) tcnt[t] = -1;
    ttouched.clear();

    // the heads of <pl> without a strategy play to region p if they can (only heads can
    // lack a strategy)
    for (const int h : heads) {
        if (owner(h) != pl || str[h] != -1) continue;
        for (auto curedge = outs(h); *curedge != -1; curedge++) {
            if (region[*curedge] == p) { str[h] = *curedge; break; }
        }
    }

    // each pulled-up vertex now sits in region p, so it is a valid <pl>-target for
    // every region of parity <pl> at or below p that it did not satisfy before;
    // re-examine those (a predecessor may now attract/close there). For a same-parity
    // move it already satisfied everything up to oldr, so only (oldr, p) is new.
    // NOTE: marking both parities is conservative but correct; restricting to
    // q&1==pl misses cases where an opponent vertex must be re-evaluated.
    // (All these ranges end at p, so their union is one range from <pull_lo>; the part
    // from <sweep_lo> is already marked, see run().)
    if (pull_lo < sweep_lo) {
        const int hi = sweep_lo < p ? sweep_lo : p;
        auto q = pull_lo == 0 ? anyregion.find_first() : anyregion.find_next(pull_lo-1);
        for (; q != bitset::npos && (int)q < hi; q = anyregion.find_next(q)) dirty[q] = true;
        sweep_lo = pull_lo;
    }

    // reset the justification closure of wrong-parity vertices pulled up
    for (unsigned i = 0; i < flipped.size(); i++) invalidateOpponent(flipped[i], p, pl);

    // is the region locally closed? (no head is open: a pl-owned head with no
    // in-region move, or an opponent head that can escape to a lower region)
    bool leaks = false;
    bool up = false; // whether an opponent head can escape to a higher region
    for (auto it = heads.rbegin(); it != heads.rend(); it++) {
        const int v = *it;
        if (owner(v) == pl) {
            if (str[v] == -1) { leaks = true; break; }
        } else {
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                const int to = *curedge;
                if (!G[to]) continue;
                if (region[to] < p) { leaks = true; break; }
                if (region[to] > p) up = true;
            }
            if (leaks) break;
        }
    }

    if (!leaks && !up && !pl_above) {
        // The region is closed in G: the opponent's heads cannot leave it, and its other
        // vertices have all their successors in region p or in a higher region of <pl>, of
        // which there is none. Every cycle in it passes a head (the highest priority), so
        // it is a dominion (as in tl, for the highest region of a player). It is closed
        // under attraction in G unless there are higher regions (of the opponent).
        dominions++;
        const bool any_above = nonempty[1-pl].find_next(p) != bitset::npos;
        if (any_above) {
            for (unsigned i = 0; i < touched.size(); i++) cnt[touched[i]] = -1;
            touched.clear();
        }
        solveDominion(Z, pl, any_above);
        Z.reset();
        mlen[p] = 0;
        restart = true;
        return;
    }
    setMembers(p);
    zlast = p;
    for (unsigned i = 0; i < touched.size(); i++) cnt[touched[i]] = -1;
    touched.clear();

    // a region that grew (gained a vertex) has new tangles to find; one that did
    // not grow but only now closed must still be extracted once
    if (new_head || pull_count != 0) extracted[p] = false;

    if (leaks) { extracted[p] = false; return; }

    if (extracted[p]) return; // already extracted in this (closed) state

    // closed region: extract tangles (and possibly dominions)
    const auto dom_before = dominions;

    pea_curidx = 1;
    for (auto it = heads.rbegin(); it != heads.rend(); it++) {
        extractTangles(*it, Z, str);
    }
    for (const int v : pea_seen) pea_vidx[v] = 0; // (all are 0 before an extraction)
    pea_seen.clear();
    extracted[p] = true;

    if (dom_before != dominions) {
        // a dominion was found (in <S>): solve it maximally and remove from G
        solveDominion(S, pl, true);
        S.reset();
        restart = true;
    }
}


/**
 * Promote the tangles found since the last promotion: link them into tin (making them
 * available for attraction) and mark the region they can now be attracted into dirty: a
 * tangle can only be attracted into a region of its player if all its escapes are in
 * regions of that player at or above it, so first into the lowest region of its escapes,
 * which attracts it as soon as it is improved (lower regions are then marked by pull).
 * New tangles are only promoted at the end of a sweep, so a region's own new tangle does
 * not immediately re-grow it (which would re-extract and explode the work).
 */
void
TLJSolver::promote()
{
    for (size_t t = promoted; t < tpr.size(); t++) {
        if (tpr[t] == -1) continue;
        const int tpl = tpr[t]&1;
        int minreg = INT_MAX;
        bool good = true; // whether all escapes are in regions of <tpl>
        for (int *ptr = &toutdata[toutpos[t]]; *ptr != -1; ptr++) {
            const int e = *ptr;
            tin[e].push_back((int)t);
            has_tin[e] = true;
            if (!G[e]) continue;
            const int r = region[e];
            if (r == -1 || (r&1) != tpl) good = false;
            else if (r < minreg) minreg = r;
        }
        if (good && minreg != INT_MAX && minreg >= tpr[t]) dirty[minreg] = true;
    }
    promoted = tpr.size();
}


void
TLJSolver::run()
{
    tin = new std::vector<int>[nodecount()];
    has_tin.resize(nodecount());
    tangled.resize(nodecount());
    str = new int[nodecount()];
    region = new int[nodecount()];

    Z.resize(nodecount());
    S.resize(nodecount());
    G = disabled;
    G.flip();

    Q.resize(nodecount());
    Dq.resize(nodecount());
    zlist.resize(nodecount());
    touched.resize(nodecount());
    flipped.resize(nodecount());

    tangleto.resize(nodecount());
    escapes.resize(nodecount());

    pea_state.resize(nodecount()*2);
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()]();
    pea_root.resize(nodecount());

    const int n = nodecount();
    cnt.assign(n, -1);
    cpos.resize(n);
    gdeg.resize(n);
    rprio.resize(n);

    const int max_prio = G.any() ? priority(G.find_last()) : -1;
    dirty.resize(max_prio+1);
    extracted.resize(max_prio+1);
    rcount.assign(max_prio+1, 0);
    nonempty[0].resize(max_prio+1);
    nonempty[1].resize(max_prio+1);
    anyregion.resize(max_prio+1);
    mbeg.assign(max_prio+1, 0);
    mlen.assign(max_prio+1, 0);
    mpool.reserve(2*n);
    gfirst.assign(max_prio+1, 0);
    gend.assign(max_prio+1, 0);

    // Compress the priorities on the fly, as tl does: a region is identified by the highest
    // priority of a run of consecutive priorities of the same parity (of the vertices in G),
    // and all vertices of that run are its heads. (The vertices are sorted by priority.)
    const bool all_in_g = disabled.none();
    int top = -1;
    for (int v = n-1; v >= 0; v--) {
        str[v] = -1;
        if (!all_in_g && disabled[v]) {
            region[v] = SOLVED;
            rprio[v] = -1;
            gdeg[v] = 0;
            continue;
        }
        region[v] = -1;
        const int pr = priority(v);
        if (top == -1 || ((pr^top)&1) != 0) {
            top = pr;
            gend[top] = v+1;
            dirty[top] = true;
        }
        rprio[v] = top;
        gfirst[top] = v;
        if (all_in_g) {
            gdeg[v] = game.outcount(v);
        } else {
            int c = 0;
            for (auto curedge = outs(v); *curedge != -1; curedge++) if (G[*curedge]) c++;
            gdeg[v] = c;
        }
    }

    while (G.any()) {
        // inner loop: improve all dirty regions using the promoted tangles. This is a sweep
        // from the highest region down: improving region p only marks regions below p dirty,
        // either directly or all nonempty regions from <sweep_lo> up to p (each such range
        // only once), except after a dominion, which restarts the sweep. (An empty region
        // only gains heads that are marked dirty directly.)
        const int nprio = (int)dirty.size();
        int cur = nprio;
        sweep_lo = INT_MAX;
        while (true) {
            const auto d = cur == nprio ? dirty.find_last() : dirty.find_prev(cur);
            if (d == bitset::npos) break;
            const int p = (int)d;
            dirty[p] = false;
            cur = p;
            if (gend[p] == 0) continue; // not the priority of a region
            iterations++;
            improveRegion(p);
            if (!G.any()) break;
            if (restart) {
                // a dominion was solved: the decomposition can change widely, so all
                // nonempty regions are re-examined (the regions that are empty only gain
                // vertices that were reset, and are marked dirty for those)
                dirty |= anyregion;
                // a new sweep (from the highest region) can use the new tangles
                promote();
                restart = false;
                cur = nprio;
                sweep_lo = INT_MAX;
            }
        }
        if (!G.any()) break;

        // the sweep is done: promote the tangles found in it
        if (promoted == tpr.size()) break; // nothing new and nothing dirty: done
        promote();
    }

    logger << "found " << dominions << " dominions and "<< tangles << " tangles." << std::endl;
    logger << "solved in " << iterations << " iterations." << std::endl;

#ifndef NDEBUG
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { THROW_ERROR("search was incomplete!"); }
    }
#endif

    delete[] tin;
    delete[] str;
    delete[] pea_vidx;
    delete[] region;
}

}
