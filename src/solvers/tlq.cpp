/*
 * Copyright 2024-2025 Tom van Dijk, University of Twente
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

#include "tlq.hpp"

namespace pg {

/**
 * Implementation of the quasi-polynomial time recursive algorithm, extended with tangles.
 * The attractor computation is extended with attracting tangles. Each step in the algorithm,
 * tangles are learned when the set of attracted vertices is locally closed.
 */

TLQSolver::TLQSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

TLQSolver::~TLQSolver()
{
}

/**
 * Returns True iff player <pl> attracts vertex <v> to region <Z> in subgame <R>.
 * (Either the vertex is owned by <pl> and can play to <Z>, or is owned by the opponent
 *  and all (and at least one) edges go to <Z>)
 */
bool
TLQSolver::attracts(const int pl, const int v, bitset &Z, bitset &R)
{
    if (owner(v) == pl) {
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            if (Z[*curedge]) return true;
        }
        return false;
    } else {
        bool hasZ = false;
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (Z[to]) { hasZ = true; continue; }
            if (R[to]) return false;
        }
        return hasZ;
    }
}


/**
 * Attract vertices in <R> to <v> in region <Z> as player <pl>.
 * The opponent may escape to vertices of <R> that are not in <Z>.
 * Attracted vertices are added to <Z> and to the queue <Q>.
 */
void
TLQSolver::attractVertices(const int pl, const int v, bitset &Z, bitset &R)
{
    auto curedge = ins(v);
    for (int from = *curedge; from != -1; from = *++curedge) {
        if (Z[from]) {
            // already in Z, set strategy if not yet set
            if (owner(from) == pl and str[from] == -1) {
#ifndef NDEBUG
                if (trace >= 3) {
                    logger << "\033[1;37msetting strategy of \033[36m" << label_vertex(from) << "\033[m to \033[1;36m" << label_vertex(v) << "\033[m" << std::endl;
                }
#endif
                str[from] = v;
            }
        } else if (R[from]) {
            // a vertex in <R> that is not yet in <Z>
            if (owner(from) != pl) {
                unsigned int e;
                // Count unprocessed successors on first touch, then consume this edge.
                // v stays unprocessed until its incoming edges have been handled, so the
                // count includes the edge to v that is consumed here.
                if (esc_epoch[from] != attract_epoch) {
                    esc_epoch[from] = attract_epoch;
                    if (count_in_G) {
                        // (no successor of <from> was processed yet, see attract)
                        e = gdeg[from];
                        if (e == 0) {
                            for (auto edge = outs(from); *edge != -1; edge++) {
                                if (G[*edge]) e++;
                            }
                            gdeg[from] = e;
                        }
#ifndef NDEBUG
                        unsigned int e2 = 0;
                        for (auto edge = outs(from); *edge != -1; edge++) {
                            if (G[*edge]) e2++;
                        }
                        if (e != e2) THROW_ERROR("wrong number of successors in G");
#endif
                        touched.push_back(from);
                    } else {
                        e = 0;
                        for (auto edge = outs(from); *edge != -1; edge++) {
                            if (R[*edge] and !Processed[*edge]) e++;
                        }
                    }
                } else {
                    e = escs[from];
                }
                e -= R[v] ? 1 : 0;
                escs[from] = e;
                if (e != 0) continue;
            }
            // vertex does not escape, so attract it
            Z[from] = true;
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
 * Drain one attraction in FIFO order, with a fresh epoch for its fixed subgame R.
 * Queued seeds are unprocessed even when already in Z. Some seeds (for example k_U)
 * are outside R: they trigger attraction but must not decrement its escape counters.
 * If crossing is given, report whether the attraction visits it.
 * With resume, continue the last attraction (to the same Z in the same R) with new seeds that
 * were not processed yet, as for the next batch of heads of a region: Z has only grown since,
 * so its escape counts and completed tangles are still valid, and need not be computed again.
 *
 * If R is the unsolved game G (the attraction of a dominion, or of the top level while it is
 * the whole game), the escape count of a vertex starts from its number of successors in G,
 * kept in <gdeg>, instead of counting them. This is the same count: the seeds in Z are
 * queued, Z holds no other vertex of R, and a successor processed earlier would have counted
 * the vertex then. The counted vertices are recorded in <touched> (see updateDegrees).
 */
bool
TLQSolver::attract(int pl, bitset &Z, bitset &R, const bitset *crossing, bool resume)
{
    if (Q.empty()) return false;
    if (!resume) {
        if (++attract_epoch == 0) {
            std::fill(esc_epoch, esc_epoch+nodecount(), 0u);
            std::fill(tdone.begin(), tdone.end(), 0u);
            ++attract_epoch;
        }
        count_in_G = &R == &G or (top_is_G and &R == Rs[0]);
        touched.clear();
        if (!count_in_G) {
            Processed = Z;
            for (unsigned i=0; i<Q.size(); i++) Processed[Q[i]] = false;
        }
#ifndef NDEBUG
        if (count_in_G and R != G) THROW_ERROR("top level is not the whole game");
#endif
    }
    bool crossed = false;
    for (unsigned qi=0; qi<Q.size(); qi++) {
        const int v = Q[qi];
        if (crossing and (*crossing)[v]) crossed = true;
        attractVertices(pl, v, Z, R);
        Processed[v] = true;
        attractTangles(pl, v, R, Z, R);
    }
    Q.clear();
    return crossed;
}


/**
 * After removing a dominion from G: if it was attracted in G (see attract), each vertex left
 * in G with a successor in it was counted, and its count of successors not in the dominion
 * (in <escs>) is its new number of successors in G. Otherwise, count them again lazily.
 */
void
TLQSolver::updateDegrees(bool counted)
{
    if (counted) {
        for (const int u : touched) if (G[u]) gdeg[u] = escs[u];
    } else {
        std::fill_n(gdeg, nodecount(), 0u);
    }
    touched.clear();
}


/**
 * Try to attract tangle <t> for player <pl> to attractor set <Z>
 * All vertices in tangle <t> must be in <Z+R> and the opponent may not escape to subgame <R\Z>
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy
 */
bool
TLQSolver::attractTangle(int t, const int pl, bitset &R, bitset &Z, bitset &G)
{
    /**
     * Check if tangle is won by player <pl> and not deleted
     */
    {
        const int tangle_pr = tpr[t];
        if (tangle_pr == -1) return false; // deleted tangle
        if (pl != (tangle_pr&1)) return false; // not of desired parity
        if (tdone[t] == attract_epoch) return false;
    }

    /**
     * Check if the tangle can escape to G\Z.
     */
    {
        const int *ptr = &toutdata[toutpos[t]];
        int v;
        while ((v=*ptr++) != -1) {
            if (Z[v]) continue;
            if (G[v]) return false; // opponent escapes
        }
    }

    /**
     * Check if tangle is contained in Z+R and if any vertices are not already in Z
     * We require at least one vertex in Z\R. Otherwise, don´t attract.
     */
    {
        bool can_attract_new = false;
        const int *ptr = &tvdata[tvpos[t]];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            ptr++; // skip strategy
            if (!this->G[v]) {
                // on-the-fly detect out-of-game tangles
                tpr[t] = -1; // delete the tangle
                return false; // is now a deleted tangle
            } else if (Z[v]) {
                continue; // already attracted
            } else if (!R[v]) {
                return false; // not contained in Z+R
            } else {
                can_attract_new = true; // has vertices not yet attracted
            }
        }
        if (!can_attract_new) {
            tdone[t] = attract_epoch;
            return false;
        }
        // R and global G stay fixed and Z only grows during this attraction.
        tdone[t] = attract_epoch;
    }

    /**
     * Attract!
     */
    {
        const int *ptr = &tvdata[tvpos[t]];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            const int s = *ptr++;
            if (Z[v]) continue; // already in <Z>
            Z[v] = true;
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
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 */
void
TLQSolver::attractTangles(int pl, int v, bitset &R, bitset &Z, bitset &G)
{
    if (has_tin.empty() or !has_tin[v]) return;
    const Backlinks &links = tin_pages[v/TIN_PAGE_SIZE][v%TIN_PAGE_SIZE];
    const int *ts = &tin_pool[links.offset];
    for (unsigned int i = 0; i < links.size; i++) attractTangle(ts[i], pl, R, Z, G);
}


void
TLQSolver::addBacklinks(const int tidx)
{
    if (has_tin.empty()) {
        has_tin.resize(nodecount());
        tin_pages.resize((nodecount()+TIN_PAGE_SIZE-1)/TIN_PAGE_SIZE);
    }
    for (const int e : tangleto) {
        auto &page = tin_pages[e/TIN_PAGE_SIZE];
        if (!page) page.reset(new Backlinks[TIN_PAGE_SIZE]);
        auto &links = page[e%TIN_PAGE_SIZE];
        if (!has_tin[e]) {
            has_tin[e] = true;
            links = {};
        }
        if (links.size == links.capacity) {
            // move the list of <e> to a new block of twice the size at the end of the pool
            const unsigned int cap = links.capacity ? 2*links.capacity : 2;
            const size_t off = tin_pool.size();
            tin_pool.resize(off+cap);
            if (links.size) std::copy_n(&tin_pool[links.offset], links.size, &tin_pool[off]);
            links.offset = off;
            links.capacity = cap;
        }
        tin_pool[links.offset+links.size++] = tidx;
    }
}


/**
 * Compute SCCs in subgraph induced by <R> and <str>.
 * Start the SCC computation at vertex <startvertex>.
 * Every SCC is then processed as a tangle.
 * If the tangle is closed, it is a dominion and added to <S> and <Q>.
 */
bool
TLQSolver::extractTangles(int startvertex, bitset &R)
{
    bool new_tangles = false;
    const int pr = priority(startvertex);
    const int pl = pr&1;

    /**
     * The following is the nonrecursive implementation of David Pearce,
     * "A space-efficient algorithm for finding strongly connected components" (IPL, 2016),
     * modified to on-the-fly restrict the graph by <R> and <str>.
     */

    if (pea_vidx[startvertex] != 0) return false; // already did?

    // beginVisiting
    unsigned int next_vertex = startvertex;
begin_visiting:
    pea_state.push_back({next_vertex, 0});
    pea_root[next_vertex] = true;
    pea_vidx[next_vertex] = pea_curidx++;
    while (!pea_state.empty()) {
        // visitLoop
        const unsigned int n = pea_state.back().vertex;
        unsigned int i = pea_state.back().edge;

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
                        pea_state.back().edge = i+1;
                        next_vertex = to;
                        goto begin_visiting;
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
                    pea_state.back().edge = 1;
                    next_vertex = s;
                    goto begin_visiting;
                }
            }
            // finishEdge
            if (pea_vidx[s] < pea_vidx[n]) {
                pea_vidx[n] = pea_vidx[s];
                pea_root[n] = false;
            }
        }
        // finishVisiting
        pea_state.pop_back();
        if (pea_root[n]) {
            pea_curidx -= 1;
            tangle.push_back(n);
            while (!pea_S.empty()) {
                const int t = pea_S.back();
                if (pea_vidx[n]>pea_vidx[t]) break;
                pea_S.pop_back();
                pea_curidx -= 1;
                pea_vidx[t] = (unsigned int)-1;
                tangle.push_back(t);
            }
            pea_vidx[n] = (unsigned int)-1;
        } else {
            pea_S.push_back(n);
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
                        tangleto.push_back(to);
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
            logger << " with " << tangleto.size() << " escapes";
            if (trace >= 2) {
                for (unsigned int x = 0; x < tangleto.size(); x++) {
                    logger << " \033[1;36m" << label_vertex(tangleto[x]) << "\033[m";
                }
            }
            logger << std::endl;
        }
#endif

        // new tangle idx
        const int tidx = tpr.size();

        addBacklinks(tidx);

        // store the tangle: its vertex-strategy pairs and its escapes (each ending with -1),
        // in flat arrays as in TL, instead of two allocations per tangle
        toutpos.push_back(toutdata.size());
        toutdata.insert(toutdata.end(), tangleto.begin(), tangleto.end());
        toutdata.push_back(-1);
        tvpos.push_back(tvdata.size());
        for (const int v : tangle) {
            tvdata.push_back(v);
            tvdata.push_back(str[v]);
        }
        tvdata.push_back(-1);

        // and set p to pr
        tpr.push_back(pr);
        tdone.push_back(0);

        tangles++;
        tangle.clear();
        tangleto.clear();
        new_tangles = true;
    }

    pea_S.clear();
    return new_tangles;
}


void
TLQSolver::solve(bitset &SG, int vtop, int pe, int po)
{
    /**
     * This is based on a universal tree of height pr/2, with parameter n associated with pe or po
     */

    iterations++;  // record the number of recursive calls (visits)

    /**
     * If precision is 0 for a player, it means we are only looking for tangles of size <= 0,
     * that is, we presume everything is won by the other player.
     */

    if (po <= 0) {
        W0 |= SG;

#ifndef NDEBUG
        if (trace>=2 and SG.any()) {
            logger << "presumed won by player 0:";
            for (int v=vtop; v>=0; v--) {
                if (SG[v]) logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }
#endif

        return;
    }

    if (pe <= 0) {
        W1 |= SG;

#ifndef NDEBUG
        if (trace>=2 and SG.any()) {
            logger << "presumed won by player 1:";
            for (int v=vtop; v>=0; v--) {
                if (SG[v]) logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }
#endif

        return;
    }

    /**
     * Move <vtop> to the true top (highest vertex that is in SG)
     */

    if (!SG[vtop]) {
        vtop = SG.find_prev(vtop);
        if (((size_t)vtop) == bitset::npos) return;
    }

    /**
     * Obtain the priority and player of the subgame
     */

    if (Ps.size() < Rs.size()) {
        Ps.push_back(priority(vtop));
        // optimisation when we enter this method for the first time!
        auto size = SG.count();
        while ((size_t)(po/2) > size) po/=2;
        while ((size_t)(pe/2) > size) pe/=2;
    }

    const int pr = Ps.back();
    const int pl = pr&1;

    /**
     * Expand the left side of the current level of the universal tree
     */

    if (pl == 0) {
        solve(SG, vtop, pe, po/2);
    } else {
        solve(SG, vtop, pe/2, po);
    }

#ifndef NDEBUG
    if (trace >= 2) logger << "in pr=" << pr << ", pe=" << pe << ", po=" << po << std::endl;
#endif

    auto &Wm = pl == 0 ? W0 : W1; // my W
    auto &Wo = pl == 0 ? W1 : W0; // opponent's W

    /**
     * Update vtop to the highest vertex in the intersection
     */

    if (!SG[vtop]) {
        vtop = SG.find_prev(vtop);
        if (((size_t)vtop) == bitset::npos) return;
    }

    /**
     * Compute H := Attr(vertices of <pr> in <R>)
     */

    bitset H(nodecount());
    Hs.push_back(&H);

    // Seed the same-parity heads together, as in TL, before expanding their attractor.
    // (the batches after the first continue its attraction, see attract)
    int next = vtop;
    while (next >= 0) {
        const bool resume = !heads.empty();
        for (; next >= 0; next--) {
            if (!SG[next] or H[next]) continue;
            if ((priority(next)&1) != pl) break;
            heads.push_back(next);
            H[next] = true;
            str[next] = -1;
            Q.push(next);
        }
        if (Q.empty()) break;
        attract(pl, H, SG, nullptr, resume);
        // Attraction can remove the opposite-parity vertex that stopped this batch.
    }

#ifndef NDEBUG
    if (trace) {
        logger << "\033[1;33mregion \033[36m" << pr << "\033[m";
        for (auto v = H.find_last(); v != bitset::npos; v = H.find_prev(v)) {
            logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
            if (trace >= 2) {
                if (owner(v) == pl) {
                    if (str[v] == -1) logger << "\033[38;5;37m\033[33m=>\033[38;5;37m-\033[m";
                    else logger << "\033[33m=>\033[38;5;37m =>" << label_vertex(str[v]) << "\033[m";
                }
            }
        }
        logger << std::endl;
    }
#endif

    /**
     * Now check if H has any escapes to R-H
     * A head of <pl> can stay in H iff it has a strategy: it got one from the first of its
     * successors in H that attracted from it (one attracting before the head was added would
     * have attracted the head). If H is the whole subgame (no vertex was left to stop the
     * heads), there is nothing to check, as every vertex of a subgame has a successor in it.
     */

    bool leaky = false;
    if (next >= 0) {
        for (auto v : heads) {
            if (owner(v) == pl ? str[v] == -1 : !attracts(pl, v, H, SG)) {
                leaky = true;
                break;
                // TODO: consider learning tangles from partially leaky regions
            }
        }
    }
#ifndef NDEBUG
    for (auto v : heads) {
        const bool closed = attracts(pl, v, H, SG);
        if (next < 0 and !closed) THROW_ERROR("head of the whole subgame is open");
        if (owner(v) == pl and closed != (str[v] != -1)) THROW_ERROR("strategy of a head does not match");
    }
#endif

    /**
     * If not leaky, then we can learn tangles from H
     */

    if (!leaky) {
        const auto T = tangles, D = dominions;
        if (!heads.empty() and Rs.size() == 1 and *Rs[0] == G) {
            /**
             * A closed region of the top level, while the top level is the whole unsolved game
             * <G>, is a dominion, as the first region in TL. It was attracted in <G> with every
             * tangle of <pl>, so the opponent can only leave it to dominions of <pl> that were
             * already removed from <G> (a vertex with an edge into a dominion of the opponent was
             * attracted to that dominion), and it is already closed under attraction in <G>.
             * So it is solved at once: without extracting its SCCs, which could only be learned
             * as tangles inside it, and without maximizing it again. Every vertex of <pl> in it
             * has its strategy: from its attraction, from its tangle, or (a head) from the first
             * successor in <H> that attracted from it. It is assigned to <Wm> here, and the
             * re-maximization below removes it from the top level.
             * This is not safe below the top level, or after the top level removed vertices
             * presumed won by the opponent (its right side): the opponent may escape to vertices
             * outside the subgame that are only presumed to be won by <pl> or by itself.
             */
#ifndef NDEBUG
            if (trace) {
                logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m (top region)";
                if (trace >= 2) {
                    for (auto v = H.find_last(); v != bitset::npos; v = H.find_prev(v)) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
                logger << std::endl;
            }
            // (check that it is closed in <G> and that <pl> has a strategy in it)
            for (auto v = H.find_first(); v != bitset::npos; v = H.find_next(v)) {
                if (owner(v) == pl) {
                    if (str[v] == -1 or !H[str[v]]) THROW_ERROR("top region dominion without strategy");
                } else {
                    for (auto curedge = outs(v); *curedge != -1; curedge++) {
                        if (G[*curedge] and !H[*curedge]) THROW_ERROR("top region dominion is not closed");
                    }
                }
            }
#endif
            dominions++;
            Wm |= H;
            Wo -= H;
            G -= H;
            // (if <H> was attracted in <G>, see attract, and the top level is <G> again)
            updateDegrees(count_in_G);
            top_is_G = true;
        } else {
            if (!pea_vidx) {
                pea_vidx.reset(new unsigned int[nodecount()]);
                pea_root.resize(nodecount());
                S.resize(nodecount());
                escapes.resize(nodecount());
            }
            std::fill_n(pea_vidx.get(), nodecount(), '\0');
            pea_curidx = 1;
            // Extract this batch against a fixed H and strategy. Maximize its dominions
            // before using the new tangles to re-maximize all regions once.
            for (const int v : heads) extractTangles(v, H);
            if (D != dominions) {
                // maximize the found dominions
                for (auto v = S.find_last(); v != bitset::npos; v = S.find_prev(v)) {
                    Q.push(v);
                }
                attract(pl, S, G);
                // remove them from the game
                Wm |= S;
                Wo -= S;
                G -= S;
                updateDegrees(true);
                S.reset();
            }
        }
        if (D != dominions or T != tangles) {
            // now we need to re-maximize all regions
            for (unsigned int k=0; k<Rs.size(); k++) {
                int k_pr = Ps[k];
                int k_pl = k_pr & 1;
                bitset& k_R = *Rs[k];
                bitset& k_H = *Hs[k];
                bitset k_U;
                // first: R[k] - H[k] >= R[k+1]
                // and compute what remains of k_U
                if (k == 0) {
                    k_R &= G;
                    k_U = G - k_R;
                } else {
                    bitset SG = *Rs[k-1] - *Hs[k-1];
                    k_R &= SG;
                    k_U = SG - k_R;
                }
                // maximize k_U
                for (auto v = k_U.find_last(); v != bitset::npos; v = k_U.find_prev(v)) {
                    Q.push(v);
                }
                attract(1-k_pl, k_U, k_R);
                k_R -= k_U;
                (k_pl == 1 ? W0 : W1) |= k_U;
                (k_pl == 1 ? W1 : W0) -= k_U;
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "\033[1;33mlevel " << k << "-R\033[m";
                    for (auto v = k_R.find_last(); v != bitset::npos; v = k_R.find_prev(v)) {
                        logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
                    }
                    logger << std::endl;
                }
#endif
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "\033[1;33mlevel " << k << "-U\033[m";
                    for (auto v = k_U.find_last(); v != bitset::npos; v = k_U.find_prev(v)) {
                        logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
                    }
                    logger << std::endl;
                }
#endif
                // maximize H
                k_H &= k_R;
                for (auto v = k_H.find_last(); v != bitset::npos; v = k_H.find_prev(v)) {
                    Q.push(v);
                }
                attract(k_pl, k_H, k_R);
                (k_pl == 1 ? W1 : W0) |= k_H;
                (k_pl == 1 ? W0 : W1) -= k_H;
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "\033[1;33mlevel " << k << "-H\033[m";
                    for (auto v = k_H.find_last(); v != bitset::npos; v = k_H.find_prev(v)) {
                        logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
                    }
                    logger << std::endl;
                }
#endif
            }
        }
    }

    heads.clear();

    // Reset Wm for vertices not in SG
    Wm -= SG;

#ifndef NDEBUG
    /**
     * If so desired, report our current status...
     */
    if (trace >= 2) {
        // Format: Subgame: [H] | [Mine in subgame] | [Theirs in subgame] | [Earlier removed from SG]
        logger << "\033[1;33mpre subgame \033[36m" << pr << "\033[m";
        for (int v=vtop; v>=0; v--) {
            if (H[v] && !Wo[v]) logger << " \033[38;5;46m" << label_vertex(v) << "\033[m";
        }
        logger << " |";
        for (int v=vtop; v>=0; v--) {
            if (SG[v] && !H[v]) logger << " \033[38;5;38m" << label_vertex(v) << "\033[m";
        }
        logger << " |";
        for (int v=vtop; v>=0; v--) {
            if (SG[v] && Wo[v]) {
                if (Wm[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m"; // orange -- freshly attracted
                else logger << " \033[38;5;196m" << label_vertex(v) << "\033[m"; // red
            }
        }
        logger << std::endl;
        // Format: Game [player 0] | [player 1] | [no player]
        logger << "\033[1;33mgame";
        for (int v=nodecount()-1; v>=0; v--) {
            if (G[v] && W0[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m";
        }
        logger << " |";
        for (int v=nodecount()-1; v>=0; v--) {
            if (G[v] && W1[v]) logger << " \033[38;5;201m" << label_vertex(v) << "\033[m"; // green
        }
        logger << " |";
        for (int v=nodecount()-1; v>=0; v--) {
            if (G[v] && !W0[v] && !W1[v]) logger << " \033[38;5;38m" << label_vertex(v) << "\033[m"; // blue
        }
        logger << std::endl;
    }
#endif

    /**
     * Go recursive!
     */

    bitset R = SG ^ H;
    if (R.any()) {
        Rs.push_back(&R);
        solve(R, vtop, pe, po);
        Rs.pop_back();
        if (Ps.size() > Rs.size()) Ps.pop_back();
    }
    Hs.pop_back();

    if (SG.none()) return;

    /**
     * Let the opponent attract from our region...
     * The intersection of H and Wo is the opponent's subgame
     */

    bool opponent_attracted_from_us = false;
    // SG only shrinks, so the last word has no subgame vertices above vtop.
    for (int k=0; k<=vtop/64; k++) {
        uint64_t seeds = SG.data()[k] & Wo.data()[k];
        while (seeds) {
            Q.push(k*64 + intrinsics::countr_zero64(seeds));
            seeds &= seeds - 1;
        }
    }
    if (Q.nonempty()) {
        opponent_attracted_from_us = attract(1-pl, Wo, SG, &Wm);
    }

#ifndef NDEBUG
    /**
     * If so desired, report our current status...
     */
    if (trace >= 2) {
        // Format: Subgame: [H] | [Mine in subgame] | [Theirs in subgame] | [Earlier removed from SG]
        logger << "\033[1;33mpost subgame \033[36m" << pr << "\033[m (pe=" << pe << " po=" << po << ")";
        for (int v=vtop; v>=0; v--) {
            if (H[v] && !Wo[v]) logger << " \033[38;5;38m" << label_vertex(v) << "\033[m"; // blueish
        }
        logger << " |";
        for (int v=vtop; v>=0; v--) {
            if (SG[v] && !H[v] && !Wo[v]) logger << " \033[38;5;46m" << label_vertex(v) << "\033[m"; // green
        }
        logger << " |";
        for (int v=vtop; v>=0; v--) {
            if (SG[v] && Wo[v]) {
                if (Wm[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m"; // orange -- freshly attracted
                else logger << " \033[38;5;196m" << label_vertex(v) << "\033[m"; // red
            }
        }
        logger << std::endl;
        // Format: Game [player 0] | [player 1] | [no player]
        logger << "\033[1;33mgame\033[m";
        for (int v=nodecount()-1; v>=0; v--) {
            if (G[v] && W0[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m"; // blueish
        }
        logger << " |";
        for (int v=nodecount()-1; v>=0; v--) {
            if (G[v] && W1[v]) logger << " \033[38;5;201m" << label_vertex(v) << "\033[m"; // green
        }
        logger << std::endl;
    }
#endif

    /**
     * Finally, if the opponent attracted from us, then we do the right side of the universal tree.
     * Also only do this if there is anything remaining in the subgame.
     */

    if (opponent_attracted_from_us) {
        Wm -= SG;  // reset Wm prior to recursion
        SG -= Wo;  // remove opponent won vertices from subgame
        if (Rs.size() == 1) top_is_G = false;

        if (SG.any()) {
            // Expand the right side of the current level of the universal tree
            if (pl == 0) solve(SG, vtop, pe, po/2);
            else solve(SG, vtop, pe/2, po);
        }
    }
}

void
TLQSolver::run()
{
    iterations = 0;

    str = new int[nodecount()];
    Processed.resize(nodecount());
    escs = new unsigned int[nodecount()];
    esc_epoch = new unsigned int[nodecount()]();
    gdeg = new unsigned int[nodecount()];
    // if no vertex is disabled, every successor is in <G>; otherwise count lazily
    if (disabled.none()) {
        for (int v=0; v<nodecount(); v++) gdeg[v] = game.outcount(v);
    } else {
        std::fill_n(gdeg, nodecount(), 0u);
    }

    Q.resize(nodecount());
    W0.resize(nodecount());
    W1.resize(nodecount());
    G = disabled;
    G.flip();

    bitset GG(G);
    Rs.push_back(&GG);
    for (;;) {
        const auto previous_dominions = dominions;
        top_is_G = true;
        solve(GG, nodecount()-1, nodecount(), nodecount());
        if (G.none() or dominions == previous_dominions) break;

        // Removing a dominion can invalidate older regions: their attractor
        // strategies may lead to heads that have now left G. The precision
        // bounds of this traversal need not suffice to revisit those regions.
        // Recheck the remaining game with fresh approximations whenever a
        // traversal removed dominions. Each retry strictly decreases G.
        GG = G;
        W0 -= G;
        W1 -= G;
        Ps.clear();
    }

    // if you want to find any games with vertices that are not in found dominions,
    // then uncomment the following line
    // if (G.any()) THROW_ERROR("found a game that is not solved via tangles");

#ifndef NDEBUG
    if (trace >= 2) {
        for (int v=0; v<nodecount(); v++) {
            if (disabled[v]) continue;
            logger << "vertex " << label_vertex(v) << " is solved by";
            if (W0[v]) {
                logger << " even";
                if (owner(v) == 0) {
                    logger << " (";
                    if (str[v] == -1) logger << "-1";
                    else logger << label_vertex(str[v]);
                    logger << ")";
                }
            }
            if (W1[v]) {
                logger << " odd";
                if (owner(v) == 1) {
                    logger << " (";
                    if (str[v] == -1) logger << "-1";
                    else logger << label_vertex(str[v]);
                    logger << ")";
                }
            }
            logger << std::endl;
        }
    }
#endif

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (W0[v]) Solver::solve(v, 0, str[v]);
        else if (W1[v]) Solver::solve(v, 1, str[v]);
        else THROW_ERROR("unsolved vertex?");
    }

    logger << "found " << dominions << " dominions and "<< tangles << " tangles." << std::endl;
    logger << "solved with " << iterations << " iterations." << std::endl;

    delete[] str;
    delete[] escs;
    delete[] esc_epoch;
    delete[] gdeg;
}

}
