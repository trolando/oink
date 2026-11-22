/*
 * Copyright 2022 Tom van Dijk
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
#include <cstring>
#include <memory>
#include <type_traits>

#include "pmtl.hpp"
#include "oink/uintqueue.hpp"

// #define CHECK_UNIQUE 1

namespace pg {


PMTLSolver::PMTLSolver(Oink& oink, Game& game, MeasureKind kind, bool enable_tangles)
    : Solver(oink, game), measure_kind(kind), use_tangles(enable_tangles)
{
}

PMTLSolver::~PMTLSolver()
{
}

/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R> from subgame <G>.
 * If <max_prio> >= 0, then only attract vertices with priority 0 <= pr <= max_prio.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
PMTLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (Z[from]) {
            // already in Z, maybe set strategy (for vertices in the original target set)
            if (owner(from) == pl and str[from] == -1) {
                str[from] = v;
            }
        } else if (R[from] and (max_prio < 0 or (priority(from) & 1) == pl or priority(from) <= max_prio)) {
            if (owner(from) != pl) {
                // check if opponent can escape
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (G[to] and !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            Zlist.push_back(from);
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
 * All vertices in tangle <t> must be in <Z+R> and the opponent may not escape to subgame <G\Z>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
bool
PMTLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    /**
     * Check if tangle is won by player <pl> and not deleted.
     */
    {
        const int tangle_pr = tpr[t];
        if (tangle_pr == -1) return false; // deleted tangle
        if (pl != -1 and pl != (tangle_pr&1)) return false; // not of desired parity
    }

    /**
     * Check if tangle is contained in Z+R and if any vertices are not already in Z
     * If no new vertices could be attracted, then don't attract...
     */
    {
        bool can_attract_new = false;
        int *ptr = tv[t];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            ptr++;
            if (!this->G[v]) {
                // on-the-fly detect out-of-game tangles
                tpr[t] = -1; // delete the tangle
                return false; // is now a deleted tangle
            } else if (Z[v]) {
                continue; // already attracted
            } else if (!R[v]) {
                return false; // not contained in Z+R
            } else if (!G[v]) {
                // the top of a higher region, which update() removes from Z again (the special
                // case), so the other vertices of the tangle would play out of the region
                return false;
            } else if (max_prio >= 0 and priority(v) > max_prio and ((priority(v) & 1) != pl)) {
                return false;
            } else {
                can_attract_new = true; // has vertices not yet attracted
            }
        }
        if (!can_attract_new) return false;
    }

    /**
     * Check if the tangle can escape to G\Z.
     */
    {
        int v, *ptr = tout[t];
        while ((v=*ptr++) != -1) {
            if (Z[v]) continue;
            if (G[v]) return false; // opponent escapes
        }
    }

    /**
     * Attract!
     */
    {
        int *ptr = tv[t];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            const int s = *ptr++;
            if (Z[v]) continue; // already in <Z>
            Z[v] = true;
            Zlist.push_back(v);
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
 * Current subgame is <G>.
 * Only tangles that are contained in <R>+<Z> and vertices of maximum priority max_prio (if >= 0)
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
inline void
PMTLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    if (has_tin.empty() or !has_tin[v]) return;
    const Backlinks &links = tin_pages[v/TIN_PAGE_SIZE][v%TIN_PAGE_SIZE];
    const int *ts = &tin_pool[links.offset];
    for (unsigned int i = 0, n = links.size; i < n; i++) attractTangle(ts[i], pl, R, Z, G, max_prio);
}


/**
 * Add the new tangle <tidx> to the backlinks of its escapes in <tangleto>, after the
 * tangles found before it.
 */
void
PMTLSolver::addBacklinks(const int tidx)
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
 * Allocate the workspace of extractTangles when it is first used.
 */
void
PMTLSolver::prepareExtraction()
{
    if (pea_vidx) return;
    pea_vidx.reset(new unsigned int[nodecount()]);
    pea_root.resize(nodecount());
    escapes.resize(nodecount());
}

/**
 * Compute SCCs in subgraph induced by <R> and <str>.
 * Start the SCC computation at vertex <startvertex>.
 * Every SCC is then processed as a tangle.
 * If the tangle is closed, it is a dominion and added to <S> and <Q>.
 */
bool
PMTLSolver::extractTangles(int startvertex, bitset &R, int *str)
{
    bool new_tangles = false;
    const int pr = priority(startvertex);
    const int pl = pr&1;

    /**
     * The following is the nonrecursive implementation of David Pearce,
     * "A space-efficient algorithm for finding strongly connected components" (IPL, 2016),
     * modified to on-the-fly restrict the graph by <R> and <str>.
     */

    // beginVisiting
    pea_vS.push_back(startvertex);
    pea_iS.push_back(0);
    pea_root[startvertex] = true;
    pea_vidx[startvertex] = pea_curidx++;
    while (!pea_vS.empty()) {
pearce_again:
        // visitLoop
        const unsigned int n = pea_vS.back();
        unsigned int i = pea_iS.back();

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
                        pea_iS.back() = i+1;
                        // beginVisiting
                        pea_vS.push_back(to);
                        pea_iS.push_back(0);
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
#ifndef NDEBUG
            if (str[n] == -1 or !R[str[n]]) LOGIC_ERROR; // our strategy must stay in the region
#endif
            const int s = str[n];
            if (i == 0) {
                // beginEdge
                if (pea_vidx[s] == 0) {
                    pea_iS.back() = 1;
                    // beginVisiting
                    pea_vS.push_back(s);
                    pea_iS.push_back(0);
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
        pea_vS.pop_back();
        pea_iS.pop_back();
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

        const bool is_tangle = (tangle.size() > 1) or
            ((unsigned int)str[n] == n) or
            (str[n] == -1 and game.has_edge(n, n));
        if (!is_tangle) {
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
                    if (G[to] and !escapes[to]) {
                        escapes[to] = true;
                        tangleto.push_back(to);
                    }
                }
            }
        }

        if (tangle.size()+tangleto.size() < escapes.size()/64) {
            for (const int v : tangle) escapes[v] = false;
            for (const int v : tangleto) escapes[v] = false;
        } else {
            escapes.reset();
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
            for (const int v : tangle) {
                // if not yet added to solve queue, mark and add it
                if (S[v] == false) {
                    S[v] = true;
                    dom_vector.push_back(v);
                    dom_vector.push_back(str[v]);
                }
            }
            dominions++;
            tangle.clear();
            new_tangles = true;
            continue;
        }

        /**
         * We're not a dominion, we're a tangle.
         */

#if CHECK_UNIQUE
        /**
         * First we sort.
         */
        std::sort(tangle.begin(), tangle.end());
        int tsize = tangle.size();

        /**
         * Then we check if it already exists.
         */
        bool already_exists = false;
        for (auto &_tv : tv) {
            for (int i=0; ; i++) {
                if (i == tsize) {
                    if (_tv[2*i] == -1) already_exists = true;
                    break;
                } else if (_tv[2*i] != tangle[i]) {
                    break;
                }
            }
            if (already_exists) break;
        }
        if (already_exists) {
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[1;38;5;198mduplicate tangle\033[m";
                if (trace >= 2) {
                    for (int &v : tangle) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
                logger << std::endl;
            }
#endif
            tangle.clear();
            tangleto.clear();
            continue;
        }
#endif

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

        // new tangle idx
        const int tidx = tpr.size();

        // add back links to all normal vertices in our [out]
        addBacklinks(tidx);

        // move tangleto into vout
        int* _tout = new int[tangleto.size()+1];
        std::copy(tangleto.begin(), tangleto.end(), _tout);
        _tout[tangleto.size()] = -1;
        tout.push_back(_tout);

        // move tangle into vv
        int* _vv = new int[tangle.size()*2+1], c=0;
        for (const int v : tangle) {
            _vv[c++] = v;
            _vv[c++] = str[v];
        }
        _vv[c] = -1;
        tv.push_back(_vv);

        // and set p to pr
        tpr.push_back(pr);
        if ((int)tangles_of_pr.size() <= pr) tangles_of_pr.resize(pr+1);
        tangles_of_pr[pr].push_back(tidx);

        tangles++;
        new_tangles = true;
        tangle.clear();
        tangleto.clear();
    }

    pea_S.clear();
    return new_tangles;
}


template<class Measure, bool GoUp, bool UseTangles>
bool
PMTLSolver::updateMeasure(Measure &pm, Measure &target_pm, const int player)
{
    // sort function to sort by value (desc) then priority (desc)
    auto max_pm_first = [&](int const & a, int const & b) {
            int pmc = pm.compare(a, b);
            // return pmc > 0;
            if (pmc > 0) return true;
            if (pmc < 0) return false;
            int pr_a = priority(a);
            int pr_b = priority(b);
            if ((pr_a&1) != player) pr_a = 0-pr_a;
            if ((pr_b&1) != player) pr_b = 0-pr_b;
            if (pr_a > pr_b) return true;
            if (pr_a < pr_b) return false;
            return a > b;
        };

    // order the unsolved vertices by value then priority. Only the vertices whose measure
    // changed since the last update of <player> are sorted and then merged into the order of
    // that update (the order is total, so this is the same as sorting all vertices)
    auto &order = this->order[player];
    if (!sorted[player]) {
        order.clear();
        for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) order.push_back(v);
        std::sort(order.begin(), order.end(), max_pm_first);
        sorted[player] = true;
    } else {
        moved.clear();
        size_t j = 0;
        for (const int v : order) {
            if (!G[v]) continue; // solved
            if (Changed[player][v]) moved.push_back(v);
            else order[j++] = v;
        }
        order.resize(j);
        if (!moved.empty()) {
            std::sort(moved.begin(), moved.end(), max_pm_first);
            merged.resize(order.size() + moved.size());
            std::merge(order.begin(), order.end(), moved.begin(), moved.end(), merged.begin(), max_pm_first);
            order.swap(merged);
        }
    }
    Changed[player].reset();
    const int norder = (int)order.size();

#ifndef NDEBUG
    if (trace >= 2) {
        // echo current order
        for (int x=0; x<norder; x++) {
            auto v = order[x];
            if (!G[v]) continue;
            logger << "vertex " << label_vertex(v) << ": ";
            pm.stream(logger, v);
            logger << std::endl;
        }
    }
#endif

    R = G; // initialize remainder (including unattracted top vertices) with G
    WasHead.reset();
    int scanned = 0; // the ties up to here were scanned for heads
    E = G; // initialize escapes (for opponent) with G
    Updated.reset();
    Deferred.reset();
    bool any_updated = false;

    // go from "best" to "worst" vertex
    for (int x=0; x<norder; x++) {
        auto top = order[x];
        if (!R[top]) continue; // not in remainder

        // Source measures and their bounds remain fixed during this update.
        // Ordered::see uses -3 internally; -4 is reserved here for this cache.
        int cached_priority = -1;
        auto see_top = [&](int p) {
            if (p != cached_priority) {
                pm.copy(top, -4);
                pm.see(-4, p);
                cached_priority = p;
            }
        };

        // found a new top vertex
        bool known_tangle = false;
        if (UseTangles) {
            // first check if <top> happens to be a the top of a tangle that we already have
            // if so, this "initializes" Z with that tangle
            // (only the tangles of this priority, in order of creation; deleted ones have tpr -1)
            const int pr = priority(top);
            if (pr < (int)tangles_of_pr.size()) {
                const auto& ts = tangles_of_pr[pr];
                for (size_t j=0; j<ts.size(); j++) {
                    const int i = ts[j];
                    if (tpr[i] != pr) continue;
                    if (attractTangle(i, player, R, Z, E, pr)) {
                        if (Z[top]) {
                            known_tangle = true;
                            break;
                        }
                    }
                }
            }
        }

        // attract with max prio pr(top), to <top> and to the next vertices in the order with
        // the same measure and priority (one region with all these heads, as in TL, rather
        // than a region per vertex)
        int max_pr = priority(top);
        heads.clear();
        heads.push_back(top);
        // A known tangle already added the top to <Z> and <Q>, with its strategy. Adding it
        // again would process it as the top of a higher region on its second visit.
        if (!known_tangle) {
            Z[top] = true; // add to <Z>
            Zlist.push_back(top);
            str[top] = -1; // don't care, but needed for correct trace reporting
            Q.push(top);
        }
        if (!Updated[top] and x >= scanned) {
            // (each run of tied vertices is scanned once: a vertex skipped here cannot become
            // a head later in this round, as it stays out of R, Updated, or was a head)
            int y = x+1;
            for (; y < norder; y++) {
                const int w = order[y];
                if (priority(w) != max_pr or !pm.eq(w, top)) break;
                if (!R[w] or Z[w] or Updated[w] or WasHead[w]) continue;
                WasHead[w] = true; // (an open head is later processed as a single top)
                heads.push_back(w);
                Z[w] = true;
                Zlist.push_back(w);
                str[w] = -1;
                Q.push(w);
            }
            scanned = y;
        }

#ifndef NDEBUG
        if (trace) {
            logger << "attracting to " << label_vertex(top) << ", max priority " << max_pr << std::endl;
        }
#endif
        while (Q.nonempty()) {
            const int v = Q.pop();
            if (E[v]) {
                E[v] = false; // remove from E
                attractVertices(player, v, R, Z, E, max_pr);
                if (UseTangles) {
                    attractTangles(player, v, R, Z, E, max_pr);
                }
            } else {
                // oops, we have a special case (aka the top vertex of a higher region)
                Z[v] = false; // remove from Z again
                R[v] = false; // and now remove from R
#ifndef NDEBUG
                if (trace) {
                    logger << "special case " << label_vertex(v) << " with U==" << Updated[v] << std::endl;
                }
#endif
                if (!Updated[v]) {
                    // in this case we might lift v too high if our top vertex was already lifted to higher measure in a higher region.
                    // An example of when this occurs is if player(v) != player, there is an edge from v to itself and priority(v) & 1
                    // != player, i.e. v has a winning self loop for !player.
                    if (pm.compare(top, v) <= 0) {
                        see_top(priority(v));
                        if (pm.compare(-4, v) > 0) {
                            Updated[v] = true;
                            any_updated = true;
                            //pm.copy(-4, v);
                            target_pm.copy_from(pm, -4, v);
#ifndef NDEBUG
                            if (trace >= 2) {
                                logger << "updating vertex " << label_vertex(v) << " to ";
                                target_pm.stream(logger, v);
                                logger << std::endl;
                            }
#endif
                            lifts++;
                        }
                    }
                }
            }
        }

        // a head is open if it can escape the region (or has no strategy into it); open heads
        // are not lifted and stay in R, and the region is closed if no head is open
        bool closed_region = true;
        for (const int h : heads) {
            bool open = false;
            if (owner(h) == player) {
                open = str[h] == -1;
            } else {
                for (auto curedge = outs(h); *curedge != -1; curedge++) {
                    if (E[*curedge]) {
                        open = true;
                        break;
                    }
                }
            }
            if (open) {
                OpenHead[h] = true;
                closed_region = false;
            }
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m ";
            logger << "\033[1;36m" << priority(top) << "\033[m";
            for (auto v = Z.find_last(); v != bitset::npos; v = Z.find_prev(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif

        for (const int v : Zlist) if (Z[v]) R[v] = false; // R -= Z
        // put open heads back into R: <top> as before (it can still be attracted to a later
        // region), the other heads as if they were never heads, so they get a region later
        for (const int h : heads) {
            if (!OpenHead[h]) continue;
            R[h] = true;
            if (h != top) E[h] = true;
        }

        // If the measure of top was already increased in a higher we might not be able to immediately use that measure
        // so we skip the region. (This is also why there is the if statement in the special case above) Another possible
        // solution is to keep track of a separate new set of measures pm2 and only replace pm with pm2 at the end. That
        // is likely faster.
        if (Updated[top]) {
            // the lifts of this region are postponed to the next round, so shortcuts() must
            // treat its vertices as changed (like lifted vertices)
            for (const int v : Zlist) if (Z[v]) Deferred[v] = true;
            for (const int v : Zlist) Z[v] = false; // Z.reset()
            Zlist.clear();
            for (const int h : heads) OpenHead[h] = false;
            continue;
        }

        // update all vertices with prog
        for (const int v : Zlist) {
            if (!Z[v]) continue;
            if (OpenHead[v]) continue;
            if (Updated[v]) continue; // already updated to a higher pm

            see_top(priority(v));
            if (!pm.eq(-4, v)) {
                Updated[v] = true;
                any_updated = true;
                //pm.copy(-4, v);
                target_pm.copy_from(pm, -4, v);
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "updating vertex " << label_vertex(v) << " to ";
                    target_pm.stream(logger, v);
                    logger << std::endl;
                }
#endif
                lifts++;
            }
        }

        {
            // Closed-region dominions do not require tangle learning.
            // As in TL: a closed region from which the opponent cannot escape in the game is a
            // dominion, without decomposing it into tangles. Its heads have our parity and a
            // higher priority than any of the opponent in it, and the attractor strategy either
            // stays in a tangle of ours or reaches a head again, so <player> wins every play.
            bool dominion = false;
            if (closed_region && (priority(top)&1) == player) {
                dominion = true;
                for (const int v : Zlist) {
                    if (!Z[v]) continue;
                    if (owner(v) == player) {
                        // (our strategy must stay in the region too: a vertex attracted via a
                        // vertex that was then removed from Z again, as the top of a higher
                        // region (the special case above), or via a tangle containing such a
                        // vertex, still plays to it)
                        if (str[v] == -1 or !Z[str[v]]) {
                            dominion = false;
                            break;
                        }
                        continue;
                    }
                    for (auto curedge = outs(v); *curedge != -1; curedge++) {
                        const int to = *curedge;
                        if (G[to] and !Z[to]) {
                            dominion = false;
                            break;
                        }
                    }
                    if (!dominion) break;
                }
                if (dominion) {
#ifndef NDEBUG
                    if (trace) logger << "\033[1;38;5;201mdominion (closed region) \033[36m" << priority(top) << "\033[m" << std::endl;
#endif
                    for (const int v : Zlist) {
                        if (Z[v] and !S[v]) {
                            S[v] = true;
                            dom_vector.push_back(v);
                            dom_vector.push_back(str[v]);
                        }
                    }
                    dominions++;
                }
            }
            // otherwise learn tangles (could be duplicate!) if the region is closed, but only if
            // correct priority; with open heads, in the part of the region that does not depend on them
            if (UseTangles && !dominion && !known_tangle && (priority(top)&1) == player) {
                bool any_closed_head = false;
                for (const int h : heads) if (!OpenHead[h]) { any_closed_head = true; break; }
                if (any_closed_head) {
                    bitset *TR = &Z;
                    if (!closed_region) {
                        // TZ := Z without the open heads and the vertices that depend on them
                        // (opponent vertices with an edge to them, and our vertices that play to them)
                        if (TZ.empty()) TZ.resize(nodecount());
                        for (const int v : Zlist) if (Z[v]) TZ[v] = true;
                        for (const int h : heads) {
                            if (OpenHead[h]) {
                                TZ[h] = false;
                                Q.push(h);
                            }
                        }
                        while (Q.nonempty()) {
                            const int v = Q.pop();
                            for (auto curedge = ins(v); *curedge != -1; curedge++) {
                                const int from = *curedge;
                                if (TZ[from] and (str[from] == -1 or str[from] == v)) {
                                    TZ[from] = false;
                                    Q.push(from);
                                }
                            }
                        }
                        TR = &TZ;
                    }
                    prepareExtraction();
                    for (const int v : Zlist) pea_vidx[v] = 0; // extractTangles only visits <Z>
                    pea_curidx=1;
                    for (const int h : heads) {
                        // (not if open, or already visited from another head)
                        if ((*TR)[h] and pea_vidx[h] == 0) extractTangles(h, *TR, str);
                    }
                    if (!closed_region) for (const int v : Zlist) TZ[v] = false;
                }
            }
        }

        if (GoUp) {
            // attract with increasing priority
            const int d = priority(nodecount()-1);
            for (max_pr=priority(top)+1; max_pr<=d; max_pr++) {
                for (const int v : Zlist) if (Z[v]) Q.push(v);

#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "attracting with max priority " << max_pr << std::endl;
                }
#endif
                while (Q.nonempty()) {
                    const int v = Q.pop();

                    attractVertices(player, v, R, Z, E, max_pr);
                    if (UseTangles) {
                        attractTangles(player, v, R, Z, E, max_pr);
                    }

                    if (!Updated[v] && priority(v) >= max_pr) {
                        see_top(priority(v));
                        if (pm.compare(-4, v) > 0) {
                            Updated[v] = true;
                            any_updated = true;
                            //pm.copy(-4, v);
                            target_pm.copy_from(pm, -4, v);
#ifndef NDEBUG
                            if (trace >= 2) {
                                logger << "updating vertex " << label_vertex(v) << " to ";
                                target_pm.stream(logger, v);
                                logger << std::endl;
                            }
#endif
                            lifts++;
                        }
                    }
                }
            }
        } else { // !GoUp
            // we now update all vertices that are attracted to any vertex in Z...
            for (const int v : Zlist) {
                if (!Z[v]) continue;
                for (auto curedge = ins(v); *curedge != -1; curedge++) {
                    int from = *curedge;
                    if (!R[from]) continue; // but only if still in R

#ifndef NDEBUG
                    if (trace >= 3) {
                        logger << "Attempting expansion where " << label_vertex(from) << " with owner " << owner(from) << " plays to " << label_vertex(v) << std::endl;
                    }
#endif
                    if (owner(from) == player) {
                        if (!Updated[from]) {
                            see_top(priority(from));
                            if (pm.compare(-4, from) > 0) {
                                Updated[from] = true;
                                any_updated = true;
                                // pm.copy(-4, from);
                                target_pm.copy_from(pm, -4, from);
#ifndef NDEBUG
                                if (trace >= 2) {
                                    logger << "updating vertex " << label_vertex(from) << " to ";
                                    target_pm.stream(logger, from);
                                    logger << std::endl;
                                }
#endif
                                lifts++;
                            }
                        }
                    } else {
                        bool attracted = true;
                        for (auto edge = outs(from); *edge != -1; edge++) {
                            int to = *edge;
                            if (E[to] && !Z[to]) {
#ifndef NDEBUG
                                if (trace >= 3) {
                                    logger << "Found escape " << label_vertex(to) << std::endl;
                                }
#endif
                                attracted = false;
                                break;
                            }
                        }
                        if (attracted &&!Updated[from]) {
                            see_top(priority(from));
                            if (pm.compare(-4, from) > 0) {
                                Updated[from] = true;
                                any_updated = true;
                                // pm.copy(-4, from);
                                target_pm.copy_from(pm, -4, from);
#ifndef NDEBUG
                                if (trace >= 2) {
                                    logger << "updating vertex " << label_vertex(from) << " to ";
                                    target_pm.stream(logger, from);
                                    logger << std::endl;
                                }
#endif
                                lifts++;
                            }
                        }
                    }
                }
            }
        }

        for (const int v : Zlist) Z[v] = false; // Z.reset()
        Zlist.clear();
        for (const int h : heads) OpenHead[h] = false;

        if (any_updated) {
#ifndef NDEBUG
            if (trace >= 3) {
                logger << "Updated:";
                for (auto v=Updated.find_last(); v!=bitset::npos; v=Updated.find_prev(v)) {
                    logger << " " << label_vertex(v);
                }
                logger << std::endl;
            }
#endif
            // break; // only run once...
        }
    }

    /**
     * Extend any dominions that were found.
     * (Any solved vertices are now in <S> and <dom_vector>.)
     */

    {
        if (!dom_vector.empty()) {
            for (unsigned i = 0; i<dom_vector.size(); i+=2) {
                int v = dom_vector[i];
                int s = dom_vector[i+1];
                str[v] = s;
                Q.push(v);
            }
            dom_vector.clear();

            // disable
            while (Q.nonempty()) {
                const int v = Q.pop();
                Solver::solve(v, player, str[v]);
                attractVertices(player, v, G, S, G, -1);
                if (UseTangles) attractTangles(player, v, G, S, G, -1);
            }

            // the solved vertices leave G, which can take away the only escape of an opponent
            // vertex, so shortcuts() must treat them as changed (like lifted vertices)
            Updated |= S;
            G -= S; // and remove from G
            S.reset();
        }
    }

    if (UseTangles and !GoUp) {
        // A tangle head can follow its strategy until the opponent escapes, or win inside
        // the tangle. Its priority dominates every opponent priority on such a path, so
        // seeing the head priority at the least escape measure is a lower bound for its
        // lift. Aura only expands ordinary predecessors; without this tangle lift, a head
        // can keep counting through internal cycles instead of reaching the escape measure.
        for (int t = 0; t < (int)tpr.size(); t++) {
            if (tpr[t] == -1 or (tpr[t]&1) != player) continue;
            bool active = true;
            for (int *ptr = tv[t]; *ptr != -1; ptr += 2) {
                if (!G[*ptr]) {
                    active = false;
                    break;
                }
            }
            if (!active) {
                tpr[t] = -1;
                continue;
            }
            pm.copy_top(-2);
            for (int *ptr = tout[t]; *ptr != -1; ptr++) {
                if (G[*ptr] and pm.compare(*ptr, -2) < 0) pm.copy(*ptr, -2);
            }
            for (int *ptr = tv[t]; *ptr != -1; ptr += 2) {
                const int v = *ptr;
                if (priority(v) < tpr[t] or (priority(v)&1) != player) continue;
                pm.copy(-2, -1);
                pm.see(-1, priority(v));
                target_pm.copy_from(pm, -1, -1);
                if (target_pm.compare(-1, v) > 0) {
                    target_pm.copy_from(pm, -1, v);
                    Updated[v] = true;
                    any_updated = true;
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "updating tangle " << t << " head " << label_vertex(v) << " to ";
                        target_pm.stream(logger, v);
                        logger << std::endl;
                    }
#endif
                    lifts++;
                }
            }
        }
    }

    Changed[player] |= Updated; // (the lifted vertices)
    return any_updated;
}


template<bool GoUp, bool UseTangles>
bool
PMTLSolver::update(Measures &pm, Measures &target_pm, const int player)
{
    // (both tables of a player have the same concrete class, see measure_code_size)
    return with_measure(measure_kind, pm, [&](auto &concrete) {
        using Measure = std::remove_reference_t<decltype(concrete)>;
        return updateMeasure<Measure, GoUp, UseTangles>(concrete, static_cast<Measure&>(target_pm), player);
    });
}


void
PMTLSolver::shortcuts(const int player, Measures &pm, Measures &target1, Measures &target2)
{
    // Source measures and G stay fixed while vertices and tangles grow L.
    // Start a fresh epoch because measures can change between calls.
    if (++escape_epoch == 0) {
        std::fill(escape_stamp.begin(), escape_stamp.end(), 0);
        escape_epoch = 1;
    }
    auto& L = Updated;

    // we want to find all potentially liftable vertices for <player>
    // add all tops
    // L&=G;
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (pm.is_top(v)) L[v] = true;
    }
    // and the vertices of the regions that update() skipped: they may be liftable already
    L |= Deferred;
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        //logger << "pushing " << label_vertex(v) << std::endl;
        Q.push(v);
    }
    while (Q.nonempty()) {
        const int v = Q.pop();
        // update attracts tangles as well as vertices: a tangle of <player> with an escape to
        // <v> can be attracted into a higher region once <v> is lifted, which lifts all its
        // vertices, so they are potentially liftable too (without this, a vertex of such a
        // tangle could be given to the opponent and then lifted to top by <player>, leaving it
        // unsolved)
        if (use_tangles and !has_tin.empty() and has_tin[v]) {
            const Backlinks &links = tin_pages[v/TIN_PAGE_SIZE][v%TIN_PAGE_SIZE];
            for (unsigned int i = 0, n = links.size; i < n; i++) {
                const int t = tin_pool[links.offset+i];
                if (tpr[t] == -1 or (tpr[t]&1) != player) continue;
                for (int *ptr = tv[t]; *ptr != -1; ptr += 2) {
                    const int u = *ptr;
                    if (L[u] || !G[u]) continue;
                    L[u] = true;
                    Q.push(u);
                }
            }
        }
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (L[from] || !G[from]) continue; // don't attract
            if (owner(from) != player) {
                int first = 0;
                if (escape_stamp[from] == escape_epoch) {
                    const int previous = escape_offset[from];
                    if (!L[outs(from)[previous]]) continue;
                    // Earlier edges cannot become escapes while L grows.
                    first = previous + 1;
                }
                bool escapes = false;
                for (auto curedge = outs(from) + first; *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!G[to] || L[to]) continue; // not a potential escape
                    // it is only a valid escape if the measure obtained by going from
                    // <from> to <to> actually is == the current measure...
                    pm.copy(to, -1);
                    pm.see(-1, priority(from));
                    if (pm.compare(-1, from) <= 0) {
                        escape_offset[from] = int(curedge - outs(from));
                        escape_stamp[from] = escape_epoch;
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            L[from] = true;
            Q.push(from);
        }

    }
    L ^= G;
    L &= G; // only the unsolved vertices that are not potentially liftable
    for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
        // logger << "\033[1;36mplayer " << (1-player) << " should win " << label_vertex(v) << "\033[m" << std::endl;
        target1.copy_top(v);
        target2.copy_top(v);
        Changed[1-player][v] = true;
    }
}



/**
 * Given progress measure of player X, win for player ~X
 */
void
PMTLSolver::solve(Measures &pm, const int player)
{
    // extract strategy
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (G[v] && !pm.is_top(v)) {
            // won by ~X
            str[v] = -1;
            if (owner(v) == player) {
                // find lowest prog...
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (G[to] && !pm.is_top(to)) {
                        if (str[v] == -1) {
                            pm.copy(to, -1);
                            str[v] = to;
                        } else if (pm.compare(to, -1) < 0) {
                            pm.copy(to, -1);
                            str[v] = to;
                        }
                    }
                }

                if (str[v] == -1) LOGIC_ERROR;
            }
        }
    }
    // mark as solved
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (G[v] && !pm.is_top(v)) {
            Solver::solve(v, player, str[v]);
            G[v] = false;
        }
    }
    Solver::flush();
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (disabled[v]) {
            G[v] = false;
        }
    }
}

bool
PMTLSolver::parseOptions(std::string& opts) {
    if (opts.empty()) return true; // keep the defaults (or the measure kind of pmtl-small etc.)
    std::string remaining = opts;

    //go_up = false;
    //use_tangles = false;
    //measure_kind = MeasureKind::Small;

    do {
        size_t comma_index = remaining.find(',');
        std::string opt = remaining.substr(0, comma_index);
        if (opt == "up" || opt == "UP") {
            go_up = true;
        } else if (opt == "noup" || opt == "aura" || opt == "NOUP" || opt == "AURA") {
            go_up = false;
        } else if (opt == "tangle" || opt == "tangles" || opt == "TANGLE" || opt == "TANGLES") {
            use_tangles = true;
        } else if (opt == "notangle" || opt == "notangles" || opt == "NOTANGLE" || opt == "NOTANGLES") {
            use_tangles = false;
        } else if (opt == "noshortcuts" || opt == "no-shortcuts" || opt == "NOSHORTCUTS" || opt == "NO-SHORTCUTS") {
            use_shortcuts = false;
        } else {
            std::optional<MeasureKind> parsed = parse_measure_kind(opt);
            if (parsed.has_value()) {
                measure_kind = *parsed;
            } else {
                logger << "unrecognised option for PMTL solver, possible options: ";
                stream_measure_kinds(logger);
                logger << ", up, noup, tangles, notangles (separated by commas and no spaces)" << std::endl;
                return false;
            }
        }
        if (comma_index == std::string::npos) {
            remaining = remaining.substr(remaining.size());
        } else {
            remaining = remaining.substr(comma_index + 1);
        }
    } while (!remaining.empty());

#ifndef NDEBUG
    if (trace) {
        logger << "selected ";
        stream_measure_kind(logger, measure_kind);
        logger << " for PMTL" << std::endl;
        if (go_up) {
            logger << "enabled attracting up for PMTL" << std::endl;
        } else {
            logger << "disabled attracting up for PMTL" << std::endl;
        }
        if (use_tangles) {
            logger << "enabled tangle attractors for PMTL" << std::endl;
        } else {
            logger << "disabled tangle attractors for PMTL" << std::endl;
        }
        if (!use_shortcuts) {
            logger << "disabled shortcuts for PMTL" << std::endl;
        }
    }
#endif

    return true;
}

/**
 * Make the vertex order of player <pl> for its first update.
 * All measures start at bottom, so the priority order alone defines the order: first the
 * vertices of priority parity <pl> from high to low priority, then the others from low to
 * high priority. Shortcuts mark any changes before the first update of <pl>.
 */
void
PMTLSolver::initOrder(const int pl)
{
    order[pl].clear();
    order[pl].reserve(G.count());
    Changed[pl].resize(nodecount());
    for (int v=nodecount()-1; v>=0; v--) {
        if (G[v] and (priority(v)&1) == pl) order[pl].push_back(v);
    }
    for (int first=0; first<nodecount();) {
        int last=first+1;
        while (last<nodecount() and priority(last) == priority(first)) last++;
        if ((priority(first)&1) != pl) {
            for (int v=last-1; v>=first; v--) if (G[v]) order[pl].push_back(v);
        }
        first=last;
    }
    sorted[pl] = true;
}

void
PMTLSolver::run()
{
    str = new int[nodecount()];
    if (use_shortcuts) {
        escape_offset.resize(nodecount());
        escape_stamp.resize(nodecount());
    }
    G = disabled;
    G.flip();
    S.resize(nodecount());
    // (the tangle storage and the workspace of tangle learning are allocated when used)
    Q.resize(nodecount());

    // Even updates first. The order and the measures of Odd are only made once an update of
    // Even leaves unsolved vertices, before the shortcuts of Even can write Odd's measures.
    initOrder(0);
    std::unique_ptr<Measures> pm0(new_measure(measure_kind, game, 0));
    std::unique_ptr<Measures> pm0b(new_measure(measure_kind, game, 0));
    std::unique_ptr<Measures> pm1, pm1b;

    G = disabled;
    G.flip();

    Z.resize(nodecount());
    R.resize(nodecount());
    E.resize(nodecount());
    Updated = disabled;
    Deferred.resize(nodecount());
    OpenHead.resize(nodecount());
    WasHead.resize(nodecount());
    bool c0 = true;
    bool c1 = true;

    std::uint64_t evenlifts = 0;
    std::uint64_t oddlifts = 0;

    while (true) {
#ifndef NDEBUG
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "\033[m" << std::endl;
#endif
        iterations++;

        if (c0) {
#ifndef NDEBUG
            if (trace) logger << "\033[1;38;5;196mfor player Even\033[m" << std::endl;
#endif
            const auto previous_lifts = lifts;
            if (go_up) {
                if (use_tangles) {
                    c0 = update<true, true>(*pm0, *pm0b, 0);
                } else  {
                    c0 = update<true, false>(*pm0, *pm0b, 0);
                }
            } else if (use_tangles) {
                c0 = update<false, true>(*pm0, *pm0b, 0);
            } else {
                c0 = update<false, false>(*pm0, *pm0b, 0);
            }
            evenlifts += lifts - previous_lifts;
            // All remaining vertices may have been solved by dominion attraction.
            if (G.none()) break;
            if (!pm1) {
                initOrder(1);
                pm1.reset(new_measure(measure_kind, game, 1));
                pm1b.reset(new_measure(measure_kind, game, 1));
            }
            pm0->set(*pm0b);
            if (use_shortcuts) shortcuts(0, *pm0, *pm1, *pm1b);
            if (!c0) {
                solve(*pm0, 1);
                logger << "done with Even at iteration " << iterations << std::endl;
                if (G.none()) break;
            }
        }

        if (c1) {
#ifndef NDEBUG
            if (trace) logger << "\033[1;38;5;196mfor player Odd\033[m" << std::endl;
#endif
            const auto previous_lifts = lifts;
            if (go_up) {
                if (use_tangles) {
                    c1 = update<true, true>(*pm1, *pm1b, 1);
                } else  {
                    c1 = update<true, false>(*pm1, *pm1b, 1);
                }
            } else if (use_tangles) {
                c1 = update<false, true>(*pm1, *pm1b, 1);
            } else {
                c1 = update<false, false>(*pm1, *pm1b, 1);
            }
            oddlifts += lifts - previous_lifts;
            // All remaining vertices may have been solved by dominion attraction.
            if (G.none()) break;
            pm1->set(*pm1b);
            if (use_shortcuts) shortcuts(1, *pm1, *pm0, *pm0b);
            if (!c1) {
                solve(*pm1, 0);
                logger << "done with Odd at iteration " << iterations << std::endl;
                if (G.none()) break;
            }
        }

        if (!c0 && !c1) break;
    }

#ifndef NDEBUG
    // Check if the whole game is now solved
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { THROW_ERROR("search was incomplete!"); }
    }
#endif

    logger << "found " << dominions << " dominions and " << tangles << " tangles." << std::endl;
    logger << "solved in " << iterations << " iterations and " << lifts << " lifts." << std::endl;
    logger << "number of lifts: " << evenlifts << " even lifts and " << oddlifts << " odd lifts." << std::endl;

    // Free all explicitly allocated memory
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] str;
}

}
