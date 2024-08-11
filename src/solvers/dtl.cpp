/*
 * Copyright 2017-2019 Tom van Dijk
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

#include <algorithm> // for sort
#include <climits> // for INT_MAX
#include <cstring> // for memset

#include "dtl.hpp"

#define CHECK_UNIQUE (1 or !NDEBUG) // detect duplicate dominions

namespace pg {

DTLSolver::DTLSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

DTLSolver::~DTLSolver()
{
}

/**
 * Returns True iff player <pl> attracts vertex <v> to region <Z> in subgame <R>.
 * (Either the vertex is owned by <pl> and can play to <Z>, or is owned by the opponent
 *  and all (and at least one) edges go to <Z>)
 */
inline bool
DTLSolver::attracts(const int pl, const int v, bitset &Z, bitset &R)
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
 * Attract as player <pl> via <v> to <Z>, vertices in <R> from subgame <G>.
 * If <max_prio> >= 0, then only attract vertices with priority 0 <= pr <= max_prio.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
DTLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (Z[from]) {
            // already in Z, maybe set strategy (for vertices in the original target set)
            if (owner(from) == pl and str[from] == -1) str[from] = v;
        } else if (R[from] and (max_prio < 0 or priority(from) <= max_prio)) {
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
DTLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio)
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
            ptr++; // skip strategy
            if (!this->G[v]) {
                // on-the-fly detect out-of-game tangles
                tpr[t] = -1; // delete the tangle
                return false; // is now a deleted tangle
            } else if (Z[v]) {
                continue; // already attracted
            } else if (!R[v]) {
                return false; // not contained in Z+R
            } else if (max_prio >= 0 and priority(v) > max_prio) {
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
DTLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    const auto &in_cur = tin[v];
    for (int from : in_cur) attractTangle(from, pl, R, Z, G, max_prio);
}

/**
 * Compute SCCs in subgraph induced by <R> and <str>.
 * Start the SCC computation at vertex <startvertex>.
 * Every SCC is then processed as a tangle.
 * If the tangle is closed, it is a dominion and added to <S> and <Q>.
 */
bool
DTLSolver::extractTangles(int startvertex, bitset &R, int *str)
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
    pea_vS.push(startvertex);
    pea_iS.push(0);
    pea_root[startvertex] = true;
    pea_vidx[startvertex] = pea_curidx++;
    while (pea_vS.nonempty()) {
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
                        pea_vS.push(to);
                        pea_iS.push(0);
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
                    pea_iS.back() = 1;
                    // beginVisiting
                    pea_vS.push(s);
                    pea_iS.push(0);
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
        pea_vS.pop();
        pea_iS.pop();
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
                        tangleto.push(to);
                    }
                }
            }
        }

        escapes.reset();

        /**
         * If there are no outgoing edges, then we have found a dominion.
         */

        if (tangleto.empty()) {
            // dominion
            if (trace) {
                logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m";
                logger << " (depth " << cur_depth <<")";
                if (cur_depth != 0) {
                    LOGIC_ERROR;
                    exit(-1);
                }
#ifndef NDEBUG
                if (trace >= 2) {
                    for (const int v : tangle) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
#endif
                logger << std::endl;
            }
            for (const int v : tangle) {
                // if not yet added to solve queue, mark and add it
                if (S[v] == false) {
                    S[v] = true;
                    SQ.push(v);
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
            if (trace >= 1) {
                logger << "\033[1;38;5;198mduplicate tangle\033[m";
#ifndef NDEBUG
                if (trace >= 2) {
                    for (int &v : tangle) {
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                    }
                }
#endif
                logger << std::endl;
            }
            tangle.clear();
            tangleto.clear();
            continue;
        }
#endif

        if (trace >= 1) {
            logger << "\033[1;38;5;198mnew tangle " << pr << "\033[m (" << tpr.size() << ")";
            logger << " (depth " << cur_depth <<")";
            if (trace >= 2) {
                for (const int v : tangle) {
                    logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                }
            }
            logger << " with " << tangleto.size() << " escape vertices.";
            logger << std::endl;
        }

        // new tangle idx
        const int tidx = tpr.size();

        // add back links to all normal vertices in our [out]
        for (unsigned int x = 0; x < tangleto.size(); x++) tin[tangleto[x]].push_back(tidx);

        // move tangleto into vout
        int* _tout = new int[tangleto.size()+1];
        std::copy(&tangleto[0], &tangleto[tangleto.size()], _tout);
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

        tangles++;
        new_tangles = true;
        tangle.clear();
        tangleto.clear();
    }

    pea_S.clear();
    return new_tangles;
}

/**
 * Find new tangles for <player> given the set of best vertices <V> in the subgame <R>.
 */
bool
DTLSolver::sptl(bitset &V, bitset &W, bitset &R, const int player)
{
    bool changes = false;

    for (int top=nodecount()-1; top>=0; top--) {
        if (!V[top] or !R[top]) continue; // meh

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            R[v] = false; // remove from <R>
#ifndef NDEBUG
            if (trace >= 2) Zvec.push(v);
#endif
            attractVertices(player, v, R, Z, R, priority(top));
            attractTangles(player, v, R, Z, R, priority(top));
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m \033[1;36m" << priority(top) << "\033[m";
            for (unsigned int i=0; i<Zvec.size(); i++) {
                int v = Zvec[i];
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
            Zvec.clear();
        }
#endif

        bool closed_region = true;

        if (owner(top) == player) {
            if (str[top] == -1) {
                closed_region = false;
            }
        } else {
            for (auto curedge = outs(top); *curedge != -1; curedge++) {
                int to = *curedge;
                if (R[to]) {
                    closed_region = false;
                    break;
                }
            }
        }

        if (closed_region and !W[top]) {
            W[top] = true;

            /**
             * Extract tangles by computing SCCs starting at each top vertex.
             * Note: each bottom SCC must contain a head vertex.
             */

            memset(pea_vidx, 0, sizeof(int[nodecount()]));
            pea_curidx = 1;

            if (extractTangles(top, Z, str)) changes = true;
        }

        Z.reset();
    }

    return changes;
}

/**
 * Computes the distances of <player>-priority vertices <V> in region <R> in subgame <G>.
 * Meaning the opponent can escape to <G>.
 */
bool
DTLSolver::computeDistanceValues(bitset &V, bitset &R, bitset &G, int* val, const int player)
{
    int dist = 1; // start distance 1

    while (V.any()) {
        // make the partition (attract lower to higher)
        for (int top=nodecount()-1; top>=0; top--) {
            if (V[top] and !Z[top]) {
                Z[top] = true;
                str[top] = -1; // avoid valgrind warnings...
                Q.push(top);
                while (Q.nonempty()) {
                    const int v = Q.pop();
                    attractVertices(player, v, R, Z, G, priority(top));
                    attractTangles(player, v, R, Z, G, priority(top));
                }
            }
        }

        // find candidate vertices to remove from V
        int n_candidates = 0;

        for (int v=nodecount()-1; v>=0; v--) {
            if (V[v] and !attracts(player, v, Z, G)) Candidates[n_candidates++] = v;
        }
        if (n_candidates == 0) {
            for (int v=0; v<nodecount(); v++) if (V[v]) val[v] = INT_MAX;
            Z.reset();
            return true;
        }

        // attract higher to lower
        int t = 0;
        for (int v=0; v<=nodecount();) {
            if (v == nodecount()) {
                if (Q.empty()) break;
            } else if (!R[v] or Z[v]) {
                v++;
                continue; // only attract vertices in <R\Z>
            } else {
                if (priority(v) > t and Q.empty()) t = priority(v);
                if (priority(v) <= t) {
                    if (attracts(player, v, Z, G)) {
                        Z[v] = true;
                        Q.push(v);
                    }
                    v++;
                    continue;
                } 
            }

            while (Q.nonempty()) {
                const int u = Q.pop();
                attractVertices(player, u, R, Z, G, t);
                attractTangles(player, u, R, Z, G, t);
            }

            for (int x=0; x<n_candidates; x++) {
                const int u = Candidates[x];
                if (u == -1) continue;
                if (priority(u) <= t) break; // stop checking once we are below the threshold
                if (attracts(player, u, Z, G)) Candidates[x] = -1;
            }
        }

        // remaining Candidates have no path from high to low, remove them from V

        bool V_changed = false;
        for (int x=0; x<n_candidates; x++) {
            const int u = Candidates[x];
            if (u == -1) continue;
            val[u] = dist;
            V[u] = false;
            V_changed = true;
        }

        Z.reset();
        dist++;

        if (!V_changed) {
            for (int v=0; v<nodecount(); v++) if (V[v]) val[v] = INT_MAX;
            return true;
        }
    }

    return false;
}



int
DTLSolver::best_vertices(bitset &V, bitset &R, const int player)
{
#ifndef NDEBUG
    // Assume: V is empty
    assert(V.none());
#endif

    bitset prevV(nodecount()); // the previous V
    bitset Z(nodecount()); // the current Z
    bitset CurR(nodecount());
    bitset CurZ(nodecount());
    int curdist = 0;

    /**
     * Initialize V := vertices in <R> with priority of <player>
     */

    bool hasV = false;
    for (int v=0; v<nodecount(); v++) {
        if (R[v] and (priority(v)&1) == player) {
            hasV = true;
            V[v] = true;
        }
    }
    if (!hasV) return 0;

    /**
     * Now we prune until no more vertices are removed from V, or V is empty.
     */

    while (true) {
        curdist++;
        prevV = V;

#ifndef NDEBUG
        if (trace >= 3) {
            logger << "\033[1;38;5;202mcurrent best vertices:\033[m";
            for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }
#endif

        // First compute the current <player> region Z

        CurR = R;
        for (int top=nodecount()-1; top>=0; top--) {
            if (V[top] and !Z[top]) {
                CurZ[top] = true; // add to <Z>
                str[top] = -1; // don't care, but needed for correct trace reporting
                Q.push(top);
                const int pr = priority(top);
                while (Q.nonempty()) {
                    const int v = Q.pop();
                    CurR[v] = false;
                    attractVertices(player, v, CurR, CurZ, CurR, pr);
                    attractTangles(player, v, CurR, CurZ, CurR, pr);
                }
                Z |= CurZ;
                CurZ.reset();
            }
        }

        /**
         * Now find all Candidates: v \in V that do not yet have a good path to Z
         */

        int n_candidates = 0;
        for (int v=nodecount()-1; v>=0; v--) {
            if (V[v] and !attracts(player, v, Z, R)) Candidates[n_candidates++] = v;
        }
        if (n_candidates == 0) return 1; // no more vertices are removed

        // Now attract from R to Z, "slowly"
        // Note: we have to be carefully increasing t. Why?
        // because we might attract tangles OR vertices at some threshold...
        const int max_prio = priority(nodecount()-1);
        for (int t=0; t<=max_prio; t++) {
            // add vertices in Z to Q
            for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) Q.push(v);

            // run tangle attractor to Z with the given threshold
            while (Q.nonempty()) {
                const int u = Q.pop();
                attractVertices(player, u, R, Z, R, t);
                attractTangles(player, u, R, Z, R, t);
            }

            // Now check if any candidate is now attracted (via the good path)
            for (int x=0; x<n_candidates; x++) {
                const int u = Candidates[x];
                if (u != -1) {
                    if (priority(u) <= t) break; // stop checking once we are below the threshold
                    if (attracts(player, u, Z, R)) Candidates[x] = -1;
                }
            }
        }

        /**
         * We have now attracted all possible vertices to <Z>.
         * Any remaining vertices in Candidates are now distractions...
         */

        bool changed_V = false;
        for (int x=0; x<n_candidates; x++) {
            const int u = Candidates[x];
            if (u == -1) continue;
            V[u] = false;
            if (trace >= 2) {
                if (!changed_V) {
                    logger << "\033[1;38;5;202mdistance " << curdist << "\033[m:";
                }
                logger << " \033[1;36m" << label_vertex(u) << "\033[m";
            }
            changed_V = true;
        }
        if (changed_V and trace >= 1) logger << std::endl;

        Z.reset();

        if (!changed_V) {
            // No vertices were removed from V
            return 1;
        } else if (!V.any()) {
            // All vertices were removed from V
            V.swap(prevV);
            return -1;
        }
    }

    return true;
}

/*
void
DTLSolver::new_algorithm(bitset &R, const int player)
{
    // COMPUTE DISTANCES

    int dodd[nodecount()], deven[nodecount()];
    bitset V(nodecount()); // the current V
    for (int v=0; v<nodecount(); v++) {
        if (R[v] and (priority(v)&1) == 0) V[v] = true;
    }
    computeDistanceValues(V, R, R, deven, 0);
    V.reset();
    for (int v=0; v<nodecount(); v++) {
        if (R[v] and (priority(v)&1) == 1) V[v] = true;
    }
    computeDistanceValues(V, R, R, dodd, 1);

    // COMPUTE PARTITION

    bitset Z(nodecount());

    for (; top!=-1; top--) {
        if (!R[top]) continue;

        const int pr = priority(top);
        const int pl = priority(top)&1;

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            attractVertices(pl, v, R, Z, R);
            attractTangles(pl, v, R, Z, R);
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
            for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif
      
        // now find tangles in the region
    }
}
*/


void
DTLSolver::search_rec(bitset &R, const int player, const int depth)
{
    bitset X(nodecount()); // vertices removed from R
    bitset V(nodecount()); // the best vertices each iteration
    bitset Z(nodecount()); // holds various subregions
    bitset Z2(nodecount()); // holds various subregions

    if (trace >= 2) {
        logger << "\033[1;38;5;202mcurrent region:\033[m";
        for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
            logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
        }
        logger << std::endl;
    }

    while (true) {
        // Compute the best vertices V
        V.reset();
        int result = best_vertices(V, R, player);
        if (result == 0) return; // done

        if (trace >= 1) {
            logger << "\033[1;38;5;202mbest vertices:\033[m";
            for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            if (result == 1) logger << " (stable)";
            logger << std::endl;
        }

        // If the best vertices are stable, run SPTL to learn local dominions and add them to X

        if (result == 1) {
            // Run SPTL for a while
            bool more = true;
            while (more) {
                Z = R;
                cur_depth = depth;
                more = sptl(V, W, Z, player);
            }
            Z.reset();

            // Add Dom(V,R) to X
            for (int top=nodecount()-1; top>=0; top--) {
                if (V[top] and !X[top]) {
                    X[top] = true; // add to <X>
                    str[top] = -1;
                    Q.push(top);
                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        R[v] = false;
                        attractVertices(player, v, R, X, R, priority(top));
                        attractTangles(player, v, R, X, R, priority(top));
                    }
                }
            }

            continue;
        }

        // If the best vertices are not stable, first try to attract to X.
        // Find the lowest threshold to attract to X.

        int t = 0; // initial threshold
        for (int v=0; v<nodecount(); v++) {
            if (R[v]) {
                if (priority(v) > t) {
                    if (Q.nonempty()) break;
                    // update threshold
                    t = priority(v);
                }
                if (attracts(player, v, X, R)) {
                    Z[v] = true;
                    str[v] = -1;
                    Q.push(v);
                }
            }
        }

        if (Q.nonempty()) {
            // Apparently, we found the lowest threshold t, they are in Z
            // Now attract with threshold t from R to Z

            while (Q.nonempty()) {
                const int u = Q.pop();
                R[u] = false;
                attractVertices(player, u, R, Z, R, t);
                attractTangles(player, u, R, Z, R, t);
            }

#ifndef NDEBUG
            if (trace >= 2) {
                logger << "X ∪= attracted with threshold " << t << ":\033[m";
                for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) {
                    logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                }
                logger << std::endl;
            }
#endif

            // Now remove the attractor to threshold vertices from Z (by moving them to X)
            for (int u=nodecount()-1; u >= 0; u--) {
                if (priority(u) < t) break; // done
                if (Z[u] and priority(u) == t) {
                    X[u] = true;
                    Q.push(u);
                }
            }

            while (Q.nonempty()) {
                const int u = Q.pop();
                Z[u] = false;
                attractVertices(1-player, u, Z, X, Z, t);
            }

            // And add all vertices in Z also to X
            X |= Z;

            // Recursively Z
            /*if (depth == 0)*/ search_rec(Z, player, depth+1);
            Z.reset();
        } else {
#ifndef NDEBUG
            if (trace >= 2) logger << "X ∪= Dom(V)" << std::endl;
#endif

            // Add Dom(V,R) to X
            for (int top=nodecount()-1; top>=0; top--) {
                if (V[top] and !X[top]) {
                    X[top] = true; // add to <X>
                    str[top] = -1;
                    Q.push(top);
                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        R[v] = false;
#ifndef NDEBUG
                        if (trace >= 3) Zvec.push(v);
#endif
                        attractVertices(player, v, R, X, R, priority(top));
                        attractTangles(player, v, R, X, R, priority(top));
                    }
#ifndef NDEBUG
                    if (trace >= 2) {
                        // report region
                        logger << "\033[1;33mregion\033[m ";
                        logger << "\033[1;36m" << priority(top) << "\033[m";
                        for (unsigned int i=0; i<Zvec.size(); i++) {
                            int v = Zvec[i];
                            logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                            if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                        }
                        logger << std::endl;
                        Zvec.clear();
                    }
#endif
                }
            }
        }
    }
}


bool
DTLSolver::search(const int player)
{
    const int T = tangles;
    const int D = dominions;

#ifndef NDEBUG
    int val1[nodecount()], val2[nodecount()];
    if (trace) {
        bitset V(nodecount()); // the current V
        for (int v=0; v<nodecount(); v++) {
            if (G[v] and (priority(v)&1) == player) V[v] = true;
        }
        computeDistanceValues(V, G, G, val1, player);
    }
#endif

    bitset CurG(G);
    search_rec(CurG, player, 0);
    W.reset();

    /**
     * Extend any dominions that were found.
     * (Any solved vertices are now in <SQ>.)
     */

    if (!dom_vector.empty()) {
        for (unsigned i = 0; i<dom_vector.size(); i+=2) {
            int v = dom_vector[i];
            int s = dom_vector[i+1];
            str[v] = s;
            Q.push(v);
        }
        dom_vector.clear();

        while (Q.nonempty()) {
            const int v = Q.pop();
            Solver::solve(v, player, str[v]);
            G[v] = false; // remove from Game
            attractVertices(player, v, G, S, G, -1);
            attractTangles(player, v, G, S, G, -1);
        }
        S.reset();
    }

#ifndef NDEBUG
    if (trace >= 3) {
        if (tangles != T and dominions == D) {
            bitset V(nodecount()); // the current V
            for (int v=0; v<nodecount(); v++) {
                if (G[v] and (priority(v)&1) == player) V[v] = true;
            }
            computeDistanceValues(V, G, G, val2, player);
            bool a_change = false;
            for (int v=0; v<nodecount(); v++) {
                if (val1[v] != val2[v]) {
                    logger << label_vertex(v) << ": " << val1[v] << " => ";
                    if (val2[v] == INT_MAX)  logger << "MAX\n";
                    else logger << val2[v] << "\n";
                    a_change = true;
                }
            }
            if (!a_change) {
                for (int v=0; v<nodecount(); v++) logger << label_vertex(v) << ": " << val1[v] << "\n";
                logger << "nothing changed!\n";
                // throw "nothing changed!";
            }
        }
    }
#endif

    return tangles != T or dominions != D;
}

void
DTLSolver::run()
{
    tin = new std::vector<int>[nodecount()];
    str = new int[nodecount()];

    dvalue = new int[nodecount()];

    Z.resize(nodecount());
    S.resize(nodecount());
    W.resize(nodecount());
    G = disabled;
    G.flip();

    Candidates = new int[nodecount()];

    Player.resize(nodecount());
    for (int v=0; v<nodecount(); v++) Player[v] = priority(v)&1;

    Q.resize(nodecount());
    SQ.resize(nodecount());

    Zvec.resize(nodecount());
    tangleto.resize(nodecount());
    escapes.resize(nodecount());

    pea_vS.resize(nodecount());
    pea_iS.resize(nodecount());
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()];
    pea_root.resize(nodecount());

    int even_iterations = 0;
    int odd_iterations = 0;

    if (!interleaved) {
        // First solve for player Odd, then solve for player Even

        while (G.any()) {
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-even\033[m\n";
            iterations++;
            even_iterations++;
            if (!search(0)) break;
        }

        while (G.any()) {
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-odd\033[m\n";
            iterations++;
            odd_iterations++;
            if (!search(1)) break;
        }
    } else {
        // Interleave solving for player Even and player Odd

        while (true) {
            iterations++;

            if (!G.any()) break;
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-even\033[m\n";
            bool found0 = search(0);

            if (!G.any()) break;
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-odd\033[m\n";
            bool found1 = search(1);

            if (!found0 and !found1) THROW_ERROR("solver stuck");
            
        }
    }

    logger << "found " << dominions << " dominions and "<< tangles << " tangles.\n";
    logger << "solved in " << iterations << " iterations and " << steps << " pruning steps.\n";
    logger << "odd iterations: " << odd_iterations << std::endl;
    logger << "even iterations: " << even_iterations << std::endl;

#ifndef NDEBUG
    // Check if the whole game is now solved
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { THROW_ERROR("search was incomplete!"); }
    }
#endif

    // Free all explicitly allocated memory
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] tin;
    delete[] str;
    delete[] pea_vidx;
    delete[] Candidates;
    delete[] dvalue;
}

}
