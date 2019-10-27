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
#include <cstring> // for memset

#include "dtl.hpp"

#define CHECK_UNIQUE !NDEBUG // detect duplicate dominions
#define FASTDOMINION 0

namespace pg {

DTLSolver::DTLSolver(Oink *oink, Game *game) : Solver(oink, game)
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
    const int *_out = outs + outa[v];
    if (owner[v] == pl) {
        for (int to = *_out; to != -1; to = *++_out) {
            if (Z[to]) return true;
        }
        return false;
    } else {
        bool hasZ = false;
        for (int to = *_out; to != -1; to = *++_out) {
            if (Z[to]) { hasZ = true; continue; }
            if (R[to]) return false;
        }
        return hasZ;
    }
}

/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R> from subgame <G>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
DTLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z, bitset &G)
{
    // attract vertices with an edge to <v>
    const int *_in = ins + ina[v];
    for (int from = *_in; from != -1; from = *++_in) {
        if (Z[from]) {
            // already in Z, maybe set strategy (for vertices in the original target set)
            if (owner[from] == pl and str[from] == -1) str[from] = v;
        } else if (R[from]) {
            if (owner[from] != pl) {
                // check if opponent can escape
                bool escapes = false;
                const int *_out = outs + outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (G[to] and !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            str[from] = owner[from] == pl ? v : -1;
            Q.push(from);
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << "\033[m by \033[1;36m" << pl << "\033[m";
                if (owner[from] == pl) logger << " (via " << label_vertex(v) << ")" << std::endl;
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
DTLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G)
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
            } else if (R[v]) {
                can_attract_new = true; // has vertices not yet attracted
            } else {
                return false; // not contained in Z+R
            }
        }
        if (!can_attract_new) return false;
    }

    /**
     * Check if the tangle can escape to R\Z.
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
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
inline int
DTLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G)
{
    int added = 0;
    const auto &in_cur = tin[v];
    for (int from : in_cur) {
        if (attractTangle(from, pl, R, Z, G)) {
            added++;
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[1;36m" << tpr[from] << "\033[m-tangle " << from << " to \033[1;36m" << pl << "\033[m" << std::endl;
            }
#endif
        }
    }
    return added;
}

/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R> from subgame <G>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
DTLSolver::attractVerticesM(const int pl, const int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    // attract vertices with an edge to <v>
    const int *_in = ins + ina[v];
    for (int from = *_in; from != -1; from = *++_in) {
        if (priority[from] > max_prio) continue;
        if (Z[from]) {
            // already in Z, maybe set strategy (for vertices in the original target set)
            if (owner[from] == pl and str[from] == -1) str[from] = v;
        } else if (R[from]) {
            if (owner[from] != pl) {
                // check if opponent can escape
                bool escapes = false;
                const int *_out = outs + outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (G[to] and !Z[to]) {
                        escapes = true;
                        break;
                    }
                }
                if (escapes) continue;
            }
            // attract
            Z[from] = true;
            str[from] = owner[from] == pl ? v : -1;
            Q.push(from);
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << "\033[m by \033[1;36m" << pl << "\033[m";
                if (owner[from] == pl) logger << " (via " << label_vertex(v) << ")" << std::endl;
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
DTLSolver::attractTangleM(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio)
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
            } else if (priority[v] > max_prio) {
                return false;
            } else if (Z[v]) {
                continue; // already attracted
            } else if (R[v]) {
                can_attract_new = true; // has vertices not yet attracted
            } else {
                return false; // not contained in Z+R
            }
        }
        if (!can_attract_new) return false;
    }

    /**
     * Check if the tangle can escape to R\Z.
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
 * Current subgame is <Z> + <R>.
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
inline int
DTLSolver::attractTanglesM(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
{
    int added = 0;
    const auto &in_cur = tin[v];
    for (int from : in_cur) {
        if (attractTangleM(from, pl, R, Z, G, max_prio)) {
            added++;
#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracted \033[1;36m" << tpr[from] << "\033[m-tangle " << from << " to \033[1;36m" << pl << "\033[m" << std::endl;
            }
#endif
        }
    }
    return added;
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
    const int pr = priority[startvertex];
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

        if (owner[n] != pl) {
            const int *_out = outs + outa[n];
            if (i>0) {
                // finishEdge
                const int w = _out[i-1];
                if (pea_vidx[w] < pea_vidx[n]) {
                    pea_vidx[n] = pea_vidx[w];
                    pea_root[n] = false;
                }
            }
            for (;;) {
                const int to = _out[i];
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
            (std::find(out[n].begin(), out[n].end(), n) != out[n].end());
        if (!is_tangle) {
            tangle.clear();
            continue;
        }

        /**
         * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
         */

        for (const int v : tangle) escapes[v] = true;

        for (const int v : tangle) {
            if (owner[v] != pl) {
                const int *_out = outs + outa[v];
                for (int to = *_out; to != -1; to = *++_out) {
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
                }
            }
            dominions++;
            new_tangles = true;
            tangle.clear();
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
#if !FASTDOMINION
            // This occasionally happens when we just learned a dominion.
            // For example game vb163: tangle 1 is attracted towards dominion 5 or to 7.
            // When dominion 5 is learned, the tangle would now attract to 7 but because we are
            // just continuing our loop, now tangle 1 is in a region of its own temporarily.
            logger << "duplicate tangle" << std::endl;
            // exit(-1);
#endif
            // tangle.clear();
            // tangleto.clear();
            // continue;
        }
#endif

        if (trace >= 1) {
            logger << "\033[1;38;5;198mnew tangle " << pr << "\033[m (" << tpr.size() << ")";
#ifndef NDEBUG
            if (trace >= 2) {
                for (const int v : tangle) {
                    logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                }
            }
            logger << " with " << tangleto.size() << " escape vertices.";
#endif
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

        // and set p to pr and current region to BOT
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
 * Single-player tangle learning.
 * Find tangles for player <player> in subgame <R>.
 * The final partition is saved to <Even> and <Odd>.
 */
bool
DTLSolver::sptl(bitset &R, int top, int player)
{
    bool changes = false;

    for (; top!=-1; top--) {
        if (!R[top]) continue; // not in the remaining game

        const int pl = priority[top]&1;

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            R[v] = false; // remove from <R>
#ifndef NDEBUG
            if (trace >= 2) Zvec.push(v);
#endif
            attractVertices(pl, v, R, Z, R);
            attractTangles(pl, v, R, Z, R);
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m \033[1;36m" << priority[top] << "\033[m";
            for (unsigned int i=0; i<Zvec.size(); i++) {
                int v = Zvec[i];
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
            Zvec.clear();
        }
#endif

        if (player == pl) {
            /**
             * Now Z is the lowest region, <top> the top vertex
             */

            bool closed_region = true;

            /**
             * Lowest region, check if globally closed.
             */

            if (owner[top] == player) {
                if (str[top] == -1) {
                    closed_region = false;
                }
            } else {
                const int *_out = outs + outa[top];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (R[to]) {
                        closed_region = false;
                        break;
                    }
                }
            }

            if (closed_region) {
                /**
                 * Extract tangles by computing SCCs starting at each top vertex.
                 * Note: each bottom SCC must contain a head vertex.
                 */

                memset(pea_vidx, 0, sizeof(int[n_nodes]));
                pea_curidx = 1;

                if (extractTangles(top, Z, str)) changes = true;

#if FASTDOMINION
                /**
                 * Extend any dominions that were found.
                 * (Any solved vertices are now in <SQ>.)
                 */

                if (SQ.nonempty()) {
                    SQ.swap(Q);
                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        oink->solve(v, player, str[v]);
                        G[v] = false; // remove from Game
                        R[v] = false; // remove from region
                        attractVertices(player, v, G, S, G);
                        attractTangles(player, v, G, S, G);
                    }

                    S.reset();
                }
#endif
            }
        }
        Z.reset();
    }

    return changes;
}

void
DTLSolver::search_rec(bitset &R, const int player)
{
    bitset Ra(n_nodes); // region for recursive call
    bitset V(n_nodes); // the current V
    bitset prevV(n_nodes); // the previous V
    bitset Z(n_nodes); // the current Z
    bitset lastZ(n_nodes); // copy of Z, for after pruning

    /**
     * This is a do-while loop that avoids some of the recursion.
     */

    do {
        /**
         * Initialize V := vertices in <R> with priority of <player>
         * We monotonically remove vertices from V that are distractions, w.r.t. the
         * tangles that we currently know.
         */

        bool hasV = false;
        for (int v=0; v<n_nodes; v++) {
            if (R[v] and (priority[v]&1) == player) {
                hasV = true;
                V[v] = true;
            }
        }
        if (!hasV) {
            if (trace) {
                for (int v=0; v<n_nodes; v++) if (R[v]) reg[v] = regidx;
                regidx++;
            }

            return; // No vertices left with <player>'s priority, return.
        }

        // we only count partition steps where we actually do work.

        if (trace >= 3) {
            logger << "\033[1;38;5;196mstep\033[m \033[1;36m" << steps << (player ? "-odd\033[m\n": "-even\033[m\n");
        }
        steps++;

        /**
         * Now we prune until no more vertices are removed from V, or V is empty.
         */

        bool lastIsEmpty = false;
        while (true) {
            prevV = V;

#ifndef NDEBUG
            if (trace >= 3) {
                logger << "\033[1;38;5;202mcurrent best vertices:\033[m";
                for (int v=n_nodes-1; v>=0; v--) {
                    if (V[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                }
                logger << std::endl;
            }
#endif

            // First compute the current <player> region Z
            for (int top=n_nodes-1; top>=0; top--) {
                if (!V[top] or Z[top]) continue;

                Z[top] = true; // add to <Z>
                str[top] = -1; // don't care, but needed for correct trace reporting
                Q.push(top);
                while (Q.nonempty()) {
                    const int v = Q.pop();
#ifndef NDEBUG
                    if (trace >= 3) Zvec.push(v);
#endif
                    attractVerticesM(player, v, R, Z, R, priority[top]);
                    attractTanglesM(player, v, R, Z, R, priority[top]);
                }

#ifndef NDEBUG
                if (trace >= 3) {
                    // report region
                    logger << "\033[1;33mregion\033[m ";
                    logger << "\033[1;36m" << priority[top] << "\033[m";
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

            lastZ = Z; // Make a copy, in case this was the last one. (Saves recomputing cost.)

            /**
             * Now find all Candidates: vertices in V that could be removed from V if the player
             * has no good path through the opponent's regions back to the player's regions.
             */

            int n_candidates = 0;
            for (int v=n_nodes-1; v>=0; v--) {
                // vertices in V that the opponent attracts out of Z
                if (V[v] and !attracts(player, v, Z, R)) {
#ifndef NDEBUG
                    if (trace >= 3) {
                        logger << "\033[1;38;5;202mcandidate distraction\033[m \033[1;36m" << label_vertex(v) << "\033[m\n";
                    }
#endif
                    Candidates[n_candidates++] = v;
                }
            }
            if (n_candidates == 0) break; // No more vertices can be removed from V.

            // Now attract from R to Z, "slowly"
            int t = 0; // initial threshold
            for (int v=0; v<=n_nodes; v++) {
                if (v == n_nodes) {
                    if (Q.empty()) break; // done (otherwise: still vertices in the queue)
                } else if (!R[v]) {
                    continue;
                } else {
                    if (priority[v] > t and Q.empty()) t = priority[v];
                    if (priority[v] <= t) {
                        if (!Z[v] and attracts(player, v, Z, R)) {
                            Z[v] = true;
                            Q.push(v);
#ifndef NDEBUG
                            if (trace >= 2) {
                                logger << "\033[1;37mattracted opponent vertex\033[m " << label_vertex(v) << " with pr=" << priority[v] << ", t=" << t << std::endl;
                            }
#endif
                        }
                        continue;
                    } 
                }
                v--;

                // The above quickly skips to the correct threshold
                // and find all vertices <= threshold that can be attracted (added to Q)
                // if we're here, we now found vertices to attract that are <= threshold

                while (Q.nonempty()) {
                    const int u = Q.pop();
                    attractVerticesM(player, u, R, Z, R, t);
                    attractTanglesM(player, u, R, Z, R, t);
                }

                // Now check if any candidate distraction stays
                // meaning we check which candidate vertices actually have a good path back to Z
                for (int x=0; x<n_candidates; x++) {
                    const int u = Candidates[x];
                    if (u == -1) continue;
                    if (priority[u] <= t) break; // stop checking once we are below the threshold
                    if (attracts(player, u, Z, R)) {
                        Candidates[x] = -1;
#ifndef NDEBUG
                        if (trace >= 3) { 
                            logger << "\033[1;38;5;202mnon-distraction\033[m \033[1;36m" << label_vertex(u) << "\033[m\n";
                        }
#endif
                    }
                }
            }

            /**
             * We have now attracted all possible vertices to <Z>.
             * Any remaining vertices in Candidates are now distractions...
             */

            bool new_distractions = false;
            for (int x=0; x<n_candidates; x++) {
                const int u = Candidates[x];
                if (u == -1) continue;
                V[u] = false;
                new_distractions = true;
#ifndef NDEBUG
                if (trace >= 2) { 
                    logger << "\033[1;38;5;202mdistraction\033[m \033[1;36m" << label_vertex(u) << "\033[m\n";
                }
#endif
            }

            Z.reset();

            if (!new_distractions) {
                // No vertices were removed from V
                break;
            } else if (!V.any()) {
                // All vertices were removed from V
                V.swap(prevV);
                lastIsEmpty = true;
                break;
            }
        }

        if (trace >= 2) {
            logger << "\033[1;38;5;202mbest vertices in current block:\033[m";
            for (int v=n_nodes-1; v>=0; v--) {
                if (V[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }

#if 1
        // Just use last computed Z
        Z = lastZ;
#else
        // Compute the current <player> region Z
        for (int top=n_nodes-1; top>=0; top--) {
            if (V[top] and !Z[top]) {
                Z[top] = true; // add to <Z>
                str[top] = -1;
                Q.push(top);
                while (Q.nonempty()) {
                    const int v = Q.pop();
                    attractVerticesM(player, v, R, Z, R, priority[top]);
                    attractTanglesM(player, v, R, Z, R, priority[top]);
                }
            }
        }
#endif

        if (!lastIsEmpty) {
            // If <Z> contains new tangles, then the fixpoint was not empty.
            // So if the fixpoint was empty, then there are no new tangles to learn.
            Ra = Z;
            sptl(Ra, n_nodes-1, player);
            Ra.reset();
        }

        if (trace) {
            for (int v=0; v<n_nodes; v++) if (Z[v]) reg[v] = regidx;
            if (!lastIsEmpty) regflag[regidx] = true;
            regidx++;
        }

        /**
         * Compute the R_p sets, of vertices attracted to <R> but via a vertex of opponent's
         * priority p, and p is then the highest priority to reach any vertex in V.
         */
        int t = 0; // initial threshold
        for (int v=0; v<=n_nodes; v++) {
            if (v == n_nodes) {
                if (Q.empty()) break; // done (otherwise: still vertices in the queue)
            } else if (!R[v]) {
                continue; // only vertices in R
            } else {
                if (priority[v] > t and Q.empty()) t = priority[v];
                if (priority[v] <= t) {
                    if (!Z[v] and attracts(player, v, Z, R)) {
                        Z[v] = true;
                        Q.push(v);
                    }
                    continue;
                } 
            }
            v--;

            while (Q.nonempty()) {
                const int u = Q.pop();
                Ra[u] = true; // also add to Ra
                attractVerticesM(player, u, R, Z, R, t);
                attractTanglesM(player, u, R, Z, R, t);
            }

            // We now have Ra := the R_t set associated with threshold t.

#ifndef NDEBUG
            if (trace >= 3) {
                logger << "\033[1;38;5;202mattracted with threshold " << t << ":\033[m";
                for (int v=n_nodes-1; v>=0; v--) {
                    if (Ra[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                }
                logger << std::endl;
            }
#endif

            // Recursively partition set R_t <Ra>.

            search_rec(Ra, player);
            Ra.reset();

#if FASTDOMINION
            R &= G; // in case there were new dominions
#endif
        }

        // Continue partitioning the remainder
        R -= Z;
        Z.reset();
        V.reset();

#ifndef NDEBUG
        if (trace >= 3) {
            logger << "\033[1;38;5;202mremainder:\033[m";
            for (int v=n_nodes-1; v>=0; v--) {
                if (R[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }
#endif
    } while (R.any());
}

bool
DTLSolver::search(const int player)
{
    if (trace) {
        memset(reg, 0, sizeof(int[n_nodes]));
        regflag.reset();
        regidx = 1; // initialize
    }

    const int T = tangles;
    const int D = dominions;

    bitset CurG(G);
    search_rec(CurG, player);

    if (trace) {
#ifndef NDEBUG
        for (int v=0; v<n_nodes; v++) {
            if (G[v] and (reg[v] <= 0 or reg[v] >= regidx)) {
                logger << "logic error\n";
                exit(-1);
            }
        }
#endif
        for (int i=1; i<regidx; i++) {
            logger << "\033[1;38;5;33mblock\033[m \033[1;38;5;231m" << i;
            if (regflag[i]) logger << "*";
            logger << "\033[m:";
            for (int v=0; v<n_nodes; v++) {
                if (reg[v] == i) {
                    if ((priority[v]&1) == player) logger << " \033[38;5;231;1m" << label_vertex(v) << "\033[m";
                    else logger << " " << label_vertex(v);
                }
            }
            logger << std::endl;
        }
    }

#if !FASTDOMINION
    /**
     * Extend any dominions that were found.
     * (Any solved vertices are now in <SQ>.)
     */

    if (SQ.nonempty()) {
        SQ.swap(Q);
        while (Q.nonempty()) {
            const int v = Q.pop();
            oink->solve(v, player, str[v]);
            G[v] = false; // remove from Game
            attractVertices(player, v, G, S, G);
            attractTangles(player, v, G, S, G);
        }
        S.reset();
    }
#endif

    return tangles != T or dominions != D;
}

void
DTLSolver::run()
{
    tin = new std::vector<int>[n_nodes];
    str = new int[n_nodes];

    reg = new int[n_nodes];
    regflag.resize(n_nodes);

    Z.resize(n_nodes);
    S.resize(n_nodes);
    G = disabled;
    G.flip();

    Candidates = new int[n_nodes];

    Q.resize(n_nodes);
    SQ.resize(n_nodes);

    Zvec.resize(n_nodes);
    tangleto.resize(n_nodes);
    escapes.resize(n_nodes);

    pea_vS.resize(n_nodes);
    pea_iS.resize(n_nodes);
    pea_S.resize(n_nodes);
    pea_vidx = new unsigned int[n_nodes];
    pea_root.resize(n_nodes);

#if 1
    // First solve for player Even, then solve for player Odd

    while (G.any()) {
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-even\033[m\n";
        iterations++;
        if (!search(0)) break;
    }

    while (G.any()) {
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-odd\033[m\n";
        iterations++;
        if (!search(1)) break;
    }
#else
    // Interleave solving for player Even and player Odd

    while (true) {
        iterations++;

        if (!G.any()) break;
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-even\033[m\n";
        bool found0 = search(0);

        if (!G.any()) break;
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-odd\033[m\n";
        bool found1 = search(1);

        if (!found0 and !found1) {
            logger << "stuck\n";
            exit(-1);
        }
    }
#endif

    logger << "found " << dominions << " dominions and "<< tangles << " tangles.\n";
    logger << "solved in " << iterations << " iterations and " << steps << " pruning steps.\n";

#ifndef NDEBUG
    // Check if the whole game is now solved
    for (int i=0; i<n_nodes; i++) {
        if (!disabled[i]) { logger << "search was incomplete!" << std::endl; exit(-1); }
    }
#endif

    // Free all explicitly allocated memory
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] tin;
    delete[] str;
    delete[] pea_vidx;
    delete[] Candidates;
    delete[] reg;
}

}
