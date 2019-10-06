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

// TODO: on-the-fly compress instead of on-the-fly inflate... (purely for performance)

#define CHECK_UNIQUE !NDEBUG // we should not see duplicate tangles anymore

namespace pg {

DTLSolver::DTLSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

DTLSolver::~DTLSolver()
{
}

/**
 * Returns True if player <pl> attracts vertex <v> to region <Z>.
 * Returns False if the vertex escapes to region <R>\<Z>.
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
        for (int to = *_out; to != -1; to = *++_out) {
            if (Z[to]) continue;
            if (R[to]) return false;
        }
        return true;
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
 * Very simple: partition subgame R into Even and Odd regions using top-down attractor computation.
 * Vertex <top> must be >= the highest vertex in R.
 */
void
DTLSolver::partition(bitset &R, int top, bitset &Even, bitset &Odd, bool check_distractions)
{
    for (; top!=-1; top--) {
        if (!R[top]) continue;
        if (check_distractions and Distractions[top]) continue;
#ifndef NDEBUG
        const int pr = priority[top];
#endif
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
            logger << "\033[1;33mregion\033[m ";
            logger << "\033[1;36m" << pr << "\033[m";
            for (unsigned int i=0; i<Zvec.size(); i++) {
                int v = Zvec[i];
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
            Zvec.clear();
        }
#endif
        if (pl == 0) Even |= Z;
        else Odd |= Z;
        Z.reset();
    }

#ifndef NDEBUG
    for (int v=0; v<n_nodes; v++) {
        if (R[v]) {
            logger << "vertex still in R: " << label_vertex(v) << std::endl;
            exit(-1);
        }
    }
#endif
}

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

        bool is_tangle = (tangle.size() > 1) or
            ((unsigned int)str[n] == n) or
            (std::find(out[n].begin(), out[n].end(), n) != out[n].end());
        if (!is_tangle) {
            tangle.clear();
            continue;
        }

        /**
         * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
         */

        for (const int v : tangle) bs_exits[v] = true;

        for (const int v : tangle) {
            if (owner[v] != pl) {
                const int *_out = outs + outa[v];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (G[to] and !bs_exits[to]) {
                        bs_exits[to] = true;
                        tangleto.push(to);
                    }
                }
            }
        }

        bs_exits.reset();

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
                    Q.push(v);
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
#if 0
            logger << "duplicate tangle" << std::endl;
            exit(-1);
#endif
            tangle.clear();
            tangleto.clear();
            continue;
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
DTLSolver::sptl(bitset &R, int top, int player, bitset &Even, bitset &Odd)
{
    bool changes = false;

    for (; top!=-1; top--) {
        if (!R[top]) continue; // not in the remaining game

        // If the vertex is a distraction, just skip it.
        // It will be attracted by the opponent. And not used in a tangle.
        if (Distractions[top]) continue;
#ifndef NDEBUG
        const int pr = priority[top];
#endif
        const int pl = priority[top]&1;

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            R[v] = false; // remove from <R>
            // if it's attracted by <player>, no longer Distractions (turn into non-distraction)
            if (player == pl) Distractions[v] = false;
#ifndef NDEBUG
            if (trace >= 2) Zvec.push(v);
#endif
            attractVertices(pl, v, R, Z, R);
            attractTangles(pl, v, R, Z, R);
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
            for (unsigned int i=0; i<Zvec.size(); i++) {
                int v = Zvec[i];
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
            Zvec.clear();
        }
#endif

        if (pl == 0) Even |= Z;
        else Odd |= Z;

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

                /**
                 * Extend any dominions that were found.
                 * (Any solved vertices are now in <Q>.)
                 */

                if (Q.nonempty()) {
                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        oink->solve(v, player, str[v]);
                        G[v] = false; // remove from Game
                        R[v] = false; // remove from region
                        // Distractions[v] = false; // remove from Distractions (doesn't matter though)
                        attractVertices(player, v, G, S, G);
                        attractTangles(player, v, G, S, G);
                    }

                    S.reset();
                }
            }
        }
        Z.reset();
    }

    return changes;
}


void
DTLSolver::go(const int player)
{
    /**
     * We currently have two variations of the algorithm.
     * The "pure" variation prunes first, then learns tangles.
     * A more "greedy" variation alternates tangle learning and pruning.
     */
#if 1
    Distractions.reset(); // start with a clean slate

    while (true) {
        if (trace) logger << "\033[1;38;5;196mprune step\033[m \033[1;36m" << steps << (player ? "-odd\033[m\n": "-even\033[m\n");
        steps++;

        Even.reset();
        Odd.reset();
        CurG = G;
        partition(CurG, n_nodes-1, Even, Odd, true);
        if (!(player?Odd:Even).any()) return; // no more regions of player <player>
        if (!prune(player)) break;
    }

    if (trace) logger << "\033[1;38;5;196msptl\033[1;36m " << (player ? "odd\033[m\n": "even\033[m\n");
    sptl_loop(player);
#else
    Distractions.reset(); // start with a clean slate

    while (true) {
        if (trace) logger << "\033[1;38;5;196mstep\033[m \033[1;36m" << steps << "\033[m\n";
        steps++;

        if (!sptl_loop(player)) break; // stop if no regions left
        if (!prune(player)) break; // stop if no new distractions
    }
#endif
}

/**
 * Run SPTL until no new tangles are found.
 * SPTL treats vertices in Distractions differently.
 * Vertices in Distractions are not top vertices for a region,
 * but simply added to opponent's region.
 */
bool
DTLSolver::sptl_loop(const int player)
{
    while (true) {
        Even.reset();
        Odd.reset();
        CurG = G;
        bool new_tangles = sptl(CurG, n_nodes-1, player, Even, Odd);
        if (!(player?Odd:Even).any()) return false; // no more regions of player <player>
        if (!new_tangles) return true; // no new tangles
    }
}

bool
DTLSolver::prune(const int player)
{
    bitset &RA = player == 0 ? Even : Odd;
    bitset &RO = player == 0 ? Odd : Even;

    /**
     * Find candidate distractions, i.e., vertices directly attracted to the opponent's region.
     * (Notice: every vertex in <RA> attracted to <RO> is a top vertex in <RA>.)
     */

    int n_candidates = 0;
    for (int v=n_nodes-1; v>=0; v--) {
        if (RA[v] and attracts(1-player, v, RO, RA)) {
            if (trace) {
                logger << "\033[1;38;5;202mcandidate distraction\033[m \033[1;36m" << label_vertex(v) << "\033[m\n";
            }
            Candidates[n_candidates++] = v;
        }
    }
    if (n_candidates == 0) return false;

    /**
     * Greedy attract from the opponent to the player, but maintain increasing threshold
     */

    int t = 0; // initial threshold

    // find lowest vertex in <RO> that can be attracted to <RA>.
    for (int v=0; v<=n_nodes; v++) {
        if (v == n_nodes) {
            if (Q.empty()) break; // done (otherwise: still vertices in the queue)
        } else {
            if (priority[v] > t and Q.empty()) t = priority[v];
            if (priority[v] <= t) {
                if (RO[v] and attracts(player, v, RA, RO)) {
                    RA[v] = true;
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

        // the above code is basically an optimization to quickly skip to the correct threshold
        // and find all vertices <= threshold that can be attracted (added to Q)
        // if we're here, we now found vertices to attract that are <= threshold
        
        bool changes = true;
        while (changes) {
            changes = false;

            // attract to RA
            while (Q.nonempty()) {
                const int u = Q.pop();
                RO[u] = false; // remove from <RO>
                attractVerticesM(player, u, RO, RA, RO, t);
                attractTanglesM(player, u, RO, RA, RO, t);
            }

            // now add all regions of player <player> in subgame <RO> to <RA>.
            partition(RO, n_nodes-1, SubEven, SubOdd, false);
            auto & SubA = player ? SubOdd : SubEven;
            for (int u=0; u<n_nodes; u++) {
                if (SubA[u]) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mattracted vertex\033[m " << label_vertex(u) << std::endl;
                    }
#endif
                    RA[u] = true;
                    Q.push(u);
                    changes = true;
                }
            }

            RO.swap(player ? SubEven : SubOdd); // quickly set RO := SubO
            SubEven.reset();
            SubOdd.reset();
        }

        // now check if any candidate distraction is now a non-distraction
        for (int x=0; x<n_candidates; x++) {
            const int u = Candidates[x];
            if (u == -1) continue;
            if (priority[u] <= t) break; // stop checking once we are below the threshold
            if (attracts(player, u, RA, G)) {
                Candidates[x] = -1;
                if (trace) { 
                    logger << "\033[1;38;5;202mnon-distraction\033[m \033[1;36m" << label_vertex(u) << "\033[m\n";
                }
            }
        }
    }

    /**
     * We're done, any remaining vertices in Candidates are now surely Distractions...
     */

    bool new_distractions = false;
    for (int x=0; x<n_candidates; x++) {
        const int u = Candidates[x];
        if (u == -1) continue;
        Distractions[u] = true;
        new_distractions = true;
    }

    if (trace) {
        logger << "\033[1;38;5;202mcurrent distractions\033[m";
        for (int u=n_nodes-1; u>=0; u--) {
            if (G[u] and Distractions[u]) logger << " \033[1;38;5;15m" << label_vertex(u) << "\033[m";
        }
        logger << std::endl;
    }

    return new_distractions;
}

void
DTLSolver::run()
{
    tin = new std::vector<int>[n_nodes];
    str = new int[n_nodes];

    Z.resize(n_nodes);
    S.resize(n_nodes);
    G = disabled;
    G.flip();

    Candidates = new int[n_nodes];
    Distractions.resize(n_nodes);

    Q.resize(n_nodes);
    Zvec.resize(n_nodes);
    tangleto.resize(n_nodes);
    bs_exits.resize(n_nodes);

    pea_vS.resize(n_nodes);
    pea_iS.resize(n_nodes);
    pea_S.resize(n_nodes);
    pea_vidx = new unsigned int[n_nodes];
    pea_root.resize(n_nodes);

    Even.resize(n_nodes);
    Odd.resize(n_nodes);
    CurG.resize(n_nodes);
    SubEven.resize(n_nodes);
    SubOdd.resize(n_nodes);

    while (true) {
        iterations++;

#ifndef NDEBUG
        int D=dominions;
        int T=tangles;
#endif

        if (!G.any()) break;
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-even\033[m\n";
        go(0);

        if (!G.any()) break;
        if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-odd\033[m\n";
        go(1);

#ifndef NDEBUG
        if (tangles == T and dominions == D) {
            logger << "stuck\n";
            exit(-1);
        }
#endif
    }

    logger << "solved with " << dominions << " dominions, "<< tangles << " tangles, in " << iterations << " iterations and " << steps << " pruning steps." << std::endl;

#ifndef NDEBUG
    // check if actually all solved
    for (int i=0; i<n_nodes; i++) {
        if (!disabled[i]) { logger << "search was incomplete!" << std::endl; exit(-1); }
    }
#endif

    // delete[] tangles
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] tin;
    delete[] str;
    delete[] pea_vidx;
    delete[] Candidates;
}

}
