/*
 * Copyright 2019 Tom van Dijk, University of Twente
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
#include <cstring> // for memset

#include "ptl.hpp"

#define CHECK_UNIQUE !NDEBUG // we should not see duplicate tangles anymore

namespace pg {

PTLSolver::PTLSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

PTLSolver::~PTLSolver()
{
}

/**
 * Attract as player <pl> via <v> to <Z> vertices in <R> from subgame <G>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
PTLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z, bitset &G)
{
    // attract vertices with an edge to <v>
    for (auto curedge = ins(v); *curedge != -1; curedge++) {
        int from = *curedge;
        if (Z[from]) {
            // already in Z, maybe set strategy (for vertices in the original target set)
            if (owner(from) == pl and str[from] == -1) str[from] = v;
        } else if (R[from]) {
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
PTLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G)
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
 * Try to attract tangles for player <pl> to attractor set <Z>.
 * We look at all tangles with an escape edge to <v>.
 * All vertices in tangle <t> must be in <Z+R> and the opponent may not escape to subgame <G\Z>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
inline int
PTLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G)
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
 * Find strongly connected components in subgame <R> with strategy <str>.
 */
bool
PTLSolver::extractTangles(int startvertex, bitset &R, int *str)
{
    bool good = false;

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

        bool is_tangle = (tangle.size() > 1) or
            ((unsigned int)str[n] == n) or
            (str[n] == -1 and game.has_edge(n, n));
        if (!is_tangle) {
            tangle.clear();
            continue;
        }

        /**
         * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
         */

        for (const int v : tangle) bs_exits[v] = true;

        for (const int v : tangle) {
            if (owner(v) != pl) {
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
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
            if (pl == 0) {
                for (const int v : tangle) {
                    // if not yet added to solve queue, mark and add it
                    if (S0[v] == false) {
                        S0[v] = true;
                        Solver::solve(v, 0, str[v]);
                        SolvedQ0.push(v);
                    }
                }

            } else {
                for (const int v : tangle) {
                    // if not yet added to solve queue, mark and add it
                    if (S1[v] == false) {
                        S1[v] = true;
                        Solver::solve(v, 1, str[v]);
                        SolvedQ1.push(v);
                    }
                }
            }
            good = true;
            dominions++;
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
            // We can currently find duplicate tangles if they lead to a dominion?
            // logger << "duplicate tangle" << std::endl;
            // exit(-1);
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

        // and set p to pr
        tpr.push_back(pr);

        good = true;
        tangles++;
        tangle.clear();
        tangleto.clear();
    }

    pea_S.clear();
    return good;
}

bool
PTLSolver::search(bitset &R, int top, int player)
{
    bool changes = false;

    bitset Z(nodecount());

    // Partition the game <R>

    for (; top!=-1; top--) {
        // find next top
        if (!R[top]) continue;

        const int pr = priority(top);
        const int pl = priority(top)&1;

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            /*/ R[v] = false; // remove from <R> */
            attractVertices(pl, v, R, Z, R);
            // attract tangles if correct player, or <player> equals -1
            if ((1-pl) != player) attractTangles(pl, v, R, Z, R);
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

        if (player == -1 or pl == player) {
            /**
             * Now Z is the current region and <top> the top vertex
             */

            bool closed_region = true;

            /**
             * Check if Z is locally closed by looking at <top>.
             */

            if (owner(top) == pl) {
                if (str[top] == -1) closed_region = false;
            } else {
                for (auto curedge = outs(top); *curedge != -1; curedge++) {
                    if (!Z[*curedge] and R[*curedge]) {
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

                memset(pea_vidx, 0, sizeof(int[nodecount()]));
                pea_curidx = 1;
                changes |= extractTangles(top, Z, str);
            } else {
                /**
                 * Not a closed region, go recursive.
                 */

               // initialize path (for trace output) to <pr>
               path = std::to_string(pr);
               changes |= search_rec(Z, top, pl, R);
            }
        }

        R -= Z;
        Z.reset();
    }

    return changes;
}

bool
PTLSolver::search_rec(bitset &Region, int vtop, int player, bitset &Subgame)
{
    bool changes = false; // whether we found tangles/dominions

    bitset Y(Region); // the remaining subgame of <R>
    bitset Z(nodecount()); // move me
    bitset Lower(Subgame);

#if 0 /**/
    // First, attract to the open "exit" for the opponent
    // That is, remove Attr(top) from <Y>.

    Z[vtop] = true; // add to <Z>
    str[vtop] = -1;
    Q.push(vtop);

    while (Q.nonempty()) {
        const int v = Q.pop();
        attractVertices(1-player, v, Y, Z, Y);
        if (multiplayer) attractTangles(1-player, v, Y, Z, Y);
    }

    Y -= Z;
    Z.reset();
#else
    Y[vtop] = false; // ensure <vtop> is never attracted...
    // Lower[vtop] = false; // and that 
#endif

    /////////////////
    // After each iteration:
    // <Y> and <Lower> are updated, removing subregions
    // if top is "good": attract to Z, remove from <Y> and <Lower>
    // if <Z> is a "bad" ... ..
    /////////////////

    for (int top=vtop-1; top!=-1; top--) {
        // find next top below <vtop> in Y
        if (!Y[top]) continue;

        /*
        logger << "Contents of Y:";
        for (auto v = Y.find_first(); v != bitset::npos; v = Y.find_next(v)) logger << " " << label_vertex(v);
        logger << std::endl;
        logger << "Contents of Lower:";
        for (auto v = Lower.find_first(); v != bitset::npos; v = Lower.find_next(v)) logger << " " << label_vertex(v);
        logger << std::endl;
        */

        const int pr = priority(top);
        const int pl = priority(top)&1;

        if (pl == player) {
            Z[top] = true; // add to <Z>
            str[top] = -1;
            Q.push(top);

            while (Q.nonempty()) {
                const int v = Q.pop();
                attractVertices(player, v, Y, Z, Lower);
                attractTangles(player, v, Y, Z, Lower);
            }

#ifndef NDEBUG
            if (trace >= 2) {
                // report region
                logger << "\033[1;33mregion\033[m \033[1;36m" << path << "-" << pr << "\033[m";
                for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) {
                    logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1) logger << "->" << label_vertex(str[v]);
                }
                logger << std::endl;
            }
#endif

            /**
             * Now Z is the lowest region, <top> the top vertex
             */

            bool closed_region = true;

            /**
             * Lowest region, check if closed.
             */

            if (owner(top) == player) {
                if (str[top] == -1) closed_region = false;
            } else {
                for (auto curedge = outs(top); *curedge != -1; curedge++) {
                    if (!Z[*curedge] and Lower[*curedge]) {
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

                memset(pea_vidx, 0, sizeof(int[nodecount()]));
                pea_curidx = 1;
                changes |= extractTangles(top, Z, str);
            } else {
                /**
                 * Not a closed region, go recursive.
                 */

                std::string px = path;
                path += "-" + std::to_string(pr);
                changes |= search_rec(Z, top, player, Lower);
                path = px;
            }

            Y -= Z;
            Lower -= Z;
        } else {
            // Now recompute the region: reattract to <vtop>,
            // but not via higher (already removed) regions or via the bad vertex <top>

            // Y := remaining region except <top>
            // Lower := remaining region plus remaining lower game

            Y[top] = false; // remove <top> from consideration

            // logger << "Contents of Y2:";
            // for (auto v = Y2.find_first(); v != bitset::npos; v = Y2.find_next(v)) logger << " " << label_vertex(v);
            // logger << std::endl;

            Z[vtop] = true; // add to <Z>
            str[vtop] = -1;
            Q.push(vtop);

            while (Q.nonempty()) {
                const int v = Q.pop();
                attractVertices(player, v, Y, Z, Lower);
                attractTangles(player, v, Y, Z, Lower);
            }

            // Now Z is everything that can go to <vtop> avoiding <top>
            // And Y\Z is everything that must go via <top>

            // Set Y := everything to <vtop> not via <top>
            // Set Z := everything only to <vtop> via <top>

            Y.swap(Z);
            Y[vtop] = false; // but remove <vtop> again
            Z -= Y;
            Z[top] = true;

#ifndef NDEBUG
            if (trace >= 2) {
                // report region
                logger << "\033[1;33mregion\033[m \033[1;36m" << path << "-" << pr << "\033[m";
                for (auto v = Z.find_last(); v != bitset::npos; v = Z.find_prev(v)) {
                    logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                    if (str[v] != -1 and Z[str[v]]) logger << "->" << label_vertex(str[v]);
                }
                logger << std::endl;
            }
#endif

            std::string px = path;
            path += "-" + std::to_string(pr);
            changes |= search_rec(Z, top, player, Lower);
            path = px;
        }

        Z.reset();
    }

    return changes;
}


void
PTLSolver::solve()
{
    bitset CurG(nodecount());

    iterations = 0;

    bool play0 = true;
    bool play1 = true;

    while (G.any()) {
        iterations++;

        if (multiplayer) {
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "\033[m\n";
            CurG = G;
            search(CurG, nodecount()-1, -1);
        } else {
            // play for player 0 until no changes, then player 1

            while (play0) {
                iterations++;
                if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "\033[m, player Even\n";
                CurG = G;
                play0 = search(CurG, nodecount()-1, 0);

                if (SolvedQ0.nonempty()) {
                    Q.swap(SolvedQ0);
                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        // TODO: should we need the check?
                        if (!game.isSolved(v)) Solver::solve(v, 0, str[v]);
                        G[v] = false; // remove from Game
                        attractVertices(0, v, G, S0, G);
                        attractTangles(0, v, G, S0, G);
                    }
                }
            }

            if (play1) {
                if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "\033[m, player Odd\n";
                CurG = G;
                play1 = search(CurG, nodecount()-1, 1);
            }
        }

        /**
         * Extend any dominions that were found.
         * (Any solved vertices are now in <Q>.)
         */

        if (SolvedQ0.nonempty()) {
            Q.swap(SolvedQ0);
            while (Q.nonempty()) {
                const int v = Q.pop();
                if (!game.isSolved(v)) Solver::solve(v, 0, str[v]);
                G[v] = false; // remove from Game
                attractVertices(0, v, G, S0, G);
                attractTangles(0, v, G, S0, G);
            }
        }

        if (SolvedQ1.nonempty()) {
            Q.swap(SolvedQ1);
            while (Q.nonempty()) {
                const int v = Q.pop();
                if (!game.isSolved(v)) Solver::solve(v, 1, str[v]);
                G[v] = false; // remove from Game
                attractVertices(1, v, G, S1, G);
                attractTangles(1, v, G, S1, G);
            }
        }

        if (!play0 and !play1) break;

#if NDEBUG and 0
        if (iterations == 1000) {
            logger << "stuck\n";
            exit(-1);
        }
#endif
    }
}

void
PTLSolver::run()
{
    iterations = 0;
    dominions = 0;
    tangles = 0;

    tin = new std::vector<int>[nodecount()];
    str = new int[nodecount()];

    H.resize(nodecount());
    S0.resize(nodecount());
    S1.resize(nodecount());
    G = disabled;
    G.flip();

    Q.resize(nodecount());
    SolvedQ0.resize(nodecount());
    SolvedQ1.resize(nodecount());
    Zvec.resize(nodecount());
    tangleto.resize(nodecount());
    bs_exits.resize(nodecount());

    pea_vS.resize(nodecount());
    pea_iS.resize(nodecount());
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()];
    pea_root.resize(nodecount());

    solve();

    logger << "found " << dominions << " dominions." << std::endl;
    logger << "solved with " << tangles << " tangles and " << iterations << " iterations." << std::endl;

    // check if actually all solved
#ifndef NDEBUG
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { logger << "search was incomplete!" << std::endl; exit(-1); }
    }
#endif

    // delete[] tangles
    for (auto &x : tv) delete[] x;
    for (auto &x : tout) delete[] x;
    delete[] tin;
    delete[] str;
    delete[] pea_vidx;
}

}
