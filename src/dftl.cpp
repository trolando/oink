/*
 * Copyright 2020-2021 Tom van Dijk
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

#include "dftl.hpp"

namespace pg {

/**
 * The "distraction-free" tangle learning algorithm.
 *
 * PROCEDURE iteration(player)
 * 1. Take <player>-priority vertices V as targets
 * 2. Tangle-attract to V from top to bottom.
 * 3. Closed regions become tangles
 * 4. Repeat 2-3 until no more closed regions
 * 5. Remove lowest region top from V
 * 6. Repeat 2-5 until V is empty
 *
 * Possible executions:
 * A. Run iteration(even) until no new tangles found; then iteration(odd) until game is solved
 * B. Run iteration(odd) until no new tangles found; then iteration(even) until game is solved
 * C. Run iteration(even) once then iteration(odd) once; repeat until game is solved
 *
 * Each iteration of this algorithm runs in polynomial time and consists of O(n²) many "steps".
 * Each "step" (2-3) is simply one iteration of single-player tangle learning.
 * After each step, either there are fewer regions, or V is smaller.
 *
 * Each DFTL iteration starts with V set to all even-priority (or odd-priority) vertices.
 * Each step then decomposes the game into regions dominated by the vertices in V.
 * That is, starting with the highest v in V, attract all vertices+tangles with priority <= pr(v).
 * Each highest vertex in v (the target of the attractor) is called the top vertex of its region.
 * Every *closed* region is then analysed to find new tangles, as in standard tangle learning.
 * If however no region is closed, then we assume that at least one vertex in V is a distraction.
 * We assume that "surely" the lowest top vertex is distracting, so we remove that vertex from V.
 *
 * Every step thus either removes a vertex from V, or learns new tangles (removing a region).
 * Hence each iteration has at most O(n²) steps and resuilts in at most O(n²) new tangles.
 * 
 * For the possible executiong A, B, C, the counter_dftl generator is a lower bound that requires
 * exponentially many steps.
 *
 * When the ``prepartition'' flag is set to True, this lower bound does not work anymore. It is
 * yet unclear whether the counterexample can be adapted to deal with the prepartition variation.
 * The prepartition variation modifies step 1 of iteration, by first computing a partition of the
 * game into even regions and odd regions, and only considering all even/odd vertices in even/odd
 * regions. In a way, this combines the "standard" mechanism of distraction removal (attract to
 * opponent) with the remove-lowest-open-top mechanism.
 */


DFTLSolver::DFTLSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}


DFTLSolver::~DFTLSolver()
{
}


/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R> from subgame <G>.
 * If <max_prio> >= 0, then only attract vertices with priority 0 <= pr <= max_prio.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
DFTLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
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
DFTLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio)
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
DFTLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio)
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
DFTLSolver::extractTangles(int startvertex, bitset &R, int *str)
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
            (str[n] == -1 and game->has_edge(n, n));
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
 * Single-player tangle learning.
 * Find new tangles for <player> given a set of vertices <V> in the subgame <R>.
 * <R> is typically the full remaining game.
 */
bool
DFTLSolver::sptl(bitset &V, bitset &R, const int player, int &lowest_top)
{
    bool changes = false;
    lowest_top = -1;

    for (auto top = V.find_last(); top != bitset::npos; top = V.find_prev(top)) {
        if (!R[top]) continue;

        lowest_top = top;
        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            R[v] = false; // remove from <R>
            attractVertices(player, v, R, Z, R, priority(top));
            attractTangles(player, v, R, Z, R, priority(top));
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            logger << "\033[1;33mregion\033[m \033[1;36m" << priority(top) << "\033[m";
            for (auto v = Z.find_last(); v != bitset::npos; v = Z.find_prev(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif

        bool closed_region = true;

        if (owner(top) == player) {
            closed_region = (str[top] != -1);
        } else {
            for (auto curedge = outs(top); *curedge != -1; curedge++) {
                if (R[*curedge]) {
                    closed_region = false;
                    break;
                }
            }
        }

        if (closed_region) {
            /**
             * Extract tangles by computing bottom SCCs starting at each top vertex.
             * Note: each bottom SCC must contain a top vertex.
             */

            std::fill(pea_vidx, pea_vidx+nodecount(), '\0');
            pea_curidx = 1;

            if (extractTangles(top, Z, str)) changes = true;
        }

        Z.reset();
    }

    return changes;
}


void
DFTLSolver::partition(bitset &R, int top, bitset &Even, bitset &Odd)
{
    for (; top!=-1 ;top--) {
        if (!R[top]) continue;

        const int pl = priority(top)&1;
        auto &Z = pl == 0 ? Even : Odd;

#ifndef NDEBUG
        if (trace >= 2) W = Z;
#endif

        Z[top] = true; // add to <Z>
        str[top] = -1;
        Q.push(top);

        while (Q.nonempty()) {
            const int v = Q.pop();
            R[v] = false; // remove from <R>
            attractVertices(pl, v, R, Z, R, priority(top));
            attractTangles(pl, v, R, Z, R, priority(top));
        }

#ifndef NDEBUG
        if (trace >= 2) {
            // report region
            W ^= Z;
            logger << "\033[1;33mregion\033[m ";
            logger << "\033[1;36m" << priority(top) << "\033[m";
            for (auto v = W.find_last(); v != bitset::npos; v = W.find_prev(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                if (str[v] != -1) logger << "->" << label_vertex(str[v]);
            }
            logger << std::endl;
        }
#endif
    }
}


bool
DFTLSolver::search(const int player)
{
    const int T = tangles;
    const int D = dominions;

    // just start with V is all <player>'s vertices
    V = G;
    if (player == 0) V -= Parity;
    else V &= Parity;

    // if we restrict V based on partition, do that now
    if (prepartition) {
#ifndef NDEBUG
        if (trace) {
            logger << "\033[1;38;5;33mcomputing initial partition\033[m " << std::endl;
        }
#endif

        CurG = G;
        Even.reset();
        Odd.reset();
        partition(CurG, nodecount()-1, Even, Odd);

        if (player == 0) V &= Even;
        else V &= Odd;

#ifndef NDEBUG
        if (trace) {
            logger << "\033[1;38;5;33mselected vertices\033[m:";
            for (auto v = V.find_last(); v != bitset::npos; v = V.find_prev(v)) {
                logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
            logger << "\033[1;38;5;33mrunning SPTL steps\033[m " << std::endl;
        }
#endif
    }

    while (V.any()) {
        steps++;

        CurG = G;
        int lowest_top;
        if (sptl(V, CurG, player, lowest_top)) {
            // Extend any dominions that were found.
            // (Any solved vertices are now in <S> and <dom_vector>.)
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
                    oink->solve(v, player, str[v]);
                    attractVertices(player, v, G, S, G, -1);
                    attractTangles(player, v, G, S, G, -1);
                }

                V -= S; // remove from V
                G -= S; // remove from G
                S.reset();
            }
        } else {
            // No new tangles, so remove the lowest top
            // We assume it is a distraction
            V[lowest_top] = false;
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[1;38;5;33mremoving \033[36m" << label_vertex(lowest_top) << "\033[m" << std::endl;
            }
#endif
        }
    }

    return tangles != T or dominions != D;
}


void
DFTLSolver::run()
{
    tin = new std::vector<int>[nodecount()];
    str = new int[nodecount()];

    Z.resize(nodecount());
    S.resize(nodecount());
    G = disabled;
    G.flip();
    CurG.resize(nodecount());
    V.resize(nodecount());
    W.resize(nodecount());

    Parity.resize(nodecount());
    for (int v=0; v<nodecount(); v++) Parity[v] = priority(v)&1;

    if (prepartition) {
        Even.resize(nodecount());
        Odd.resize(nodecount());
    }

    Q.resize(nodecount());

    tangleto.resize(nodecount());
    escapes.resize(nodecount());

    pea_vS.resize(nodecount());
    pea_iS.resize(nodecount());
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()];
    pea_root.resize(nodecount());

    if (interleaved) {
        // Interleave solving for player Even and player Odd

        int first_dominion = -1;
        bool continue_even = true;
        bool continue_odd = true;
        while (true) {
            iterations++;

            if (continue_odd) {
                if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-odd\033[m\n";
                continue_odd = search(1);
                if (dominions != 0 && first_dominion == -1) first_dominion = iterations;
                if (!G.any()) break;
            }

            if (continue_even) {
                if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations-1 << "-even\033[m\n";
                continue_even = search(0);
                if (dominions != 0 && first_dominion == -1) first_dominion = iterations;
                if (!G.any()) break;
            }

            if (!continue_even and !continue_odd) THROW_ERROR("solver stuck");
        }

        logger << "first dominion: " << first_dominion << std::endl;
    } else {
        // First solve for player Odd, then solve for player Even

        int even_iterations = 0;
        int odd_iterations = 0;
        int first_dominion = -1;

        while (G.any()) {
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-odd\033[m\n";
            iterations++;
            odd_iterations++;
            if (!search(1)) break;
            if (dominions != 0 && first_dominion == -1) first_dominion = odd_iterations;
        }

        while (G.any()) {
            if (trace) logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "-even\033[m\n";
            iterations++;
            even_iterations++;
            if (!search(0)) break;
            if (dominions != 0 && first_dominion == -1) first_dominion = even_iterations;
        }

        logger << "odd iterations: " << odd_iterations << std::endl;
        logger << "even iterations: " << even_iterations << std::endl;
        logger << "first dominion: " << first_dominion << std::endl;
    }

    logger << "found " << dominions << " dominions and "<< tangles << " tangles.\n";
    logger << "solved in " << iterations << " iterations and " << steps << " SPTL steps.\n";

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
}

}
