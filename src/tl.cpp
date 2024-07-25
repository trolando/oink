/*
 * Copyright 2021 Tom van Dijk
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

#include "tl.hpp"


#define PARTIALLY_CLOSED 0 // find tangles in partially locally closed regions, but costs time on practical games

namespace pg {

TLSolver::TLSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}


TLSolver::~TLSolver()
{
}


/**
 * Attract as player <pl> via <v> to <Z>, vertices in <R>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
void
TLSolver::attractVertices(const int pl, const int v, bitset &R, bitset &Z)
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
                unsigned int e = escs[from];
                if (e == 0) {
                    for (auto curedge = outs(from); *curedge != -1; curedge++) {
                        if (G[*curedge]) e++;
                    }
                }
                escs[from] = --e;
                if (e > 0) continue; // escapes
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
 * All vertices in tangle <t> must be in <Z+R> and the opponent may not escape to subgame <R\Z>.
 * Add attracted vertices to <Z> and to queue <Q>.
 * Update <str> with the obtained attractor strategy.
 */
bool
TLSolver::attractTangle(const int t, const int pl, bitset &R, bitset &Z, int maxpr)
{
    /**
     * Check if tangle is won by player <pl> and not deleted.
     */
    {
        const int tangle_pr = tpr[t];
        if (tangle_pr == -1) return false; // deleted tangle
        if (tangle_pr > maxpr) return false; // too high...
        if (pl != (tangle_pr&1)) return false; // not of desired parity
    }

    /**
     * Check if the tangle can escape to R\Z.
     */
    {
        int e = tescs[t];
        if (e == 0) {
            int v, *ptr = tout[t];
            while ((v=*ptr++) != -1) {
                if (G[v]) e++;
            }
        }
        tescs[t] = --e;
        if (e > 0) return false;
    }

    /**
     * On-the-fly detect out-of-game tangles
     */
    {
        int *ptr = tv[t];
        for (;;) {
            const int v = *ptr++;
            if (v == -1) break;
            ptr++; // skip strategy
            if (!G[v]) {
                // on-the-fly detect out-of-game tangles
                tpr[t] = -1; // delete the tangle
                return false; // is now a deleted tangle
            }
        }
    }

    if (maxpr == INT_MAX) {
        // attracted to a dominion, so delete the tangle
        tpr[t] = -1;
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
            if (!R[v]) continue; // not in <R> (might be higher region)
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
 * Current subgame is <R>.
 * Only tangles that are contained in <R>+<Z>
 * Write strategy to <str>.
 * Current attracting vertex is <v>.
 * Attracting for player <pl>.
 * (for trace) Attracting to region with priority <pr>.
 */
__attribute__((always_inline)) void
TLSolver::attractTangles(const int pl, int v, bitset &R, bitset &Z, int maxpr)
{
    const auto &in_cur = tin[v];
    for (int from : in_cur) attractTangle(from, pl, R, Z, maxpr);
}


/**
 * Compute SCCs in subgraph induced by <R> and <str>.
 * Start the SCC computation at vertex <startvertex>.
 * Every SCC is then processed as a tangle.
 * If the tangle is closed, it is a dominion and added to <S> and <Q>.
 */
bool
TLSolver::extractTangles(int startvertex, bitset &R, int *str)
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
        tescs.push_back(0);

        tangles++;
        tangle.clear();
        tangleto.clear();
        new_tangles = true;
    }

    pea_S.clear();
    return new_tangles;
}


bool
TLSolver::tl(void)
{
    bool new_tangles = false;

    std::fill(escs, escs+nodecount(), '\0');
    std::fill(tescs.begin(), tescs.end(), '\0');

    bool highest0 = true;
    bool highest1 = true;

    R = G;
    auto top = R.find_last();

    while (top != bitset::npos) {
        const int pr = priority(top);
        const int pl = priority(top)&1;
        unsigned int vcnt = 0;

        // attract from all heads with priority <pr> that are in <R>
        for (; top != bitset::npos; top = R.find_prev(top)) {
            if ((priority(top)&1) != (pr&1)) break; // otf compress
            // if (priority(top) != pr) break; // disable otf compress

            vcnt++;
            V[top] = true; // heads
            Z[top] = true; // add to <Z>
            str[top] = -1;
            Q.push(top);

#ifndef NDEBUG
            // maybe report event
            if (trace >= 3) {
                logger << "\033[1;37mattracting towards \033[36m" << label_vertex(top) << "\033[m for \033[1;36m" << pl << "\033[m" << std::endl;
            }
#endif

            while (Q.nonempty()) {
                const int v = Q.pop();
                R[v] = false; // remove from <R>
                attractVertices(pl, v, R, Z);
                attractTangles(pl, v, R, Z, pr);
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

                std::fill(pea_vidx, pea_vidx+nodecount(), '\0');
                pea_curidx = 1;

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

                    // start anew with escs
                    std::fill(escs, escs+nodecount(), '\0');
                    std::fill(tescs.begin(), tescs.end(), '\0');

                    for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) Q.push(v);

                    while (Q.nonempty()) {
                        const int v = Q.pop();
                        attractVertices(pl, v, G, Z);
                        attractTangles(pl, v, G, Z, INT_MAX);
                    }

                    for (auto v = Z.find_first(); v != bitset::npos; v = Z.find_next(v)) {
                        Solver::solve(v, pl, str[v]);
                    }

                    G -= Z; // remove from G
                    V.reset();
                    Z.reset();
                    return true; // end tl so that we can restart...
                } else {
                    const auto D = dominions;
                    for (auto v = V.find_last(); v != bitset::npos; v = V.find_prev(v)) {
                        if (extractTangles(v, Z, str)) new_tangles = true;
                    }

                    // Extend any dominions that were found.
                    // (Any solved vertices are now in <S>.)
                    if (D != dominions) {
                        // start anew with escs
                        std::fill(escs, escs+nodecount(), '\0');
                        std::fill(tescs.begin(), tescs.end(), '\0');

                        for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) Q.push(v);

                        while (Q.nonempty()) {
                            const int v = Q.pop();
                            attractVertices(pl, v, G, S);
                            attractTangles(pl, v, G, S, INT_MAX);
                        }

                        for (auto v = S.find_first(); v != bitset::npos; v = S.find_next(v)) {
                            Solver::solve(v, pl, str[v]);
                        }

                        G -= S; // remove from G
                        S.reset();
                        V.reset();
                        Z.reset();
                        return true; // end tl so that we can restart...
                    }
                }
            }
        }

        V.reset();
        Z.reset();

        if (pr&1) highest1 = false;
        else highest0 = false;
    }

    return new_tangles;
}


void
TLSolver::run()
{
    tin = new std::vector<int>[nodecount()];
    str = new int[nodecount()];
    escs = new unsigned int[nodecount()];

    V.resize(nodecount());
    W.resize(nodecount());
    R.resize(nodecount());
    Z.resize(nodecount());
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
    delete[] escs;
}

}
