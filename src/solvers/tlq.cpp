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
 * Attract vertices in <R> from subgame <Y> to <v> in region <Z> as player <pl>.
 * We attract only vertices in <R>, but the opponent may escape to <R> and <Y>.
 * Attracted vertices are added to <Z> and to the queue <Q>.
 */
void
TLQSolver::attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y)
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
                // check each exit
                bool escapes = false;
                auto curedge = outs(from);
                for (int to = *curedge; !escapes and to != -1; to = *++curedge) {
                    // check if escapes to a vertex in <R> or <Y> that is not in <Z>
                    if (!Z[to] and (R[to] or Y[to])) escapes = true;
                }
                if (escapes) continue;
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
    }

    /**
     * Check if tangle is contained in Z+R and if any vertices are not already in Z
     * We require at least one vertex in Z\R. Otherwise, don´t attract.
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
            } else {
                can_attract_new = true; // has vertices not yet attracted
            }
        }
        if (!can_attract_new) return false; // either no vertices in R\Z or strategy leaves R+Z
    }

    /**
     * Check if the tangle can escape to G\Z.
     */
    {
        int *ptr = tout[t];
        int v;
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
 */
void
TLQSolver::attractTangles(int pl, int v, bitset &R, bitset &Z, bitset &G)
{
    const auto &in_cur = tin[v];
    for (int from : in_cur) attractTangle(from, pl, R, Z, G);
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

    for (int v=vtop; v!=-1; v--) {
        if (!SG[v] || H[v]) continue;
        // if (priority(v) != pr) break; // stop attracting when we reach the next priority
        if ((priority(v)&1) != pl) break; // on-the-fly compression
        heads.push_back(v);
        H[v] = true;
        str[v] = -1;
        Q.push(v);
        while (Q.nonempty()) {
            const int u = Q.pop();
            attractVertices(pl, u, H, SG, SG);
            attractTangles(pl, u, SG, H, SG);
        }
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
     */

    bool leaky = false;
    for (auto v : heads) {
        if (!attracts(pl, v, H, SG)) {
            leaky = true;
            break;
            // TODO: consider learning tangles from partially leaky regions
        }
    }

    /**
     * If not leaky, then we can learn tangles from H
     */

    if (!leaky) {
        std::fill(pea_vidx, pea_vidx+nodecount(), '\0');
        pea_curidx = 1;
        for (auto v : heads) {
            auto T = tangles, D = dominions;
            extractTangles(v, H);
            if (D != dominions) {
                // maximize the found dominions
                for (auto v = S.find_last(); v != bitset::npos; v = S.find_prev(v)) {
                    Q.push(v);
                }
                while (Q.nonempty()) {
                    const int u = Q.pop();
                    attractVertices(pl, u, S, G, G);
                    attractTangles(pl, u, G, S, G);
                }
                // remove them from the game
                Wm |= S;
                Wo -= S;
                G -= S;
                S.reset();
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
                    while (Q.nonempty()) {
                        const int u = Q.pop();
                        attractVertices(1-k_pl, u, k_U, k_R, k_R);
                        attractTangles(1-k_pl, u, k_R, k_U, k_R);
                    }
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
                    while (Q.nonempty()) {
                        const int u = Q.pop();
                        attractVertices(k_pl, u, k_H, k_R, k_R);
                        attractTangles(k_pl, u, k_R, k_H, k_R);
                    }
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
    for (int v=0; v<=vtop; v++) {
        if (SG[v] and Wo[v]) {
            Q.push(v);
        }
    }
    if (Q.nonempty()) {
        while (Q.nonempty()) {
            const int v = Q.pop();
            if (Wm[v]) opponent_attracted_from_us = true;
            attractVertices(1-pl, v, Wo, SG, SG);
            attractTangles(1-pl, v, SG, Wo, SG);
        }
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

    tin = new std::vector<int>[nodecount()];
    str = new int[nodecount()];

    Q.resize(nodecount());
    W0.resize(nodecount());
    W1.resize(nodecount());
    S.resize(nodecount());
    G = disabled;
    G.flip();

    tangleto.resize(nodecount());
    escapes.resize(nodecount());

    pea_state.resize(nodecount()*2);
    pea_S.resize(nodecount());
    pea_vidx = new unsigned int[nodecount()];
    pea_root.resize(nodecount());

    bitset GG(G);
    Rs.push_back(&GG);
    solve(GG, nodecount()-1, nodecount(), nodecount());

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
    delete[] tin;
    delete[] pea_vidx;
}

}
