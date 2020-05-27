/*
 * Copyright 2019-2020 Tom van Dijk, University of Twente
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

#include "zlkq.hpp"

namespace pg {

/**
 * This is an implementation of the universal recursive algorithm run on a quasi-polynomial tree
 * The original version was based on https://www.mimuw.edu.pl/~parys/publications/2018-parity-algorithm.pdf
 * I modified the algorithm to produce a strategy as well.
 * I changed it with the optimizations from Lehtinen et al (arXiv)
 * Furthermore I added some optimizations, i.e., shortcuts in the tree
 */

ZLKQSolver::ZLKQSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

ZLKQSolver::~ZLKQSolver()
{
}

/**
 * Attract vertices in <R> from subgame <Y> to <v> in region <Z> as player <pl>.
 * We attract only vertices in <R>, but the opponent may escape to <R> and <Y>.
 * Attracted vertices are added to <Z> and to the queue <Q>.
 */
void
ZLKQSolver::attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y)
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

void
ZLKQSolver::solve(bitset &SG, int vtop, int pr, const int pe, const int po)
{
    /**
     * This is based on a universal tree of height pr/2, with parameter n associated with pe or po
     * We use several implementation tricks to make optimal use of memory and allocate as few as possible
     * For example, we record partial winning sets for Even and Odd in W0 and W1, because all partial
     * winning sets between calls are disjoint (so, safe to share the same bit array)
     */

    /**
     * Record the number of recursive calls (visits)
     */

    iterations++;

    /**
     * Check if we have run out of precision / nodes! (parameter n)
     * If our opponent's precision is gone, assume it's all won by us.
     * If our own precision is gone, assume it's all won by the opponent.
     */

    if (po <= 0) {
        W0 |= SG;
        return;
    }

    if (pe <= 0) {
        W1 |= SG;
        return;
    }

    /**
     * Move <vtop> to the true top (highest vertex that is in SG)
     */

    while (vtop >= 0 and !SG[vtop]) vtop--;

    /**
     * Check if the game is empty
     */
    if (vtop == -1) return; // empty game, bye

    /**
     * Optimization: set pr/pl to top priority/player
     * This is like taking a shortcut in the Tree
     */

    pr = priority(vtop);
    const int pl = pr&1;

#ifndef NDEBUG
    if (trace >= 2) logger << "IN: " << pr << " (" << pl << ") E=" << pe << " O=" << po << "\n";
#endif

    /**
     * Expand the left side of the current level of the universal tree
     */

    if (pl == 0) solve(SG, vtop, pr, pe, po/2);
    else solve(SG, vtop, pr, pe/2, po);

    auto &Wm = pl == 0 ? W0 : W1; // my W
    auto &Wo = pl == 0 ? W1 : W0; // opponent's W

    /**
     * Now compute <R> := vertices in subgame that we won in the left side
     * That is, the intersection of SG and Wm.
     *
     * Consider R to be the game that remains after the steps on the left in the universal tree.
     * Fmore, SG is immutable, R is a local bitset
     */

    bitset R(SG);
    R &= Wm;

    /**
     * Update vtop to the highest vertex in the intersection
     */
    while (vtop>=0 and !R[vtop]) vtop--;

    /**
     * Compute H := Attr(vertices of <pr> in <R>)
     */
    bitset H(nodecount());

    for (int v=vtop; v!=-1; v--) {
        if (priority(v) != pr) break; // stop attracting ... ;) [todo: add otf compression ?]
        if (R[v]) {
            H[v] = true;
            str[v] = -1;
            Q.push(v);
            while (Q.nonempty()) {
                attractVertices(pl, Q.pop(), H, R, R);
            }
        }
    }

    /**
     * Reset Wo and Wm for vertices in <R> - <H>;
     * Set H to the remaining subgame (G' in the universal paper)
     */

    // reset opponent and our region
    Wo -= R;
    Wm -= R;
    // but set everything in H to ours
    Wm |= H;
    // and set H := R - H
    H ^= R;

#ifndef NDEBUG
    /**
     * If so desired, report our current status...
     */
    if (trace) {
        logger << "Pre-rec subgame of (" << pr << " pe=" << pe << " po=" << po << "):";
        for (int v=vtop; v>=0; v--) {
            if (SG[v]) {
                if (R[v]) {
                    if (Wm[v]) {
                        if (Wo[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m"; 
                        else if (!H[v]) logger << " \033[1;38;5;46m" << label_vertex(v) << "\033[m";
                        else logger << " \033[38;5;38m" << label_vertex(v) << "\033[m"; 
                    }
                    else logger << " \033[1;38;5;196m" << label_vertex(v) << "\033[m"; 
                }
                else logger << " \033[38;5;160m" << label_vertex(v) << "\033[m"; 
            }
        }
        logger << std::endl;
    }
#endif

    /**
     * Go recursive!
     */
    solve(H, vtop, pr-1, pe, po);
 
    /**
     * Let the opponent attract from our region...
     * The intersection of H and Wo is the opponent's subgame
     */
    bool opponent_attracted_from_us = false;
    for (int v=0; v<=vtop; v++) {
        if (H[v] and Wo[v]) Q.push(v);
    }
    if (Q.nonempty()) {
        while (Q.nonempty()) {
            const int v = Q.pop();
            if (Wm[v]) {
#ifndef NDEBUG
                if (trace and priority(v) == pr) {
                    logger << "\033[1;37mfound distraction\033[m: " << label_vertex(v) << std::endl;
                }
#endif
                opponent_attracted_from_us = true;
            }
            attractVertices(1-pl, v, Wo, R, R);
        }
    }

#ifndef NDEBUG
    /**
     * If so desired, report our current status...
     */
    if (trace) {
        logger << "Subgame of (" << pr << " pe=" << pe << " po=" << po << "):";
        for (int v=vtop; v>=0; v--) {
            if (SG[v]) {
                if (R[v]) {
                    if (Wm[v]) {
                        if (Wo[v]) logger << " \033[38;5;202m" << label_vertex(v) << "\033[m"; 
                        else if (!H[v]) logger << " \033[1;38;5;46m" << label_vertex(v) << "\033[m";
                        else logger << " \033[38;5;38m" << label_vertex(v) << "\033[m"; 
                    }
                    else logger << " \033[1;38;5;196m" << label_vertex(v) << "\033[m"; 
                }
                else logger << " \033[38;5;160m" << label_vertex(v) << "\033[m"; 
            }
        }
        logger << std::endl;
    }
#endif

    /**
     * Finally, if the opponent attracted from us, then we do the right side
     * of the universal tree. Otherwise: update strategies and return.
     * Notice this is also an optimization.
     */
    if (opponent_attracted_from_us) {
        // Reset my winning region <Wm>
        Wm -= R;

        // Reduce subgame by removing opponent's winning region <Wo>
        R -= Wo;

        // Expand the right side of the current level of the universal tree
        if (pl == 0) solve(R, vtop, pr, pe, po/2);
        else solve(R, vtop, pr, pe/2, po);
    } else {
        // Set strategy for vertices that do not yet have a strategy
        for (int v=vtop; v>=0; v--) {
            if (R[v] and Wm[v] and owner(v) == pl and str[v] == -1) {
                auto curedge = outs(v);
                for (int to = *curedge; to != -1; to = *++curedge) {
                    if (Wm[to]) {
                        str[v] = to;
                        break;
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    if (trace >= 2) logger << "OUT: " << pr << " (" << pl << ") E=" << pe << " O=" << po << "\n";
#endif
}

void
ZLKQSolver::run()
{
    iterations = 0;

    str = new int[nodecount()];

    Q.resize(nodecount());
    W0.resize(nodecount());
    W1.resize(nodecount());
    
    bitset G(nodecount());
    G = disabled;
    G.flip();

    solve(G, nodecount()-1, priority(nodecount()-1), nodecount(), nodecount());

#ifndef NDEBUG
    if (trace) {
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
        if (W0[v]) oink->solve(v, 0, str[v]);
        if (W1[v]) oink->solve(v, 1, str[v]);
    }

    logger << "solved with " << iterations << " iterations." << std::endl;

    // check if actually all solved
#ifndef NDEBUG
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { logger << "search was incomplete!" << std::endl; exit(-1); }
    }
#endif

    delete[] str;
}

}
