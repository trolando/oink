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

#include "zlkq.hpp"

namespace pg {

/**
 * This is an implementation of Pawel Parys's Zielonka variant
 * as presented by himself on his website.
 * https://www.mimuw.edu.pl/~parys/publications/2018-parity-algorithm.pdf
 *
 * I modified the algorithm to produce a strategy as well.
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
ZLKQSolver::solve(bitset &SG, int vtop, const int pr, const int pe, const int po)
{
#ifndef NDEBUG
    if (trace >= 3) {
        logger << "solve(" << pr << " pE=" << pe << " pO=" << po << ")" << std::endl;
    }
#endif

    // find actual top
    while (true) {
        if (vtop == -1) return; // empty game, bye
        if (SG[vtop]) break;  
        vtop--;
    }

    // uncomment the following line to use the priority of the next top
    // pr = priority(vtop);
    const int pl = pr&1;

    if (pl == 1 and po <= 1) return; // out of <po>
    if (pl == 0 and pe <= 1) return; // out of <pe>

    iterations++;

    auto &Wm = pl == 0 ? W0 : W1; // my W
    auto &Wo = pl == 0 ? W1 : W0; // opponent's W

    bitset R(nodecount()); // remaining subgame (initially SG)
    bitset Z(nodecount()); // current region (initially empty)
    R = SG;

    int phase = 1;
    int top = -1;

#ifndef NDEBUG
    int debug_counter = 0;
#endif

    while (true) {
#ifndef NDEBUG
        if (debug_counter++ > nodecount()+1) {
            logger << "logic error, debug counter too large\n";
            exit(-1);
        }
#endif
        // compute the top region Z (attractor to nodes of priority <pr>)
        for (top=vtop; top != -1; top--) {
            if (!R[top]) continue;
            if (priority(top) != pr) break;
            Z[top] = true;
            str[top] = -1;
            Q.push(top);
            while (Q.nonempty()) {
                const int v = Q.pop();
                R[v] = false;
                attractVertices(pl, v, Z, R, R);
            }
        }
#ifndef NDEBUG
        // report
        if (trace >= 2) {
            logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
            for (int v=vtop; v>=0; v--) {
                if (Z[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
        }
#endif
        bool subgame_has_opponent = false;
        if (R.any()) {
            // ensure W0 and W1 are reset for the remaining subgame
            W0 -= R;
            W1 -= R;
            // solve the remaining subgame with the correct precision
            if (phase == 2) {
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1;33mrunning full precision below \033[1;36m" << pr << "\033[m" << std::endl;
                }
#endif
                solve(R, top, pr-1, pe, po); // solve once with full precision
                phase = 3; // only once
            } else if (pl == 0) {
                solve(R, top, pr-1, pe, po/2); // solve with half precision
            } else /* pl == 1 */ {
                solve(R, top, pr-1, pe/2, po); // solve with half precision
            }
            subgame_has_opponent = SG.intersects(Wo); // some vertices in SG are now won by opponent
        }
        if (subgame_has_opponent) {
            // add all subgame vertices in Wo to the queue for attracting
            for (int v=vtop; v>=0; v--) {
                if (SG[v] and Wo[v]) Q.push(v);
            }
            // attract from subgame to Wo and remove from subgame
            while (Q.nonempty()) {
                const int v = Q.pop();
                SG[v] = false; // remove from the current subgame
                Wm[v] = false; // because anything now won by opponent is not won by me
                attractVertices(1-pl, v, Wo, SG, SG);
            }
            // reset R and Z for recomputing in next part
            R = SG;
            Z.reset();
            continue;
        }
        // if we are here, then the recursive call did not have opponent's won vertices
        if (phase == 1) {
            // ok we are done with the first phase, now do it once with full precision
            phase = 2;
            R = SG;
            Z.reset();
            continue;
        }
        // if we are here, then whatever remains is won by us; set strategy correctly
#ifndef NDEBUG
        if (trace) {
            if (SG.any()) {
                logger << "marking remaining game <= " << pr << " as won by " << pl << " (" << SG.count() << " vertices)";
                if (trace >= 2) {
                    for (int v=vtop; v>=0; v--) {
                        if (SG[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
                    }
                }
                logger << std::endl;
            }
        }
#endif
        // any top vertices in Z might have no strategy, so set if needed
        if (Z.any()) {
            for (int v=vtop; v>top; v--) {
                if (!Z[v]) continue;
                if (owner(v) == pl and str[v] == -1) {
                    auto curedge = outs(v);
                    for (int to = *curedge; to != -1; to = *++curedge) {
                        if (SG[to]) {
                            str[v] = to;
                            break;
                        }
                    }
                }
            }
        }
        Wm |= SG;
        return;
    }
}

void
ZLKQSolver::solve2(bitset &SG, int vtop, int pl, int pr, int pe, int po)
{
    /**
     * Record number of calls
     */
    iterations++;

    /**
     * Check if we have run out of precision/depth!
     * If our opponent's precision is gone, assume it's all won by us.
     * If our own precision is gone, assume it's all won by the opponent.
     */

    if (po <= 1) {
        W0 |= SG;
        return;
    }
    if (pe <= 1) {
        W1 |= SG;
        return;
    }

    /**
     * Check if the game is empty and move vtop to true top
     */

    while (true) {
        if (vtop == -1) return; // empty game, bye
        if (SG[vtop]) break; // found a vertex in the subgame 
        vtop--;
    }

    /**
     * Optimization: set pr/pl to top priority/player
     */

    pr = priority(vtop);
    pl = pr&1;

    if (trace) logger << "IN: " << pr << " (" << pl << ") E=" << pe << " O=" << po << "\n";

    /**
     * Expand the left side of the current level of the universal tree
     */

    if (pl == 0) solve2(SG, vtop, pl, pr, pe, po/2);
    else solve2(SG, vtop, pl, pr, pe/2, po);

    /**
     * Now compute <R> := vertices in subgame that I won in the left side
     * Simultaneously, reduce vtop and compute size
     */

    auto &Wm = pl == 0 ? W0 : W1; // my W
    auto &Wo = pl == 0 ? W1 : W0; // opponent's W

    bitset R(nodecount());

    bool first = true;
    unsigned long Rsize = 0;
    for (int v=vtop; v>=0; v--) {
        if (SG[v] and Wm[v]) {
            if (first) {
                vtop=v;
                first = false;
            }
            R[v] = 1;
            Rsize++;
        }
    }

    /**
     * Compute H := Attr(vertices of <pr> in <R>)
     */
    bitset H(nodecount());

    for (int v=vtop; v != -1; v--) {
        if (!R[v]) continue;
        if (priority(v) != pr) break; // stop attracting ... ;) [todo: add compression otf]
        Hs.push(v);
        H[v] = true;
        str[v] = -1;
        Q.push(v);
        while (Q.nonempty()) {
            attractVertices(pl, Q.pop(), H, R, R);
        }
    }

    /**
     * Check if the attractor is closed
     */
    if (Hs.nonempty()) {
        bool open = false;
        while (Hs.nonempty()) {
            const int v = Hs.pop();
            if (owner(v) == pl) {
                if (str[v] == -1) open = true;
            } else {
                /* find ONE successor */
                auto curedge = outs(v);
                for (int to = *curedge; to != -1; to = *++curedge) {
                    if (R[to] and !H[to]) {
                        open = true;
                        break;
                    }
                }
            }
            if (open) {
                Hs.clear();
                break;
            }
        }
        if (!open) {
            //logger << "Region is closed!!" << std::endl;
            /*for (int v=0; v<=vtop; v++) {
                if (R[v]) {
                    Wo[v] = 0;
                    Wm[v] = 1;
                }
            }
            return;*/
        }
    }
                

    /**
     * Reset Wo and Wm for vertices in <R> - <H>;
     * Set <H> := <R> - <H>
     */
    for (int v=0; v<=vtop; v++) {
        if (R[v]) {
            Wo[v] = 0;
            Wm[v] = H[v];
            H[v] = !H[v];
        }
    }

    // if (pe > vtop) pe = vtop;
    // if (po > vtop) po = vtop;

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

    /**
     * Go recursive!
     */
    solve2(H, vtop, 1-pl, pr-1, pe, po);
 
    /**
     * Let the opponent attract from our region...
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
                if (trace and priority(v) == pr)
                    logger << "\033[1;37mfound distraction\033[m: " << label_vertex(v) << std::endl;
#endif
                opponent_attracted_from_us = true;
            }
            attractVertices(1-pl, v, Wo, R, R);
        }
    }

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

    /**
     * Finally, if the opponent attracted from us, then we do the right side
     * of the universal tree. Otherwise: update strategies and return.
     */
    if (opponent_attracted_from_us) {
        /**
         * Reduce <R> to the remaining region.
         * Unset Wm in <R>
         */
        for (int v=vtop; v>=0; v--) {
            if (R[v]) {
                if (Wo[v]) R[v] = 0;
                if (Wm[v]) Wm[v] = 0;
            }
        }

        /**
         * Expand the right side of the current level of the universal tree
         */

        if (pl == 0) solve2(R, vtop, pl, pr, pe, po/2);
        else solve2(R, vtop, pl, pr, pe/2, po);
    } else {
        // set strategy for vertices that do not yet have a strategy
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
    if (trace) logger << "OUT: " << pr << " (" << pl << ") E=" << pe << " O=" << po << "\n";
}

void
ZLKQSolver::run()
{
    iterations = 0;

    str = new int[nodecount()];

    Q.resize(nodecount());
    Hs.resize(nodecount());
    W0.resize(nodecount());
    W1.resize(nodecount());
    
    bitset G(nodecount());
    G = disabled;
    G.flip();

#if 1
    solve(G, nodecount()-1, priority(nodecount()-1), nodecount(), nodecount());
#else
    solve2(G, nodecount()-1, priority(nodecount()-1)&1, priority(nodecount()-1), nodecount(), nodecount());
#endif

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
