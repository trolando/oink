/*
 * Copyright 2019 Tom van Dijk, Johannes Kepler University Linz
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
    const int *_in = ins + ina[v];
    for (int from = *_in; from != -1; from = *++_in) {
        if (Z[from]) {
            // already in Z, set strategy if not yet set
            if (owner[from] == pl and str[from] == -1) {
#ifndef NDEBUG
                if (trace >= 3) {
                    logger << "\033[1;37msetting strategy of \033[36m" << label_vertex(from) << "\033[m to \033[1;36m" << label_vertex(v) << "\033[m" << std::endl;
                }
#endif
                str[from] = v;
            }
        } else if (R[from]) {
            // a vertex in <R> that is not yet in <Z>
            if (owner[from] != pl) {
                // check each exit
                bool escapes = false;
                const int *_out = outs + outa[from];
                for (int to = *_out; !escapes and to != -1; to = *++_out) {
                    // check if escapes to a vertex in <R> or <Y> that is not in <Z>
                    if (!Z[to] and (R[to] or Y[to])) escapes = true;
                }
                if (escapes) continue;
            }
            // vertex does not escape, so attract it
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

void
ZLKQSolver::solve(bitset &SG, int vtop, int pr, int pe, int po)
{
    // find actual top
    while (true) {
        if (vtop == -1) return; // empty game, bye
        if (SG[vtop]) break;  
        vtop--;
    }

    // uncomment the following line to use the priority of the next top
    // pr = priority[vtop];
    const int pl = pr&1;

    if (pl == 1 and po <= 0) return; // out of <po>
    if (pl == 0 and pe <= 0) return; // out of <pe>

    iterations++;

    auto &Wm = pl == 0 ? W0 : W1;
    auto &Wo = pl == 0 ? W1 : W0;

    int _pe = pl == 0 ? pe/2 : pe;
    int _po = pl == 1 ? po/2 : po;

    bitset R(n_nodes); // remaining subgame (initially SG)
    bitset Z(n_nodes); // current region (initially empty)
    R = SG;

    int phase = 1;
    int top = -1;

    while (true) {
        if (phase != 2) {
            // compute the top region Z (attractor to nodes of priority <pr>)
            for (top=vtop; top != -1; top--) {
                if (!R[top]) continue;
                if (priority[top] != pr) break;
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
            logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
            for (int v=vtop; v>=0; v--) {
                if (Z[v]) logger << " \033[1;38;5;15m" << label_vertex(v) << "\033[m";
            }
            logger << std::endl;
#endif
        }
        bool subgame_has_opponent = false;
        if (R.any()) {
            // ensure W0 and W1 are reset for the remaining subgame
            W0 -= R;
            W1 -= R;
            // solve the remaining subgame with the correct precision
            if (phase == 2) {
                if (trace) {
                    logger << "\033[1;33mrunning full precision below \033[1;36m" << pr << "\033[m" << std::endl;
                }
                solve(R, top, pr-1, pe, po); // solve once with full precision
            } else {
                solve(R, top, pr-1, _pe, _po); // otherwise solve with half precision
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
            if (phase == 2) phase = 3; // only do full precision once
            R = SG;
            Z.reset();
            continue;
        }
        // if we are here, then the recursive call did not have opponent's won vertices
        if (phase == 1) {
            // ok we are done with the first phase, now do it once with full precision
            // we don't actually need to recompute R/Z for the full precision one
            phase = 2;
            continue;
        } // else if (phase == 2) {
            // phase = 3;
            // continue;
            // uncomment the above lines in order to not skip the (redundant) phase 3 call
            // the idea is that if opponent does not win anything in the full precision call, then
            // it will certainly not win anything in the half precision call
            // indeed, in Pawel's algorithm, this half precision call is skipped too
        //}
        // if we are here, then whatever remains is won by us; set strategy correctly
        if (Z.any()) {
            if (trace) {
                logger << "marking region " << pr << " as won by " << pl << " (" << Z.count() << " vertices)\n";
            }
            for (int v=vtop; v>top; v--) {
                if (!Z[v]) continue;
                if (owner[v] == pl and str[v] == -1) {
                    for (const int to : out[v]) {
                        if (SG[to]) {
                            str[v] = to;
                            break;
                        }
                    }
                }
            }
            Wm |= Z;
        }
        return;
    }
}

void
ZLKQSolver::run()
{
    iterations = 0;

    str = new int[n_nodes];

    Q.resize(n_nodes);
    W0.resize(n_nodes);
    W1.resize(n_nodes);
    
    bitset G(n_nodes);
    G = disabled;
    G.flip();

    solve(G, n_nodes-1, priority[n_nodes-1], n_nodes, n_nodes);

    if (trace) {
        for (int v=0; v<n_nodes; v++) {
            if (disabled[v]) continue;
            logger << "vertex " << label_vertex(v) << " is solved by";
            if (W0[v]) {
                logger << " even";
                if (owner[v] == 0) {
                    logger << " (";
                    if (str[v] == -1) logger << "-1";
                    else logger << label_vertex(str[v]);
                    logger << ")";
                }
            }
            if (W1[v]) {
                logger << " odd";
                if (owner[v] == 1) {
                    logger << " (";
                    if (str[v] == -1) logger << "-1";
                    else logger << label_vertex(str[v]);
                    logger << ")";
                }
            }
            logger << std::endl;
        }
    }

    for (int v=0; v<n_nodes; v++) {
        if (disabled[v]) continue;
        if (W0[v]) oink->solve(v, 0, str[v]);
        if (W1[v]) oink->solve(v, 1, str[v]);
    }

    logger << "solved with " << iterations << " iterations." << std::endl;

    // check if actually all solved
#ifndef NDEBUG
    for (int i=0; i<n_nodes; i++) {
        if (!disabled[i]) { logger << "search was incomplete!" << std::endl; exit(-1); }
    }
#endif

    delete[] str;
}

}
