/*
 * Copyright 2021 Tom van Dijk, University of Twente
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
#include <climits>

#include "ppq.hpp"

namespace pg {

/**
 * This implementation is my take on the idea to use PP to accelerate the ZLKQ family of algorithms.
 *
 * As a base, I take my version of ZLKQ of Lehtinen et al. (2021) see also the ZLKQ solver.
 * I however modified this to use the "r" and "u" arrays of the paper, to make it simpler.
 * In <solve>, the left/right sides of the universal tree follow from Lehtinen et al's version.
 */

PPQSolver::PPQSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

PPQSolver::~PPQSolver()
{
}

/**
 * Attract vertices in <R> from subgame <Y> to <v> in region <Z> as player <pl>.
 * We attract only vertices in <R>, but the opponent may escape to <R> and <Y>.
 * Attracted vertices are added to <Z> and to the queue <Q>.
 */
void
PPQSolver::attractVertices(const int pl, const int v, bitset &Z, bitset &R, bitset &Y)
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
 * Given current region R and lower area L, check if R is closed in L
 * Simultaneously find lowest escape out of R (for the opponent)
 * If closed, promote R and maximize the result
 */
bool
PPQSolver::maybePromote(int pr, int pl, bitset &R, bitset &L)
{
    bool is_closed = true;
    int best_r = INT_MAX; // the best escape towards r
    int best_u = INT_MAX; // the best escape towards u

    // if <R> is empty, there's nothing to do here, just quit
    if (!R.any()) return false;

    // go through all vertices in <R> and check outgoing edges
    for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
        if (owner(v) == pl) {
            // we own it, so just check if we can stay in R
            // this assumes HIGHER regions are all maximized
            is_closed = false;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int to = *curedge;
                if (R[to]) {
                    is_closed = true;
                    break;
                }
            }
            if (!is_closed) return false; // quickly abort
        } else {
            // opponent owns it, so check if it can 'escape' to L
            // again, this assumes HIGHER regions are all maximized
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int to = *curedge;
                if (L[to] && !R[to]) {
                    return false; // quickly abort
                }
                // opponent does not (yet) escape, so update best_r / best_u
                if (!R[to]) {
                    if (r[to] != -1) {
                        best_r = std::min(best_r, r[to]);
                    } else if (u[to] != -1) {
                        best_u = std::min(best_u, u[to]);
                    } else {
#ifndef NDEBUG
                        if (G[to]) LOGIC_ERROR; // this is bad.
#endif
                    }
                }
            }
        }
    }

    // if we are here, it is closed

    if (best_r == INT_MAX && best_u == INT_MAX) {
        // dominion
        if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m" << std::endl;

        // the idea is to add R to the winning region, then attract more from G to the winning region
        // and remove won vertices from G, and from r, and from u.

        auto& W = pl ? W1 : W0;
        W |= R;

        for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
            Q.push(v);
        }

        while (Q.nonempty()) {
            int v = Q.pop();
            G[v] = false;
            r[v] = -1;
            u[v] = -1;
            attractVertices(pl, v, W, G, G);
        }
    } else if (best_u < best_r) {
        if (trace) logger << "\033[1;33mpromoted \033[36m" << pr << " \033[37mto \033[36mu " << best_u << "\033[m" << std::endl;

        // attract R to u, then continue attracting...
        for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
            u[v] = best_u;
            r[v] = -1;
            Q.push(v);
        }

        // just set L to everything that could be attracted now from all lower u/r (below best_u)
        L.reset();
        for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
            if ((r[v] != -1 && r[v] < best_u) || (u[v] != -1 && u[v] < best_u)) L[v] = true;
            if (priority(v) > best_u) break; // only with priorities <= best_u
        }

        /**
         * Since R is everything that was attracted, we can just attract to R from L...
         */
        while (Q.nonempty()) {
            int v = Q.pop();
            u[v] = best_u;
            r[v] = -1;
            attractVertices(pl, v, R, L, L);
        }
    } else {
        if (trace) logger << "\033[1;33mpromoted \033[36m" << pr << " \033[37mto \033[36mr " << best_r << "\033[m" << std::endl;

        // attract R to r, then continue attracting...
        for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
            r[v] = best_r;
            Q.push(v);
        }

        // just set L to everything that could be attracted now from all lower u/r (below best_r)
        L.reset();
        for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
            if ((r[v] != -1 && r[v] < best_r) || (u[v] != -1 && u[v] < best_r)) L[v] = true;
            if (priority(v) > best_r) break; // only with priorities <= best_r
        }

        /**
         * Since R is everything that was attracted, we can just attract to R from L...
         */
        while (Q.nonempty()) {
            int v = Q.pop();
            r[v] = best_r;
            u[v] = -1;
            attractVertices(pl, v, R, L, L);
        }
    }

    // At this point, because of the promotions, it is possible that some vertices are attracted to higher regions...
    // So we perform the "maximize" operation now; it's not precise, but it's sufficient.

    // The maximize operation here considers the entire game.
    // We find the next highest vtop, which is the highest u or r, s.t. u 8 > r 8 > u 7 > r 7
    // Then simply set R to all vertices currently in that region, and attract from the remainder of the game
    // Continue this until the game is empty or the highest r equals pr.

    L = G;
    while (L.any()) {
        // Find highest vtop to attract to
        int vtop = L.find_last();
        for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
            if (r[vtop] != -1) {
                if (r[v] > r[vtop] or u[v] >= r[vtop]) vtop = v;
            } else {
                if (r[v] > u[vtop] or u[v] > u[vtop]) vtop = v;
            }
        }
        if (r[vtop] != -1) {
            // it's an r
            int apr = r[vtop];
            if (apr <= pr) break; // stop attracting below pr, because that results in errors
#ifndef NDEBUG
            if (trace >= 2) logger << "maximizing to r " << r[vtop] << std::endl;
#endif
            R.reset();
            for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
                if (r[v] == apr) {
                    R[v] = true;
                    Q.push(v);
                }
            }
            while (Q.nonempty()) {
                int v = Q.pop();
                r[v] = apr;
                u[v] = -1;
                attractVertices(apr&1, v, R, L, L);
            }
        } else {
            // it's a u
            int apr = u[vtop];
            if (apr < pr) break; // no need to investigate further (but this also probably never occurs)
#ifndef NDEBUG
            if (trace >= 2) logger << "maximizing to u " << u[vtop] << std::endl;
#endif
            R.reset();
            for (auto v = L.find_first(); v != bitset::npos; v = L.find_next(v)) {
                if (u[v] == apr) {
                    R[v] = true;
                    Q.push(v);
                }
            }
            while (Q.nonempty()) {
                int v = Q.pop();
                u[v] = apr;
                r[v] = -1;
                attractVertices(1-(apr&1), v, R, L, L);
            }
        }

        L -= R;
    }

#ifndef NDEBUG
    if (trace >= 3) {
        for (auto v = G.find_last(); v != bitset::npos; v = G.find_prev(v)) {
            logger << label_vertex(v) << " r=" << r[v] << " u=" << u[v] << " str=" << label_vertex(str[v]) << std::endl;
        }
    }
#endif

    return true;
}

void
PPQSolver::solve(int vtop, const int pr, const int prc, const int pe, const int po)
{
    // track the number of visits
    iterations++;

    // no state change
    if (po == 0 || pe == 0) return;

    // move vtop to the true top, then check for game emptiness
    while (vtop >= 0 and (!G[vtop] or r[vtop] > pr or u[vtop] > pr)) vtop--;
    if (vtop == -1) return; // empty game, bye

    // const int pr = r[vtop]; // set to highest priority
    const int pl = pr&1;

    // Expand the left side of the current level of the universal tree
    // This is where we use Lehtinen's version, not Pawel's.
    if (pl == 0) solve(vtop, pr, prc, pe, po/2);
    else solve(vtop, pr, prc, pe/2, po);

#ifndef NDEBUG
    if (trace >= 2) logger << "in pr=" << pr << ", pe=" << pe << ", po=" << po << " (caller pr=" << prc << ")" << std::endl;
#endif

    // after the left side, a lot can be changed, find our new vtop

    while (vtop >= 0 and (!G[vtop] or r[vtop] > pr or u[vtop] > pr)) vtop--;
    if (vtop == -1) return; // empty game, bye

    // TODO handle some special cases: entire region in u; nothing left in r; L is smaller than pe/po

    // set R := all vertices with r[v] == pr
    // set L := all vertices with r[v] <= pr and r[v] != -1
    // also add each vertex in R to the attract queue Q

    R.reset();
    L.reset();
    for (int v = vtop; v >= 0; v--) {
        // assumption: if v not in G, then r[v] == -1
        if (r[v] != -1) {
            if (r[v] == pr) {
                R[v] = true; // add to R
                Q.push(v);   // add to Q
            } 
            if (r[v] <= pr) {
                L[v] = true; // add to L
            }
        }
    } 

    // Attract for pl from L to R

    while (Q.nonempty()) {
        int v = Q.pop();
        r[v] = pr;    // set r ; we already know/assume that u[v] == -1
        // L[v] = false; // remove from L
        attractVertices(pl, v, R, L, L);
    }

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1;33mregion \033[36m" << pr << "\033[m";
        for (auto v = R.find_last(); v != bitset::npos; v = R.find_prev(v)) {
            logger << " \033[38;5;38m" << label_vertex(v) << "\033[m"; 
        }
        logger << std::endl;
    }
#endif

    // Now that we have our R, check if it is closed.
    // Simultaneously find lowest escape out of R (for the opponent)

    bool promoted = maybePromote(pr, pl, R, L);

    // if promoted, we don't need to consider the recursive game (according to the paper)
    if (promoted) return;

    // if (promoted) R.reset(); // alternative...

    /** /
    // optimization to not continue if we can be sure that the recursion already took all
    // possible opponent's small dominions
    //
    // sadly, this does not work well with the PP optimization
    L -= R;
    unsigned long count = L.count(); // what remains ?
    if (pl) {
        if (count <= (po/2)) return;
    } else {
        if (count <= (pe/2)) return;
    }
    // */

    // make a local copy of R, because R/L can be destroyed in the recursion
    bitset Z(R);

    // compute next vtop (for subgame) and the next priority
    int next_vtop = vtop;
    while (next_vtop >= 0 and (!G[next_vtop] or r[next_vtop] >= pr or u[next_vtop] >= pr)) next_vtop--;
    if (next_vtop == -1) return; // subgame is actually empty
    const int subpr = r[next_vtop];

    // run recursively
    solve(vtop, subpr, pr, pe, po);

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "back from recursive " << subpr << ", we are pr " << pr << ", pe=" << pe << ", po=" << po << std::endl;
    }
#endif

    // Now we must first perform Und
    // There may be a lower u. If same pl pr, then merge with our u. If different pl pr, reset.
    // There may be a lower r. If other pl, then merge to our u

    for (int v = vtop; v >= 0; v--) {
        if (u[v] != -1 and u[v] < pr) {
            // lower u
            if ((u[v]&1) == pl) {
                // same pl pr, so merge
                u[v] = pr;
            } else {
                // other pl pr, so reset
                u[v] = -1;
                r[v] = priority(v);
                str[v] = -1;
            }
        } else if ((subpr&1)!=pl and r[v] != -1 and r[v] < pr) {
            // it's a lower other pl r, merge to our u
            r[v] = -1;
            u[v] = pr;
        }
    }

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1mstate after Und\033[m" << std::endl;
        for (auto v = G.find_last(); v != bitset::npos; v = G.find_prev(v)) {
            logger << label_vertex(v) << " r=" << r[v] << " u=" << u[v] << " str=" << label_vertex(str[v]) << std::endl;
        }
    }
#endif

    // after the recursion, R/U/L could be anything...
    // set R := all vertices with r[v] == pr
    // set U := all vertices with u[v] == pr
    // set L := all vertices with r[v] <= pr
    // also add each vertex in U to the attract queue Q

    R.reset();  
    L.reset();
    U.reset();
    for (int v = vtop; v >= 0; v--) {
        if (r[v] == pr) {
            R[v] = true; // add to R
        }
        if (r[v] != -1 && r[v] <= pr) {
            L[v] = true; // add to L
        }
        if (u[v] == pr) {
            U[v] = true; // add to U
            Q.push(v);   // add to Q for attracting from R/L to U
        }
    } 

    // Attract to U from L
    while (Q.nonempty()) {
        int v = Q.pop();
        u[v] = pr;
        r[v] = -1;
        attractVertices(1-pl, v, U, L, L);
    }

    // Update R and L
    R -= U;
    L -= U;

    // Add original R to current R in order to check if we need to reset R
    R |= Z;

    // check if we need to reset the region
    bool reset_R = false;
    for (auto v = R.find_first(); v != bitset::npos; v = R.find_next(v)) {
        // check if a vertex is now in a dominion or in u
        if (r[v] == -1) {
            reset_R = true; // either dominion or in u
            break;
        }
        // check if a vertex is in a higher opponent r
        if (r[v] > pr && (r[v]&1) != pl) {
            reset_R = true;
            break;
        }
        // check if a strategy leaves R
        if (str[v] != -1 && r[str[v]] != pr) {
            reset_R = true;
            break;
        }
    }

    if (reset_R) {
        // reset everything with r <= pr
        for (auto v = vtop; v >= 0; v--) {
            if (r[v] != -1 && r[v] <= pr) {
#ifndef NDEBUG
                if (trace >= 2) logger << "reset " << label_vertex(v) << std::endl;
#endif
                r[v] = priority(v);
                u[v] = -1;
                str[v] = -1;
            }
        }
    } else {
        // we did not reset, try to promote again
        maybePromote(pr, pl, R, L);
        // regardless of the outcome, if nothing reset, no need to go recursive again
        return;
    }

    // move vtop then expand right side of the universal tree
    while (vtop >= 0 and (!G[vtop] or r[vtop] > pr or u[vtop] > pr)) vtop--;
    if (vtop == -1) return; // and we're done, there's nothing left

    // Expand the right side of the current level of the universal tree
    if (pl == 0) solve(vtop, pr, prc, pe, po/2);
    else solve(vtop, pr, prc, pe/2, po);
}

void
PPQSolver::run()
{
    iterations = 0;

    str = new int[nodecount()];
    std::fill(str, str+nodecount(), -1);

    r = new int[nodecount()];
    for (int v=0; v<nodecount(); v++) {
        if (!disabled[v]) r[v] = priority(v);
        else r[v] = -1;
    }

    u = new int[nodecount()];
    std::fill(u, u+nodecount(), -1);

    Q.resize(nodecount());
    L.resize(nodecount());
    R.resize(nodecount());
    U.resize(nodecount());

    W0.resize(nodecount());
    W1.resize(nodecount());
    
    G.resize(nodecount()); 
    G = disabled;
    G.flip();

    if (G.any()) {
        solve(nodecount()-1, r[G.find_last()], INT_MAX, nodecount(), nodecount()); // caller has INT_MAX priority...
    }

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
        if (W1[v]) Solver::solve(v, 1, str[v]);
    }

    logger << "solved with " << iterations << " iterations." << std::endl;

    // check if actually all solved
#ifndef NDEBUG
    for (int i=0; i<nodecount(); i++) {
        if (!disabled[i]) { THROW_ERROR("search was incomplete!"); }
    }
#endif

    delete[] str;
}

}
