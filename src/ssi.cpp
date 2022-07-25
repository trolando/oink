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

#include "ssi.hpp"

namespace pg {

SSISolver::SSISolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

SSISolver::~SSISolver()
{
}

/**
 * Returns true if strategy valuation of "a" is less than of "b" for the Even player
 */
bool
SSISolver::si_val_less(int a, int b)
{
    // if a == b, then of course return false
    if (a == b) return false;
    const int *a_val = val + k*a;
    const int *b_val = val + k*b;
    // find highest priority where they differ
    for (int i=k-1; i>=0; i--) {
        const int a_i = a == -1 ? 0 : a_val[i];
        const int b_i = b == -1 ? 0 : b_val[i];
        if (a_i == b_i) continue;
        if (i&1) return a_i > b_i; // for odd priorities
        else return a_i < b_i;     // for even priorities
    }
    // equal
    return false;
}

void
SSISolver::compute_vals_ll(int side)
{
    auto & H = side == 0 ? halt0 : halt1;
    auto & S = side == 0 ? str0 : str1;
    auto & W = side == 0 ? W0 : W1;
    auto & Wopp = side == 0 ? W1 : W0;

    // set all vertices as on cycle
    C = G;

    // find all halting vertices
    std::fill(first_in, first_in+nodecount(), '\xff');
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (W[v]) continue;
        // check if the successor of v is halted
        int s = S[v];
        if (H[s]) {
            Q.push(v);
        } else {
            next_in[v] = first_in[s];
            first_in[s] = v;
        }
    }

    // walk the finite paths
    while (Q.nonempty()) {
        int v = Q.pop();

        // mark as not on a cycle
        C[v] = false;

        // set valuation
        int *val_v = val + k*v;
        int s = S[v];
        if (H[s]) {
            std::fill(val_v, val_v+k, '\x00');
        } else {
            std::copy(val+k*s, val+k*(s+1), val_v);
        }
        val_v[priority(v)]++;

        // add predecessors
        int from = first_in[v];
        while (from != -1) {
            Q.push(from);
            from = next_in[from];
        }
    }
}

void
SSISolver::compute_vals_bw(int side)
{
    auto & H = side == 0 ? halt0 : halt1;
    auto & S = side == 0 ? str0 : str1;
    auto & W = side == 0 ? W0 : W1;

    // set all vertices as on cycle
    C = G;

    // find all halting vertices
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (H[v]) {
            // add all predecessors with the strategy to v
            for (auto curedge = ins(v); *curedge != -1; curedge++) {
                auto from = *curedge;
                if (!G[from]) continue;
                if (S[from] == v) {
                    Q.push(from);
                    C[from] = false;
                }
            }
        }
    }

    // walk the finite paths
    while (Q.nonempty()) {
        int v = Q.pop();

        // set valuation
        int *val_v = val + k*v;
        int s = S[v];
        if (s == -1) LOGIC_ERROR;
        if (H[s]) {
            std::fill(val_v, val_v+k, '\x00');
        } else {
            std::copy(val+k*s, val+k*(s+1), val_v);
        }
        val_v[priority(v)]++;

        // add predecessors
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (!G[from]) continue;
            if (!C[from]) continue; // already seen
            if (S[from] == v) {
                Q.push(from);
                C[from] = false;
            }
        }
    }
}

int
SSISolver::mark_solved(int side)
{
#ifndef NDEBUG
    if (trace >= 2) {
        V = C;
        V -= (side == 0 ? W0 : W1);
        for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
            logger << "\033[1;37mvertex " << label_vertex(v) << " is now won by " << (side == 0 ? "Even" : "Odd") << "!\033[m" << std::endl;
        }
    }
#endif
    if (side == 0) {
        V = C;
        V -= W0;
        W0 = C;
        return V.count();
    } else {
        V = C;
        V -= W1;
        W1 = C;
        return V.count();
    }
}

int
SSISolver::switch_opp_strategy(int side)
{
    int count = 0;

    auto & H = side == 0 ? halt0 : halt1; // side's halt
    auto & S = side == 0 ? str0 : str1;   // side's str

    // V is precomputed all non-won vertices in G, controlled by the opponent
    // for won vertices, we are not interested in computing a better strategy

    for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
        int cur_strat = S[v];
        bool cur_is_C = !H[cur_strat] and C[cur_strat];

        bool changed = false;
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (to == cur_strat) continue; // only check different successors
            if (!G[to]) continue; // not in G
            if (!H[to] and C[to]) continue; // never play to a cycle!
            if (side == 1) {
                // improving for player Even
                if (cur_is_C or si_val_less(H[cur_strat] ? -1 : cur_strat, H[to] ? -1 : to)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdating Even's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
#endif
                    cur_strat = to;
                    cur_is_C = false;
                    changed = true;
                }
            } else {
                // improving for player Odd
                if (cur_is_C or si_val_less(H[to] ? -1 : to, H[cur_strat] ? -1 : cur_strat)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdating Odd's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
#endif
                    cur_strat = to;
                    cur_is_C = false;
                    changed = true;
                }
            }
        }
        if (changed) {
            S[v] = cur_strat;
            count++;
        }
    }

    return count;
}

int
SSISolver::switch_sym_strategy(int side)
{
    int count = 0; // number of changed vertices

    auto & H = side == 0 ? halt0 : halt1; // side's halt
    auto & S = side == 0 ? str0 : str1;   // side's str
    auto & Sbr = side == 0 ? str1 : str0; // best response
    auto & W = side == 0 ? W0 : W1;       // side's winning region

    V = G;
    V -= W0;
    V -= W1;
    V -= side == 0 ? V1 : V0;

    for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
        int cur_strat = S[v];
        if (C[cur_strat]) continue; // do not change a winning strategy

        int br_strat = Sbr[v];
        bool better = false;

        // now check if it is better
        if (cur_strat != br_strat) {
            if (H[br_strat]) {
                if (!H[cur_strat]) {
                    if (side == 0) {
                        if (si_val_less(cur_strat, -1)) better = true;
                    } else {
                        if (si_val_less(-1, cur_strat)) better = true;
                    }
                }
            } else {
                if (C[br_strat]) {
                    better = true;
                } else if (side == 0) {
                    if (si_val_less(H[cur_strat] ? -1 : cur_strat, br_strat)) better = true;
                } else {
                    if (si_val_less(br_strat, H[cur_strat] ? -1 : cur_strat)) better = true;
                }
            }
        }

        if (better) {
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[1;37mupdating " << (side == 0 ? "Even" : "Odd") << "'s strategy\033[m: " << label_vertex(v) << " => " << label_vertex(br_strat) << std::endl;
            }
#endif
            S[v] = br_strat;
            count++;
        } else {
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int to = *curedge;
                if (to == cur_strat) continue; // only check different successors
                if (to == br_strat) continue; // already checked
                if (!G[to]) continue; // not in G
                if (H[to]) {
                    if (H[cur_strat]) continue; // both halting...
                    if (side == 0) {
                        if (si_val_less(cur_strat, -1)) better = true;
                    } else {
                        if (si_val_less(-1, cur_strat)) better = true;
                    }
                } else {
                    if (C[to]) {
                        better = true;
                    } else if (side == 0) {
                        if (si_val_less(H[cur_strat] ? -1 : cur_strat, to)) better = true;
                    } else {
                        if (si_val_less(to, H[cur_strat] ? -1 : cur_strat)) better = true;
                    }
                }
                if (better) {
                    count++;
                    break;
                }
            }
        }
    }

    // update halting
    for (auto v = H.find_first(); v != bitset::npos; v = H.find_next(v)) {
        if (!G[v]) continue; // only vertices in G anyway
        bool good = false;
        if (W[v]) {
            good = true;
        } else if (side == 0) {
            if (si_val_less(-1, v)) good = true;
        } else {
            if (si_val_less(v, -1)) good = true;
        }
        if (good) {
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[1;37m" << (side == 0 ? "Even" : "Odd") << " no longer halts\033[m before vertex " << label_vertex(v) << std::endl;
            }
#endif
            H[v] = false; // stop halting
            count++;
        }
    }

    return count;
}


void
SSISolver::run()
{
    // set G to the game
    G = disabled;
    G.flip();
    if (G.none()) return; // nothing to do

    // assume the game is ordered
    k = 1 + priority(G.find_last());

    Q.resize(nodecount());
    C.resize(nodecount());
    V.resize(nodecount());
    W0.resize(nodecount());
    W1.resize(nodecount());
    halt0.resize(nodecount());
    halt1.resize(nodecount());

    V0 = G;
    V1 = G;
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (owner(v) == 0) V1[v] = false;
        else V0[v] = false;
    }

    val = new int[k*nodecount()];
    str0 = new int[nodecount()];
    str1 = new int[nodecount()];

    first_in = new int[nodecount()];
    next_in = new int[nodecount()];

    // initialize the datastructure
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        // select first available edge
        int to = -1;
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int u = *curedge;
            if (G[u]) {
                to = u;
                break;
            }
        }
        if (to == -1) LOGIC_ERROR;
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "initial strategy: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
        }
#endif
        // set initial strategy of <v> to <to>.
        str0[v] = to;
        str1[v] = to;
        // halt for the right player
        if (priority(v)&1) {
            // odd priority, so halt for even
            halt0[v] = true;
        } else {
            // even priority, so halt for odd
            halt1[v] = true;
        }
    }

    long remaining = G.count();

    for (;;) {
        ++major;
        if (trace) {
            logger << "\033[1;38;5;208mMajor iteration " << major << "\033[m" << std::endl;
        }
        {
            // First compute Odd's best response for Even's strategy
            V = V1;
            V -= W0;
            if (trace) logger << "Computing Odd's best response..." << std::endl;
            for (;;) {
                ++minor;
                compute_vals_ll(0);
                int count = switch_opp_strategy(0);
                if (count == 0) break; // if nothing left, done
            }
            remaining -= mark_solved(0); // mark nodes won by Even after Odd's best response
        }
        {
            // Now compute Even's best response for Odd's strategy
            V = V0;
            V -= W1;
            if (trace) logger << "Computing Even's best response..." << std::endl;
            for (;;) {
                ++minor;
                compute_vals_ll(1);
                int count = switch_opp_strategy(1);
                if (count == 0) break; // if nothing left, done
            }
            remaining -= mark_solved(1); // mark nodes won by Odd after Even's best response
        }
        if (remaining <= 0) break;
        if (trace) logger << "Improving Odd's strategy..." << std::endl;
        if (switch_sym_strategy(1) == 0) {
            // Now all vertices not in a winning region are won by Even with strategy str1
            V = G;
            V -= W0;
            V -= W1;
            W0 |= V;
            for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
                str0[v] = str1[v];
            }
            if (trace) logger << "No improvements - Even wins remainder." << std::endl;
            break;
        }
        if (trace) logger << "Improving Even's strategy..." << std::endl;
        compute_vals_ll(0);
        if (switch_sym_strategy(0) == 0) {
            // Now all vertices not in a winning region are won by Odd with strategy str0
            V = G;
            V -= W0;
            V -= W1;
            W1 |= V;
            for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
                str1[v] = str0[v];
            }
            if (trace) logger << "No improvements - Odd wins remainder." << std::endl;
            break;
        }
    }

    // Now set dominions and derive strategy for odd.
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
#ifndef NDEBUG
        if (W0[v] == W1[v]) LOGIC_ERROR; // in both winning regions?
#endif
        if (W0[v]) {
            oink->solve(v, 0, owner(v) == 0 ? str0[v] : -1);
        } else {
            oink->solve(v, 1, owner(v) == 1 ? str1[v] : -1);
        }
    }

    delete[] val;
    delete[] first_in;
    delete[] next_in;

    logger << "solved with " << major << " major iterations, " << minor << " minor iterations." << std::endl;
}

}
