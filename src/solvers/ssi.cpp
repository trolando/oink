/*
 * Copyright 2021-2022 Tom van Dijk
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

SSISolver::SSISolver(Oink& oink, Game& game) : Solver(oink, game)
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

    // set all vertices in G to "on cycle"
    C = G;

    // find all halting vertices and set the linked list of in edges
    std::fill(first_in, first_in+nodecount(), '\xff');
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        int s = S[v];
        if (H[s]) {
            // successor of v is halted
            Q.push(v);
        } else {
            // successor of v is not halted
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

int
SSISolver::mark_solved(int side)
{
    if (C.none()) return 0;

    auto & W = side == 0 ? W0 : W1;
    auto & S = side == 0 ? str0 : str1;

    // the vertices in C are the new winning region (dominion)
    W |= C;

    // attract to the winning region
    int count = 0;
    for (auto v = C.find_first(); v != bitset::npos; v = C.find_next(v)) {
        Q.push(v);
    }
    while (Q.nonempty()) {
        auto v = Q.pop();
        count++;
        for (auto curedge = ins(v); *curedge != -1; curedge++) {
            auto from = *curedge;
            if (G[from] and !W[from]) {
                if (owner(from) != side) {
                    // check if escapes
                    bool escapes = false;
                    for (auto curedge = outs(from); *curedge != -1; curedge++) {
                        auto to = *curedge;
                        if (G[to] && !W[to]) {
                            escapes = true;
                            break;
                        }
                    }
                    if (escapes) continue;
                } else {
                    S[from] = v;
                }
                C[from] = true;
                W[from] = true;
                Q.push(from);
            }
        }
    }

#ifndef NDEBUG
    if (trace >= 2) {
        for (auto v = C.find_first(); v != bitset::npos; v = C.find_next(v)) {
            if (owner(v) == side) {
                logger << "\033[1;37mvertex " << label_vertex(v) << " is now won by " << (side == 0 ? "Even" : "Odd") << " with strategy " << label_vertex(S[v]) << "!\033[m" << std::endl;
            } else {
                logger << "\033[1;37mvertex " << label_vertex(v) << " is now won by " << (side == 0 ? "Even" : "Odd") << "!\033[m" << std::endl;
            }
        }
    }
#endif

    // remove maximized winning region from remaining game
    G -= W;
    V0 &= G;
    V1 &= G;

    // partially reset the strategies and reset halting
    halt0.reset();
    halt1.reset();
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        if (!G[str0[v]] || !G[str1[v]]) {
            // select first available edge
            int to = -1;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int u = *curedge;
                if (G[u]) {
                    to = u;
                    break;
                }
            }
            if (!G[str0[v]]) str0[v] = to;
            if (!G[str1[v]]) str1[v] = to;
        }

        if (priority(v)&1) {
            halt0[v] = true;
        } else {
            halt1[v] = true;
        }
    }

    return count;
}

int
SSISolver::switch_opp_strategy(int side)
{
    int count = 0;

    auto & H = side == 0 ? halt0 : halt1; // side's halt
    auto & S = side == 0 ? str0 : str1;   // side's str
    auto & V = side == 0 ? V1 : V0;       // opponent vertices in G

    for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
        int cur_strat = S[v];
        bool cur_is_G = G[cur_strat];
        bool cur_is_C = !H[cur_strat] and C[cur_strat];
        bool changed = false;

        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            if (to == cur_strat) {
                // only check different successors
            } else if (!G[to]) {
                // only check successors in G
            } else if (!cur_is_G) {
                // if current strategy is outside G, take anything
#ifndef NDEBUG
                if (trace >= 2) {
                    if (side == 1) {
                        logger << "\033[1;37mupdating Even's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    } else {
                        logger << "\033[1;37mupdating Odd's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
                }
#endif
                cur_strat = to;
                cur_is_G = true;
                cur_is_C = !H[cur_strat] and C[cur_strat];
                changed = true;
            } else if (!H[to] and C[to]) {
                // never play to a cycle!
            } else if (cur_is_C) {
                // if current strategy is to a cycle, take anything
#ifndef NDEBUG
                if (trace >= 2) {
                    if (side == 1) {
                        logger << "\033[1;37mupdating Even's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    } else {
                        logger << "\033[1;37mupdating Odd's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
                }
#endif
                cur_strat = to;
                cur_is_C = false;
                changed = true;
            } else if (side == 1) {
                // improving for player Even?
                if (si_val_less(H[cur_strat] ? -1 : cur_strat, H[to] ? -1 : to)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdating Even's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
#endif
                    cur_strat = to;
                    changed = true;
                }
            } else {
                // improving for player Odd?
                if (si_val_less(H[to] ? -1 : to, H[cur_strat] ? -1 : cur_strat)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdating Odd's best response\033[m: " << label_vertex(v) << " => " << label_vertex(to) << std::endl;
                    }
#endif
                    cur_strat = to;
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
    auto & Sbr = side == 0 ? str1 : str0; // side's best response
    auto & V = side == 0 ? V0 : V1;       // side's controlled vertices in G

    for (auto v = V.find_first(); v != bitset::npos; v = V.find_next(v)) {
        int cur_strat = S[v];
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
            // best response is not improving, check if there exists an improving edge
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                int to = *curedge;
                if (to == cur_strat) {
                    // only check different successors
                } else if (to == br_strat) {
                    // already checked
                } else if (!G[to]) {
                    // only check successors in G
                } else if (H[to]) {
                    if (H[cur_strat]) continue; // both halting...
                    if (side == 0) {
                        if (si_val_less(cur_strat, -1)) count++;
                    } else {
                        if (si_val_less(-1, cur_strat)) count++;
                    }
                } else if (C[to]) {
                    count++;
                } else if (side == 0) {
                    if (si_val_less(H[cur_strat] ? -1 : cur_strat, to)) count++;
                } else {
                    if (si_val_less(to, H[cur_strat] ? -1 : cur_strat)) count++;
                }
            }
        }
    }

    // update halting
    for (auto v = H.find_first(); v != bitset::npos; v = H.find_next(v)) {
        if (!G[v]) continue; // only vertices in G anyway
        bool good = false;
        if (side == 0) {
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
    V0.resize(nodecount());
    V1.resize(nodecount());

    val = new int[k*nodecount()];
    str0 = new int[nodecount()];
    str1 = new int[nodecount()];

    first_in = new int[nodecount()];
    next_in = new int[nodecount()];

    // initialize the datastructure
    for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
        // set V0 or V1
        (owner(v) == 0 ? V0 : V1)[v] = true;
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

    int minor = 0, major = 0;
    auto remaining = G.count();

    for (;;) {
        ++major;
        if (trace) {
            logger << "\033[1;38;5;208mMajor iteration " << major << "\033[m" << std::endl;
        }
        {
            // First compute Odd's best response for Even's strategy
            if (trace) logger << "Computing Odd's best response..." << std::endl;
            for (;;) {
                ++minor;
                compute_vals_ll(0);
                int count = switch_opp_strategy(0);
                if (count == 0) break; // if nothing left, done
            }
            int solved = mark_solved(0); // mark nodes won by Even after Odd's best response
            if (solved != 0) {
                remaining -= solved;
                if (remaining <= 0) break;
                continue; // restart major loop
            }
        }
        {
            // Now compute Even's best response for Odd's strategy
            if (trace) logger << "Computing Even's best response..." << std::endl;
            for (;;) {
                ++minor;
                compute_vals_ll(1);
                int count = switch_opp_strategy(1);
                if (count == 0) break; // if nothing left, done
            }
            int solved = mark_solved(1); // mark nodes won by Odd after Even's best response
            if (solved != 0) {
                remaining -= solved;
                if (remaining <= 0) break;
                continue; // restart major loop
            }
        }
        if (trace) logger << "Improving Odd's strategy..." << std::endl;
        if (switch_sym_strategy(1) == 0) {
            // Now all vertices not in a winning region are won by Even with strategy str1
            if (trace) logger << "No (possible) improvements - Even wins remainder." << std::endl;
            W0 |= G;
            for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
                str0[v] = str1[v];
            }
            break;
        }
        if (trace) logger << "Improving Even's strategy..." << std::endl;
        compute_vals_ll(0);
        if (switch_sym_strategy(0) == 0) {
            // Now all vertices not in a winning region are won by Odd with strategy str0
            if (trace) logger << "No (possible) improvements - Odd wins remainder." << std::endl;
            W1 |= G;
            for (auto v = G.find_first(); v != bitset::npos; v = G.find_next(v)) {
                str1[v] = str0[v];
            }
            break;
        }
    }

    for (auto v = W0.find_first(); v != bitset::npos; v = W0.find_next(v)) {
        Solver::solve(v, 0, owner(v) == 0 ? str0[v] : -1);
    }

    for (auto v = W1.find_first(); v != bitset::npos; v = W1.find_next(v)) {
        Solver::solve(v, 1, owner(v) == 1 ? str1[v] : -1);
    }

    delete[] val;
    delete[] first_in;
    delete[] next_in;
    delete[] str0;
    delete[] str1;

    logger << "solved with " << major << " major iterations, " << minor << " minor iterations." << std::endl;
}

}
