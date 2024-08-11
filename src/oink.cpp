/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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
#include <queue>
#include <stack>
#include <iomanip>
#include <iostream>
#include <chrono>

#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "oink/solver.hpp"
#include "lace.h"

namespace pg {

Oink::Oink(Game &game, std::ostream &out) : game(&game), logger(out), todo(game.vertexcount()), disabled(game.getSolved())
{
    // ensure the vertices are ordered properly
    game.ensure_sorted();
    // ensure arrays are built, but don't rebuild
    game.build_in_array(false);

    // initialize outcount (for flush)
    outcount = new int[game.vertexcount()];
    for (int i=0; i<game.vertexcount(); i++) {
        outcount[i] = 0;
        const int *ptr = game.outedges() + game.firstout(i);
        for (int to = *ptr; to != -1; to = *++ptr) {
            if (!disabled[to]) outcount[i]++; // only count the non-disabled subgame
        }
    }
}

Oink::~Oink()
{
    delete[] outcount;
}

/**
 * Find all SCCs at priority p, game limited to nodes with priority <= p
 * Then for every SCC that contains edges and a node with priority p, can win
 * with priority p and is therefore a counterexample to the strategy.
 */
int
Oink::solveTrivialCycles()
{
    // Record number of trivial cycles
    int count = 0;

    // Allocate and initialize datastructures
    const int n_nodes = game->vertexcount();
    int *done = new int[n_nodes];
    int64_t *low = new int64_t[n_nodes];
    for (int i=0; i<n_nodes; i++) done[i] = disabled[i] ? -2 : -1;
    for (int i=0; i<n_nodes; i++) low[i] = 0;

    std::vector<int> res;
    std::vector<int> scc;
    std::stack<int> st;
    std::queue<int> q;

    int64_t pre = 0;

    for (int i=n_nodes-1; i>=0; i--) {
        if (disabled[i]) continue;
        if (done[i] == -2) continue; // already know to skip this

        /**
         * We're going to search all winner-controlled SCCs reachable from node <i>
         */
        const int pr = game->priority(i);
        const int pl = pr & 1;

        /**
         * Only start at unsolved winner-controlled and not yet seen for this priority
         */
        if (game->owner(i) != pl) {
            done[i] = -2; // highest priority but not winner-controlled, never check again
            continue;
        }

        if (done[i] == pr) continue; // already visited in this search

        /**
         * Set bot to pre...
         */
        int64_t bot = pre;

        /**
         * Start forward DFS at <i>
         */
        st.push(i);
        while (!st.empty()) {
            int idx = st.top();

            /**
             * When we see it for the first time, we assign the next number to it and add it to <res>.
             */
            if (low[idx] <= bot) {
                low[idx] = ++pre;
                if (pre < 0) LOGIC_ERROR; // overflow on a 64-bit integer...
                res.push_back(idx);
            }

            /**
             * Now we check all outgoing (allowed) edges.
             * If seen earlier, then update "min"
             * If new, then 'recurse'
             */
            int min = low[idx];
            bool pushed = false;

            auto ptr = game->outedges() + game->firstout(idx);
            for (int to = *ptr; to != -1; to = *++ptr) {
                /**
                 * Only go to lower priority nodes, controlled by <pl>, that are not excluded or seen this round.
                 */
                if (disabled[i]) continue;
                if (to > i or done[to] == -2 or done[to] == pr or game->owner(to) != pl) continue;
                if (low[to] <= bot) {
                    // not visited, add to <st> and break!
                    st.push(to);
                    pushed = true;
                    break;
                } else {
                    // visited, update min
                    if (low[to] < min) min = low[to];
                }
            }
            if (pushed) continue; // we pushed...

            /**
             * If we're here, then there was no new edge and we check if we're the root of an SCC
             */
            if (min < low[idx]) {
                // not the root
                low[idx] = min;
                st.pop();
                continue;
            }

            /**
             * We're the root of an scc. Move the scc from <res> to <scc>.
             * Record highest prio, highest prio of good parity, highest node of good parity.
             */
            int max_pr = -1, max_pr_pl = -1, max_pr_n = -1;
            for (;;) {
                if (res.empty()) LOGIC_ERROR;
                int n = res.back();
                res.pop_back();
                scc.push_back(n);
                done[n] = pr; // dont check again this round
                if (low[n] != min) low[n] = min; // set it [for strategy]
                int d = game->priority(n);
                if (d > max_pr) max_pr = d;
                if ((d&1) == pl and d > max_pr_pl) {
                    max_pr_pl = d;
                    max_pr_n = n;
                }
                if (n == idx) break; // end when done
            }

            /**
             * Check if a single-node SCC without a self-loop
             */
            if (scc.size() == 1) {
                if (!game->has_edge(idx, idx)) {
                    // no self-loop
                    done[idx] = -2; // never check again
                    scc.clear();
                    st.pop();
                    continue;
                }
            }

            /**
             * Check if the highest priority in the SCC is actually won by <pl>.
             */

            if ((max_pr & 1) != pl) {
                for (auto n : scc) if (game->priority(n) > max_pr_pl) done[n] = -2; // never check again
                scc.clear();
                st.pop();
                continue;
                // Note that this SCC will be found again in lower runs, but without offending nodes.
            }

            /**
             * OK, got a winner!
             */

            if (trace) {
                logger << "winner-controlled scc with win priority \033[1;34m" << max_pr << "\033[m" << std::endl;
            }

            // Set strategies for all nodes in the SCC via backward search
            q.push(max_pr_n);
            while (!q.empty()) {
                int cur = q.front();
                q.pop();
                auto ptr = game->inedges() + game->firstin(cur);
                for (int from = *ptr; from != -1; from = *++ptr) {
                    if (low[from] != min or disabled[from]) continue;
                    solve(from, pl, cur); // also sets "disabled"
                    q.push(from);
                }
            }
            flush(); // this triggers attraction of every winner-controlled node in the stack...

            /**
             * For obvious reasons, all nodes on the stack are now solved.
             */
            while (!st.empty()) st.pop();
            res.clear();
            scc.clear();
            count++;
        }
    }

    delete[] done;
    delete[] low;
    return count;
}

int
Oink::solveSelfloops()
{
    int res = 0;
    for (int v=0; v<game->vertexcount(); v++) {
        if (disabled[v]) continue; // skip vertices that are hidden

        if (game->has_edge(v, v)) {
            if (game->owner(v) == (game->priority(v)&1)) {
                // a winning selfloop
                if (trace) logger << "winning self-loop with priority \033[1;34m" << game->priority(v) << "\033[m" << std::endl;
                solve(v, game->owner(v), v);
                res++;
            } else {
                // a losing selfloop
                if (game->outcount(v) == 1) {
                    // it is a losing dominion
                    solve(v, 1 - game->owner(v), -1);
                    res++;
                }
            }
        }
    }

    if (res != 0) flush();
    return res;
}

bool
Oink::solveSingleParity()
{
    int parity = -1;
    for (int v=0; v<game->vertexcount(); v++) {
        if (disabled[v]) continue;
        if (parity == -1) {
            parity = game->priority(v)&1;
        } else if (parity == (game->priority(v)&1)) {
            continue;
        } else {
            return false;
        }
    }
    if (parity == 0 or parity == 1) {
        // solved with random strategy
        logger << "parity game only has parity " << (parity ? "odd" : "even") << std::endl;
        for (int v=0; v<game->vertexcount(); v++) {
            if (disabled[v]) continue;
            if (game->owner(v) == parity) {
                // set random strategy for winner
                for (auto curedge = game->outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (!disabled[to]) {
                        solve(v, parity, to);
                        break;
                    }
                }
            } else {
                solve(v, parity, -1);
            }
        }
        flush();
        return true;
    } else {
        // all disabled
        return false;
    }
}

void
Oink::solve(int node, int win, int strategy)
{
    /*
    if (trace) {
        logger << "\033[1;32msolved " << (win ? "(odd)" : "(even)") << "\033[m " << game->label_vertex(node);
        if (strategy != -1) logger << " to " << game->label_vertex(strategy);
        logger << std::endl;
    }
    // */

#ifndef NDEBUG
    if (game->isSolved(node) or disabled[node]) LOGIC_ERROR;
#endif

    game->solve(node, win, strategy);
    disabled[node] = true; // disable
    todo.push(node);
}

void
Oink::flush()
{
    // the <todo> queue contains vertex that are marked as solved

    while (todo.nonempty()) {
        int v = todo.pop();
        bool winner = game->getWinner(v);

        for (auto curedge = game->ins(v); *curedge != -1; curedge++) {
            int from = *curedge;
            if (!game->isSolved(from) and !disabled[from]) {
                if (game->owner(from) == winner) {
                    // node of winner
                    solve(from, winner, v);
                } else {
                    // node of loser
                    if (--outcount[from] == 0) solve(from, winner, -1);
                }
            }
        }
    }
}

void
Oink::setSolver(std::string solver)
{
    this->solver = solver;
}

void _solve_loop(Oink* s)
{
    s->solveLoop();
}

VOID_TASK_1(solve_loop, Oink*, s)
{
    _solve_loop(s);
}

void
Oink::solveLoop()
{
    /**
     * Report chosen solver.
     */
    if (bottomSCC) {
        do {
            // disable all solved vertices
            disabled = game->getSolved();  // copy assignment

            // solve bottom SCC
            std::vector<int> sel;
            getBottomSCC(sel);
            assert(sel.size() != 0);
            disabled.set();
            for (int i : sel) disabled[i] = false;

            logger << "solving bottom SCC of " << sel.size() << " nodes (";
            logger << game->count_unsolved() << " nodes left)" << std::endl;

            // solve current subgame
            auto s = Solvers::construct(*solver, *this, *game);
            if (!s->parseOptions(options)) {
                logger << "error parsing options: " << options << std::endl;
                exit(-1);
            }
            s->run();

            // flush the todo buffer
            flush();
        } while (!game->game_solved());
    } else {
        do {
            // disable all solved vertices
            disabled = game->getSolved();

            // solve current subgame
            auto s = Solvers::construct(*solver, *this, *game);
            bool fullSolver = true; // every solver is actually a full solver
            if (!s->parseOptions(options)) {
                logger << "error parsing options: " << options << std::endl;
                exit(-1);
            }
            s->run();

            if (fullSolver) {
                // trash the todo buffer
                todo.clear();
                return;
            } else {
                // flush the todo buffer
                flush();
                auto c = game->count_unsolved();
                logger << c << " nodes left." << std::endl;
                if (c == 0) return;
            }
        } while (true);
    }
}

void
Oink::run()
{
    using namespace std::chrono;
    auto time_before = high_resolution_clock::now();

    /**
     * Now inflate / compress / renumber...
     */
    if (inflate) {
        int d = game->inflate();
        logger << "parity game inflated (" << d << " priorities)" << std::endl;
    } else if (compress) {
        int d = game->compress();
        logger << "parity game compressed (" << d << " priorities)" << std::endl;
    } else if (renumber) {
        int d = game->renumber();
        logger << "parity game renumbered (" << d << " priorities)" << std::endl;
    }

    /**
     * Deal with partial solutions
     * TODO: test this code, or maybe disable partial solutions and only accept full solutions for verification??
     */
    if (game->getSolved().any()) {
        for (int v=0; v<game->vertexcount(); v++) {
            if (game->isSolved(v)) todo.push(v);
        }
        flush();
    }

    if (solveSingle and solveSingleParity()) {
        // already reported in solveSingleParity
        auto time_after = high_resolution_clock::now();
        double diff = duration_cast<duration<double>>(time_after - time_before).count();
        logger << "preprocessing took " << std::fixed << std::setprecision(6) << diff << " sec." << std::endl;
        logger << "solved by preprocessor." << std::endl;
        return;
    }

    if (removeLoops) {
        int count = solveSelfloops();
        if (count == 0) logger << "no self-loops removed." << std::endl;
        else if (count == 1) logger << "1 self-loops removed." << std::endl;
        else logger << count << " self-loops removed." << std::endl;
    }

    if (removeWCWC) {
        int count = solveTrivialCycles();
        if (count == 0) logger << "no trivial cycles removed." << std::endl;
        else if (count == 1) logger << "1 trivial cycle removed." << std::endl;
        else logger << count << " trivial cycles removed." << std::endl;
    }

    auto time_mid = high_resolution_clock::now();

    if (game->game_solved()) {
        double preprocess_time = duration_cast<duration<double>>(time_mid - time_before).count();
        logger << "preprocessing took " << std::fixed << std::setprecision(6) << preprocess_time << " sec." << std::endl;
        logger << "solved by preprocessor." << std::endl;
        return;
    }
        
    if (solver == std::nullopt) {
        logger << "no solver selected!" << std::endl;
        return;
    }

    // building these arrays is not really a preprocessing step in my opinion
    time_mid = high_resolution_clock::now();

    /**
     * Start Lace if we are parallel
     * - if parallel solver, -w [0..N] and Lace is not running, start Lace
     * - if parallel solver, -w -1, run sequential anyway
     * - if sequential solver, run sequantial
     */

    logger << "solving using " << Solvers::desc(*solver) << std::endl;

    if (Solvers::isParallel(*solver)) {
        if (workers >= 0) {
            if (lace_workers() == 0) {
                lace_start(workers, 0);
                logger << "initialized Lace with " << lace_workers() << " workers" << std::endl;
                RUN(solve_loop, this);
                lace_stop();
            } else {
                logger << "running parallel (Lace already initialized)" << std::endl;
                solveLoop();
            }
        } else {
            logger << "running sequentially" << std::endl;
            solveLoop();
        }
    } else {
        solveLoop();
    }

    auto time_after = high_resolution_clock::now();
    double preprocess_time = duration_cast<duration<double>>(time_mid - time_before).count();
    logger << "preprocessing took " << std::fixed << std::setprecision(6) << preprocess_time << " sec." << std::endl;
    double solving_time = duration_cast<duration<double>>(time_after - time_mid).count();
    logger << "solving took " << std::fixed << std::setprecision(6) << solving_time << " sec." << std::endl;
}

}
