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
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <queue>
#include <deque>
#include <stack>

#include "verifier.hpp"

using namespace std;

namespace pg {

void
Verifier::verify(bool fullgame, bool even, bool odd)
{
    // The caller must pass a game that is already sorted by priority; the SCC
    // search below relies on vertices being ordered by priority. Verification
    // only reads outgoing edges, so no incoming-edge array is needed.
    const int n_vertices = game.vertexcount();

    /**
     * The first loop removes all edges from "won" vertices that are not the strategy.
     * This turns each dominion into a single player game.
     * Also some trivial checks are performed.
     */
    for (int v=0; v < n_vertices; v++) {
        // (for full solutions) check whether every vertex is won
        if (!game.isSolved(v)) {
            if (fullgame) throw std::runtime_error("not every vertex is won");
            else continue;
        }

        const bool winner = game.getWinner(v) ;

        if (winner == 0 and !even) continue; // whatever
        if (winner == 1 and !odd) continue; // whatever

        if (winner == game.owner(v)) {
            if (game.hasMultiStrategyFor(v)) {
                // multi-strategy: check that *every* recorded strategy edge
                // stays in the dominion (each edge is a real move by construction)
                const int* outbase = game.outedges();
                for (auto curedge = game.outs(v); *curedge != -1; curedge++) {
                    if (!game.isStrategyEdgeIndex(curedge - outbase)) continue;
                    const int to = *curedge;
                    if (!game.isSolved(to) or game.getWinner(to) != winner) {
                        throw std::runtime_error("strategy leaves dominion");
                    }
                }
                n_strategies++; // number of checked strategies (>=1 by hasMultiStrategyFor)
            } else {
                // if winner, check whether the strategy stays in the dominion
                int str = game.getStrategy(v);
                if (str == -1) {
                    throw std::runtime_error("winning vertex has no strategy");
                } else if (!game.has_edge(v, str)) {
                    throw std::runtime_error("strategy is not a valid move");
                } else if (!game.isSolved(str) or game.getWinner(str) != winner) {
                    throw std::runtime_error("strategy leaves dominion");
                }
                n_strategies++; // number of checked strategies
            }
        } else {
            // if loser, check whether the loser can escape
            for (int to : game.out_edges(v)) {
                if (!game.isSolved(to) or game.getWinner(to) != winner) {
                    logger << "escape edge from " << game.label_vertex(v) << " to " << game.label_vertex(to) << std::endl;
                    throw std::runtime_error("loser can escape");
                }
            }
            // and of course check that no strategy is set
            if (game.getStrategy(v) != -1) throw std::runtime_error("losing vertex has strategy");
        }
    }

    // Allocate datastructures for Tarjan search
    int *done = new int[n_vertices];
    int64_t *low = new int64_t[n_vertices];
    for (int i=0; i<n_vertices; i++) done[i] = -1;
    for (int i=0; i<n_vertices; i++) low[i] = 0;

    std::vector<int> res;
    std::stack<int> st;

    int64_t pre = 0;

    for (int top = n_vertices - 1; top >= 0; top--) {
        // only if a dominion
        if (!game.isSolved(top)) continue;

        int prio = game.priority(top);
        int winner = game.getWinner(top);

        // only compute SCC for a (probably) top vertex
        if (winner == 0 and !even) continue; // don't check even dominions
        if (winner == 1 and !odd) continue; // don't check odd dominions

        // only try to find an SCC where the loser wins
        if (winner == (prio&1)) continue;

        // only run the check if not yet done at priority <prio>
        if (done[top] == prio) continue;

        // set <bot> (in tarjan search) to current pre
        int64_t bot = pre;

        // start the tarjan search at vertex <top>
        st.push(top);

        while (!st.empty()) {
            int v = st.top();

            /**
             * When we see it for the first item, we assign the next number to it and add it to <res>.
             */
            if (low[v] <= bot) {
                low[v] = ++pre;
                res.push_back(v);
            }

            /**
             * Now we check all outgoing (allowed) edges.
             * If seen earlier, then update "min"
             * If new, then 'recurse'
             */
            int min = low[v];
            bool pushed = false;
            if (game.hasMultiStrategyFor(v)) {
                // multi-strategy: the winner may use *any* recorded strategy edge,
                // so follow all of them (like the loser's all-edges branch below,
                // but restricted to the strategy edges)
                const int* outbase = game.outedges();
                for (auto curedge = game.outs(v); *curedge != -1; curedge++) {
                    if (!game.isStrategyEdgeIndex(curedge - outbase)) continue;
                    int to = *curedge;
                    // skip if to higher priority
                    if (game.priority(to) > prio) continue;
                    // skip if already found scc (done[to] set to prio)
                    if (done[to] == prio) continue;
                    // check if visited in this search
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
            } else if (game.getStrategy(v) != -1) {
                int to = game.getStrategy(v);
                if (game.priority(to) > prio) {
                    // skip if to higher priority
                } else if (done[to] == prio) {
                    // skip if already found scc (done[to] set to prio)
                } else if (low[to] <= bot) {
                    // not visited, add to <st> and break!
                    st.push(to);
                    pushed = true;
                } else {
                    // visited, update min
                    if (low[to] < min) min = low[to];
                }
            } else {
                for (auto curedge = game.outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    // skip if to higher priority
                    if (game.priority(to) > prio) continue;
                    // skip if already found scc (done[to] set to prio)
                    if (done[to] == prio) continue;
                    // check if visited in this search
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
            }
            if (pushed) continue; // we pushed a new vertex to <st>...

            // there was no edge to a new vertex
            // check if we are the root of a SCC
            if (min < low[v]) {
                // not the root (outgoing edge to lower ranked vertex)
                low[v] = min;
                st.pop();
                continue;
            }

            /**
             * We're the root!
             * Now we need to figure out if we have cycles with p...
             * Also, mark every vertex in the SCC as "done @ search p"
             */
            int max_prio = -1;
            int scc_size = 0;

            for (auto it=res.rbegin(); it!=res.rend(); it++) {
                int n = *it;
                if (game.priority(n) > max_prio) max_prio = game.priority(n);
                scc_size++;
                done[n] = prio; // mark as done at prio
                if (n == v) break;
            }

            bool self_strategy;
            if (game.hasMultiStrategyFor(v)) {
                // winner uses a self-loop only if v->v is a recorded strategy edge
                self_strategy = game.hasStrategyEdgeTo(v, v);
            } else if (game.getStrategy(v) != -1) {
                self_strategy = (game.getStrategy(v) == v);
            } else {
                // loser: any self-edge keeps v in its own SCC
                self_strategy = game.has_edge(v, v);
            }
            bool cycles = scc_size > 1 or self_strategy;

            if (cycles && (max_prio&1) == (prio&1)) {
                /**
                 * Found! Report.
                 */
                logger << "\033[1;31mscc where loser wins\033[m with priority \033[1;34m" << max_prio << "\033[m";
                for (auto it=res.rbegin(); it!=res.rend(); it++) {
                    int n = *it;
                    logger << " " << game.label_vertex(n);
                    if (n == v) break;
                }
                logger << std::endl;
                delete[] done;
                delete[] low;

                throw std::runtime_error("loser can win");
            }

            /**
             * Not found! Continue.
             */
            for (;;) {
                int n = res.back();
                res.pop_back();
                if (n == v) break;
            }
            st.pop();
        }
    }

    delete[] done;
    delete[] low;
}

}
