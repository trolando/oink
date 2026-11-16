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
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "verifier.hpp"

using namespace std;

namespace pg {

inline char
Verifier::successorMode(int v) const
{
    if (solution.strategy(v) != -1) return 0;
    // without a multi-strategy, the first part ensures that v is won by its opponent
    if (solution.has_multi() and getWinner(v) == game.owner(v)) return 1;
    return 2;
}

inline long
Verifier::firstPosition(int v, char mode) const
{
    return mode == 0 ? 0 : game.firstout(v);
}

inline int
Verifier::nextSuccessor(int v, char mode, long& pos) const
{
    if (mode == 0) {
        // a single strategy: the only successor
        if (pos != 0) return -1;
        pos = 1;
        return solution.strategy(v);
    }

    // the moves of the multi-strategy (mode 1) or all moves (mode 2)
    const int* edges = game.outedges();
    while (edges[pos] != -1) {
        const long e = pos++;
        if (mode == 2 or solution.has_edge_index(e)) return edges[e];
    }
    return -1;
}

void
Verifier::verify(bool fullgame, bool even, bool odd)
{
    // The result does not depend on the order of the vertices, but the searches below
    // are faster if the vertices are sorted by priority. Verification only reads
    // outgoing edges, so no incoming-edge array is needed.
    const int n_vertices = game.vertexcount();

    /**
     * The first loop removes all edges from "won" vertices that are not the strategy.
     * This turns each dominion into a single player game.
     * Also some trivial checks are performed.
     */
    for (int v=0; v < n_vertices; v++) {
        // (for full solutions) check whether every vertex is won
        if (!solution.is_solved(v)) {
            if (fullgame) throw std::runtime_error("not every vertex is won");
            else continue;
        }

        const bool winner = getWinner(v) ;

        if (winner == 0 and !even) continue; // whatever
        if (winner == 1 and !odd) continue; // whatever

        if (winner == game.owner(v)) {
            int str = solution.strategy(v);
            if (str != -1) {
                // single strategy: check it stays in the dominion
                if (!game.has_edge(v, str)) {
                    throw std::runtime_error("strategy is not a valid move");
                } else if (!solution.is_solved(str) or getWinner(str) != winner) {
                    throw std::runtime_error("strategy leaves dominion");
                }
                n_strategies++; // number of checked strategies
            } else {
                // sentinel: -1 with winner==owner means the moves are in the
                // multi-strategy (if one exists); a winner-owned vertex with
                // neither a single strategy nor a multi-strategy is an error.
                if (!solution.has_multi()) {
                    throw std::runtime_error("winning vertex has no strategy");
                }
                // check that *every* recorded edge stays in the dominion
                const int* outbase = game.outedges();
                int count = 0;
                for (auto curedge = game.outs(v); *curedge != -1; curedge++) {
                    if (!solution.has_edge_index(curedge - outbase)) continue;
                    const int to = *curedge;
                    if (!solution.is_solved(to) or getWinner(to) != winner) {
                        throw std::runtime_error("strategy leaves dominion");
                    }
                    count++;
                }
                if (count == 0) throw std::runtime_error("winning vertex has no strategy");
                n_strategies++; // number of checked strategies
            }
        } else {
            // if loser, check whether the loser can escape
            for (int to : game.out_edges(v)) {
                if (!solution.is_solved(to) or getWinner(to) != winner) {
                    logger << "escape edge from " << game.label_vertex(v) << " to " << game.label_vertex(to) << std::endl;
                    throw std::runtime_error("loser can escape");
                }
            }
            // and of course check that no strategy is set
            if (solution.strategy(v) != -1) throw std::runtime_error("losing vertex has strategy");
        }
    }

    /**
     * The second part searches, for each priority <prio> of the loser's parity, the strongly
     * connected components (SCCs) of the solution graph restricted to the vertices with a
     * priority of at most <prio>, from the vertices with priority <prio>. The solution graph
     * has the moves of the strategy from each vertex won by its owner and all moves from the
     * other vertices. An SCC with a cycle and a highest priority of the loser's parity has a
     * cycle that the loser wins.
     *
     * Each vertex on the stack of the search keeps its position in its successors, so that
     * the depth-first search visits each edge once (as in Tarjan's algorithm). If the
     * vertices are sorted by priority, then the searches at one priority follow each other,
     * and together they visit each vertex at most once.
     */

    // <done> has the priority of the last search that found the SCC of a vertex. A vertex
    // gets a number in <low> when the search visits it; the numbers increase over all
    // searches, so a vertex is visited in the current search if its number is above <bot>.
    // Then <low> decreases when the vertex reaches a lower number on the stack.
    std::vector<int> done(n_vertices, -1);
    std::vector<int64_t> low(n_vertices, 0);

    std::vector<int> res; // the vertices of the SCCs that are not yet complete

    // the search stack: a vertex, its position in its successors, whether it is the root of
    // its SCC so far (its number in <low> did not decrease), and how to find its successors
    struct Frame {
        long pos;
        int v;
        bool root;
        char mode;
    };
    std::vector<Frame> st;

    int64_t pre = 0;

    for (int top = n_vertices - 1; top >= 0; top--) {
        // only if a dominion
        if (!solution.is_solved(top)) continue;

        int prio = game.priority(top);
        int winner = getWinner(top);

        // only compute SCC for a (probably) top vertex
        if (winner == 0 and !even) continue; // don't check even dominions
        if (winner == 1 and !odd) continue; // don't check odd dominions

        // only try to find an SCC where the loser wins
        if (winner == (prio&1)) continue;

        // only run the check if not yet done at priority <prio>
        if (done[top] == prio) continue;

        // set <bot> (in tarjan search) to current pre
        const int64_t bot = pre;

        auto visit = [&](int v) {
            low[v] = ++pre;
            res.push_back(v);
            const char mode = successorMode(v);
            st.push_back({firstPosition(v, mode), v, true, mode});
        };

        // start the tarjan search at vertex <top>
        visit(top);

        while (!st.empty()) {
            Frame& f = st.back();
            const int v = f.v;

            /**
             * Check the next outgoing (allowed) edge.
             * If new, then 'recurse'; if on the stack, then update low.
             */
            const int to = nextSuccessor(v, f.mode, f.pos);
            if (to != -1) {
                // skip if to higher priority
                if (game.priority(to) > prio) continue;
                // skip if already found scc (done[to] set to prio)
                if (done[to] == prio) continue;
                if (low[to] <= bot) {
                    // not visited, add to <st>
                    visit(to);
                } else if (low[to] < low[v]) {
                    // visited and on the stack, update low
                    low[v] = low[to];
                    f.root = false;
                }
                continue;
            }

            // there are no more edges; pass low to the parent
            const bool root = st.back().root;
            st.pop_back();
            if (!st.empty()) {
                Frame& parent = st.back();
                if (low[v] < low[parent.v]) {
                    low[parent.v] = low[v];
                    parent.root = false;
                }
            }

            // check if we are the root of a SCC
            if (!root) continue;

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

            // a single vertex has a cycle if it is its own successor
            bool cycles = scc_size > 1;
            if (!cycles) {
                const char mode = successorMode(v);
                long pos = firstPosition(v, mode);
                for (int next; (next = nextSuccessor(v, mode, pos)) != -1; ) {
                    if (next == v) {
                        cycles = true;
                        break;
                    }
                }
            }

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
        }
    }
}

}
