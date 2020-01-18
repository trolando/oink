/*
 * Copyright 2017-2020 Tom van Dijk, Johannes Kepler University Linz
 * Copyright 2019-2020 Ruben Lapauw, KU Leuven
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
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "zlkj.hpp"
#include "lace.h"
#include "printf.hpp"

namespace pg {

static const int DIS = 0x80000000; // permanently disabled vertex
static const int BOT = 0x80000001; // bottom state for vertex


ZLKJSolver::ZLKJSolver(Oink *oink, Game *game) : Solver(oink, game), Q(game->n_nodes), D(game->n_nodes)
{
}

ZLKJSolver::~ZLKJSolver()
{
    delete[] inverse;
}

int
ZLKJSolver::attractExt(int i, int r, std::vector<int> *R)
{
    const int pr = priority[i];
    const int pl = pr & 1;

    /**
     * Starting at <i>, attract head nodes until "inversion"
     */
    for (; i>=0; i--) {
        if (region[i] == DIS or region[i] >= 0) continue;         // cannot be attracted

        // uncomment the next line to attract until lower priority instead of until inversion
        //if (priority[i] != pr) break; // until other priority
        if ((priority[i]&1) != pl) break;         // until parity inversion (Maks Verver optimization)

        region[i] = r;
        winning[i] = pl;
        strategy[i] = -2;         // head nodes do not have a strategy yet!
        Q.push(i);
        R->push_back(i);
        while (!Q.empty()) {
            int cur = Q.pop();

            // attract to <cur>
            auto end = just[cur].end();
            auto it = just[cur].begin();
            for (; it != end; it++) {
                int from = *it;
                if (from > i or region[from] == DIS or region[from] >= 0) continue;                 // cannot be attracted

                if (owner[from] == pl) {
                    // owned by same parity
                    region[from] = r;
                    winning[from] = pl;
                    strategy[from] = cur;
                    Q.push(from);

                    // Invariant: strategy[from] = cur <=> {x: from in just[x]} = {cur};
                    const int *_out = outs + outa[from];
                    for (int to = *_out; to != -1; to = *++_out) {
                        if (to == cur) continue;
                        if (region[to] == DIS) continue;
                        just[to].erase(from);
                    }
                } else {
                    // owned by other parity
                    int count = region[from];
                    if (count == BOT) {
                        // compute count (to negative)
                        count = 1;
                        const int *_out = outs + outa[from];
                        for (int to = *_out; to != -1; to = *++_out) {
                            if (region[to] == DIS) continue;
                            if (region[to] >= 0 and region[to] < r) continue;
                            count--;
                        }
                    } else {
                        count++;
                    }
                    if (count == 0) {
                        region[from] = r;
                        winning[from] = pl;
                        strategy[from] = -1;
                        Q.push(from);
                        const int *_out = outs + outa[from];
                        for (int to = *_out; to != -1; to = *++_out) {
                            if (to == cur) continue;
                            if (region[to] == DIS) continue;
                            just[to].insert(from);
                        }
                    } else {
                        region[from] = count;
                    }
                }
            }
        }
    }

    return i;
}

/**
 * Find all nodes in S (and in the game with region >= r) that attract to the other player.
 */
std::pair<int, int> ZLKJSolver::attractLosing(int i, int r, std::vector<int> *S, int next_r)
{
    int count = 0;
    const int pr = priority[i];
    const int pl = pr & 1;
    int new_i = -1;

    // NOTE: this algorithm could be improved using an "out counter"

#ifndef NDEBUG
    for (int i : *S) if (winning[i] != pl) LOGIC_ERROR;
#endif

    /**
     * First check all region nodes...
     * In reality, we just want to check the "head nodes", because all other nodes are attracted
     * to the head nodes and cannot be attracted to the opponent directly.
     * But we do not record which nodes are head nodes.
     * TODO! DONE!
     */
    S->erase(std::remove_if(S->begin(), S->end(), [this](int i) {
            return this->strategy[i] != -2;
        }), S->end());
    for (int i : *S) {
        // check if the node is attracted
        if (owner[i] == pl) {
            // "loser" attraction
            bool can_escape = false;
            const int *_out = outs + outa[i];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] < r) continue;                 // not in subgame, or -1/-2, would otherwise be attracted to subgame previously
                // if not previously attracted, then either DIS or winning[to] != pl;
                if (winning[to] != pl) continue;                 // not an escape
                can_escape = true;
                //just[from].insert(to);
                break;
            }
            if (!can_escape) {
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", i, priority[i], 1-pl);
                region[i] = r;
                winning[i] = 1-pl;
                strategy[i] = -1;
                Q.push(i);
                // strategy[from] = cur <=> {x: from in just[x]} = {cur};
                const int *_out = outs + outa[i];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (to == i) continue;                     // no concurrent modification error
                    if (region[to] == DIS) continue;
                    just[to].insert(i);
                }
            }
        } else {
            // "winner" attraction
            const int *_out = outs + outa[i];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] < r) continue;                 // not in subgame, or -1/-2
                if (winning[to] == pl) continue;                 // not attracting
                region[i] = r;
                winning[i] = 1-pl;
                strategy[i] = to;
                Q.push(i);
                const int *second = outs + outa[i];
                for (int other = *second; other != -1; other = *++second) {
                    if (other == to) continue;                     // no concurrent modification error
                    if (region[other] == DIS) continue;
                    just[other].erase(i);
                }
                break;
            }
        }
    }

    /**
     * Now attract anything in this region/subregions of <pl> to 1-<pl>
     */

    while (!Q.empty()) {
        int cur = Q.pop();
        ++count;
        // attract to <cur>
        auto end = just[cur].end();
        auto it = just[cur].begin();
        for (; it != end; it++) {
            int from = *it;
            // if (region[from] == -1) LOGIC_ERROR;
            if (region[from] < r && region[from] >= 0) {
                continue;                 // not in subgame, or disabled
            }
            if (winning[from] == 1-pl) {
                continue;                 // already lost
            }

            if (owner[from] != pl) {
                // owned by other
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = cur;
                Q.push(from);
                const int *_out = outs + outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (to == cur) continue;                     // no concurrent modification error
                    if (region[to] == DIS) continue;
                    just[to].erase(from);
                }
            } else {
                // owned by us
                bool can_escape = false;
                {
                    const int *_out = outs + outa[from];
                    for (int to = *_out; to != -1; to = *++_out) {
                        // if (region[to] == -1) LOGIC_ERROR;
                        if (region[to] >= 0 and region[to] < r) continue;                         // not in subgame, or disabled
                        if (winning[to] == 1 - pl) continue;                         // not an escape
                        can_escape = true;
                        //switch strategy
                        //just[strategy[from]].erase(from);
                        just[to].insert(from);
                        strategy[from] = to;
                        break;
                    }
                }
                if (can_escape) {
                    if(winning[from] != -1) {
                        winning[from] = -1;
                        D.push(from);
                    }
                    continue;
                }
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = -1;
                Q.push(from);
                // strategy[from] = cur <=> {x: from in just[x]} = {cur};
                const int *_out = outs + outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (to == cur) continue;                     // no concurrent modification error
                    if (region[to] == DIS) continue;
                    just[to].insert(from);
                }
            }
        }
    }
    for(unsigned int idx = 0; idx < D.size(); idx++) {
        int cur = D[idx];
        // destroy dependents on <cur>
        if(winning[cur] != -1) {
            continue;
        }
        region[cur] = BOT;
        auto end = just[cur].end();
        auto it = just[cur].begin();
        for (; it != end; it++) {
            int to = *it;
            if (region[to] < r) {
                continue;                 // not in subgame, or disabled
            }
            if (winning[to] != pl) {
                continue;
            }
            region[to] = BOT;
            winning[to] = -1;
            D.push(to);
        }
    }

    for(unsigned int idx = 0; idx < D.size(); idx++) {
        int cur = D[idx];
        if(region[cur] != BOT) {
            continue;
        }
        strategy[cur] = -3;
        if (owner[cur] == pl) {
            // "winner" attraction
            const int *_out = outs + outa[cur];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] < r) continue;                 // not in subgame, or -1/-2
                if (winning[to] == 1 - pl) continue;                 // not attracting
                region[cur] = next_r;
                winning[cur] = pl;
                strategy[cur] = to;
                Q.push(cur);
                const int *second = outs + outa[cur];
                for (int other = *second; other != -1; other = *++second) {
                    if (other == to) continue;                     // no concurrent modification error
                    if (region[other] == DIS) continue;
                    just[other].erase(cur);
                }
                just[to].insert(cur);
                break;
            }
            if(strategy[cur] == -3) {
                const int *second = outs + outa[cur];
                for (int other = *second; other != -1; other = *++second) {
                    if (region[other] == DIS) continue;
                    just[other].insert(cur);
                }
            }
        } else {
            // "loser" attraction
            bool can_escape = false;
            const int *_out = outs + outa[cur];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] >= 0 and region[to] < r) continue;                 // not in subgame, or -1/-2, would otherwise be attracted to subgame previously
                // if not previously attracted, then either DIS or winning[to] != pl;
                if (winning[to] == 1 - pl) continue;                 // not an escape
                can_escape = true;
                break;
            }
            if (!can_escape) {
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", i, priority[i], 1-pl);
                region[cur] = next_r;
                winning[cur] = pl;
                strategy[cur] = -1;
                Q.push(cur);
                // strategy[from] = cur <=> {x: from in just[x]} = {cur};
                const int *_out = outs + outa[cur];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (to == i) continue;                     // no concurrent modification error
                    if (region[to] == DIS) continue;
                    just[to].insert(cur);
                }
            }
        }
        while (!Q.empty()) {
            int cur = Q.pop();

            // attract to <cur>
            auto end = just[cur].end();
            auto it = just[cur].begin();
            for (; it != end; it++) {
                int from = *it;
                if (from > i or region[from] == DIS or region[from] >= 0) continue;                 // cannot be attracted

                if (owner[from] == pl) {
                    // owned by same parity
                    region[from] = next_r;
                    winning[from] = pl;
                    strategy[from] = cur;
                    Q.push(from);

                    // strategy[from] = cur <=> {x: from in just[x]} = {cur};
                    const int *_out = outs + outa[from];
                    for (int to = *_out; to != -1; to = *++_out) {
                        if (to == cur) continue;
                        if (region[to] == DIS) continue;
                        just[to].erase(from);
                    }
                } else {
                    // owned by other parity
                    int count = region[from];
                    if (count == BOT) {
                        // compute count (to negative)
                        count = 1;
                        const int *_out = outs + outa[from];
                        for (int to = *_out; to != -1; to = *++_out) {
                            if (region[to] == DIS) continue;
                            if (region[to] >= 0 and region[to] < r) continue;
                            count--;
                        }
                    } else {
                        count++;
                    }
                    if (count == 0) {
                        region[from] = next_r;
                        winning[from] = pl;
                        strategy[from] = -1;
                        Q.push(from);
                        const int *_out = outs + outa[from];
                        for (int to = *_out; to != -1; to = *++_out) {
                            if (to == cur) continue;
                            if (region[to] == DIS) continue;
                            just[to].insert(from);
                        }
                    } else {
                        //TODO : Make count useful, now basic can_escape.
                        // Note reassign the count to region[from] lead to bugs.
                    }
                }
            }
        }
    }
    int lost = 0;
    while(!D.empty()) {
        int cur = D.pop();

        if(region[cur] != BOT) {
            continue;
        }
        strategy[cur] = -3;
        lost++;
        if(new_i < cur) {
            new_i = cur;
        }
    }
    return std::pair<int, int>(count, new_i);
}

void
ZLKJSolver::run()
{
    iterations = 0;

    // allocate and initialize data structures
    region = new int[n_nodes];
    winning = new int[n_nodes];
    strategy = new int[n_nodes];

    std::vector<int> history;
    std::vector<std::vector<int> > levels;

    // get number of nodes and create and initialize inverse array
    max_prio = -1;
    for (int n=n_nodes-1; n>=0; n--) {
        strategy[n] = -4;
        winning[n] = -1;
        region[n] = disabled[n] ? DIS : BOT;
        if (disabled[n]) continue;
        strategy[n] = -3;
        const int pr = game->priority[n];
        if (max_prio == -1) {
            max_prio = pr;
            inverse = new int[max_prio+1];
            memset(inverse, -1, sizeof(int[max_prio+1]));
        }
        if (inverse[pr] == -1) inverse[pr] = n;
    }
    just.reserve(n_nodes);
    for (int n=0; n<n_nodes; n++) {
        std::unordered_set<int> incoming;
        just.push_back(incoming);
        if(disabled[n]) continue;
        just[n].reserve(game->in[n].size());
        for(int from : game->in[n]) {
            if (disabled[from]) continue;
            just[n].insert(from);
        }
    }
    if (max_prio == -1) LOGIC_ERROR;     // unexpected empty game

    // pre-allocate some space (similar to NPP)
    {
        int space = max_prio / 20;
        if (space >= 500) space = 500;
        history.reserve(3*space);
        levels.reserve(space);
    }

    // start the loop at the last node (highest priority) and at depth 0
    int i = inverse[max_prio];
    int next_r = 0;

    // initialize first level (i, r=0, phase=0)
    levels.push_back(std::vector<int>());
    history.push_back(i);
    history.push_back(next_r++);
    history.push_back(0);
    while (true) {
        // obtain current frame
        const int hsize = history.size();
        if (hsize == 0) break;         // no frame on the stack

        std::vector<int> *A = &(*levels.rbegin());
        const int i = history[hsize-3];
        const int r = history[hsize-2];
        const int phase = history[hsize-1];

#ifndef NDEBUG
        if (i < 0) LOGIC_ERROR;         // just a sanity check
        if (region[i] == DIS) LOGIC_ERROR;         // just a sanity check
        if (h*3 != hsize or (int)levels.size() != h) LOGIC_ERROR;         // just a sanity check
#endif


        /**
         * Get priority and player
         */
        const int pr = priority[i];
        const int pl = pr&1;

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "\n\033[1mDepth %d phase %d\033[m: node %d priority %d\n", h-0, phase, i, pr);
#endif

        if (phase == 0) {
            /**
             * We are in the first phase.
             * Compute extended attractor (until inversion). Then recursive step.
             */

            // attract until inversion and add to A
            int j = attractExt(i, r, A);
            // j is now the next i (subgame), or -1 if the subgame is empty
            if(trace) std::cout << "attracted " << A->size() << " nodes" << " for " << pl << std::endl;

            // count number of iterations
            ++iterations;

            // update phase to 1
            history.back() = 1;

            if (j != -1) {
                // go recursive (subgame not empty)
                levels.push_back(std::vector<int>());
                history.push_back(j);
                history.push_back(next_r++);
                history.push_back(0);
            }
        } else if (phase == 1) {
            /**
             * After first recursion step
             */

#ifndef NDEBUG
            if (trace) fmt::printf(logger, "current level contains %zu nodes\n", A->size());
#endif

            /**
             * Compute part of region that attracts to the other player...
             */

            int count = -1;
            int new_i = -1;
            std::tie(count, new_i) = attractLosing(i, r, A, next_r++);
            if(trace) std::cout << "Losing " << count << " " << new_i << std::endl;
            if (trace) {
                if (count > 0) fmt::printf(logger, "%d nodes are attracted to losing region\n", count);
                else if (count == 0) fmt::printf(logger, "no nodes are attracted to losing region\n");
                else fmt::printf(logger, "no losing region\n");
            }

            if (count <= 0) {
                if(new_i >= 0) {                 // count =< 0 => new_i >= 0
                    LOGIC_ERROR;
                }
                /**
                 * Nothing attracted to opponent, add A to W0/W1, fix strategies, go up.
                 */
                for (int v : *A) {
                    /**
                     * For nodes that are won and controlled by <pl>, check if their strategy needs to be fixed.
                     */
                    if (strategy[v] != -2) {
                        if(strategy[v] != -1 && winning[strategy[v]] != pl) {
                            LOGIC_ERROR;
                        }
                        continue;                         // good strategy
                    }
                    if (owner[v] != pl) {
                        strategy[v] = -1;
                        continue;                         // not controlled by <pl>
                    }
                    /**
                     * Strategy of vertex <v> needs to be updated!
                     * We search for a successor of <v> in the subgame won by <pl>
                     */
                    //strategy[v] = -1;
                    const int *_out = outs + outa[v];
                    for (int to = *_out; to != -1; to = *++_out) {
                        just[to].erase(to);
                        if (region[to] < r) continue;                         // not in subgame
                        if (winning[to] != pl) continue;                         // not winning
                        if (strategy[v] != -2) continue;
                        strategy[v] = to;
                    }
                    just[strategy[v]].insert(v);
                    if (strategy[v] < 0) LOGIC_ERROR;
                }

                /**
                 * The strategy has been updated. Nodes are now in W0/W1.
                 * Finally, pop the stack and go up...
                 */
                levels.pop_back();
                history.pop_back();
                history.pop_back();
                history.pop_back();
            } else {
                if (new_i == -1) {
                    /**
                     * And pop the stack to go up
                     */
                    history.back() = 2;                     // set current phase to 2
                    /*levels.pop_back();
                       history.pop_back();
                       history.pop_back();
                       history.pop_back();*/
                } else {
                    /**
                     * And push the stack to go down
                     */
                    history.back() = 2;                     // set current phase to 2
                    levels.push_back(std::vector<int>());
                    history.push_back(new_i);
                    history.push_back(next_r++);
                    history.push_back(0);
                }
            }
        } else if (phase == 2) {
            /**
             * After second recursion step
             * Choose a strategy
             */
            for (int v : *A) {
                /**
                 * For nodes that are won and controlled by <pl>, check if their strategy needs to be fixed.
                 */
                if (strategy[v] != -2) {
                    if(strategy[v] != -1 && winning[strategy[v]] != winning[v]) {
                        LOGIC_ERROR;
                    }
                    continue;                     // good strategy
                }
                if (owner[v] != pl) {
                    strategy[v] = -1;
                    continue;                     // not controlled by <pl>
                }
                /**
                 * Strategy of vertex <v> needs to be updated!
                 * We search for a successor of <v> in the subgame won by <pl>
                 */
                //strategy[v] = -1;
                const int *_out = outs + outa[v];
                for (int to = *_out; to != -1; to = *++_out) {
                    just[to].erase(v);
                    if (region[to] < r) continue;                     // not in subgame
                    if (winning[to] != pl) continue;                     // not winning
                    if (strategy[v] != -2) continue;
                    strategy[v] = to;
                }
                just[strategy[v]].insert(v);
                if (strategy[v] < 0) LOGIC_ERROR;
            }

            /**
             * The strategy has been updated. Nodes are now in W0/W1.
             * Finally, pop the stack to go up...
             */
            levels.pop_back();
            history.pop_back();
            history.pop_back();
            history.pop_back();
        }
    }

    // done
    for (int i=0; i<n_nodes; i++) {
        if (region[i] == DIS) continue;
#ifndef NDEBUG
        if (winning[i] == -1) LOGIC_ERROR;
#endif
        oink->solve(i, winning[i], strategy[i]);
    }

    delete[] region;
    delete[] winning;
    delete[] strategy;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
