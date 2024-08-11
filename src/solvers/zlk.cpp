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
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "zlk.hpp"
#include "lace.h"
#include "printf.hpp"

namespace pg {

static const int DIS = 0x80000000; // permanently disabled vertex
static const int BOT = 0x80000001; // bottom state for vertex

#define KC "\033[36;1m"

ZLKSolver::ZLKSolver(Oink& oink, Game& game) : Solver(oink, game), Q(game.nodecount())
{
}

ZLKSolver::~ZLKSolver()
{
    delete[] inverse;
}

typedef struct
{
    int count;
    int items[];
} par_helper;

par_helper** pvec;

VOID_TASK_4(attractParT, int, pl, int, cur, int, r, ZLKSolver*, s)
{
    s->attractParT(__lace_worker, __lace_dq_head, pl, cur, r);
}

void
ZLKSolver::attractParT(WorkerP* __lace_worker, Task* __lace_dq_head, int pl, int cur, int r)
{
    int c = 0;
    par_helper* ours = pvec[LACE_WORKER_ID];

    // attract to <cur>
    for (auto curedge = ins(cur); *curedge != -1; curedge++) {
        int from = *curedge;
        int _r = region[from];
        if (_r == DIS or _r >= 0) continue; // not in subgame, or attracted

        if (owner(from) == pl) {
            // owned by same parity, use CAS to claim it
            while (true) {
                if (__sync_bool_compare_and_swap(&region[from], _r, r)) {
                    winning[from] = pl;
                    strategy[from] = cur;
                    ours->items[ours->count++] = from;
                    SPAWN(attractParT, pl, from, r, this);
                    c++;
                    break;
                }
                _r = *(volatile int*)&region[from];
                if (_r >= 0) break;
            }
        } else {
            // owned by other parity
            volatile int* ptr = &region[from];
            bool attracted = false;

            _r = __sync_add_and_fetch(ptr, 1); // update _r
            if (_r == (BOT+1)) {
                // we are the first, do add_and_fetch with the count
                int count = 0;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == DIS) continue; // do not count disabled
                    if (region[to] >= 0 and region[to] < r) continue; // do not count supgame
                    count--; // count to negative
                }
                // now set count (in a CAS loop)
                int new_r = count; // +1 -1 (count negative to -1)
                while (true) {
                    if (new_r == -1) {
                        // we're the last, so set to r
                        if (__sync_bool_compare_and_swap(ptr, _r, r)) attracted = true;
                        break;
                    }
                    if (__sync_bool_compare_and_swap(ptr, _r, new_r)) break;
                    _r = *ptr;
                    if (_r >= 0) break; // someone else moved to r!
                    // someone else did add and fetch, recompute and try again
                    new_r = count - (BOT - _r) - 1;
                }
            } else if (_r == -1) {
                // another CAS because we may be competing with attractPar
                if (__sync_bool_compare_and_swap(ptr, -1, r)) attracted = true;
            }
            if (attracted) {
                winning[from] = pl;
                strategy[from] = -1;
                ours->items[ours->count++] = from;
                SPAWN(attractParT, pl, from, r, this);
                c++;
            }
        }
    }

    while (c) { SYNC(attractParT); c--; }
}

TASK_4(int, attractPar, int, i, int, r, std::vector<int>*, R, ZLKSolver*, s)
{
    return s->attractPar(__lace_worker, __lace_dq_head, i, r, R);
}

int
ZLKSolver::attractPar(WorkerP* __lace_worker, Task* __lace_dq_head, int i, int r, std::vector<int>* R)
{
    const int pr = priority(i);
    const int pl = pr & 1;

    // initialize pvec (set count to 0) for all workers
    const int W = lace_workers();
    for (int j=0; j<W; j++) pvec[j]->count = 0;

    par_helper* ours = pvec[LACE_WORKER_ID];
    int spawn_count = 0;

    for (; i>=0; i--) {
        int _r = region[i];
        if (_r == DIS or _r >= 0) continue; // not in subgame or attracted
        if (!to_inversion and priority(i) != pr) break;
        if ((priority(i)&1) != pl) { // search until parity inversion
            // first SYNC on all children, who knows this node may be attracted
            while (spawn_count) { SYNC(attractParT); spawn_count--; }
            // after SYNC, check if node <i> is now attracted.
            if (_r < 0) break; // not attracted, so we're done!
            else continue; // already done
        }

        // if c != 0, then we compete with attractParT and must use compare and swap
        if (spawn_count == 0) {
            region[i] = r; // just set, no competing threads
        } else {
            // competing threads! use compare and swap [in a loop]
            while (true) {
                if (__sync_bool_compare_and_swap(&region[i], _r, r)) {
                    _r = r;
                    break;
                }
                _r = *(volatile int*)&region[i];
                if (_r < 0) continue;
                _r = BOT;
                break;
            }
            if (_r == BOT) continue; // someone else claimed!
        }

        winning[i] = pl;
        strategy[i] = -1; // head nodes have no strategy (for now)
        ours->items[ours->count++] = i;
        SPAWN(attractParT, pl, i, r, this);
        spawn_count++;
    }

    // first SYNC on all children (if any)
    while (spawn_count) { SYNC(attractParT); spawn_count--; }

    // update R
    size_t to_reserve = R->size();
    for (int j=0; j<W; j++) to_reserve += pvec[j]->count;
    R->reserve(to_reserve);

    for (int j=0; j<W; j++) {
        par_helper* x = pvec[j];
        for (int k=0; k<x->count; k++) {
#ifndef NDEBUG
            if (trace >= 2) logger << "attracted " << x->items[k] << " (" << priority(x->items[k]) << ")" << std::endl;
#endif
            R->push_back(x->items[k]);
        }
        x->count = 0;
    }

    return i;
}

int
ZLKSolver::attractExt(int i, int r, std::vector<int> *R)
{
    const int pr = priority(i);
    const int pl = pr & 1;

    /**
     * Starting at <i>, attract head nodes until "inversion"
     */

    for (; i>=0; i--) {
        if (region[i] == DIS or region[i] >= 0) continue; // cannot be attracted

        // uncomment the next line to attract until lower priority instead of until inversion
        if (!to_inversion and priority(i) != pr) break; // until other priority
        if ((priority(i)&1) != pl) break; // until parity inversion (Maks Verver optimization)

        region[i] = r;
        winning[i] = pl;
        strategy[i] = -1; // head nodes do not have a strategy yet!
        Q.push(i);

#ifndef NDEBUG
        if (trace >= 2) logger << KC"head\033[m " << label_vertex(i) << std::endl;
#endif

        while (!Q.empty()) {
            int cur = Q.pop();
            R->push_back(cur);

            // attract to <cur>
            for (auto curedge = ins(cur); *curedge != -1; curedge++) {
                int from = *curedge;
                if (from >= i or region[from] == DIS or region[from] >= 0) continue; // cannot be attracted

                if (owner(from) == pl) {
                    // owned by same parity
                    region[from] = r;
                    winning[from] = pl;
                    strategy[from] = cur;
                    Q.push(from);
#ifndef NDEBUG
                    if (trace >= 2) logger << KC"attracted\033[m " << label_vertex(from) << std::endl;
#endif
                } else {
                    // owned by other parity
                    int count = region[from];
                    if (count == BOT) {
                        // compute count (to negative)
                        count = 1;
                        auto curedge = outs(from);
                        for (int to = *curedge; to != -1; to = *++curedge) {
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
#ifndef NDEBUG
                        if (trace >= 2) logger << KC"forced\033[m " << label_vertex(from) << std::endl;
#endif
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
int
ZLKSolver::attractLosing(const int i, const int r, std::vector<int> *S, std::vector<int> *R)
{
    int count = 0;

    const int pr = priority(i);
    const int pl = pr & 1;

    // NOTE: this algorithm could be improved using an "out counter"

#ifndef NDEBUG
    for (int i : *S) if (winning[i] != pl) LOGIC_ERROR;
#endif

    /**
     * First check all region nodes...
     * In reality, we just want to check the "head nodes", because all other nodes are attracted
     * to the head nodes and cannot be attracted to the opponent directly.
     * But we do not record which nodes are head nodes.
     * TODO!
     */
    for (int i : *S) {
        // check if the node is attracted
        if (owner(i) == pl) {
            // "loser" attraction
            bool can_escape = false;
            auto curedge = outs(i);
            for (int to = *curedge; to != -1; to = *++curedge) {
                if (region[to] < r) continue; // not in subgame, or -1/-2
                if (winning[to] != pl) continue; // not an escape
                can_escape = true;
                break;
            }
            if (!can_escape) {
#ifndef NDEBUG
                if (trace >= 2) logger << KC"forced distraction\033[m " << label_vertex(i) << std::endl;
#endif
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", i, priority(i), 1-pl);
                strategy[i] = -1;
                Q.push(i);
            }
        } else {
            // "winner" attraction
            auto curedge = outs(i);
            for (int to = *curedge; to != -1; to = *++curedge) {
                if (region[to] < r) continue; // not in subgame, or -1/-2
                if (winning[to] == pl) continue; // not attracting
#ifndef NDEBUG
                if (trace >= 2) logger << KC"attracted distraction\033[m " << label_vertex(i) << std::endl;
#endif
                // if (trace) fmt::printf(logger, "attracted %d (%d) to W_%d\n", i, priority(i), 1-pl);
                strategy[i] = to;
                Q.push(i);
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

        R->push_back(cur);
        region[cur] = r;
        winning[cur] = 1-pl;

        // attract to <cur>
        auto curedge = ins(cur);
        for (int from = *curedge; from != -1; from = *++curedge) {
            // if (region[from] == -1) LOGIC_ERROR;
            if (region[from] < r) continue; // not in subgame, or disabled
            if (winning[from] != pl) continue; // already lost

            if (owner(from) != pl) {
                // owned by other
#ifndef NDEBUG
                if (trace >= 2) logger << KC"attracted\033[m " << label_vertex(from) << std::endl;
#endif
                // if (trace) fmt::printf(logger, "attracted %d (%d) to W_%d\n", from, priority(from), 1-pl);
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = cur;
                Q.push(from);
            } else {
                // owned by us
                bool can_escape = false;
                auto curedge = outs(from);
                for (int to = *curedge; to != -1; to = *++curedge) {
                    // if (region[to] == -1) LOGIC_ERROR;
                    if (region[to] < r) continue; // not in subgame, or disabled
                    if (winning[to] != pl) continue; // not an escape
                    can_escape = true;
                    break;
                }
                if (can_escape) continue;
#ifndef NDEBUG
                if (trace >= 2) logger << KC"forced\033[m " << label_vertex(from) << std::endl;
#endif
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", from, priority(from), 1-pl);
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = -1;
                Q.push(from);
            }
        }
    }

    return count;
}

void
ZLKSolver::run()
{
    iterations = 0;

    // allocate and initialize data structures
    region = new int[nodecount()];
    winning = new int[nodecount()];
    strategy = new int[nodecount()];

    std::vector<int> history;
    std::vector<int> W0, W1;
    std::vector<std::vector<int>> levels;

    // initialize arrays
    memset(winning, -1, sizeof(int[nodecount()]));
    memset(strategy, -1, sizeof(int[nodecount()]));

    // get number of nodes and create and initialize inverse array
    max_prio = -1;
    for (int n=nodecount()-1; n>=0; n--) {
        region[n] = disabled[n] ? DIS : BOT;
        if (disabled[n]) continue;
        const int pr = priority(n);
        if (max_prio == -1) {
            max_prio = pr;
            inverse = new int[max_prio+1];
            memset(inverse, -1, sizeof(int[max_prio+1]));
        }
        if (inverse[pr] == -1) inverse[pr] = n;
    }
    if (max_prio == -1) LOGIC_ERROR; // unexpected empty game

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

    bool usePar = lace_workers() != 0;
    // WorkerP* __lace_worker = NULL;
    // Task* __lace_dq_head = NULL;

    if (usePar) {
        // initialize Lace and also allocate space for pvec for each worker
        const int W = lace_workers();
        // __lace_worker = lace_get_worker();
        // __lace_dq_head = lace_get_head(__lace_worker);
        pvec = (par_helper**)malloc(sizeof(par_helper*[W]));
        for (int i=0; i<W; i++) pvec[i] = (par_helper*)malloc(sizeof(par_helper) + sizeof(int[nodecount()]));
    }

    // initialize first level (i, r=0, phase=0)
    levels.push_back(std::vector<int>());
    history.push_back(i);
    history.push_back(next_r++);
    history.push_back(0);

    while (true) {
        // obtain current frame
        const int hsize = history.size();
        if (hsize == 0) break; // no frame on the stack

        std::vector<int> *A = &(*levels.rbegin());
        const int i = history[hsize-3];
        const int r = history[hsize-2];
        const int phase = history[hsize-1];

#ifndef NDEBUG
        const int h = hsize / 3;
        if (i < 0) LOGIC_ERROR; // just a sanity check
        if (region[i] == DIS) LOGIC_ERROR; // just a sanity check
        if (h*3 != hsize or (int)levels.size() != h) LOGIC_ERROR; // just a sanity check
#endif

        /**
         * Get priority and player
         */
        const int pr = priority(i);
        const int pl = pr&1;

#ifndef NDEBUG
        if (trace) fmt::printf(logger, "\n\033[1;37mDepth %d phase %d\033[m: node %d priority %d\n", h-0, phase, i, pr);
#endif

        if (phase == 0) {
            /**
             * We are in the first phase.
             * Compute extended attractor (until inversion). Then recursive step.
             */

#ifndef NDEBUG
            // sanity checks
            if (!W0.empty() or !W1.empty()) LOGIC_ERROR;
#endif

            // attract until inversion and add to A
            int j = usePar ? RUN(attractPar, i, r, A, this) : attractExt(i, r, A);
            // j is now the next i (subgame), or -1 if the subgame is empty

#ifndef NDEBUG
            if (trace) fmt::printf(logger, "attracted %zu nodes\n", A->size());
#endif

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

            int count = 0;
            if (pl == 0) {
                if (!W1.empty()) count = attractLosing(i, r, A, &W1);
                else count = -1;
            } else {
                if (!W0.empty()) count = attractLosing(i, r, A, &W0);
                else count = -1;
            }

#ifndef DEBUG
            if (trace) {
                if (count > 0) fmt::printf(logger, "%d nodes are attracted to losing region\n", count);
                else if (count == 0) fmt::printf(logger, "no nodes are attracted to losing region\n");
                else fmt::printf(logger, "no losing region\n");
            }
#endif

            if (count < 0 or (only_recompute_when_attracted and count == 0)) {
                /**
                 * Nothing attracted to opponent, add A to W0/W1, fix strategies, go up.
                 */

#ifndef NDEBUG
                // Sanity check
                for (int v : *A) { if (winning[v] != pl) LOGIC_ERROR; }
#endif

                auto &Wm = pl == 0 ? W0 : W1;
                Wm.reserve(Wm.size() + A->size());

                for (int v : *A) {
                    /**
                     * Add each element in current region to W0 or W1
                     */
                    Wm.push_back(v);

                    /**
                     * For nodes that are won and controlled by <pl>, check if their strategy needs to be fixed.
                     */
                    if (owner(v) != pl) continue; // not controlled by <pl>
                    if (strategy[v] != -1 && winning[strategy[v]] == pl) continue; // good strategy

                    /**
                     * Strategy of vertex <v> needs to be updated!
                     * We search for a successor of <v> in the subgame won by <pl>
                     */
                    strategy[v] = -1;
                    auto curedge = outs(v);
                    for (int to = *curedge; to != -1; to = *++curedge) {
                        if (region[to] < r) continue; // not in subgame
                        if (winning[to] != pl) continue; // not winning
                        strategy[v] = to;
                        break;
                    }
#ifndef DEBUG
                    if (strategy[v] == -1) LOGIC_ERROR;
#endif
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
                /**
                 * Some nodes are attracted to the opponent.
                 * Reset everything in the winning set of <pl> / region of <pr>
                 */
                int new_i = -1; // will hold lowest node index in *A and W_me
                auto &Wm = pl == 0 ? W0 : W1; // me
                auto &Wo = pl == 0 ? W1 : W0; // other
                for (int v : *A) {
                    if (winning[v] != pl) continue; // only reset for <pl>
                    if (v > new_i) new_i = v;
                    region[v] = BOT;
                }
                for (int v : Wm) {
                    if (winning[v] != pl) continue;
                    if (v > new_i) new_i = v;
                    region[v] = BOT;
                }
                if (new_i == -1) {
                    /**
                     * The remainder is empty.
                     * Go up... keep opponent W, clear our W
                     */
                    Wm.clear();

                    /**
                     * And pop the stack to go up
                     */
                    levels.pop_back();
                    history.pop_back();
                    history.pop_back();
                    history.pop_back();
                } else {
                    /**
                     * The remainder is not empty.
                     * Go down... move opponent W to region, clear W0 and W1
                     */
                    A->swap(Wo);
                    W0.clear();
                    W1.clear();

                    /**
                     * And push the stack to go down
                     */
                    history.back() = 2; // set current phase to 2
                    levels.push_back(std::vector<int>());
                    history.push_back(new_i);
                    history.push_back(next_r++);
                    history.push_back(0);
                }
            }
        } else if (phase == 2) {
            /**
             * After second recursion step
             */

#ifndef NDEBUG
            /**
             * Sanity check
             */
            for (int v : *A) { if (winning[v] == pl) LOGIC_ERROR; }
#endif

            /**
             * Add each element in current region to W0 or W1
             */
            if (pl == 0) {
                W1.reserve(A->size() + W1.size());
                W1.insert(W1.end(), A->begin(), A->end());
            } else {
                W0.reserve(A->size() + W0.size());
                W0.insert(W0.end(), A->begin(), A->end());
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

    if (usePar) {
        const int W = lace_workers();
        for (int i=0; i<W; i++) free(pvec[i]);
        free(pvec);
    }

    // done
    for (int i=0; i<nodecount(); i++) {
        if (region[i] == DIS) continue;
#ifndef NDEBUG
        if (winning[i] == -1) LOGIC_ERROR;
#endif
        Solver::solve(i, winning[i], strategy[i]);
    }

    delete[] region;
    delete[] winning;
    delete[] strategy;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
