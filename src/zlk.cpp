#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <stack>
#include <cassert>

#include "zlk.hpp"
#include "lace.h"
#include "printf.hpp"

namespace pg {

ZLKSolver::ZLKSolver(Oink *oink, Game *game, std::ostream &lgr) : Solver(oink, game, lgr)
{
    // sanity check if the game is properly sorted
    for (int i=1; i<n_nodes; i++) assert(priority[i-1] <= priority[i]);

    // get number of nodes and create and initialize inverse array
    max_prio = game->priority[n_nodes-1];
    inverse = new int[max_prio+1];
    for (int i=0; i<n_nodes; i++) if (!disabled[i]) inverse[priority[i]] = i;
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
    int c = 0;
    par_helper* ours = pvec[LACE_WORKER_ID];

    // attract to <cur>
    const int *_in = s->ins + s->ina[cur];
    for (int from = *_in; from != -1; from = *++_in) {
        if (s->region[from] != -1) continue; // not in subgame, or attracted

        if (s->owner[from] == pl) {
            // owned by same parity, use CAS to claim it
            if (__sync_bool_compare_and_swap(&s->region[from], -1, r)) {
                s->winning[from] = pl;
                s->strategy[from] = cur;
                ours->items[ours->count++] = from;
                SPAWN(attractParT, pl, from, r, s);
                c++;
            }
        } else {
            // owned by other parity
            int* ptr = &s->outcount[from];

            int count = __sync_add_and_fetch(ptr, -1);
            if (count == -2) {
                // we are the first, do add_and_fetch with the count
                count = 1; // compensate for -1
                const int *_out = s->outs + s->outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    int _r = s->region[to];
                    if (_r == -1 or _r == r) count++;
                }
                count = __sync_add_and_fetch(ptr, count);
            }
            if (count == 0) {
                // ok count is now 0, attract...
                // another CAS because we may be competing with attractPar
                if (__sync_bool_compare_and_swap(&s->region[from], -1, r)) {
                    s->winning[from] = pl;
                    s->strategy[from] = -1;
                    ours->items[ours->count++] = from;
                    SPAWN(attractParT, pl, from, r, s);
                    c++;
                }
            }
        }
    }

    while (c) { SYNC(attractParT); c--; }
}

TASK_4(int, attractPar, int, i, int, r, std::vector<int>*, R, ZLKSolver*, s)
{
    const int pr = s->priority[i];
    const int pl = pr & 1;

    // initialize pvec (set count to 0) for all workers
    const int W = lace_workers();
    for (int j=0; j<W; j++) pvec[j]->count = 0;

    par_helper* ours = pvec[LACE_WORKER_ID];
    int spawn_count = 0;

    for (; i>=0; i--) {
        if (s->region[i] != -1) continue; // not in subgame, or attracted
        if ((s->priority[i]&1) != pl) { // search until parity inversion
            // first SYNC on all children, who knows this node may be attracted
            while (spawn_count) { SYNC(attractParT); spawn_count--; }
            // after SYNC, check if node <i> is now attracted.
            if (s->region[i] == -1) break; // not attracted, so we're done!
            else continue; // already done
        }

        // if c != 0, then we compete with attractParT and must use compare and swap
        if (spawn_count == 0) s->region[i] = r;
        else if (!__sync_bool_compare_and_swap(&s->region[i], -1, r)) continue;

        s->winning[i] = pl;
        s->strategy[i] = -1;
        ours->items[ours->count++] = i;
        SPAWN(attractParT, pl, i, r, s);
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
            if (s->trace >= 2) s->logger << "attracted " << x->items[k] << " (" << s->priority[x->items[k]] << ")" << std::endl;
            R->push_back(x->items[k]);
        }
        x->count = 0;
    }

    return i;
}

int
ZLKSolver::attractExt(int i, int r, std::vector<int> *R)
{
    const int pr = priority[i];
    const int pl = pr & 1;

    std::stack<int> q;

    /**
     * Starting at <i>, attract head nodes until "inversion"
     */

    for (; i>=0; i--) {
        if (region[i] != -1) continue; // not in subgame, or attracted
        if ((priority[i]&1) != pl) break; // until parity inversion (Maks Verver optimization)

        // uncomment next line to attract until lower priority instead of until inversion
        // if (priority[i] != pr) break; // until other priority

        region[i] = r;
        winning[i] = pl;
        strategy[i] = -1; // head nodes do not have a strategy yet!
        q.push(i);

        if (trace >= 2) fmt::printf(logger, "head node %d (%d)\n", i, priority[i]);

        while (!q.empty()) {
            int cur = q.top();
            q.pop();

            R->push_back(cur);

            // attract to <cur>
            const int *_in = ins + ina[cur];
            for (int from = *_in; from != -1; from = *++_in) {
                if (from >= i or region[from] != -1) continue; // not in subgame, or attracted

                if (owner[from] == pl) {
                    // owned by same parity
                    region[from] = r;
                    winning[from] = pl;
                    strategy[from] = cur;
                    q.push(from);
                    if (trace >= 2) fmt::printf(logger, "attracted %d (%d)\n", from, priority[from]);
                } else {
                    // owned by other parity
                    int count = outcount[from];
                    if (count == -1) {
                        const int *_out = outs + outa[from];
                        for (int to = *_out; to != -1; to = *++_out) {
                            // also count if region[to] == r!
                            if (region[to] == -1 or region[to] == r) count++;
                        }
                        // no need to --count, already started at -1
                    } else {
                        count--;
                    }
                    outcount[from] = count;
                    if (count == 0) {
                        region[from] = r;
                        winning[from] = pl;
                        strategy[from] = -1;
                        q.push(from);
                        if (trace >= 2) fmt::printf(logger, "forced %d (%d)\n", from, priority[from]);
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
ZLKSolver::attractLosing(int i, int r, std::vector<int> *S, std::vector<int> *R)
{
    int count = 0;

    const int pr = priority[i];
    const int pl = pr & 1;

    // NOTE: this algorithm could be improved using an "out counter"

    std::queue<int> q;

#ifndef NDEBUG
    for (int i : *S) if (winning[i] != pl) LOGIC_ERROR;
#endif

    /**
     * First check all region nodes...
     * In reality, we just want to check the "head nodes", because all other nodes are attracted
     * to the head nodes and cannot be attracted to the opponent directly.
     * But we do not record which nodes are head nodes.
     */
    for (int i : *S) {
        // check if the node is attracted
        if (owner[i] == pl) {
            // "loser" attraction
            bool can_escape = false;
            const int *_out = outs + outa[i];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] < r) continue; // not in subgame, or -1/-2
                if (winning[to] != pl) continue; // not an escape
                can_escape = true;
                break;
            }
            if (!can_escape) {
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", i, priority[i], 1-pl);
                region[i] = r;
                winning[i] = 1-pl;
                strategy[i] = -1;
                q.push(i);
            }
        } else {
            // "winner" attraction
            const int *_out = outs + outa[i];
            for (int to = *_out; to != -1; to = *++_out) {
                if (region[to] < r) continue; // not in subgame, or -1/-2
                if (winning[to] == pl) continue; // not attracting
                // if (trace) fmt::printf(logger, "attracted %d (%d) to W_%d\n", i, priority[i], 1-pl);
                region[i] = r;
                winning[i] = 1-pl;
                strategy[i] = to;
                q.push(i);
                break;
            }
        }
    }

    /**
     * Now attract anything in this region/subregions of <pl> to 1-<pl>
     */

    while (!q.empty()) {
        int cur = q.front();
        q.pop();
        ++count;

        R->push_back(cur);

        // attract to <cur>
        const int *_in = ins + ina[cur];
        for (int from = *_in; from != -1; from = *++_in) {
            // if (region[from] == -1) LOGIC_ERROR;
            if (region[from] < r) continue; // not in subgame, or disabled
            if (winning[from] != pl) continue; // already lost

            if (owner[from] != pl) {
                // owned by other
                // if (trace) fmt::printf(logger, "attracted %d (%d) to W_%d\n", from, priority[from], 1-pl);
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = cur;
                q.push(from);
            } else {
                // owned by us
                bool can_escape = false;
                const int *_out = outs + outa[from];
                for (int to = *_out; to != -1; to = *++_out) {
                    // if (region[to] == -1) LOGIC_ERROR;
                    if (region[to] < r) continue; // not in subgame, or disabled
                    if (winning[to] != pl) continue; // not an escape
                    can_escape = true;
                    break;
                }
                if (can_escape) continue;
                // if (trace) fmt::printf(logger, "forced %d (%d) to W_%d\n", from, priority[from], 1-pl);
                region[from] = r;
                winning[from] = 1-pl;
                strategy[from] = -1;
                q.push(from);
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
    region = new int[n_nodes];
    winning = new int[n_nodes];
    strategy = new int[n_nodes];
    outcount = new int[n_nodes];

    std::vector<int> history;
    std::vector<int> W0, W1;
    std::vector<std::vector<int>> levels;

    // initialize arrays
    for (int i=0; i<n_nodes; i++) region[i] = disabled[i] ? -2 : -1;
    for (int i=0; i<n_nodes; i++) winning[i] = -1;
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;
    for (int i=0; i<n_nodes; i++) outcount[i] = -1;

    // start the loop at the last node (highest priority) and at depth 0
    int i = n_nodes - 1;
    int next_r = 0;

    // skip all disabled nodes
    while (i >= 0 && region[i] == -2) i--;
    if (i == -1) LOGIC_ERROR; // wut? empty game. pls dont call.

    bool usePar = lace_workers() != 0;
    WorkerP* __lace_worker = NULL;
    Task* __lace_dq_head = NULL;

    if (usePar) {
        const int W = lace_workers();
        __lace_worker = lace_get_worker();
        __lace_dq_head = lace_get_head(__lace_worker);
        pvec = (par_helper**)malloc(sizeof(par_helper*[W]));
        for (int i=0; i<W; i++) pvec[i] = (par_helper*)malloc(sizeof(par_helper)+sizeof(int[n_nodes]));
    }

    levels.push_back(std::vector<int>());
    history.push_back(i); history.push_back(next_r++); history.push_back(0); // i, r=0, phase=0

    while (true) {
        const int hsize = history.size();
        if (hsize == 0) break;

        std::vector<int> *A = &(*levels.rbegin());
        const int i = history[hsize-3];
        const int r = history[hsize-2];
        const int phase = history[hsize-1];
        const int h = hsize / 3; // history.size() / 2;

#ifndef NDEBUG
        if (i < 0) LOGIC_ERROR; // just a sanity check
        if (region[i] == -2) LOGIC_ERROR; // just a sanity check
        if (h*3 != hsize or (int)levels.size() != h) LOGIC_ERROR; // just a sanity check
#endif

        /**
         * Get priority and player
         */
        const int pr = priority[i];
        const int pl = pr&1;

        if (trace) fmt::printf(logger, "\n\033[1mDepth %d phase %d\033[m: node %d priority %d\n", h-0, phase, i, pr);

        if (phase == 0) {
            /**
             * We are in the first phase.
             * Compute extended attractor (until inversion). Then recursive step.
             */

            // sanity checks
            if (!W0.empty() or !W1.empty()) LOGIC_ERROR;

            // attract until inversion and add to A
            int j = usePar ? CALL(attractPar, i, r, A, this) : attractExt(i, r, A);
            if (trace) fmt::printf(logger, "attracted %zu nodes\n", A->size());

            // count number of iterations
            ++iterations;

            // update phase to 1
            history.back() = 1;

            if (j != -1) {
                // go recursive (not empty)
                levels.push_back(std::vector<int>());
                history.push_back(j); history.push_back(next_r++); history.push_back(0);
            }
        } else if (phase == 1) {
            /**
             * After first recursion step
             */

            if (trace) fmt::printf(logger, "current level contains %zu nodes\n", A->size());

            /**
             * Compute part of region that attracts to the other player...
             */

            int count = 0;
            if (pl == 0) {
                if (!W1.empty()) count = attractLosing(i, r, A, &W1);
            } else {
                if (!W0.empty()) count = attractLosing(i, r, A, &W0);
            }

            if (trace) {
                if (count != 0) fmt::printf(logger, "%d nodes are attracted to losing region\n", count);
                else if (pl == 0 ? !W1.empty() : !W0.empty()) fmt::printf(logger, "no nodes are attracted to losing region\n");
                else fmt::printf(logger, "no losing region\n");
            }

            if (count == 0) {
                /**
                 * Nothing attracted to opponent, add A to W0/W1, fix strategies, go up.
                 */

#ifndef NDEBUG
                /**
                 * Sanity check
                 */
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
                     * For nodes that are won and controlled by <pl>,
                     * check if their strategy needs to be fixed.
                     */
                    if (owner[v] != pl) continue; // not controlled by <pl>
                    if (strategy[v] != -1 && winning[strategy[v]] == pl) continue; // good strategy

                    /**
                     * Strategy needs to be updated!
                     * We do this by searching for an edge to a node in the subgame won by <pl>
                     */
                    strategy[v] = -1;
                    const int *_out = outs + outa[v];
                    for (int to = *_out; to != -1; to = *++_out) {
                        if (region[to] < r) continue; // not in subgame
                        if (winning[to] != pl) continue; // not winning
                        strategy[v] = to;
                        break;
                    }
                    if (strategy[v] == -1) LOGIC_ERROR;
                }

                /**
                 * The strategy has been updated. Nodes are now in W0/W1.
                 * Finally, pop the stack and go up...
                 */
                levels.pop_back();
                history.pop_back(); history.pop_back(); history.pop_back();
            } else {
                /**
                 * Some nodes attracted to opponent.
                 * Reset everything in winning set of <pl> / region of <pr>
                 */
                int new_i = -1; // will hold lowest node index in *A and W_me
                auto &Wm = pl == 0 ? W0 : W1; // me
                auto &Wo = pl == 0 ? W1 : W0; // other
                for (int v : *A) {
                    if (winning[v] != pl) continue; // only reset for <pl>
                    if (v > new_i) new_i = v;
                    region[v] = -1;
                    outcount[v] = -1;
                }
                for (int v : Wm) {
                    if (winning[v] != pl) continue;
                    if (v > new_i) new_i = v;
                    region[v] = -1;
                    outcount[v] = -1;
                }
                if (new_i == -1) {
                    /**
                     * The remainder is empty.
                     * Going up... keep opponent W, clear our W
                     */
                    Wm.clear();

                    /**
                     * And pop the stack to go up
                     */
                    levels.pop_back();
                    history.pop_back(); history.pop_back(); history.pop_back();
                } else {
                    /**
                     * The remainder is not empty.
                     * Going down... move opponent W to region, clear W0 and W1
                     */
                    A->swap(Wo);
                    W0.clear();
                    W1.clear();

                    /**
                     * And push the stack to go down
                     */
                    history.back() = 2;
                    levels.push_back(std::vector<int>());
                    history.push_back(new_i); history.push_back(next_r++); history.push_back(0);
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
             * Finally, pop the stack and go up...
             */
            levels.pop_back();
            history.pop_back(); history.pop_back(); history.pop_back();
        }
    }

    if (usePar) {
        const int W = lace_workers();
        for (int i=0; i<W; i++) free(pvec[i]);
        free(pvec);
    }

    // done
    for (int i=0; i<n_nodes; i++) {
        if (region[i] == -2) continue;
        if (winning[i] == -1) LOGIC_ERROR;
        int winner = game->dominion[i] = winning[i];
        if (winner == game->owner[i]) {
            if (strategy[i] == -1) LOGIC_ERROR;
            if (winning[strategy[i]] != winner) LOGIC_ERROR;
            game->strategy[i] = strategy[i];
        }
    }

    delete[] region;
    delete[] winning;
    delete[] strategy;
    delete[] outcount;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
