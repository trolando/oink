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

#include <cassert>

#include "psi.hpp"
#include "lace.h"
#include "printf.hpp"

namespace pg {

// Yeah, globals.
// But I don't want to copy the pointers in all Lace tasks.
// And we're not using multiple PSISolver objects anyway.

static int k;         // k := 1+pr(G)
static int *str;      // strategy (player/opponent)
static int *halt;     // vertices that Even (Odd) halts _in front of_
static int *val;      // current value, k-tuple for each vertex
static int *done;     // 0=unvisited 1=visited 2=cycle 3=disabled/won
static int *won;      // if won by primary player
static int *first_in; // helper for linked-list style in-edges
static int *next_in;  // helper for linked-list style in-edges

PSISolver::PSISolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

PSISolver::~PSISolver()
{
}

/**
 * Returns true if strategy valuation of "a" is less than of "b" for the Even player
 * Also ensures that won or even_cycle are Top
 */
bool
PSISolver::si_val_less(int a, int b)
{
    // if a == b, then of course return false
    if (a == b) return false;
    // if a is won or on even cycle, then "a >= b"
    if (a != -1 && (won[a] || done[a] == 2)) return false;
    // if b is won or on even cycle and a is not, then "a < b"
    if (b != -1 && (won[b] || done[b] == 2)) return true;
    // neither a/b are won/on even cycle
    int *a_val = a == -1 ? NULL : val + k*a;
    int *b_val = b == -1 ? NULL : val + k*b;
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

int
PSISolver::si_top_val(int a)
{
    // find highest priority
    for (int i=k-1; i>=0; i--) {
        if (val[k*a+i] != 0) return i;
    }
    return -1;
}

/**
 * Recursively update the value of <v> and its ancestors according to the
 * currently chosen strategy
 */
VOID_TASK_2(compute_val, int, v, PSISolver *, s)
{
    s->compute_val(__lace_worker, __lace_dq_head, v);
}

void
PSISolver::compute_val(WorkerP* __lace_worker, Task* __lace_dq_head, int v)
{
    // mark node as visited
    done[v] = 1;

    // compute valuation of current node
    int st = str[v];
    int *dst = val + k*v;
    if (st == -1 or halt[st]) memset(dst, 0, sizeof(int[k]));
    else memcpy(dst, val + k*st, sizeof(int[k])); // copy from successor
    dst[priority(v)]++;

    // recursively update predecessor positions
    int count = 0;
    int from = first_in[v];
    while (from != -1) {
        int next = next_in[from];
        if (next != -1) {
            SPAWN(compute_val, from, this);
            count++;
        } else {
            CALL(compute_val, from, this);
            break;
        }
        from = next;
    }

    while (count--) {
        SYNC(compute_val);
    }
} 

/**
 * Fill first_in and next_in based on strategies
 */
VOID_TASK_2(set_in, int, begin, int, count)
{
    // some cut-off point...
    if (count <= 64) {
        for (int i=0; i<count; i++) {
            int n = begin+i;
            if (done[n] != 3) {
                int s = str[n];
                if (s != -1 and !halt[s]) {
                    // atomically set first_in[s] to n, and return the previous value
                    next_in[n] = __sync_lock_test_and_set(first_in+s, n);
                }
            }
        }
    } else {
        SPAWN(set_in, begin+count/2, count-count/2);
        CALL(set_in, begin, count/2);
        SYNC(set_in);
    }
}

/**
 * Resets "done" array before recomputing valuations. Only resets if done equals 1 or 2
 */
VOID_TASK_2(reset_done, int, begin, int, count)
{
    // some cut-off point...
    if (count <= 64) {
        memset(first_in+begin, -1, sizeof(int[count]));
        for (int i=0; i<count; i++) {
            int n = begin+i;
            if (done[n] != 3) {
                int s = str[n];
                if (s != -1 and !halt[s] and done[n] != 2) done[n] = 2;
            }
            // TODO printf("%d: %d\n", n, done[n]);
        }
    } else {
        SPAWN(reset_done, begin+count/2, count-count/2);
        CALL(reset_done, begin, count/2);
        SYNC(reset_done);
    }
}

/**
 * Master function for (re)computing the valuations of all unsolved nodes
 */
VOID_TASK_1(compute_all_val, PSISolver*, s)
{
    s->compute_all_val(__lace_worker, __lace_dq_head);
}

void
PSISolver::compute_all_val(WorkerP* __lace_worker, Task* __lace_dq_head)
{
    // reset "done" (for nodes that are not disabled or won)
    CALL(reset_done, 0, nodecount());
    CALL(set_in, 0, nodecount());
    // for all unsolved enabled nodes that go to sink, run compute val
    int count = 0;
    for (int i=0; i<nodecount(); i++) {
        if (done[i] == 3) continue;
        if (str[i] == -1 or halt[str[i]]) { // str[i] is -2 for disabled and not -1 for won
            SPAWN(compute_val, i, this);
            count++;
        }
    }
    while (count--) {
        SYNC(compute_val);
    }
}

void
PSISolver::compute_vals_seq(void)
{
    std::vector<int> q;

    memset(first_in, -1, sizeof(int[nodecount()]));

    for (int n=0; n<nodecount(); n++) {
        if (done[n] == 3) continue; // disabled
        int s = str[n];
        if (s == -1 or halt[s]) { // either we halt or the dummy halts...
            // assuming that valuation is properly initialized...
            q.push_back(n);
        } else {
            next_in[n] = first_in[s];
            first_in[s] = n;
            if (done[n] != 2) done[n] = 2; // mark "2", if not updated, then it's on a cycle
        }
    }

    while (!q.empty()) {
        int v = q.back();
        q.pop_back();

        // set valuation
        int *val_v = val + k*v;
        int s = str[v];
        if (s == -1 or halt[s]) memset(val_v, 0, sizeof(int[k]));
        else memcpy(val_v, val+k*s, sizeof(int[k]));
        val_v[priority(v)]++;
        done[v] = 1; // mark "1" as processed and not on a cycle

        // add predecessors
        int from = first_in[v];
        while (from != -1) {
            q.push_back(from);
            from = next_in[from];
        }
    }
}

/**
 * Mark all nodes with even_cycle (done==2) as won (done==3, won==1)
 */
TASK_3(int, mark_solved_rec, PSISolver*, s, int, begin, int, count)
{
    // some cut-off point...
    if (count < 64) {
        int res = 0;
        for (int i=0; i<count; i++) {
            int n = begin+i;
            if (done[n] == 2) { // done[n] == 3 proxies disabled and won
                won[n] = 1;
                done[n] = 3; // mark as won
                res++;
            }
        }
        return res;
    } else {
        SPAWN(mark_solved_rec, s, begin+count/2, count-count/2);
        int res = CALL(mark_solved_rec, s, begin, count/2);
        return res + SYNC(mark_solved_rec);
    }
}

int
PSISolver::mark_solved_seq(void)
{
    int res = 0;
    for (int n=0; n<nodecount(); n++) {
        if (done[n] == 2) { // done[n] == 3 proxies disabled and won
            won[n] = 1;
            done[n] = 3; // mark as won
            res++;
        }
    }
    return res;
}

/**
 * Check if Even still wants to halt the path...
 */
TASK_3(int, switch_halting, PSISolver*, s, int, begin, int, count)
{
    // some cut-off point...
    if (count <= 64) {
        int res = 0;

        for (int i=0; i<count; i++) {
            int n = begin+i;

            if (halt[n] and s->si_val_less(-1, n)) {
                halt[n] = 0; // stop halting
                res++;
            }
        }

        return res;
    } else {
        SPAWN(switch_halting, s, begin+count/2, count-count/2);
        int res = CALL(switch_halting, s, begin, count/2);
        res += SYNC(switch_halting);
        return res;
    }
}

/**
 * Implementation of greedy-all-switches strategy
 */
TASK_4(int, switch_strategy, PSISolver*, s, int, pl, int, begin, int, count)
{
    return s->switch_strategy(__lace_worker, __lace_dq_head, pl, begin, count);
}

int
PSISolver::switch_strategy(WorkerP* __lace_worker, Task* __lace_dq_head, int pl, int begin, int count)
{
    // some cut-off point...
    if (count < 64) {
        int res = 0;
        for (int i=0; i<count; i++) {
            int n = begin+i;

            if (done[n] == 3) continue; // skip "disabled or won"
            if (done[n] == 0) LOGIC_ERROR; // expecting done==1 or done==2
            if (owner(n) != pl) continue; // only change strategy if owner

            int cur_strat = str[n];
            for (auto curedge = outs(n); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue; // skip strategy to disabled
                if (to == cur_strat) continue; // skip strategy to same
                if (pl == 0) {
                    // improving for player Even
                    if (si_val_less(halt[cur_strat] ? -1 : cur_strat, halt[to] ? -1 : to)) {
                        str[n] = cur_strat = to;
                        res++;
                    }
                } else {
                    // improving for player Odd
                    if (si_val_less(halt[to] ? -1 : to, halt[cur_strat] ? -1 : cur_strat)) {
                        str[n] = cur_strat = to;
                        res++;
                    }
                }
            }
        }
        return res;
    } else {
        SPAWN(switch_strategy, this, pl, begin+count/2, count-count/2);
        int res = CALL(switch_strategy, this, pl, begin, count/2);
        res += SYNC(switch_strategy);
        return res;
    }
}

int
PSISolver::switch_strategy_seq(int pl)
{
    int res = 0;

    for (int n=0; n<nodecount(); n++) {
        if (done[n] == 3) continue; // skip "disabled or won"
        if (done[n] == 0) LOGIC_ERROR; // expecting done==1 or done==2
        if (owner(n) != pl) continue; // only change strategy for vertices of player <pl>

        bool changed = false;
        int cur_strat = str[n];
        for (auto curedge = outs(n); *curedge != -1; curedge++) {
            int to = *curedge;
            if (disabled[to]) continue; // skip strategy to disabled
            if (to == cur_strat) continue; // skip strategy to same
            if (pl == 0) {
                // improving for player Even
                if (si_val_less(halt[cur_strat] ? -1 : cur_strat, halt[to] ? -1 : to)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdate even strategy\033[m \033[1;33m" << label_vertex(n) << "\033[m from \033[1;33m" << label_vertex(cur_strat) << "\033[m (" << si_top_val(cur_strat) << (halt[cur_strat]?")H":")") << " to \033[1;33m" << label_vertex(to) << "\033[m (" << si_top_val(to) << (halt[to] ? ")H" : ")") << std::endl;
                    }
#endif
                    str[n] = cur_strat = to;
                    changed = true;
                }
            } else {
                // improving for player Odd
                if (si_val_less(halt[to] ? -1 : to, halt[cur_strat] ? -1 : cur_strat)) {
#ifndef NDEBUG
                    if (trace >= 2) {
                        logger << "\033[1;37mupdate odd strategy\033[m \033[1;33m" << label_vertex(n) << "\033[m from \033[1;33m" << label_vertex(cur_strat) << "\033[m (" << si_top_val(cur_strat) << (halt[cur_strat]?")H":")") << " to \033[1;33m" << label_vertex(to) << "\033[m (" << si_top_val(to) << (halt[to] ? ")H" : ")") << std::endl;
                    }
#endif
                    str[n] = cur_strat = to;
                    changed = true;
                }
            }
        }
        if (changed) res++;
    }

    if (pl == 0) {
        for (int n=0; n<nodecount(); n++) {
            if (halt[n] and si_val_less(-1, n)) {
                halt[n] = 0; // stop halting
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "\033[1;37mno longer halt\033[m before vertex " << label_vertex(n) << std::endl;
                }
#endif
                res++;
            }
        }
    }

    return res;
}

/**
 * Print current state (for debugging purposes) i.e. all valuations and strategies...
 */
void
PSISolver::print_debug()
{
    fmt::printf(logger, "\033[1mState\033[m\n");
    for (int i=0; i<nodecount(); i++) {
        if (done[i] == 3) continue; // disabled or won
        if (disabled[i]) LOGIC_ERROR;
        logger << "vertex " << label_vertex(i) << ": [";
        for (int j=0; j<k; j++) logger << (j?" ":"") << val[k*i+j];
        logger << "] ";
        if (done[i] == 2) logger << "c ";
        if (halt[i]) logger << "h ";
        if (owner(i)) logger << "\033[1;35;160m=>\033[m ";
        else logger << "\033[1m=>\033[m ";
        if (str[i] == -1) logger << "-";
        else logger << label_vertex(str[i]);
        logger << std::endl;
    }
}

VOID_TASK_1(psi_run_par, PSISolver*, _this)
{
    return _this->run_par(__lace_worker, __lace_dq_head);
}

void
PSISolver::run_par(WorkerP* __lace_worker, Task* __lace_dq_head)
{
    for (;;) {
        ++this->major;
        if (this->trace) fmt::printf(this->logger, "\033[1;38;5;208mMajor iteration %d\033[m\n", this->major);
        for (;;) {
            ++this->minor;
            CALL(compute_all_val, this);                            // update valuation
#ifndef NDEBUG
            if (this->trace >= 3) this->print_debug();
#endif
            int count = CALL(switch_strategy, this, 1, 0, this->nodecount()); // switch strategies
            if (this->trace) fmt::printf(this->logger, "%d changed strategies for Odd\n", count);
            if (count == 0) break;                                  // if nothing left, done
        }
        /* print selected strategies for the odd player */
        if (this->trace) {
            for (int n=0; n<this->nodecount(); n++) {
                if (this->disabled[n]) continue;
                if (this->owner(n) == 0) continue;
                this->logger << "Odd plays from \033[1;33m" << this->label_vertex(n) << "\033[m to \033[1;33m" << this->label_vertex(str[n]) << "\033[m (";
                if (halt[str[n]]) this->logger << "H";
                else this->logger << this->si_top_val(str[n]);
                this->logger << ")" << std::endl;
            }
        }
        int solved = CALL(mark_solved_rec, this, 0, this->nodecount());       // mark nodes won by Even
        if (this->trace) fmt::printf(this->logger, "%d nodes marked as won by Even\n", solved);
        int count = CALL(switch_strategy, this, 0, 0, this->nodecount());     // switch strategies
        count += CALL(switch_halting, this, 0, this->nodecount());     // switch halting strategies
        if (this->trace) fmt::printf(this->logger, "%d changed strategies for Even\n", count);
        if (count == 0) break;                                      // if nothing left, done
    }
}


/**
 * Run the parallel strategy improvement solver
 */
void
PSISolver::run()
{
    // determine k as highest priority + 1
    k = 0;
    for (int i=0; i<nodecount(); i++) if (!disabled[i] && priority(i)>k) k = priority(i);
    k++;

    // now create the data structure
    val = new int[k*nodecount()];
    str = new int[nodecount()];
    halt = new int[nodecount()];
    done = new int[nodecount()];
    won = new int[nodecount()];

    first_in = new int[nodecount()];
    next_in = new int[nodecount()];

    // initialize the datastructure
    for (int i=0; i<nodecount(); i++) {
        if (disabled[i]) {
            str[i] = -2; // set to "disabled sink"
            done[i] = 3; // set to "disabled/won"
            continue;
        } else {
            // select first available edge
            int to = -1;
            for (auto curedge = outs(i); *curedge != -1; curedge++) {
                int x = *curedge;
                if (disabled[x]) continue;
                to = x;
                break;
            }
            if (to == -1) LOGIC_ERROR;
            str[i] = to;
            halt[i] = 1;
            won[i] = 0;
            done[i] = 0;
        }
    }

    if (lace_workers() == 0) {
        for (;;) {
            ++major;
            if (trace) fmt::printf(logger, "\033[1;38;5;208mMajor iteration %d\033[m\n", major);
            for (;;) {
                ++minor;
                compute_vals_seq();
#ifndef NDEBUG
                if (trace >= 3) print_debug();
#endif
                int count = switch_strategy_seq(1);
                if (trace) fmt::printf(logger, "%d changed strategies for Odd\n", count);
                if (count == 0) break;                                  // if nothing left, done
            }
            /* print selected strategies for the odd player */
            if (trace) {
                for (int n=0; n<nodecount(); n++) {
                    if (disabled[n]) continue;
                    if (owner(n) == 0) continue;
                    logger << "Odd plays from \033[1;33m" << label_vertex(n) << "\033[m to \033[1;33m" << label_vertex(str[n]) << "\033[m (";
                    if (halt[str[n]]) logger << "H";
                    else logger << si_top_val(str[n]);
                    logger << ")" << std::endl;
                }
            }
            int solved = mark_solved_seq(); // mark nodes won by Even
            if (trace) fmt::printf(logger, "%d nodes marked as won by Even\n", solved);
            int count = switch_strategy_seq(0);
            if (trace) fmt::printf(logger, "%d changed strategies for Even\n", count);
            if (count == 0) break;                                      // if nothing left, done
        }
    } else {
        RUN(psi_run_par, this);
    }

    // Now set dominions and derive strategy for odd.
    for (int i=0; i<nodecount(); i++) {
        if (disabled[i]) continue;
        bool winner = won[i] ? 0 : 1;
        Solver::solve(i, winner, owner(i) == winner ? str[i] : -1);
    }

    delete[] val;
    delete[] str;
    delete[] halt;
    delete[] done;
    delete[] won;
    delete[] first_in;
    delete[] next_in;

    logger << "solved with " << major << " major iterations, " << minor << " minor iterations." << std::endl;
}

}
