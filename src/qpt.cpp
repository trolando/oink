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

#include <iomanip>

#include "qpt.hpp"

namespace pg {

QPTSolver::QPTSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

QPTSolver::~QPTSolver()
{
}

/**
 * Returns true if a progress measure "a" is less than "b",
 * given progress measures of <k> components for player <pl>
 */
static bool
pm_less(int* a, int* b, int k, int pl)
{
    // remember: _ < 5 < 3 < 1 < 0 < 2 < 4 < 6 (even measures, solving for odd)
    // remember: _ < 6 < 4 < 2 < 0 < 1 < 3 < 5 (odd measures, solving for even)

    for (int i=0; i<k; i++) {
        // if equal, then continue with next
        if (a[i] == b[i]) continue;
        // cases where a[i] is _ or b[i] is _
        if (a[i] == -1) return true;
        if (b[i] == -1) return false;
        // case distinction based on parity
        const bool a_is_pl = (a[i]&1) == pl;
        const bool b_is_pl = (b[i]&1) == pl;
        if (a_is_pl) {
            if (b_is_pl) return a[i] < b[i];
            else return false;
        } else {
            if (b_is_pl) return true;
            else return a[i] > b[i];
        }
    }
    // they are equal, so not less...
    return false;
}


/**
 * Compute the value of a progress measure.
 */
static unsigned long
pm_val(int* pm, int k, int pl)
{
    unsigned long res=0;
    for (int i=0; i<k; i++) {
        res <<= 1;
        if (pm[i] != -1 and (pm[i]&1)==pl) res++;
    }
    return res;
}


/**
 * Print out a progress measure <pm> (length <k>) to the given stream
 */
static void
pm_stream(std::ostream &out, int *pm, int k)
{
    for (int i=0; i<k; i++) {
        if (pm[i] != -1) out << " " << std::setfill('0') << std::setw(2) << pm[i];
        else out << " __";
    }
}


/**
 * Compute updated progress measure <dst> after seeing priority "d",
 * given progress measure <src> with length <k> recording <pl> measures.
 * Returns -1 if dst == src, or the updated component index.
 */
static int
up(int *dst, int *src, int d, int k, int pl)
{
    // check if this is a "won" witness
    // a pm is won if the highest component records a value of parity <pl>
    if (src[0] != -1 and (src[0]&1) == pl) {
        for (int i=0; i<k; i++) dst[i] = src[i];
        return -1;
    }

    // first compute j for Lemma 3.3
    int j33 = -1;
    for (int i=k-1; i>=0; i--) {
        // start from the end to find longest sequence of <pl> parity
        if (src[i] != -1 and (src[i]&1) == pl) continue;
        // ok, src[i] is not of <pl> parity, test the assumption on the other half
        // (all higher are _ or higher than d)
        j33 = i;
        for (int j=i-1; j>=0; j--) {
            if (src[j] != -1 and src[j] < d) {
                // not good!
                j33 = -1;
                break;
            }
        }
        break;
    }

    // then compute j for Lemma 3.4
    int j34 = -1;
    for (int i=0; i<k; i++) {
        // start from the begin, find first index where the Lemma applies
        if (src[i] != -1 and src[i] < d) {
            j34 = i;
            break;
        }
    }

    // determine best index min(j33, j34) or -1
    int j = j33 == -1 ? j34 : (j34 == -1 ? j33 : (j33 < j34 ? j33 : j34));

    // compute resulting progress measure
    for (int i=0; i<j; i++) dst[i] = src[i]; // same before j
    dst[j] = d;
    for (int i=j+1; i<k; i++) dst[i] = -1; // _ after j

    return j;
}


/**
 * Compute smallest larger progress measure with dst[0] set to _
 */
static bool
bump(int *dst, int *src, int k, int max, int max_opp, int pl)
{
    // reminder: _ 5 3 1 0 2 4 6...
    // reminder: _ 6 4 2 0 1 3 5...

    for (int i=k-2; i>=0; i--) {
        if (src[i] == -1) {
            // if the current is _, then the smallest increase is max_opp
            dst[i] = max_opp;
            for (int z=0; z<i; z++) dst[z] = src[z];
            for (int z=i+1; z<k; z++) dst[z] = -1;
            return true;
        } else if ((src[i]&1) != pl) {
            // if the current is of opponent parity, then subtract 2 OR change parity
            if (src[i] > (1-pl)) {
                // if odd != 1 / even != 0 then we can just subtract 2 and we are done
                dst[i] = src[i]-2;
                for (int z=0; z<i; z++) dst[z] = src[z];
                for (int z=i+1; z<k; z++) dst[z] = -1;
                return true;
            } else if (i > 0 and src[i-1] >= pl) {
                // if is 1/0, then can set to 0/1 if the previous is set and we're not highest...
                dst[i] = pl;
                for (int z=0; z<i; z++) dst[z] = src[z];
                for (int z=i+1; z<k; z++) dst[z] = -1;
                return true;
            }
        } else {
            // if it's of parity <pl>, and the next one is higher, then increase and done
            // probably the only one that is relevant, since this algorithm is used
            // when an "even" sequence must be bumped
            if ((src[i]+2) <= max and i > 0 and src[i-1] >= (src[i]+2)) {
                dst[i] = src[i]+2;
                for (int z=0; z<i; z++) dst[z] = src[z];
                for (int z=i+1; z<k; z++) dst[z] = -1;
                return true;
            }
        }
    }
    // did not find a higher progress measure
    return false;
}


/**
 * Perform antagonistic update.
 * Return "true" if dst != src, "false" otherwise.
 */
static bool
au(int *dst, int *src, int d, int k, int max, int maxopp, int pl)
{
    // first compute normal update
    int j = up(dst, src, d, k, pl);
    if (j == -1) return false; // no difference, already smallest...

    // then compute antagonistic update and update the antagonistic measure
    int tmp[k];
    if (!bump(tmp, src, k, max, maxopp, pl)) return true; // no bump possible; but dst != src so return true
    up(tmp, tmp, d, k, pl);

    // use smallest updated measure
    if (pm_less(tmp, dst, k, pl)) {
        for (int i=0; i<k; i++) dst[i] = tmp[i];
    }

    // finally, check if different
    for (int i=0; i<k; i++) if (dst[i] != src[i]) return true;
    return false;
}


bool
QPTSolver::lift(int v, int target, int pl)
{
    // check if already Top
    int * const pm = pm_nodes + k*v; // obtain ptr to current progress measure
    if (pm[0] != -1 and (pm[0]&1) == pl) return false;

    const int pr = priority[v];
    int tmp[k], res[k];

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1mupdating vertex " << label_vertex(v) << (owner[v]?" (odd)":" (even)") << "\033[m with current measure";
        pm_stream(logger, pm, k);
        logger << std::endl;
    }
#endif

    // First special cases if we are lifting triggered by a given <target>.

    if (1 and target != -1) {
        if (owner[v] == pl) {
            // if owned by <pl>, just check if target is better.
            au(tmp, pm_nodes + k*target, pr, k, max, pl ? max_even : max_odd, pl);
            if (pm_val(tmp, k, pl) > goal) tmp[0] = max; // lol
            if (pm_less(pm, tmp, k, pl)) {
                // target is better, update, report, return
                for (int i=0; i<k; i++) pm[i] = tmp[i];
#ifndef NDEBUG
                if (trace) {
                    logger << "\033[1;32mnew measure\033[m of " << label_vertex(v) << ":";
                    pm_stream(logger, pm, k);
                    logger << std::endl;
                }
#endif
                return true;
            } else {
                // target not better, not checking other directions
                return false;
            }
        } else {
            // if owned by ~<pl>, test if strategy is target
            if (strategy[v] != -1 and target != strategy[v]) {
                return false; // no need to do anything... strategy is still same
            }
        }
    }

    // Compute best measure by going over all successors

    bool first = true;
    for (int to : out[v]) {
        if (disabled[to]) continue;
        au(tmp, pm_nodes + k*to, pr, k, max, pl ? max_even : max_odd, pl);
        if (pm_val(tmp, k, pl) > goal) tmp[0] = max; // goal reached, lift to Top
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "to successor " << label_vertex(to) << ":";
            pm_stream(logger, tmp, k);
            logger << std::endl;
        }
#endif
        if (first) {
            first = false;
            for (int i=0; i<k; i++) res[i] = tmp[i];
            strategy[v] = to;
        } else if (owner[v] == pl) {
            // we want the max
            if (pm_less(res, tmp, k, pl)) {
                for (int i=0; i<k; i++) res[i] = tmp[i];
                strategy[v] = to;
            }
        } else {
            // we want the min
            if (pm_less(tmp, res, k, pl)) {
                for (int i=0; i<k; i++) res[i] = tmp[i];
                strategy[v] = to;
            }
        }
    }

    // now "res" contains the best au
    if (pm_less(pm, res, k, pl)) {
        // ok, it's an update
        for (int i=0; i<k; i++) pm[i] = res[i];
#ifndef NDEBUG
        if (trace) {
            logger << "\033[1;32mnew measure\033[m of " << label_vertex(v) << ":";
            pm_stream(logger, pm, k);
            logger << std::endl;
        }
#endif
        return true;
    } else{
        return false;
    }
}


void
QPTSolver::liftloop(int pl)
{
    /**
     * Initialize/reset progress measures / strategy
     */
    k = pl ? k1 : k0;
    max = pl ? max_odd : max_even;
    goal = pl ? goal1 : goal0;
    for (int i=0; i<k*n_nodes; i++) pm_nodes[i] = -1; // initialize all to _
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    /**
     * Run first lifting loop
     */
    for (int n=n_nodes-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        if (lift(n, -1, pl)) {
            lift_count++;
            for (int from : in[n]) {
                if (disabled[from]) continue;
                lift_attempt++;
                if (lift(from, n, pl)) {
                    lift_count++;
                    todo_push(from);
                }
            }
        }
    }

    /**
     * Lift until fixed point
     */
    while (!todo.empty()) {
        int n = todo_pop();
        for (int from : in[n]) {
            if (disabled[from]) continue;
            lift_attempt++;
            if (lift(from, n, pl)) {
                lift_count++;
                todo_push(from);
            }
        }
    }

    /**
     * Report final state.
     */
    if (trace) {
        for (int v=0; v<n_nodes; v++) {
            if (disabled[v]) continue;
            int *pm = pm_nodes + v*k;

            logger << "\033[1m" << label_vertex(v) << (owner[v]?" (odd)":" (even)") << "\033[m:";
            pm_stream(logger, pm, k);

            if (pm[0] == -1 or (pm[0]&1) != pl) {
                if (owner[v] != pl) {
                    if (strategy[v] == -1) logger << " no strategy!";
                    else logger << " => " << label_vertex(strategy[v]);
                }
            }
            
            logger << std::endl;
        }
    }

    /**
     * Derive strategies, mark as solved (only for opponent).
     */
    for (int v=0; v<n_nodes; v++) {
        if (disabled[v]) continue;
        int *pm = pm_nodes + v*k;

        if (pm[0] == -1 or (pm[0]&1) != pl) {
            if (owner[v] != pl and strategy[v] == -1) LOGIC_ERROR;
            oink->solve(v, 1-pl, owner[v] != pl ? strategy[v] : -1);
        }
    }

    oink->flush();
}


void
QPTSolver::run()
{
    // determine number of even/odd vertices and highest even/odd priority
    unsigned long even_vertices = 0, odd_vertices = 0;
    max_even = 0; max_odd = 0;
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;

        if ((priority[i]&1) == 0) {
            if (priority[i] > max_even) max_even = priority[i];
            even_vertices++;
        } else {
            if (priority[i] > max_odd) max_odd = priority[i];
            odd_vertices++;
        }
    }

    // find k s.t. 2^{k-1} >= even_vertices+1
    goal0 = even_vertices;
    even_vertices++;
    k0=1;
    while (even_vertices != 0) {
        k0++;
        even_vertices >>= 1;
    }

    goal1 = odd_vertices;
    odd_vertices++;
    k1=1;
    while (odd_vertices != 0) {
        k1++;
        odd_vertices >>= 1;
    }

    logger << "for odd with even measures: k=" << k0 << std::endl;
    logger << "for even with odd measures: k=" << k1 << std::endl;

    // for allocation, just use biggest k
    int big_k = k0 > k1 ? k0 : k1;

    // now create the data structure, for each vertex a PM
    pm_nodes = new int[big_k * n_nodes];
    strategy = new int[n_nodes];

    // initialize todo/dirty queues
    todo.resize(n_nodes);
    dirty.resize(n_nodes);

    lift_count = lift_attempt = 0;

    liftloop(0);

    logger << "after running for even, " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
    logger << game->countUnsolved() << " unsolved vertices left." << std::endl;

    liftloop(1);

    logger << "after running for odd, " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
    logger << game->countUnsolved() << " unsolved vertices left." << std::endl;

    delete[] pm_nodes;
    delete[] strategy;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
