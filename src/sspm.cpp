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

#include <cstring> // for memset

#include "sspm.hpp"

namespace pg {

SSPMSolver::SSPMSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

SSPMSolver::~SSPMSolver()
{
}

void
SSPMSolver::to_tmp_0(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) tmp_b[i] = pm_0_b[base+i];
    for (int i=0; i<l; i++) tmp_d[i] = pm_0_d[base+i];
}

void
SSPMSolver::from_tmp_0(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_0_b[base+i] = tmp_b[i];
    for (int i=0; i<l; i++) pm_0_d[base+i] = tmp_d[i];
}

void
SSPMSolver::to_best_0(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) best_b[i] = pm_0_b[base+i];
    for (int i=0; i<l; i++) best_d[i] = pm_0_d[base+i];
}

void
SSPMSolver::from_best_0(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_0_b[base+i] = best_b[i];
    for (int i=0; i<l; i++) pm_0_d[base+i] = best_d[i];
}

void
SSPMSolver::to_tmp_1(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) tmp_b[i] = pm_1_b[base+i];
    for (int i=0; i<l; i++) tmp_d[i] = pm_1_d[base+i];
}

void
SSPMSolver::from_tmp_1(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_1_b[base+i] = tmp_b[i];
    for (int i=0; i<l; i++) pm_1_d[base+i] = tmp_d[i];
}

void
SSPMSolver::to_best_1(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) best_b[i] = pm_1_b[base+i];
    for (int i=0; i<l; i++) best_d[i] = pm_1_d[base+i];
}

void
SSPMSolver::from_best_1(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_1_b[base+i] = best_b[i];
    for (int i=0; i<l; i++) pm_1_d[base+i] = best_d[i];
}

void
SSPMSolver::tmp_to_best()
{
    best_b = tmp_b;
    memcpy(best_d, tmp_d, sizeof(int[l]));
}

void
SSPMSolver::tmp_to_test()
{
    test_b = tmp_b;
    memcpy(test_d, tmp_d, sizeof(int[l]));
}

/**
 * Set tmp := min { m | m ==_p tmp }
 */
void
SSPMSolver::trunc_tmp(int pindex)
{
    if (tmp_d[0] == -1) return; // already Top
    // compute the lowest pindex >= p
    // [pindex],.,...,.. => [pindex],000 
    // if pindex is the bottom, then this simply "buries" the remainder
    for (int i=l-1; i>=0 and tmp_d[i] > pindex; i--) {
        tmp_b[i] = 0;
        tmp_d[i] = pindex+1;
    }
}

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
SSPMSolver::prog_tmp(int pindex)
{
    // Simple case 1: Top >_p Top
    if (tmp_d[0] == -1) return; // already Top

    // Simple case 2: Some bits below [pindex], ergo [pindex] can go from ..e to ..10*
    if (tmp_d[l-1] > pindex) {
        int i;
        for (i=l-1; i>=0 and tmp_d[i] > pindex; i--) {
            tmp_b[i] = 0;
            tmp_d[i] = pindex;
        }
        tmp_b[i+1] = 1;
        return;
    }

    // Case 3: no bits below [pindex], so analyze lowest nonempty level
    // * If lowest contains 0: 3a
    // * Else if lowest level is root: 3b
    // * Else append 100000000... to next higher level (3c, 3d, 3e)
    // 
    // 3a: ,..011*, => ,..,000*
    // 3b: 1111111  => Top      (if root contains only 1s)
    // 3c: ,1111111 => 100*     (if non-root contains only 1s)
    // 3d: ..,111*  => ..100*
    // 3e: ,e,111*  => ,100*

    // no bits below pindex, so just find smallest increase
    for (int i=l-1; i>=0; i--) {
        if (tmp_b[i] == 0) {
            // 3a: we found a 0, increase to [eps] and 0s on leaf
            // ..011 => ..,000
            tmp_b[i] = 0;
            int new_d = tmp_d[i]+1;
            for (int k=i; k<l; k++) tmp_d[k] = new_d;
            break;
        } else if (i == 0) {
            // we have only seen 1s and only 1 element
            if (tmp_d[0] == 0) {
                // 3b: we are already the highest, so go to top
                tmp_b[0] = 0;
                tmp_d[0] = -1;
                break;
            } else {
                // 3c: increase 1 higher...
                // ,e,111   => ,100
                int new_d = tmp_d[0]-1;
                tmp_b[0] = 1;
                for (int k=0; k<l; k++) tmp_d[k] = new_d;
                break;
            }
        } else if (tmp_d[i-1] != tmp_d[i]) {
            // 3d, 3e: next is different
            // two cases
            //  .,111   => .100
            // ,e,111   => ,100
            int new_d = tmp_d[i]-1;
            tmp_b[i] = 1;
            for (int k=i; k<l; k++) tmp_d[k] = new_d;
            break;
        } else {
            // next is same
            tmp_b[i] = 0;
        }
    }
}

/**
 * Write pm to ostream.
 */
void
SSPMSolver::stream_pm_0(std::ostream &out, int idx)
{
    int base = idx*l;
    if (pm_0_d[base] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h0; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and pm_0_d[base+j] == i) {
                c++;
                out << pm_0_b[base+j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Write pm1 to ostream.
 */
void
SSPMSolver::stream_pm_1(std::ostream &out, int idx)
{
    int base = idx*l;
    if (pm_1_d[base] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h1; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and pm_1_d[base+j] == i) {
                c++;
                out << pm_1_b[base+j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Write tmp to ostream.
 */
void
SSPMSolver::stream_tmp(std::ostream &out, int h)
{
    if (tmp_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and tmp_d[j] == i) {
                c++;
                out << tmp_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Write best to ostream.
 */
void
SSPMSolver::stream_best(std::ostream &out, int h)
{
    if (best_d[0] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and best_d[j] == i) {
                c++;
                out << best_b[j];
                j++;
            }
            if (c == 0) out << "ε";
        }
        out << " }";
    }
}

/**
 * Compare tmp and best
 * res := -1 :: tmp < best
 * res := 0  :: tmp = best
 * res := 1  :: tmp > best
 */
int
SSPMSolver::compare(int pindex)
{
    // cases involving Top
    if (tmp_d[0] == -1 and best_d[1] == -1) return 0;
    if (tmp_d[0] == -1) return 1;
    if (best_d[0] == -1) return -1;
    for (int i=0; i<l; i++) {
        if (tmp_d[i] > pindex and best_d[i] > pindex) {
            // equal until pindex, return 0
            return 0;
        } else if (tmp_d[i] < best_d[i]) {
            // equal until best has [eps]
            if (tmp_b[i] == 0) return -1;
            else return 1;
        } else if (tmp_d[i] > best_d[i]) {
            // equal until tmp has [eps]
            if (best_b[i] == 0) return 1;
            else return -1;
        } else if (tmp_b[i] < best_b[i]) {
            // equal until tmp<best
            return -1;
        } else if (tmp_b[i] > best_b[i]) {
            // equal until tmp>best
            return 1;
        }
    }
    return 0;
}

/**
 * Compare tmp and test
 * res := -1 :: tmp < test
 * res := 0  :: tmp = test
 * res := 1  :: tmp > test
 */
int
SSPMSolver::compare_test(int pindex)
{
    // cases involving Top
    if (tmp_d[0] == -1 and test_d[1] == -1) return 0;
    if (tmp_d[0] == -1) return 1;
    if (test_d[0] == -1) return -1;
    for (int i=0; i<l; i++) {
        if (tmp_d[i] > pindex and test_d[i] > pindex) {
            // equal until pindex, return 0
            return 0;
        } else if (tmp_d[i] < test_d[i]) {
            // equal until test has [eps]
            if (tmp_b[i] == 0) return -1;
            else return 1;
        } else if (tmp_d[i] > test_d[i]) {
            // equal until tmp has [eps]
            if (test_b[i] == 0) return 1;
            else return -1;
        } else if (tmp_b[i] < test_b[i]) {
            // equal until tmp<test
            return -1;
        } else if (tmp_b[i] > test_b[i]) {
            // equal until tmp>test
            return 1;
        }
    }
    return 0;
}

bool
SSPMSolver::lift_0(int v, int target, int &str)
{
    // check if already Top
    if (pm_0_d[l*v] == -1) return false; // already Top

    const int pr = priority[v];
    const int pindex = h0-(pr+1)/2-1;

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1mupdating vertex " << v << "/" << priority[v] << (owner[v]?" (odd)":" (even)") << "\033[m with current measure";
        stream_pm_0(logger, v);
        logger << std::endl;
    }
#endif

    // if even owns and target is set, just check if specific target is better
    if (owner[v] == 0 and target != -1) {
        to_tmp_0(target);
        if (pr&1) trunc_tmp(pindex);
        else prog_tmp(pindex);
        to_best_0(v);
        if (compare(pindex) > 0) {
            from_tmp_0(v);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "new measure";
                stream_tmp(logger, h0);
                logger << std::endl;
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // compute best measure
    bool first = true;
    for (int to : out[v]) {
        if (disabled[to]) continue;
        to_tmp_0(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "successor " << to << " from";
            stream_tmp(logger, h0);
        }
        // DEBUG
        tmp_to_test();
#endif
        if (pr&1) trunc_tmp(pindex);
        else prog_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << " to";
            stream_tmp(logger, h0);
            logger << std::endl;
        }
        // DEBUG test is progressive!
        if (test_d[0] != -1) {
            if (pr&1) {
                if (compare_test(pindex) != 0) LOGIC_ERROR;
            } else {
                if (compare_test(pindex) != 1) LOGIC_ERROR;
            }
        } else {
            if (tmp_d[0] != -1) LOGIC_ERROR;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
        } else if (owner[v] == 0) {
            // we want the max!
            if (compare(pindex) > 0) {
                tmp_to_best();
                str = to;
            }
        } else {
            // we want the min!
            if (compare(pindex) < 0) {
                tmp_to_best();
                str = to;
            }
        }
        first = false;
    }

    // set best to pm if higher
    to_tmp_0(v);
    if (compare(pindex) < 0) {
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "new measure";
            stream_best(logger, h0);
            logger << std::endl;
        }
#endif
        from_best_0(v);
        return true;
    } else {
        return false;
    }
}

bool
SSPMSolver::lift_1(int v, int target, int &str)
{
    // check if already Top
    if (pm_1_d[l*v] == -1) return false; // already Top

    const int pr = priority[v];
    const int pindex = h1-pr/2-1;

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1mupdating vertex " << v << "/" << priority[v] << (owner[v]?" (odd)":" (even)") << "\033[m with current measure";
        stream_pm_1(logger, v);
        logger << std::endl;
    }
#endif

    // if odd owns and target is set, just check if specific target is better
    if (owner[v] == 1 and target != -1) {
        to_tmp_1(target);
        if (pr&1) prog_tmp(pindex);
        else trunc_tmp(pindex);
        to_best_1(v);
        if (compare(pindex) > 0) {
            from_tmp_1(v);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "new measure";
                stream_tmp(logger, h1);
                logger << std::endl;
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // compute best measure
    bool first = true;
    for (int to : out[v]) {
        if (disabled[to]) continue;
        to_tmp_1(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "successor " << to << " from";
            stream_tmp(logger, h1);
        }
        // DEBUG
        tmp_to_test();
#endif
        if (pr&1) prog_tmp(pindex);
        else trunc_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << " to";
            stream_tmp(logger, h1);
            logger << std::endl;
        }
        // DEBUG test is progressive!
        if (test_d[0] != -1) {
            if (pr&1) {
                if (compare_test(pindex) != 1) LOGIC_ERROR;
            } else {
                if (compare_test(pindex) != 0) LOGIC_ERROR;
            }
        } else {
            if (tmp_d[0] != -1) LOGIC_ERROR;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
        } else if (owner[v] == 1) {
            // we want the max!
            if (compare(pindex) > 0) {
                tmp_to_best();
                str = to;
            }
        } else {
            // we want the min!
            if (compare(pindex) < 0) {
                tmp_to_best();
                str = to;
            }
        }
        first = false;
    }

    // set best to pm if higher
    to_tmp_1(v);
    if (compare(pindex) < 0) {
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "new measure";
            stream_best(logger, h1);
            logger << std::endl;
        }
#endif
        from_best_1(v);
        return true;
    } else {
        return false;
    }
}

int
ceil_log2(unsigned long long x)
{
    static const unsigned long long t[6] = {
        0xFFFFFFFF00000000ull,
        0x00000000FFFF0000ull,
        0x000000000000FF00ull,
        0x00000000000000F0ull,
        0x000000000000000Cull,
        0x0000000000000002ull
    };

    int y = (((x & (x - 1)) == 0) ? 0 : 1);
    int j = 32;
    int i;

    for (i = 0; i < 6; i++) {
        int k = (((x & t[i]) == 0) ? 0 : j);
        y += k;
        x >>= k;
        j >>= 1;
    }

    return y;
}

void
SSPMSolver::run()
{
    // determine l and h
    l = ceil_log2(n_nodes);

    int max_prio = priority[n_nodes-1];
    h0 = (max_prio/2)+1; // h for even measures
    h1 = (max_prio+1)/2; // h for odd measures

    logger << "even: " << l << "-bounded adaptive " << h0 << "-counters." << std::endl;
    logger << "odd: "  << l << "-bounded adaptive " << h1 << "-counters." << std::endl;

    // create datastructures
    Q.resize(n_nodes);
    dirty.resize(n_nodes);
    unstable.resize(n_nodes);

    pm_0_b.resize(l*n_nodes);
    pm_0_d = new int[l*n_nodes];

    pm_1_b.resize(l*n_nodes);
    pm_1_d = new int[l*n_nodes];

    tmp_b.resize(l);
    tmp_d = new int[l];

    best_b.resize(l);
    best_d = new int[l];

    test_b.resize(l);
    test_d = new int[l];

    // initialize progress measures
    // pm_0_b.reset(); // standard set to 0 already
    // pm_1_b.reset(); // standard set to 0 already
    memset(pm_0_d, 0, sizeof(int[l*n_nodes])); // every bit in the top ( = min )
    memset(pm_1_d, 0, sizeof(int[l*n_nodes])); // every bit in the top ( = min )

    // set number of lifts and lift attempts to 0
    int lift_count = 0;
    int lift_attempt = 0;

    /**
     * Compute even measures.
     */

    int s;

    for (int n=n_nodes-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        if (lift_0(n, -1, s)) {
            lift_count++;
            for (int from : in[n]) {
                if (disabled[from]) continue;
                lift_attempt++;
                if (lift_0(from, n, s)) {
                    lift_count++;
                    todo_push(from);
                }
            }
        }
    }
    
    while (!Q.empty()) {
        int n = todo_pop();
        for (int from : in[n]) {
            if (disabled[from]) continue;
            lift_attempt++;
            if (lift_0(from, n, s)) {
                lift_count++;
                todo_push(from);
            }
        }
    }

    /**
     * Initialize odd measures if they are won by Odd
     */

    for (int n=n_nodes-1; n>=0; n--) {
        if (disabled[n]) continue;
        if (pm_0_d[l*n] != -1) pm_1_d[l*n] = -1;
    }

    /**
     * Compute odd measures.
     */

    for (int n=n_nodes-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        if (lift_1(n, -1, s)) {
            lift_count++;
            for (int from : in[n]) {
                if (disabled[from]) continue;
                lift_attempt++;
                if (lift_1(from, n, s)) {
                    lift_count++;
                    todo_push(from);
                }
            }
        }
    }
    
    while (!Q.empty()) {
        int n = todo_pop();
        for (int from : in[n]) {
            if (disabled[from]) continue;
            lift_attempt++;
            if (lift_1(from, n, s)) {
                lift_count++;
                todo_push(from);
            }
        }
    }

    /**
     * Derive strategies (and run tests if debugging).
     */

    for (int v=0; v<n_nodes; v++) {
        if (disabled[v]) continue;
        // compute strategy if won and owned by Even
        if (pm_0_d[l*v] == -1 and owner[v] == 0) {
            if (lift_1(v, -1, game->strategy[v])) logger << "error: " << v << " is not progressive!" << std::endl;
        }
        // compute strategy if won and owned by Odd
        if (pm_0_d[l*v] != -1 and owner[v] == 1) {
            if (lift_0(v, -1, game->strategy[v])) logger << "error: " << v << " is not progressive!" << std::endl;
        }
#ifndef NDEBUG
        // test if exactly one of the measures is Top
        if ((pm_0_d[l*v] == -1) == (pm_1_d[l*v] == -1)) LOGIC_ERROR;
#endif
    }

    // write final measures to logger
    if (trace) {
        for (int v=0; v<n_nodes; v++) {
            logger << "\033[1m" << v << "/" << priority[v] << (owner[v]?" (odd)":" (even)") << "\033[m:";
            stream_pm_0(logger, v);
            stream_pm_1(logger, v);
            if ((pm_0_d[l*v] == -1) == (owner[v] == 0)) logger << " " << v << "=>" << game->strategy[v];
            logger << std::endl;
        }
    }

    /**
     * Finally, solve the game.
     */

    for (int v=0; v<n_nodes; v++) {
        if (disabled[v]) continue;
        oink->solve(v, pm_0_d[l*v] == -1 ? 0 : 1, game->strategy[v]);
    }

    delete[] pm_0_d;
    delete[] pm_1_d;
    delete[] tmp_d;
    delete[] best_d;
    delete[] test_d;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
