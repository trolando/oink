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
#include <iomanip>

#include "sspm.hpp"

#define ODDFIRST 1

namespace pg {

SSPMSolver::SSPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

SSPMSolver::~SSPMSolver()
{
}

void
SSPMSolver::to_tmp(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) tmp_b[i] = pm_b[base+i];
    for (int i=0; i<l; i++) tmp_d[i] = pm_d[base+i];
}

void
SSPMSolver::from_tmp(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_b[base+i] = tmp_b[i];
    for (int i=0; i<l; i++) pm_d[base+i] = tmp_d[i];
}

void
SSPMSolver::to_best(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) best_b[i] = pm_b[base+i];
    for (int i=0; i<l; i++) best_d[i] = pm_d[base+i];
}

void
SSPMSolver::from_best(int idx)
{
    int base = idx*l;
    for (int i=0; i<l; i++) pm_b[base+i] = best_b[i];
    for (int i=0; i<l; i++) pm_d[base+i] = best_d[i];
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

/*
void
SSPMSolver::prog_cap_tmp(int pindex)
{
    // Start at the bottom, keep checking caps while going up
    if (tmp_d[0] == -1) return; // nothing to do

    // compute value
    int i=0;
    for (int d=0; d<=pindex and i<l; d++) {
        int val = 0, _i = i;

        for (; i<l; i++) {
            if (tmp_d[i] != d) {
                // e found
                val |= ((1 << (l-i)) - 1);
                break;
            }

            if (tmp_b[i]) val |= (1 << (l-i));
        }

        if (val > cap[d]) {
            if (d == 0) {
                // go to Top
                tmp_d[0] = -1;
                tmp_b.reset();
            } else {
                for (int j=_i; j<l; j++) {
                    tmp_d[j] = d-1;
                    tmp_b[j] = 0;
                }
                tmp_b[_i] = 1;
                prog_cap_tmp(pindex); // again
            }
            return;
        }
    }
}
*/

/**
 * Set tmp := min { m | m >_p tmp }
 */
void
SSPMSolver::prog_tmp(int pindex, int h)
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
    // * If lowest contains 0: 3a or 3b
    // * Else if lowest level is root: 3c
    // * Else append 100000000... to next higher level (3d, 3e, 3f)
    //
    // 3a: ,..011*  => ,..100*  (if lowest nonempty is the bottom)
    // 3b: ,..011*, => ,..,000* (if lowest nonempty is not the bottom)
    // 3c: 1111111  => Top      (if root contains only 1s)
    // 3d: ,1111111 => 100*     (if non-root contains only 1s)
    // 3e: ..,111*  => ..100*
    // 3f: ,e,111*  => ,100*

    // no bits below pindex, so just find smallest increase
    for (int i=l-1; i>=0; i--) {
        if (tmp_b[i] == 0) {
            if (tmp_d[i] == h) {
                // 3a: we found a 0 on the bottom, increase to 100...
                // ..011 => 00100
                tmp_b[i] = 1;
                break;
            } else {
                // 3b: we found a 0 (not bottom), increase to [eps] and 0s on leaf
                // ..011 => ..,000
                tmp_b[i] = 0;
                int new_d = tmp_d[i]+1;
                for (int k=i; k<l; k++) tmp_d[k] = new_d;
                break;
            }
        } else if (i == 0) {
            // we have only seen 1s and only 1 element
            if (tmp_d[0] == 0) {
                // 3c: we are already the highest, so go to top
                tmp_b[0] = 0;
                tmp_d[0] = -1;
                break;
            } else {
                // 3d: increase 1 higher...
                // ,e,111   => ,100
                int new_d = tmp_d[0]-1;
                tmp_b[0] = 1;
                for (int k=0; k<l; k++) tmp_d[k] = new_d;
                break;
            }
        } else if (tmp_d[i-1] != tmp_d[i]) {
            // 3e, 3f: next is different
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
SSPMSolver::stream_pm(std::ostream &out, int idx)
{
    int base = idx*l;
    if (pm_d[base] == -1) {
        out << " \033[1;33mTop\033[m";
    } else {
        out << " { ";
        int j=0;
        for (int i=0; i<h; i++) {
            if (i>0) out << ",";
            int c=0;
            while (j<l and pm_d[base+j] == i) {
                c++;
                out << pm_b[base+j];
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

        out << " {";

        // compute value
        int i=0;
        for (int d=0; d<h; d++) {
            int val = 0;

            for (; i<l; i++) {
                if (tmp_d[i] != d) {
                    // e found
                    val |= ((1 << (l-i)) - 1);
                    break;
                }

                if (tmp_b[i]) val |= (1 << (l-i));
            }

            logger << " " << val;
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
    if (tmp_d[0] == -1 and best_d[0] == -1) return 0;
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
    if (tmp_d[0] == -1 and test_d[0] == -1) return 0;
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
SSPMSolver::lift(int v, int target, int &str, int pl)
{
    // check if already Top
    if (pm_d[l*v] == -1) return false; // already Top

    const int pr = priority(v);
    const int pindex = pl == 0 ? h-(pr+1)/2-1 : h-pr/2-1;

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[37;1mupdating vertex " << label_vertex(v) << " (" << pr << " => " << pindex << ")" << (owner(v)?" (odd)":" (even)") << "\033[m with current measure";
        stream_pm(logger, v);
        logger << std::endl;
    }
#endif

    // if even owns and target is set, just check if specific target is better
    if (owner(v) == pl and target != -1) {
        to_tmp(target);
        if (pl == (pr&1)) prog_tmp(pindex, h);
        else trunc_tmp(pindex);
        to_best(v);
        if (compare(pindex) > 0) {
            from_tmp(v);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "to successor " << label_vertex(target) << ":";
                stream_tmp(logger, h);
                logger << " =>";
            }
#endif
            if (pl == (pr&1)) prog_tmp(pindex, h);
            else trunc_tmp(pindex);
#ifndef NDEBUG
            if (trace >= 2) {
                stream_tmp(logger, h);
                logger << std::endl;
            }
#endif
#ifndef NDEBUG
            if (trace >= 1) {
                logger << "\033[32;1mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
                stream_tmp(logger, h);
                logger << " (to " << label_vertex(target) << ")\n";
            }
#endif
            return true;
        } else {
            return false;
        }
    }

    // compute best measure
    bool first = true;
    for (auto curedge = outs(v); *curedge != -1; curedge++) {
        int to = *curedge;
        if (disabled[to]) continue;
        to_tmp(to);
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "to successor " << label_vertex(to) << " from";
            stream_tmp(logger, h);
            logger << " =>";
        }
#endif
        if (pl == (pr&1)) prog_tmp(pindex, h);
        else trunc_tmp(pindex);
#ifndef NDEBUG
        if (trace >= 2) {
            stream_tmp(logger, h);
            logger << std::endl;
        }
#endif
        if (first) {
            tmp_to_best();
            str = to;
        } else if (owner(v) == pl) {
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
    to_tmp(v);
    if (compare(pindex) < 0) {
#ifndef NDEBUG
        if (trace >= 1) {
            logger << "\033[1;32mnew measure\033[m of \033[36;1m" << label_vertex(v) << "\033[m:";
            stream_best(logger, h);
            logger << " (to " << label_vertex(str) << ")\n";
        }
#endif
        from_best(v);
        return true;
    } else {
        return false;
    }
}

static int
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
SSPMSolver::run(int n_bits, int depth, int player)
{
    l = n_bits;
    h = depth;

    pm_b.resize(l*nodecount());
    pm_d = new int[l*nodecount()];

    tmp_b.resize(l);
    tmp_d = new int[l];

    best_b.resize(l);
    best_d = new int[l];

    test_b.resize(l);
    test_d = new int[l];

    // initialize progress measures
    // pm_b.reset(); // standard set to 0 already
    memset(pm_d, 0, sizeof(int[l*nodecount()])); // every bit in the top ( = min )

    // lift_counters = new uint64_t[nodecount()];
    // memset(lift_counters, 0, sizeof(uint64_t[nodecount()]));

    for (int n=nodecount()-1; n>=0; n--) {
        if (disabled[n]) continue;
        lift_attempt++;
        int s;
        if (lift(n, -1, s, player)) {
            lift_count++;
            // lift_counters[n]++;
            for (auto curedge = ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                if (disabled[from]) continue;
                lift_attempt++;
                int s;
                if (lift(from, n, s, player)) {
                    lift_count++;
                    // lift_counters[from]++;
                    todo_push(from);
                }
            }
        }
    }

    while (!Q.empty()) {
        int n = todo_pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (disabled[from]) continue;
            lift_attempt++;
            int s;
            if (lift(from, n, s, player)) {
                lift_count++;
                // lift_counters[from]++;
                todo_push(from);
            }
        }
    }

    /**
     * Derive strategies.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_d[l*v] != -1) {
            if (owner(v) != player) {
                // TODO: don't rely on the strategy array in the Game class
                if (lift(v, -1, game.getStrategy()[v], player)) logger << "error: " << v << " is not progressive!" << std::endl;
            }
        }
    }

    if (trace) {
        for (int v=0; v<nodecount(); v++) {
            if (disabled[v]) continue;

            logger << "\033[1m" << label_vertex(v) << (owner(v)?" (odd)":" (even)") << "\033[m:";
            stream_pm(logger, v);

            if (pm_d[l*v] != -1) {
                if (owner(v) != player) {
                    logger << " => " << label_vertex(game.getStrategy(v));
                }
            }

            logger << std::endl;
        }
    }

    /**
     * Mark solved.
     */

    for (int v=0; v<nodecount(); v++) {
        if (disabled[v]) continue;
        if (pm_d[l*v] != -1) Solver::solve(v, 1-player, game.getStrategy(v));
    }

    Solver::flush();

    /*
    {
        uint64_t hi_c = 0;
        int hi_n = -1;
        for (int i=0; i<nodecount(); i++) {
            if (hi_c < lift_counters[i]) {
                hi_c = lift_counters[i];
                hi_n = i;
            }
        }
        logger << "highest counter: " << hi_c << " lifts for vertex " << label_vertex(hi_n) << std::endl;
    }
    */

    delete[] pm_d;
    delete[] tmp_d;
    delete[] best_d;
    delete[] test_d;
}

void
SSPMSolver::run()
{
    int max_prio = priority(nodecount()-1);

    // compute ml (max l) and the h for even/odd
    int ml = ceil_log2(nodecount());
    int h0 = (max_prio/2)+1;
    int h1 = (max_prio+1)/2;

    // create datastructures
    Q.resize(nodecount());
    dirty.resize(nodecount());
    unstable.resize(nodecount());

    logger << "even wants " << ml << "-bounded adaptive " << h0 << "-counters." << std::endl;
    logger << "odd wants " << ml << "-bounded adaptive " << h1 << "-counters." << std::endl;

    // if running bounded sspm, start with 1-bounded adaptive counters
    int i = bounded ? 1 : ml;

    for (; i<=ml; i++) {
        int _l = lift_count, _a = lift_attempt;
        uint64_t _c = game.count_unsolved(), c;

        if (ODDFIRST) {
            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        } else {
            // run even counters
            run(i, h0, 0);
            c = game.count_unsolved();
            logger << "after even with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;

            // if now solved, no need to run odd counters
            if (c == 0) break;

            // run odd counters
            run(i, h1, 1);
            c = game.count_unsolved();
            logger << "after odd  with k=" << i << ", " << std::setw(9) << lift_count-_l << " lifts, " << std::setw(9) << lift_attempt-_a << " lift attempts, " << c << " unsolved left." << std::endl;
        }

        if (i != ml) {
            // if i == ml then we are guaranteed to be done
            // otherwise check if done
            if (c == 0) break;
            if (_c != c) i--; // do not increase i if we solved vertices with current i
        } else {
            break; // do not count higher pls
        }
    }

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts, max l " << i << "." << std::endl;
}

}
