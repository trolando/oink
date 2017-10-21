#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stack>
#include "spm.hpp"

/**
 * For easy debugging, remove the comment from one of below defines to
 * restrict the solver to even or odd measures only.
 */

// #define SINGLE(d) d == 0 // only compute even measures
// #define SINGLE(d) d == 1 // only compute odd measures

#ifndef SINGLE
#define SINGLE(d) true
#endif

namespace pg {

SPMSolver::SPMSolver(Oink *oink, Game *game, std::ostream &lgr) : Solver(oink, game, lgr)
{
}

SPMSolver::~SPMSolver()
{
}

/**
 * Returns true if a progress measure "a" is less than "b"
 * up to and including priority <d>, for player <pl>.
 */
bool
SPMSolver::pm_less(int* a, int* b, int d, int pl)
{
    // cases where a or b is Top
    if (a[pl] == -1) return false;
    if (b[pl] == -1) return true;
    // normal comparison, start with highest priority
    const int start = ((k&1) == pl) ? k-2 : k-1;
    for (int i=start; i>=d; i-=2) {
        if (a[i] == b[i]) continue;
        if (a[i] > counts[i] and b[i] > counts[i]) return false;
        return a[i] < b[i];
    }
    return false;
}

/**
 * Copy for player <pl>.
 */
void
SPMSolver::pm_copy(int *dst, int *src, int pl)
{
    for (int i=pl; i<k; i+=2) dst[i] = src[i];
}

/**
 * Obtain the highest priority for player <pl> for which this measure is cyclic, or -1 if none.
 */
int
SPMSolver::pm_cycles(int *a, int pl)
{
    int m = k-1;
    if ((k&1)==pl) m--;
    for (int i=m; i>=0; i-=2) {
        if (a[i] > counts[i]) return i;
    }
    return -1;
}

/**
 * Write pm to ostream.
 */
void
SPMSolver::pm_stream(std::ostream &out, int *pm)
{
    bool top_e = pm[0] == -1;
    bool top_o = pm[1] == -1;
    out << " {";
    if (top_e) out << " \033[1;33mTe\033[m";
    else out << " " << pm[0];
    if (top_o) out << " \033[1;33mTo\033[m";
    else out << " " << pm[1];
    for (int i=2; i<k; i++) {
        if (i&1) out << " " << (top_o ? 0 : pm[i]);
        else     out << " " << (top_e ? 0 : pm[i]);
    }
    out << " } ";
}

/**
 * Perform update for player <pl>, node with priority <d>. Does not carry-over.
 */
void
SPMSolver::Prog(int *dst, int *src, int d, int pl)
{
    // check if top
    if (src[pl] == -1) {
        dst[pl] = -1;
        return;
    }

    // set every value lower than <d> to 0.
    int i = pl;
    for (; i<d; i+=2) dst[i] = 0;
    if (d == i) {
        dst[i] = src[i] + 1;
        i+=2;
    }
    for (; i<k; i+=2) dst[i] = src[i];
}

bool
SPMSolver::canlift(int node, int pl)
{
    // obtain ptr to current progress measure
    int *pm = pms + k*node;

    // check if already Top
    if (pm[pl] == -1) return false;

    const int d = priority[node];

    if (owner[node] == pl) {
        // do max
        for (int to : out[node]) {
            if (disabled[to]) continue;
            Prog(tmp, pms + k*to, d, pl);
            if (pm_less(pm, tmp, d, pl)) return true;
        }
        return false;
    } else {
        // do min
        int best_to = -1;
        for (int to : out[node]) {
            if (disabled[to]) continue;
            Prog(tmp, pms + k*to, d, pl);
            if (best_to == -1 or pm_less(tmp, best, d, pl)) {
                for (int i=0; i<k; i++) best[i] = tmp[i];
                best_to = to;
            }
        }
        if (best_to == -1) return false;
        return pm_less(pm, best, d, pl);
    }
}

bool
SPMSolver::lift(int node, int target)
{
    // obtain ptr to current progress measure
    int *pm = pms + k*node;

    // check if already Top for both players
    if (pm[0] == -1 and pm[1] == -1) return false;

    lift_attempt++;

    // initialize stuff
    const int pl_max = owner[node];
    const int pl_min = 1 - pl_max;
    const int d = priority[node];

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1mupdating node " << node << "/" << d << (owner[node]?" (odd)":" (even)") << "\033[m with current progress measure";
        pm_stream(logger, pm);
        logger << std::endl;
    }
#endif

    int best_ch0 = -1, best_ch1 = -1;

    // do max for player <pl_max>
    if (SINGLE(pl_max) and pm[pl_max] != -1) {
#ifndef NDEBUG
        if (trace >= 2) {
            logger << "computing max" << std::endl;
            pm_copy(tmp, pm, 1-pl_max);
        }
#endif
        if (target != -1) {
            // just look at target
            Prog(tmp, pms + k*target, d, pl_max);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "successor node " << target << "/" << priority[target] << " results in";
                pm_stream(logger, tmp);
                logger << std::endl;
            }
#endif
            if (pm_less(pm, tmp, d, pl_max)) {
                pm_copy(pm, tmp, pl_max);
                if (pl_max) best_ch1 = target;
                else best_ch0 = target;
            }
        } else for (int to : out[node]) {
            if (disabled[to]) continue;
            Prog(tmp, pms + k*to, d, pl_max);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "successor node " << to << "/" << priority[to] << " results in";
                pm_stream(logger, tmp);
                logger << std::endl;
            }
#endif
            if (pm_less(pm, tmp, d, pl_max)) {
                pm_copy(pm, tmp, pl_max);
                if (pl_max) best_ch1 = to;
                else best_ch0 = to;
            }
        }
    }

    // do min for player <pl_min>
    if (SINGLE(pl_min) and pm[pl_min] != -1 and (target == -1 or target == strategy[node])) {
#ifndef NDEBUG
        if (trace >= 2) logger << "computing min" << std::endl;
        if (trace >= 2) pm_copy(tmp, pm, 1-pl_min);
#endif
        int best_to = -1;
        for (int to : out[node]) {
            if (disabled[to]) continue;
            Prog(tmp, pms + k*to, d, pl_min);
#ifndef NDEBUG
            if (trace >= 2) {
                logger << "successor node " << to << "/" << priority[to] << " results in";
                pm_stream(logger, tmp);
                logger << std::endl;
            }
#endif
            if (best_to == -1 or pm_less(tmp, best, d, pl_min)) {
                for (int i=0; i<k; i++) best[i] = tmp[i];
                best_to = to;
            }
        }

        strategy[node] = best_to;
        // note: sometimes only the strategy changes, but the lowest pm stays the same
        // now "best" contains the smallest Prog, which may be higher than the current min
        if (pm_less(pm, best, d, pl_min)) {
            pm_copy(pm, best, pl_min);
            if (pl_min) best_ch1 = best_to;
            else best_ch0 = best_to;
        }
    }

    if (best_ch0 != -1 or best_ch1 != -1) {
        if (trace) {
            logger << "\033[1;32mupdated node " << node << "/" << d << (owner[node]?" (odd)":" (even)") << "\033[m to";
            pm_stream(logger, pm);
            logger << std::endl;
        }
        // increase count and return true
        lift_count++;
        return true;
    } else {
        return false;
    }
}

void
SPMSolver::update(int pl)
{
    std::queue<int> q;

    // find unstable nodes (for measure <pl>)
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        unstable[i] = 0; // first mark as stable
        if (pms[k*i + pl] == -1 or canlift(i, pl) or pm_cycles(pms+k*i, pl) != -1) {
            unstable[i] = 1;
            q.push(i);
        }
    }

    while (!q.empty()) {
        int n = q.front();
        q.pop();
        for (int m : in[n]) {
            if (disabled[m] or unstable[m]) continue;
            if (owner[m] != pl) {
                int best_to = -1;
                const int d = priority[m];
                for (int to : out[m]) {
                    if (disabled[to]) continue;
                    if (unstable[to]) continue;
                    Prog(tmp, pms + k*to, d, pl);
                    if (best_to == -1 or pm_less(tmp, best, d, pl)) {
                        for (int i=0; i<k; i++) best[i] = tmp[i];
                        best_to = to;
                    }
                }
                if (best_to != -1 and !pm_less(pms+k*m, best, d, pl)) continue;
            }
            unstable[m] = 1;
            q.push(m);
        }
    }

    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        if (unstable[i] == 0 and pms[k*i + 1-pl] != -1) {
            if ((priority[i]&1) != pl) counts[priority[i]]--;
            pms[k*i + 1-pl] = -1;
            todo_push(i);

            if (trace) {
                logger << "\033[1;33mupdated node " << i << "/" << priority[i] << (owner[i]?" (odd)":" (even)") << "\033[m to";
                pm_stream(logger, pms + i*k);
                logger << std::endl;
            }
        }
    }
}

void
SPMSolver::run()
{
    // determine k = highest priority + 1
    k = priority[n_nodes-1]+1;
    if (k < 2) k = 2;

    // now create the data structure, for each node
    pms = new int[(size_t)k*n_nodes];
    strategy = new int[n_nodes];
    counts = new int[k];
    tmp = new int[k];
    best = new int[k];
    dirty = new int[n_nodes];
    unstable = new int[n_nodes];

    int max0 = -1, max1 = -1;

    // initialize all measures to 0
    for (int i=0; i<k*n_nodes; i++) pms[i] = 0;

    // initialize strategy to -1
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    // initialize counts for each priority
    for (int i=0; i<k; i++) counts[i] = 0;
    for (int i=0; i<n_nodes; i++) if (disabled[i] == 0) counts[priority[i]]++;
    for (int i=k-1; i>=0; i--) {
        if (counts[i] == 0) continue;
        if (i&1) { if (max1 == -1) max1 = i; }
        else { if (max0 == -1) max0 = i; }
        if (max0 != -1 and max0 != -1) break;
    }

    // initialize all nodes as not dirty
    for (int n=0; n<n_nodes; n++) dirty[n] = 0;

    // allocate and initialize additional array for cycle measures
    int *cm = new int[n_nodes];
    for (int n=0; n<n_nodes; n++) cm[n] = 0;

    // a queue and a vector for cycle measure analysis
    std::queue<int> cm_queue;
    std::vector<int> cycles;

    // set number of lifts and lift attempts to 0
    lift_count = lift_attempt = 0;

    /**
     * Strategy that updates predecessors then marks updated predecessors for processing.
     * Uses a queue to store the dirty vertices.
     */

    /**
     * Initialization loop.
     */

    for (int n=n_nodes-1; n>=0; n--) {
        if (!disabled[n] and lift(n, -1)) {
            for (int from : in[n]) if (!disabled[from] and lift(from, n)) todo_push(from);
        }
    }
    
    /**
     * The main loop.
     */

    int64_t last_update = 0;

    while (true) {
        while (!todo.empty()) {
            int n = todo_pop();
            for (int from : in[n]) if (!disabled[from] and lift(from, n)) todo_push(from);
            if (last_update + 10*n_nodes < lift_count) {
                last_update = lift_count;
                update(0);
                update(1);
            }
        }

#ifndef NDEBUG
        if (trace >= 2) {
            for (int n=0; n<n_nodes; n++) {
                if (disabled[n]) continue;
                logger << "\033[35m**\033[m \033[1mnode " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m is";
                pm_stream(logger, pms + k*n);
                logger << std::endl;
            }
        }
#endif

        /**
         * End of the main loop.
         * When no more lifting occurs, we must analyze nodes with cyclic measures.
         * We analyze for player 0 first, then for player 1.
         */

        for (int pl=0; pl<2; pl++) {
            /**
             * Put all nodes with cycle measures in "cycles"
             * Also find the best (lowest) exit for the loser
             */
            int best_from = -1, best_to = -1;
            for (int i=0; i<k; i++) best[i] = 0;

            int max = pl == 0 ? max0 : max1;

            for (int n=0; n<n_nodes; n++) {
                if (disabled[n]) continue;
                int *pm = pms + k*n;
                if (pm[pl] == -1) continue; // already won
                int c = pm_cycles(pm, pl);
                if (c == -1) continue; // not a cycle measure
#ifndef NDEBUG
                if (trace >= 2) {
                    logger << "\033[31;1m" << (pl ? "odd" : "even") << (pl == owner[n] ? " winner " : " loser ");
                    logger << "loops\033[m: " << n << "/" << priority[n] << " (cm = " << pm_cycles(pm, pl) << ")";
                    pm_stream(logger, pm);
                    logger << std::endl;
                }
#endif
                if (c == max) {
                    pms[n*k + pl] = -1;
                    todo_push(n);

                    const int d = priority[n];
                    if ((d&1) == pl) counts[d]--;

                    if (trace) {
                        logger << "\033[1;36mupdated node " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m to";
                        pm_stream(logger, pms + n*k);
                        logger << std::endl;
                    }

                    continue;
                }

                cycles.push_back(n);
                if (owner[n] != pl) {
                    const int d = priority[n];
                    for (int to : out[n]) {
                        if (disabled[to]) continue;
                        int *pm_to = pms + k*to;
                        if (pm_to[pl] == -1) continue; // already won
                        Prog(tmp, pm_to, d, pl);
                        if (!pm_less(pms + k*n, tmp, d, pl)) continue; // candidate must increase measure
                        if (best_to == -1 or pm_less(tmp, best, 0, pl))  {
                            for (int i=pl; i<k; i+=2) best[i] = tmp[i];
                            best_from = n;
                            best_to = to;
                        }
                    }
                }
            }

            if (best_to != -1) {
                /**
                 * OK we found the best exit from <best_from> to <best_to>.
                 */
                strategy[best_from] = best_to;
                pm_copy(pms + best_from*k, best, pl);
                todo_push(best_from);

                if (trace) {
                    logger << "\033[1;36mupdated node " << best_from << "/" << priority[best_from] << (owner[best_from]?" (odd)":" (even)") << "\033[m to";
                    pm_stream(logger, pms + best_from*k);
                    logger << std::endl;
                }

                cycles.clear();
            } else {
                /**
                 * No exit found, everything that remains a cycle after removing escaping nodes is now won.
                 */
                for (int n : cycles) cm[n] = 1; // mark everything in 'cycles'

                for (int n : cycles) {
                    bool escapes;
                    if (owner[n] == pl) {
                        // check if it can stay in cm
                        escapes = true;
                        for (int m : out[n]) {
                            if (disabled[m] == 0 and cm[m]) {
                                escapes = false;
                                break;
                            }
                        }
                    } else {
                        // check if it can move out of cm
                        escapes = false;
                        for (int m : out[n]) {
                            if (disabled[m] == 0 and pms[k*m+pl] != -1 and cm[m] == 0) {
                                escapes = true;
                                break;
                            }
                        }
                    }
                    if (escapes) {
                        cm[n] = 0;
                        cm_queue.push(n);
                    }
                }

                /**
                 * Run backwards search to remove nodes that can escape to non-CM
                 */

                while (!cm_queue.empty()) {
                    int n = cm_queue.front();
                    cm_queue.pop();
                    for (int m : in[n]) {
                        if (disabled[m] != 0 or cm[m] == 0) continue;
                        if (owner[m] == pl) {
                            bool escapes = false;
                            for (int to : out[m]) {
                                if (disabled[to] == 0 and cm[to]) {
                                    escapes = true;
                                    break;
                                }
                            }
                            if (escapes) continue;
                        }
                        cm[m] = 0;
                        cm_queue.push(m);
                    }
                }

                /**
                 * Finally, mark remaining nodes as won
                 */

                for (int n : cycles) {
                    if (cm[n]) {
                        cm[n] = 0;
                        pms[n*k + pl] = -1;
                        todo_push(n);

                        const int d = priority[n];
                        if ((d&1) == pl) counts[d]--;

                        if (trace) {
                            logger << "\033[1;36mupdated node " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m to";
                            pm_stream(logger, pms + n*k);
                            logger << std::endl;
                        }
                    }
                }

                cycles.clear();
            }

            if (counts[max] == 0) {
                if (trace >= 2) logger << "\033[1mlowering max\033[m for player " << (pl == 0 ? "even" : "odd") << std::endl;
                if (pl == 0) { while (max0 > 0 and counts[max0] == 0) max0 -= 2; }
                else { while (max1 > 0 and counts[max1] == 0) max1 -= 2; }
            }
        }

        if (todo.empty()) break;
    }

#ifndef NDEBUG
    if (trace >= 2) {
        for (int n=0; n<n_nodes; n++) {
            if (disabled[n]) continue;
            logger << "\033[35m**\033[m \033[1mnode " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m is";
            pm_stream(logger, pms + k*n);
            logger << std::endl;
        }
    }
#endif

    // Now set dominions and derive strategy for even.
    for (int n=0; n<n_nodes; n++) {
        if (disabled[n]) continue;
        int *pm = pms + k*n;
        if ((pm[0] == -1) == (pm[1] == -1)) LOGIC_ERROR;
        const int winner = pm[0] == -1 ? 0 : 1;
        oink->solve(n, winner, game->owner[n] == winner ? strategy[n] : -1);
    }

    delete[] pms;
    delete[] strategy;
    delete[] counts;
    delete[] tmp;
    delete[] best;
    delete[] dirty;
    delete[] cm;
    delete[] unstable;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
