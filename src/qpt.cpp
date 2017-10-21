#include <iostream>
#include <iomanip>
#include <stack>
#include "qpt.hpp"

namespace pg {

QPTSolver::QPTSolver(Oink *oink, Game *game, std::ostream &lgr) : Solver(oink, game, lgr)
{
}

QPTSolver::~QPTSolver()
{
}

/**
 * Returns true if a progress measure "a" is less than "b"
 */
static bool
pm_less(int* a, int* b, int k)
{
    // remember: _ < 5 < 3 < 1 < 0 < 2 < 4 < 6
    for (int i=0; i<k; i++) {
        // if equal, then continue with next
        if (a[i] == b[i]) continue;
        // case where a[i] is _ and b[i] is not
        if (a[i] == -1) return true;
        // case where b[i] is _ and a[i] is not
        if (b[i] == -1) return false;
        // cases where a[i] is odd
        if (a[i] & 1) {
            // .. and b[i] is also odd
            if (b[i] & 1) return a[i] > b[i];
            // .. and b[i] is even
            else return true;
        }
        // case where b[i] is odd and a[i] is even
        if (b[i] & 1) return false;
        // case where a[i] and b[i] are even
        return a[i] < b[i];
    }
    // case where a == b
    return false;
}

/**
 * Returns true if a progress measure "a" is less than "b" (odd)
 */
/*static bool
pm_less_o(int* a, int* b, int k)
{
    // remember: _ < 6 < 4 < 2 < 0 < 1 < 3 < 5
    for (int i=0; i<k; i++) {
        // if equal, then continue with next
        if (a[i] == b[i]) continue;
        // case where a[i] is _ and b[i] is not
        if (a[i] == -1) return true;
        // case where b[i] is _ and a[i] is not
        if (b[i] == -1) return false;
        // cases where a[i] is odd
        if (a[i] & 1) {
            // .. and b[i] is also odd
            if (b[i] & 1) return a[i] < b[i];
            // .. and b[i] is even
            else return false;
        }
        // case where b[i] is odd and a[i] is even
        if (b[i] & 1) return true;
        // case where a[i] and b[i] are even
        return a[i] > b[i];
    }
    // case where a == b
    return false;
}*/

/**
 * Compute updated witness with new vertex with priority "d"
 */
static int
up(int *dst, int *src, int d, int k)
{
    // check if this is a "won" witness
    if ((src[0]&1) == 0) {
        for (int i=0; i<k; i++) dst[i] = src[i];
        return -1;
    }

    // first compute j for Lemma 3.3
    int j33 = -1;
    for (int i=k-1; i>=0; i--) {
        // start from the end to find longest even sequence
        if (src[i] == -1 || src[i]&1) {
            // src[i] is _ or odd, try here!
            j33 = i;
            // and verify other half of the property (all higher are _ or bigger than d)
            for (int j=i-1; j>=0; j--) {
                if (src[j] != -1 && src[j] < d) {
                    // not good!
                    j33 = -1;
                    break;
                }
            }
            break;
        }
    }

    // then compute j for Lemma 3.4
    int j34 = -1;
    for (int i=0; i<k; i++) {
        // start from the begin, find first index where the Lemma applies
        if (src[i] != -1 && src[i] < d) {
            j34 = i;
            break;
        }
    }

    // determine best index min(j33, j34) or -1
    int j = j33 == -1 ? j34 : (j34 == -1 ? j33 : (j33 < j34 ? j33 : j34));

    // compute resulting progress measure
    for (int i=0; i<k; i++) {
        if (i == j) {
            dst[i] = d;
            while (++i != k) dst[i] = -1;
            break;
        } else {
            dst[i] = src[i];
        }
    }

    return j;
}

/**
 * Compute smallest larger progress measure with dst[0] set to _
 */
static bool
bump(int *dst, int *src, int k)
{
    // reminder: _ 5 3 1 0 2 4 6...
    dst[k-1] = -1;
    for (int i=k-2; i>=0; i--) {
        // if odd != 1 then we can just subtract 2 and done
        if ((src[i] & 1) && src[i] > 1) {
            dst[i] = src[i]-2;
            while (--i >= 0) dst[i] = src[i];
            return true;
        }
        // if it's 1, and the next one is higher, then set 2 and done
        if (src[i] == 1 && (i == 0 || src[i-1] > 1)) {
            dst[i] = 0;
            while (--i >= 0) dst[i] = src[i];
            return true;
        }
        // if it's even, and the next one is higher, then increase and done
        // probably the only one that is relevant, since this algorithm is used
        // when an "even" sequence must be bumped
        if ((src[i]&1) == 0 && (i == 0 || src[i-1] >= (src[i]+2))) {
            dst[i] = src[i]+2;
            while (--i >= 0) dst[i] = src[i];
            return true;
        }
        // no rule applies, set _
        dst[i] = -1;
    }
    // no rule applied, dst is trash, return false
    return false;
}

/**
 * Print out a progress measure <pm> to the given stream
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
 * Perform antagonistic update.
 * Return "true" if dst <> src, "false" otherwise.
 */
static bool
au(int *dst, int *src, int d, int k)
{
    int j = up(dst, src, d, k);
    if (j == -1) return false; // no difference, already smallest...

    // try antagonistic update
    int tmp[k];
    // compute a tmp ] src and tmp_0 equals _
    if (!bump(tmp, src, k)) return true; // no bump possible; but dst != src so return true
    // now update it
    up(tmp, tmp, d, k);

    // if smaller, use small version
    if (pm_less(tmp, dst, k)) {
        // printf("\033[1;32mAntagonistic update used bumped PM\033[m\n");
        for (int i=0; i<k; i++) dst[i] = tmp[i];
    }
    // and finally, check if different
    for (int i=0; i<k; i++) if (dst[i] != src[i]) return true;
    return false;
}

bool
QPTSolver::try_lift(int n, std::vector<int> &vec)
{
    int tmp[k], res[k];

    vec.clear();

    if (owner[n] == 0) {
        // in Ve so max!
        bool first = true;
        for (int to : out[n]) {
            if (disabled[to]) continue; // not looking at this
            au(tmp, pm_nodes + k*to, priority[n], k);

            if (first or pm_less(res, tmp, k)) {
                first = false;
                for (int i=0; i<k; i++) res[i] = tmp[i];
                vec.clear();
                vec.push_back(to);
            } else if (!pm_less(tmp, res, k)) {
                vec.push_back(to);
            }
        }
        int *pm = pm_nodes + k*n;
        return !vec.empty() && pm_less(pm, res, k);
    } else {
        // in Vo so min!
        bool first = true;
        for (int to : out[n]) {
            if (disabled[to]) continue; // not looking at this
            au(tmp, pm_nodes + k*to, priority[n], k);

            if (first || pm_less(tmp, res, k)) {
                first = false;
                for (int i=0; i<k; i++) res[i] = tmp[i];
                vec.clear();
                vec.push_back(to);
            } else if (!pm_less(res, tmp, k)) {
                vec.push_back(to);
            }
        }
        int *pm = pm_nodes + k*n;
        return !vec.empty() && pm_less(pm, res, k);
    }
}

bool
QPTSolver::lift(int n)
{
    int *pm = pm_nodes + k*n; // obtain ptr to current progress measure
    bool changed_this = false; // whether the progress measure was updated

    int tmp[k], res[k];

    if (trace >= 2) {
        logger << "\033[1mupdating node " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m with current progress measure";
        pm_stream(logger, pm, k);
        logger << std::endl;
    }

    if (owner[n] == 0) {
        // in Ve so max!
        // compute au for each successor and if it's larger, take it
        for (int to : out[n]) {
            if (disabled[to]) continue; // not looking at this
            au(tmp, pm_nodes + k*to, priority[n], k);

            if (trace >= 2) {
                logger << "successor node " << to << "/" << priority[to] << " results in";
                pm_stream(logger, tmp, k);
                logger << std::endl;
            }

            if (pm_less(pm, tmp, k)) {
                for (int i=0; i<k; i++) pm[i] = tmp[i];
                changed_this = true;
            }
        }
    } else {
        // in Vo so min!
        // compute au for each successor and if it's smaller, take it
        bool first = true;
        int best_to = -1;
        for (int to : out[n]) {
            if (disabled[to]) continue; // not looking at this
            au(tmp, pm_nodes + k*to, priority[n], k);

            if (trace >= 2) {
                logger << "successor node " << to << "/" << priority[to] << " results in";
                pm_stream(logger, tmp, k);
                logger << std::endl;
            }

            if (first || pm_less(tmp, res, k)) {
                first = false;
                for (int i=0; i<k; i++) res[i] = tmp[i];
                best_to = to;
            }
        }
        // now "res" contains the smallest au
        if (pm_less(pm, res, k)) {
            for (int i=0; i<k; i++) pm[i] = res[i];
            changed_this = true;
        }
        strategy[n] = best_to; // because maybe just different strategy
    }

    if (trace >= 2) {
        if (strategy[n] != -1) logger << "updated pm for node " << n << " to node " << strategy[n] << " to";
        else logger << "updated pm for node " << n << " to";
        pm_stream(logger, pm_nodes+k*n, k);
        logger << std::endl;
    }

    return changed_this;
}

bool
QPTSolver::liftR(int node, int target)
{
    int *pm = pm_nodes + k*node; // obtain ptr to current progress measure

    if (owner[node] == 0) {
        // in Ve so max!
        int tmp[k];
        au(tmp, pm_nodes + k*target, priority[node], k);
        if (!pm_less(pm, tmp, k)) return false;
        for (int i=0; i<k; i++) pm[i] = tmp[i];
        return true;
    } else {
        // in Vo so min!
        if (strategy[node] != target) return false; // current strategy is unchanged...
        bool first = true;
        int best_to = -1;
        int tmp[k], res[k];
        for (int to : out[node]) {
            if (disabled[to]) continue; // not looking at this
            au(tmp, pm_nodes + k*to, priority[node], k);
            if (first || pm_less(tmp, res, k)) {
                first = false;
                for (int i=0; i<k; i++) res[i] = tmp[i];
                best_to = to;
            }
        }
        strategy[node] = best_to; // because maybe just different strategy
        // now "res" contains the smallest au
        if (!pm_less(pm, res, k)) return false;
        for (int i=0; i<k; i++) pm[i] = res[i];
        return true;
    }
}

void
QPTSolver::print_state(std::vector<int>* choices)
{
    std::vector<int> vec;
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        int *pm = pm_nodes + k*i;
        logger << std::setfill('0') << std::setw(2) << i << " : " << std::setfill('0') << std::setw(2) << priority[i] << " :";
        pm_stream(logger, pm, k);
        logger << " : ";
        if (try_lift(i, vec)) logger << "\033[1;33m";
        logger << "=>";
        for (int to : vec) logger << " " << to;
        logger << "\033[m";
        logger << " " << "\033[1;35mcur:";
        for (auto to : choices[i]) logger << " " << to;
        logger << "\033[m\n";
    }
}

void
QPTSolver::run()
{
    // determine k
    int even_vertices = 0;
    for (int i=0; i<n_nodes; i++) if (!disabled[i] && (priority[i]&1) == 0) even_vertices++;
    // find k s.t. 2^{k-1} >= even_vertices+1
    even_vertices++;
    k=1;
    while (even_vertices != 0) {
        k++;
        even_vertices >>= 1;
    }

    // now create the data structure, for each node
    pm_nodes = new int[k*n_nodes];
    for (int i=0; i<k*n_nodes; i++) pm_nodes[i] = -1; // initialize all to _

    strategy = new int[n_nodes];
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    /**
     * Interactive lifting (fun!)
     */
    /*
    std::vector<int> *choices = new std::vector<int>[n_nodes];
    while (true) {
        print_state(choices);
        int n;
        std::cin >> n;
        if (n < 0 || n >= n_nodes) continue;
        try_lift(n, choices[n]);
        lift(n);
    }
    delete[] choices;
    */

    lift_count = lift_attempt = 0;

    const bool use_r = true; // use liftR instead of lift (look at predecessors)
    const bool use_queue = false; // use a queue
    const bool use_stack = true; // use a stack (if both are set, queue has priority)

    if (use_r) {
        if (use_queue or use_stack) {
            /**
             * Strategy that updates predecessors then marks updated predecessors for processing.
             * Uses a queue/stack to store the dirty vertices.
             */

            int *dirty = new int[n_nodes];
            for (int n=0; n<n_nodes; n++) dirty[n] = 0;

            std::deque<int> todo;
            for (int n=0; n<n_nodes; n++) {
                if (disabled[n]) continue;
                lift_attempt++;
                if (lift(n)) {
                    if (dirty[n]) dirty[n] = 0;
                    lift_count++;
                    for (int from : in[n]) {
                        if (disabled[from]) continue;
                        lift_attempt++;
                        if (!liftR(from, n)) continue;
                        lift_count++;
                        if (dirty[from]) continue;
                        dirty[from] = 1;
                        todo.push_back(from);
                    }
                }
            }

            while (!todo.empty()) {
                int n;
                if (use_queue) {
                    n = todo.front();
                    todo.pop_front();
                } else {
                    n = todo.back();
                    todo.pop_back();
                }
                dirty[n] = 0;
                if ((lift_count % n_nodes) == 0) {
                    if (trace) logger << "\033[1;38;5;208mIteration " << (lift_count / n_nodes) << "\033[m" << std::endl;
                }
                for (int from : in[n]) {
                    if (disabled[from]) continue;
                    lift_attempt++;
                    if (!liftR(from, n)) continue;
                    lift_count++;
                    if (trace) {
                        logger << "\033[1mUpdated node " << from << "/" << priority[from] << (owner[from]?" (odd)":" (even)") << "\033[m due to node " << n << "/" << priority[n] << " to";
                        pm_stream(logger, pm_nodes + from*k, k);
                        logger << std::endl;
                    }
                    if (dirty[from]) continue;
                    dirty[from] = 1;
                    todo.push_back(from);
                }
            }

            delete[] dirty;
        } else {
            /**
             * Strategy that updates predecessors then marks updated predecessors for processing.
             * But does not use a queue or a stack.
             */
            bool changed;
            int *dirty = new int[n_nodes]; // dirty is initialized in the first round

            ++iterations;
            if (trace) logger << "\033[1;38;5;208mIteration " << iterations << "\033[m" << std::endl;
            for (int n=0; n<n_nodes; n++) {
                if (disabled[n]) continue;
                lift_attempt++;
                if (lift(n)) {
                    if (dirty[n]) dirty[n] = 0;
                    lift_count++;
                    for (int from : in[n]) {
                        if (disabled[from]) continue;
                        lift_attempt++;
                        if (!liftR(from, n)) continue;
                        lift_count++;
                        if (dirty[from] == 0) dirty[from] = 1;
                    }
                }
            }

            do {
                ++iterations;
                if (trace) logger << "\033[1;38;5;208mIteration " << iterations << "\033[m" << std::endl;

                changed = false;
                for (int n=0; n<n_nodes; n++) {//n_nodes-1; n>=0; n--) {
                    if (disabled[n]) continue; // not looking at this
                    if (!dirty[n]) continue; // not dirty, continue
                    dirty[n] = 0; // unmark

                    for (int from : in[n]) {
                        if (disabled[from]) continue;
                        lift_attempt++;
                        if (!liftR(from, n)) continue;
                        lift_count++;
                        changed = true;
                        if (dirty[from] == 0) dirty[from] = 1;
                    }
                }
            } while (changed);

            delete[] dirty;
        }
    } else if (use_queue or use_stack) {
        /**
         * Strategy that updates nodes then marks predecessors for processing.
         * Uses a queue/stack to store the dirty vertices.
         */

        int *dirty = new int[n_nodes];
        for (int n=0; n<n_nodes; n++) dirty[n] = 0;

        std::deque<int> todo;
        for (int n=0; n<n_nodes; n++) {
            if (disabled[n]) continue;
            lift_attempt++;
            if (lift(n)) {
                lift_count++;
                for (int from : in[n]) {
                    if (disabled[from] || dirty[from]) continue;
                    dirty[from] = 1;
                    todo.push_back(from);
                }
            }
        }

        while (!todo.empty()) {
            int n;
            if (use_queue) {
                n = todo.front();
                todo.pop_front();
            } else {
                n = todo.back();
                todo.pop_back();
            }
            dirty[n] = 0;
            lift_attempt++;
            if (lift(n)) {
                if (trace) {
                    logger << "\033[1mUpdated node " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m to";
                    pm_stream(logger, pm_nodes + n*k, k);
                    logger << std::endl;
                }

                lift_count++;
                for (int from : in[n]) {
                    if (disabled[from] || dirty[from]) continue;
                    dirty[from] = 1;
                    todo.push_back(from);
                }
            }
        }

        delete[] dirty;
    } else {
        /**
         * Strategy that updates nodes then marks predecessors for processing.
         * But does not use a queue or a stack.
         */
        bool first_round = true;
        bool changed;
        int *dirty = new int[n_nodes]; // dirty is initialized in the first round
        do {
            ++iterations;
            if (trace) logger << "\033[1;38;5;208mIteration " << iterations << "\033[m" << std::endl;

            changed = false;
            for (int n=0; n<n_nodes; n++) {//n_nodes-1; n>=0; n--) {
                if (disabled[n]) continue; // not looking at this
                if (!first_round && !dirty[n]) continue; // not dirty or first round, continue
                dirty[n] = 0; // unmark

                lift_attempt++;
                if (lift(n)) {
                    if (trace) {
                        logger << "\033[1mUpdated node " << n << "/" << priority[n] << (owner[n]?" (odd)":" (even)") << "\033[m to";
                        pm_stream(logger, pm_nodes + n*k, k);
                        logger << std::endl;
                    }
                    lift_count++;
                    changed = true;
                    for (int from : in[n]) dirty[from] = 1;
                }
            }
            first_round = false;
        } while (changed);
        delete[] dirty;
    }

    // Now set dominions and derive strategy for odd.
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        int *pm = pm_nodes + i*k;
        int winner = pm[0]&1;
        game->dominion[i] = winner;
        if (game->owner[i] == 1 && winner == 1) game->strategy[i] = strategy[i];
    }

    delete[] pm_nodes;
    delete[] strategy;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
