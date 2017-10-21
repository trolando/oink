#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stack>
#include "tspm.hpp"

namespace pg {

TSPMSolver::TSPMSolver(Oink *oink, Game *game, std::ostream &lgr) : Solver(oink, game, lgr)
{
}

TSPMSolver::~TSPMSolver()
{
}

/**
 * Returns true if a progress measure "a" is less than "b"
 * up to and including priority <d>, for player <pl>.
 */
bool
TSPMSolver::pm_less(int* a, int* b, int d, int pl)
{
    // cases where a or b is Top
    if (b[pl] == -1) return a[pl] != -1;
    else if (a[pl] == -1) return false;
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
TSPMSolver::pm_copy(int *dst, int *src, int pl)
{
    for (int i=pl; i<k; i+=2) dst[i] = src[i];
}

/**
 * Write pm to ostream.
 */
void
TSPMSolver::pm_stream(std::ostream &out, int *pm)
{
    bool tope = pm[0] == -1;
    bool topo = pm[1] == -1;
    out << " {";
    if (tope) out << " \033[1;33mTe\033[m";
    else out << " " << pm[0];
    if (topo) out << " \033[1;33mTo\033[m";
    else out << " " << pm[1];
    for (int i=2; i<k; i++) {
        if (i&1) out << " " << (topo ? 0 : pm[i]);
        else     out << " " << (tope ? 0 : pm[i]);
    }
    out << " } ";
}

/**
 * Perform update for player <pl>, node with priority <d>.
 */
void
TSPMSolver::Prog(int *dst, int *src, int d, int pl)
{
    // check if top
    if (src[pl] == -1) {
        dst[pl] = -1;
        return;
    }

    // set every value lower than <d> to 0.
    int i = pl;
    for (; i<d; i+=2) dst[i] = 0;

    int carry = (d == i) ? 1 : 0; // only increase if <d> has parity <pl>

    for (; i<k; i+=2) {
        // increase or copy for same parity
        int v = src[i] + carry;
        if (v > counts[i]) {
            dst[i] = 0;
            carry = 1;
        } else {
            dst[i] = v;
            carry = 0;
        }
    }

    // check if top
    if (carry) dst[pl] = -1;
}

bool
TSPMSolver::canlift(int node, int pl)
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
TSPMSolver::lift(int node, int target)
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
    if (pm[pl_max] != -1) {
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
    if (pm[pl_min] != -1 and (target == -1 or target == strategy[node])) {
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

    bool ch0 = best_ch0 != -1;
    bool ch1 = best_ch1 != -1;

    if (ch0 or ch1) {
        if (trace) {
            logger << "\033[1;32mupdated node " << node << "/" << d << (owner[node]?" (odd)":" (even)") << "\033[m to";
            pm_stream(logger, pm);
            logger << std::endl;
        }
        // update counts if a changed pm is now Top, but only if priority has same parity as winner
        if (ch0 and pm[0] == -1 and (d&1) == 0) counts[d]--;
        if (ch1 and pm[1] == -1 and (d&1) == 1) counts[d]--;
        // increase count and return true
        lift_count++;
        return true;
    } else {
        return false;
    }
}

void
TSPMSolver::update(int pl)
{
    std::queue<int> q;

    // find unstable nodes (for measure <pl>)
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        unstable[i] = 0; // first mark as stable
        if (pms[k*i + pl] == -1 or canlift(i, pl)) {
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
TSPMSolver::run()
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

    // initialize all measures to 0
    for (int i=0; i<k*n_nodes; i++) pms[i] = 0;

    // initialize strategy to -1
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    // initialize counts for each priority
    for (int i=0; i<k; i++) counts[i] = 0;
    for (int i=0; i<n_nodes; i++) if (disabled[i] == 0) counts[priority[i]]++;

    // initialize all nodes as not dirty
    for (int n=0; n<n_nodes; n++) dirty[n] = 0;

    // set number of lifts and lift attempts to 0
    lift_count = lift_attempt = 0;

    /**
     * Strategy that updates predecessors then marks updated predecessors for processing.
     * Uses a queue/stack to store the dirty vertices.
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
    delete[] unstable;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
