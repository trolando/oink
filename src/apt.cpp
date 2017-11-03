#include <iostream>
#include <cassert>
#include <cstring>

#include "apt.hpp"

namespace pg {

APTSolver::APTSolver(Oink *oink, Game *game, std::ostream &lgr) : Solver(oink, game, lgr)
{
    // sanity check if the game is properly sorted
    for (int i=1; i<n_nodes; i++) assert(priority[i-1] <= priority[i]);
}

APTSolver::~APTSolver()
{
}

void
APTSolver::run()
{
    unsigned int iterations = 0;

    // allocate and initialize data structures
    int *F_in_Y = new int[n_nodes];
    memset(F_in_Y, 0, sizeof(int[n_nodes]));

    for (;;) {
        iterations++;

        /**
         * We update F_in_Y from low to high.
         * We process per block of priorities.
         * When a block changes, all lower F_in_Y are reset.
         *
         * The algorithm refines a partition into two winning regions.
         * APT terminates when the winning regions are inductive, i.e.,
         * no player escapes in one step.
         * If Y (set of vertices of priority pr not won by player pr&1) for some
         * priority is updated, then all lower priorities are reset.
         */
        bool changed = false;
        int cur_pr = priority[0]; // priority of current block
        int first_idx = 0;        // first idx of current block

        for (int i=0; i<n_nodes;i++) {
            if (disabled[i]) continue;
            // like Verver's Zielonka trick, not per priority change,
            // but per parity change (mimics compression)
            if ((priority[i]&1) != (cur_pr&1)) {
                // new block!
                if (changed) break;
                // the previous block did not change, continue with new block
                cur_pr = priority[i];
                first_idx = i;
            }
            if (F_in_Y[i]) continue; // Y increases monotonically...
            // compute Y for current vertex
            int y;
            if (owner[i] == 0) {
                y = 1;
                const int *_out = outs + outa[i];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    if ((priority[to] & 1) != F_in_Y[to]) continue;
                    y = 0;
                    break;
                }
            } else {
                y = 0;
                const int *_out = outs + outa[i];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    if ((priority[to] & 1) == F_in_Y[to]) continue;
                    y = 1;
                    break;
                }
            }
            // check if now in Y
            if ((cur_pr&1) == y) continue; // not in Y
            changed = true;
            F_in_Y[i] = 1;
            if (trace >= 2) logger << "vertex " << i << " is now in Y^" << cur_pr << std::endl;
        }

        if (!changed) break; // nothing changed, we're done

        // the previous block changed, reset all F_in_Y from 0 to first_idx
        memset(F_in_Y, 0, sizeof(int[first_idx]));
        if (trace) logger << "restarting after update to Y of prio " << cur_pr << std::endl;
    }

    // done -- unfortunately, APT does not compute a strategy!
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        oink->solve(i, F_in_Y[i] ? 1-(priority[i]&1) : (priority[i]&1), -1); 
    }

    // free allocated data structures
    delete[] F_in_Y;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
