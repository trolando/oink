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
    iterations = 0;

    // allocate and initialize data structures
    int *F_in_Y = new int[n_nodes];
    int *Y = new int[n_nodes];
    memset(F_in_Y, 0, sizeof(int[n_nodes]));

    for (;;) {
        iterations++;
        // first compute Y
        for (int n=0; n<n_nodes; n++) {
            if (disabled[n]) continue;
            int y;
            if (owner[n] == 0) {
                y = 1;
                const int *_out = outs + outa[n];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    if ((priority[to] & 1) != F_in_Y[to]) continue;
                    y = 0;
                    break;
                }
            } else {
                y = 0;
                const int *_out = outs + outa[n];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    if ((priority[to] & 1) == F_in_Y[to]) continue;
                    y = 1;
                    break;
                }
            }
            Y[n] = y;
        }
        // report V and Y if trace >= 2
        if (trace >= 2) {
            logger << "V for 0:";
            for (int n=0; n<n_nodes; n++) {
                if ((priority[n] & 1) != F_in_Y[n]) continue;
                logger << " " << n << "/" << priority[n];
            }
            logger << std::endl;
            logger << "Y for 0:";
            for (int n=0; n<n_nodes; n++) {
                if (Y[n]) continue;
                logger << " " << n << "/" << priority[n];
            }
            logger << std::endl;
        }
        // now go up, steps of priority, until something changed
        bool changed = false;
        int tpr = -1; // for trace
        int i = -1;
        for (;;) {
            // find next j
            int j = i + 1;
            while (j != n_nodes and disabled[j]) j++;
            if (j == n_nodes) break;
            // update F_in_Y for next F
            const int pr = tpr = priority[j];
            const int pl = pr & 1;
            for (;j!=n_nodes;j++) {
                if (disabled[j]) continue;
                if (priority[j] != pr) break;
                int in_Y = pl != Y[j] ? 1 : 0;
                if (F_in_Y[j] != in_Y) {
                    changed = true;
                    F_in_Y[j] = in_Y;
                }
            }
            // update state
            if (changed) break;
            if (j == n_nodes) break;
            // it didn't change and we're not done, go up
            i = j - 1;
        }
        // if nothing changed, then we're done
        if (!changed) break;
        // something changed, reset all F_in_Y from 0 to i
        memset(F_in_Y, 0, sizeof(int[i+1]));
        if (trace) logger << "update to Y of prio " << tpr << std::endl;
    }

    // done -- unfortunately, APT does not compute a strategy!
    for (int i=0; i<n_nodes; i++) game->dominion[i] = Y[i];

    // free allocated data structures
    delete[] F_in_Y;
    delete[] Y;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
