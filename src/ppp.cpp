#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "ppp.hpp"

namespace pg {

PPPSolver::PPPSolver(Oink *oink, Game *game, std::ostream &lgr) : PPSolver(oink, game, lgr)
{
}

void
PPPSolver::run()
{
    // obtain highest priority and allocate arrays
    max_prio = game->priority[n_nodes-1];
    regions = new std::vector<int>[max_prio+1];
    region = new int[n_nodes];
    strategy = new int[n_nodes];
    inverse = new int[max_prio+1];

    // initialize arrays
    for (int i=0; i<n_nodes; i++) region[i] = disabled[i] ? -2 : priority[i];
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    // initialize reset values
    reset0 = -1;
    reset1 = -1;

    // start loop at last node (highest priority)
    int i = n_nodes - 1;

    /**
     * Two loops: the outer (normal) loop, and the inner (promotion-chain) loop.
     * The outer loop does region setup and attractor on the full region.
     * The inner loop only attracts from the promoted region.
     */

    while (i >= 0) {
        // get current priority and skip all disabled/attracted nodes
        int p = priority[i];
        while (i >= 0 and priority[i] == p and (disabled[i] or region[i] > p)) i--;
        if (i < 0) break;

        // if empty, possibly reset and continue with next
        if (priority[i] != p) {
            if (!regions[p].empty()) resetRegion(p);
            continue;
        }

        inverse[p] = i;

        // PPP: reset if lower than value
        bool reset = false;
        if (p&1) {
            if (p <= reset1) {
                reset = true;
                reset1 = p-2;
            }
        } else {
            if (p <= reset0) {
                reset = true;
                reset0 = p-2;
            }
        }
        if (setupRegion(i, p, reset)) {
            // region not empty, maybe promote
            while (true) {
                if (trace >= 2) reportRegion(p);
                int res = getRegionStatus(i, p);
                if (res == -2) {
                    // not closed, skip to next priority and break inner loop
                    while (i >= 0 and priority[i] >= p) i--;
                    break;
                } else if (res == -1) {
                    // found dominion, return
                    setDominion(p);
                    // restart algorithm and break inner loop
                    i = n_nodes - 1;
                    break;
                } else {
                    // found promotion, promote
                    promote(p, res);
                    // increase reset value if needed
                    if (res&1) {
                        if (res > reset0) reset0 = res-1;
                    } else {
                        if (res > reset1) reset1 = res-1;
                    }
                    // continue inner loop with higher priority
                    i = inverse[res];
                    p = res;
                }
            }
        } else {
            // skip to next priority
            while (i >= 0 and priority[i] == p) i--;
        }
    }

    delete[] regions;
    delete[] region;
    delete[] strategy;
    delete[] inverse;

    logger << "solved with " << promotions << " promotions." << std::endl;
}

}
