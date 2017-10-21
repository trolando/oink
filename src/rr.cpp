#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "rr.hpp"

namespace pg {

RRSolver::RRSolver(Oink *oink, Game *game, std::ostream &lgr) : PPSolver(oink, game, lgr)
{
}

/**
 * RR: only reset the region if what remains of the region can escape or has a bad strategy
 */
bool
RRSolver::checkRegion(int p)
{
    // check if the region escape is gone (thus require reset)
    if (regions[p].empty()) return true; // technically good

    // remove nodes that are no longer in the region (higher measure)
    // check for the remaining nodes that their strategy still stays
    // in the region and that losing nodes (except top nodes) cannot escape lower

    // first remove nodes no longer in the region
    auto &Rp = regions[p];
    Rp.erase(std::remove_if(Rp.begin(), Rp.end(),
        [&](const int n) {return region[n] > p;}), Rp.end());

    for (auto j : Rp) {
        // assert(priority[j] <= p && region[j] == p);
        if (disabled[j]) {
            // now disabled, requires a reset...
            return false;
        } else if (priority[j] == p) {
            // an escape node; if its strategy leaves the region, reset it
            if (strategy[j] != -1 && region[strategy[j]] != p) {
                strategy[j] = -1;
            }
        } else if (owner[j] == (p&1)) {
            // not-top winner
            // check if the strategy stays in the region
            // in very rare cases, strategy[j] == -1 (when fields are reset)
            if (strategy[j] == -1) return false;
            // requires a reset...
            if (region[strategy[j]] != p) return false;
        } else /*if (priority[j] != p)*/ {
            // not-top loser
            // check if it can escape in the subgame
            for (auto to : out[j]) {
                // it may be able to escape to a lower region if there have been resets
                if (region[to] != -2 && region[to] < p) return false;
            }
        }
    }

    return true;
}

void
RRSolver::run()
{
    // obtain highest priority and allocate arrays
    int max_prio = priority[n_nodes-1];
    regions = new std::vector<int>[max_prio+1];
    region = new int[n_nodes];
    strategy = new int[n_nodes];
    inverse = new int[max_prio+1];

    // initialize arrays
    for (int i=0; i<n_nodes; i++) region[i] = disabled[i] ? -2 : priority[i];
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

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

        // RR: only reset the region if:
        // - current node is promoted or attracted
        // - or region is empty
        // - or region does not fulfill conditions
        // This is checked by checkRegion()
        if (setupRegion(i, p, !checkRegion(p))) {
            // region not empty, maybe promote
            while (true) {
                if (trace >= 2) reportRegion(p);
                int res = getRegionStatus(i, p);
                if (res == -2) {
                    // not closed, skip to next priority and break inner loop
                    while (i >= 0 and priority[i] == p) i--;
                    break;
                } else if (res == -1) {
                    // found dominion, return
                    setDominion(p);
                    // restart algorithm and break inner loop
                    i = n_nodes - 1;
                    break;
                } else {
                    // found promotion, promote
                    if (trace >= 2) printState();
                    promote(p, res);
                    if (trace >= 2) printState();
                    // continue inner loop with the higher priority
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
