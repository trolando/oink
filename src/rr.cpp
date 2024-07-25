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

#include <algorithm>
#include <cassert>

#include "rr.hpp"

namespace pg {

RRSolver::RRSolver(Oink& oink, Game& game) : PPSolver(oink, game)
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
        // assert(priority(j) <= p && region[j] == p);
        if (disabled[j]) {
            // now disabled, requires a reset...
            return false;
        } else if (priority(j) == p) {
            // an escape node; if its strategy leaves the region, reset it
            if (strategy[j] != -1 && region[strategy[j]] != p) {
                strategy[j] = -1;
            }
        } else if (owner(j) == (p&1)) {
            // not-top winner
            // check if the strategy stays in the region
            // in very rare cases, strategy[j] == -1 (when fields are reset)
            if (strategy[j] == -1) return false;
            // requires a reset...
            if (region[strategy[j]] != p) return false;
        } else /*if (priority(j) != p)*/ {
            // not-top loser
            // check if it can escape in the subgame
            for (auto curedge = outs(j); *curedge != -1; curedge++) {
                int to = *curedge;
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
    int max_prio = priority(nodecount()-1);
    regions = new std::vector<int>[max_prio+1];
    region = new int[nodecount()];
    strategy = new int[nodecount()];
    inverse = new int[max_prio+1];

    // initialize arrays
    for (int i=0; i<nodecount(); i++) region[i] = disabled[i] ? -2 : priority(i);
    for (int i=0; i<nodecount(); i++) strategy[i] = -1;

    // start loop at last node (highest priority)
    int i = nodecount() - 1;

    // reset statistics
    promotions = 0;

    /**
     * Two loops: the outer (normal) loop, and the inner (promotion-chain) loop.
     * The outer loop does region setup and attractor on the full region.
     * The inner loop only attracts from the promoted region.
     */

    while (i >= 0) {
        // get current priority and skip all disabled/attracted nodes
        int p = priority(i);
        while (i >= 0 and priority(i) == p and (disabled[i] or region[i] > p)) i--;
        if (i < 0) break;

        // if empty, possibly reset and continue with next
        if (priority(i) != p) {
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
                    while (i >= 0 and priority(i) == p) i--;
                    break;
                } else if (res == -1) {
                    // found dominion, return
                    setDominion(p);
                    // restart algorithm and break inner loop
                    i = nodecount() - 1;
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
            while (i >= 0 and priority(i) == p) i--;
        }
    }

    delete[] regions;
    delete[] region;
    delete[] strategy;
    delete[] inverse;

    logger << "solved with " << promotions << " promotions." << std::endl;
}

}
