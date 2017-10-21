#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <deque>
#include <cassert>

#include "dp.hpp"

namespace pg {

DPSolver::DPSolver(Oink *oink, Game *game, std::ostream &lgr) : PPPSolver(oink, game, lgr)
{
}

int
DPSolver::getRegionStatus(int i, int p)
{
    const int pl = p&1;

    // check if the region is closed in the subgame
    for (int j=i; j>=0 && priority[j]==p; j--) {
        if (disabled[j] or region[j] > p) {
            // escape not in region
            continue;
        } else if (owner[j] == pl) {
            // region won by escape owner
            // it is therefore closed if it has a strategy
            if (strategy[j] == -1) return -2;
        } else {
            // region lost by escape owner
            // it is therefore closed if there are no edges to lower regions
            for (auto to : out[j]) {
                if (disabled[to]) continue;
                int r = region[to];
                if (region_[to] > r) r = region_[to]; // update from region_
                if (r < p) return -2; // open
            }
        }
    }
    // closed in the subgame, so find lowest higher region for possible promotion
    int lowest = -1;
    for (auto j : regions[p]) {
        if (owner[j] != pl) {
            // losing node, find lowest higher ...
            for (auto to : out[j]) {
                /** HOT SPOT of the iterator and minor hot spot of obtaining the region **/
                if (disabled[to]) continue;
                int r = region[to];
                if (region_[to] > r) r = region_[to]; // update from region_
                if (r > p && (r < lowest || lowest == -1)) lowest = r;
                else if (r < p) LOGIC_ERROR;
            }
        }
    }
    return lowest; // -1 if dominion, <region> otherwise
}

void
DPSolver::run()
{
    // obtain highest priority and allocate arrays
    max_prio = priority[n_nodes-1];
    regions = new std::vector<int>[max_prio+1];
    region = new int[n_nodes];
    region_ = new int[n_nodes];
    strategy = new int[n_nodes];
    inverse = new int[max_prio+1];

    // initialize arrays
    for (int i=0; i<n_nodes; i++) region[i] = disabled[i] ? -2 : priority[i];
    for (int i=0; i<n_nodes; i++) region_[i] = -1;
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;

    // initialize reset values
    reset0 = -1;
    reset1 = -1;

    // start loop at last node (highest priority)
    int i = n_nodes - 1;

    // set to track delayed promotions
    std::vector<int> P;

    /**
     * Two loops: the outer (normal) loop, and the inner (promotion-chain) loop.
     * The outer loop does region setup and attractor on the full region.
     * The inner loop only attracts from the promoted region.
     */

    while (true) {
        // get current priority and skip all disabled/attracted nodes
        int p = i < 0 ? -1 : priority[i];
        while (i >= 0 and priority[i] == p and (disabled[i] or region[i] > p)) i--;

        if (i < 0) {
            int max = -1;
            for (int i=0; i<n_nodes; i++) if (max < region_[i]) max = region_[i];
            if (max == -1) break; // done

            // perform delayed promotions of highest player
            if (trace) logger << "performing delayed promotions of player " << (max&1) << std::endl;
            for (int i=0; i<n_nodes; i++) {
                if (region[i] == -2) continue;
                if (region_[i] != -1) {
                    if ((region_[i]&1) == (max&1)) {
                        region[i] = region_[i];
                        regions[region[i]].push_back(i);
                    }
                    region_[i] = -1;
                }
            }
            P.clear();
            // increase reset value if needed
            if (max&1) {
                if (max > reset0) reset0 = max-1;
            } else {
                if (max > reset1) reset1 = max-1;
            }
            if (trace) logger << "finished performing delayed promotions" << std::endl;
            i = inverse[max];
            continue;
        }

        // if empty, possibly reset and continue with next
        if (priority[i] != p) {
            if (!regions[p].empty()) resetRegion(p);
            continue;
        }

        inverse[p] = i;

        // PPP-DP: reset if lower than value
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
                    while (i >= 0 and priority[i] == p) i--;
                    break;
                } else if (res == -1) {
                    // found dominion
                    setDominion(p);
                    // restart algorithm and break inner loop
                    i = n_nodes - 1;
                    // reset everything... (sadly)
                    for (int j=0; j<n_nodes; j++) region[j] = disabled[j] ? -2 : priority[j];
                    for (int j=0; j<n_nodes; j++) strategy[j] = -1;
                    for (int j=0; j<n_nodes; j++) region_[j] = -1;
                    reset0 = reset1 = -1;
                    P.clear();
                    break;
                } else {
                    // check if maybe already delayed
                    if (res == region_[i]) break;
                    // check if we are locked or not
                    bool locked = false;
                    for (auto l : P) {
                        if ((l&1) == (res&1) && l < res) {
                            // locked for reason a
                            locked = true;
                            break;
                        }
                    }
                    if (!locked) {
                        for (int i=0; i<n_nodes; i++) {
                            // this loop is expensive, by the way
                            if (disabled[i]) continue;
                            if (region[i] < res && res < region_[i]) {
                                // locked for reason b
                                locked = true;
                                break;
                            }
                        }
                    }
                    if (locked) {
                        // found delayed promotion
                        if (trace) logger << "\033[1;33mdelayed \033[36m" << p << " \033[37mto \033[36m" << res << "\033[m" << std::endl;
                        delayed++;
                        // emulate delayed promotion in region_
                        for (int i : regions[p]) region_[i] = res;
                        // skip to next priority and break inner loop
                        while (i >= 0 and priority[i] == p) i--;
                        break;
                    } else {
                        // found instant promotion, perform it
                        promote(p, res);
                        // increase reset value if needed
                        if (res&1) {
                            if (res > reset0) reset0 = res-1;
                        } else {
                            if (res > reset1) reset1 = res-1;
                        }
                        // update the set of promotion targets P
                        // remove "from" from promotion targets
                        P.erase(std::remove(P.begin(), P.end(), p), P.end());
                        // add "to" to promotion targets
                        P.push_back(res);
                        // remove from region_ below res
                        for (int i=0; i<n_nodes; i++) if (region[i] != -2 && region_[i] <= res) region_[i] = -1;
                        // continue loop with the promotion target
                        i = inverse[res];
                        p = res;
                    }
                }
            }
        } else {
            // skip to next priority
            while (i >= 0 and priority[i] == p) i--;
        }
    }

    delete[] regions;
    delete[] region;
    delete[] region_;
    delete[] strategy;
    delete[] inverse;

    logger << "solved with " << promotions << "+" << delayed << " promotions." << std::endl;
}

}
