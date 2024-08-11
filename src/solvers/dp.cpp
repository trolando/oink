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

#include "dp.hpp"

namespace pg {

DPSolver::DPSolver(Oink& oink, Game& game) : PPPSolver(oink, game)
{
}

int
DPSolver::getRegionStatus(int i, int p)
{
    const int pl = p&1;

    // check if the region is closed in the subgame
    for (int j=i; j>=0 && priority(j)==p; j--) {
        if (disabled[j] or region[j] > p) {
            // escape not in region
            continue;
        } else if (owner(j) == pl) {
            // region won by escape owner
            // it is therefore closed if it has a strategy
            if (strategy[j] == -1) return -2;
        } else {
            // region lost by escape owner
            // it is therefore closed if there are no edges to lower regions
            for (auto curedge = outs(j); *curedge != -1; curedge++) {
                int to = *curedge;
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
        if (owner(j) != pl) {
            // losing node, find lowest higher ...
            for (auto curedge = outs(j); *curedge != -1; curedge++) {
                int to = *curedge;
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
    max_prio = priority(nodecount()-1);
    regions = new std::vector<int>[max_prio+1];
    region = new int[nodecount()];
    region_ = new int[nodecount()];
    strategy = new int[nodecount()];
    inverse = new int[max_prio+1];

    // initialize arrays
    for (int i=0; i<nodecount(); i++) region[i] = disabled[i] ? -2 : priority(i);
    for (int i=0; i<nodecount(); i++) region_[i] = -1;
    for (int i=0; i<nodecount(); i++) strategy[i] = -1;

    // initialize reset values
    reset0 = -1;
    reset1 = -1;

    // start loop at last node (highest priority)
    int i = nodecount() - 1;

    // reset statistics
    promotions = 0;
    int performances = 0;

    // set to track delayed promotions
    std::vector<int> P;
    int del0=0, del1=0, discarded=0;

    /**
     * Two loops: the outer (normal) loop, and the inner (promotion-chain) loop.
     * The outer loop does region setup and attractor on the full region.
     * The inner loop only attracts from the promoted region.
     */

    while (true) {
        // get current priority and skip all disabled/attracted nodes
        int p = i < 0 ? -1 : priority(i);
        while (i >= 0 and priority(i) == p and (disabled[i] or region[i] > p)) i--;

        if (i < 0) {
            int max = -1;
            for (int i=0; i<nodecount(); i++) if (max < region_[i]) max = region_[i];
            if (max == -1) break; // done

            // perform delayed promotions of highest player
            if (trace) logger << "performing delayed promotions of player " << (max&1) << std::endl;
            performances++;
            for (int i=0; i<nodecount(); i++) {
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
                promotions += del1;
                discarded += del0;
            } else {
                if (max > reset1) reset1 = max-1;
                promotions += del0;
                discarded += del1;
            }
            if (trace) logger << "finished performing delayed promotions" << std::endl;
            i = inverse[max];
            del1 = del0 = 0;
            continue;
        }

        // if empty, possibly reset and continue with next
        if (priority(i) != p) {
            if (!regions[p].empty()) {
                resetRegion(p);
                // but then we must also reset everything lower...
                if (p&1) reset1 = p-2;
                else reset0 = p-2;
            }
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
                    while (i >= 0 and priority(i) == p) i--;
                    break;
                } else if (res == -1) {
                    // found dominion
                    setDominion(p);
                    // restart algorithm and break inner loop
                    i = nodecount() - 1;
                    // reset everything... (sadly)
                    // for (int j=0; j<nodecount(); j++) region[j] = disabled[j] ? -2 : priority(j);
                    // for (int j=0; j<nodecount(); j++) strategy[j] = -1;
                    for (int j=0; j<nodecount(); j++) region_[j] = -1;
                    reset0 = priority(nodecount()-1);
                    reset1 = priority(nodecount()-1);
                    if (reset0&1) reset0--;
                    else reset1--;
                    P.clear();
                    break;
                } else {
                    // check if maybe already delayed
                    if (res == region_[i]) break;
                    // check if we are locked or not
                    bool locked = false;
                    for (auto l : P) {
                        if ((l&1) != (res&1) && l < res) {
                            // locked for reason a
                            locked = true;
                            break;
                        }
                    }
                    if (!locked) {
                        for (int i=0; i<nodecount(); i++) {
                            // this loop is expensive, by the way
                            if (disabled[i]) continue;
                            if (region[i] < res && res <= region_[i]) {
                                // locked for reason b
                                locked = true;
                                break;
                            }
                        }
                    }
                    if (locked) {
                        // found delayed promotion
                        if (trace) {
                            if (trace >= 2) {
                                logger << "\033[1;33mdelayed \033[36m" << p << " \033[37mto \033[36m" << res << "\033[m:";
                                for (int n : regions[p]) {
                                    logger << " \033[37m" << label_vertex(n) << "\033[m";
                                }
                                logger << std::endl;
                            } else {
                                logger << "\033[1;33mdelayed \033[36m" << p << " \033[37mto \033[36m" << res << "\033[m" << std::endl;
                            }
                        }
                        delayed++;
                        if (p&1) del1++;
                        else del0++;
                        // emulate delayed promotion in region_
                        for (int i : regions[p]) region_[i] = res;
                        // skip to next priority and break inner loop
                        while (i >= 0 and priority(i) == p) i--;
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
                        for (int i=0; i<nodecount(); i++) if (region[i] != -2 && region_[i] <= res) region_[i] = -1;
                        // continue loop with the promotion target
                        i = inverse[res];
                        p = res;
                    }
                }
            }
        } else {
            // skip to next priority
            while (i >= 0 and priority(i) == p) i--;
        }
    }

    delete[] regions;
    delete[] region;
    delete[] region_;
    delete[] strategy;
    delete[] inverse;

    logger << "solved with " << promotions << " promotions, " << performances << "x performing delayed promotions (delayed " << delayed << ", discarded " << discarded << ", total " << promotions+discarded << ")" << std::endl;
}

}
