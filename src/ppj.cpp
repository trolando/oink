/*
 * Copyright 2017-2020 Tom van Dijk, Johannes Kepler University Linz
 * Copyright 2019-2020 Ruben Lapauw, KU Leuven
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
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "ppj.hpp"

namespace pg {

struct escape {
    int region;
    int play;
    int to;
};

int pm_cmp(int pl, int p_left, int p_right) {
    p_left += pl;
    p_right += pl;
    return p_left == p_right ? 0 :
           (p_right & 1) != (p_left & 1) ? (p_right & 1) - (p_left & 1) :
           ((p_left & 1) ^ (p_left < p_right)) ? -1 : 1;
}

PPJSolver::PPJSolver(Oink *oink, Game *game) : PPSolver(oink, game), waitingPriority(0), waiting(NULL)
{
}

void PPJSolver::unattracted(int node) {
    if (is_lost[node]) return;
    if(trace >= 3) logger << "Node may be lost " << node << std::endl;
    is_lost[node] = true;
    lost.push_back(node);
}
escape PPJSolver::bestEscape(int node) {
    int pl = owner[node];
    int r_node = region[node];
    //A higher priority than max_prio of opposite parity of pl;
    int highest_escape = max_prio + 2 - ((1 + pl + max_prio) & 1);
    int play = -1;
    const int *_out = outs + outa[node];
    for (int to = *_out; to != -1; to = *++_out) {
        int r_to = region[to];
        if (r_to == -2) continue;
        //logger << node << " " << r_node << " " << to << " @" << r_to << std::endl;
        if(pm_cmp(pl, highest_escape, r_to) < 0) {
            highest_escape = r_to;
            play = to;
        }
    }
    int strategy = (highest_escape & 1) != pl ? -1 : play;
    //strategy = (highest_escape < priority[node]) ? -1 : strategy;
    escape out = escape{highest_escape, strategy, play};
    return out;
}
void PPJSolver::endAttracting(int prio) {
    int pl = prio & 1;
    std::vector<int> killQueue;
    // if (trace >= 3) logger << "Cleanup after attracting " << prio << std::endl;
    for(auto it = lost.begin(); it != lost.end(); it++) {
        int node = *it;
        is_lost[node] = false;
        if(region[node] == prio) continue;         // has been attracted
        if(owner[node] == pl) LOGIC_ERROR;
        if(strategy[node] != -1 && (region[strategy[node]] & 1) != pl) {
            // the node is attracted to an agreeing edge that was not flipped.
            // the attraction stays valid. Thus ignore.
            // if (trace >= 3) logger << "Not in justification " << node << std::endl;
            continue;
        }
        escape highest_escape = bestEscape(node);
        int r_escape = highest_escape.region;
        if (((r_escape & 1) == pl) != (highest_escape.play == -1)) LOGIC_ERROR;
        if (r_escape > prio) {
            // At least one non-disabled edge should exist.
            // At least one lower than prio should exist otherwise should have been attracted.
            logger << "Node " << node << " should have been attracted to " << highest_escape.to << std::endl;
            LOGIC_ERROR;
        } else if (r_escape >= region[node]) {
            if ((r_escape & 1) != pl) {
                // if (trace >= 2) logger << "Node " << node << " with owner " << owner[node] << " may switch play to " << highest_escape.to << " ... rebuild" << std::endl;
                killQueue.push_back(node);
            } else if (r_escape == region[node]) {
                // if (trace >= 3) logger << "Unmovable " << node << " to " << r_escape << std::endl;
            } else {
                // all nodes are disagreeing with owner, the node is forced.
                // but it was not attracted to prio. Thus waiting.
                // if(trace >= 2) logger << "Delayed force " << node << " to @" << r_escape << std::endl;
                setWaiting(node, r_escape);
            }
        } else if (r_escape > priority[node]) {
            // The region can flip from attracted to a lower forced node.
            if((region[node] & 1) == pl) {
                // The region cannot both drop and be attracted while edge-nodes increase in priority (of the same parity).
                logger << "Drop halfway " << node << " to @" << r_escape << " is " << region[node] << " over " << priority[node] << std::endl;
                LOGIC_ERROR;
            }
            // Kill dependents and reset waiting.
            // TODO? don't recalculate highest_escape when rebuilding but only kill dependents of node?
            // if (trace >= 2) logger << "Flip-halfway " << node << " to " << r_escape << std::endl;
            killQueue.push_back(node);
        } else {         //if (r_escape <= priority[node]) <= region[node]
            if (region[node] == priority[node]) {
                if (trace >= 3) logger << "Don't drop heads " << node << " of " << priority[node] << " at " << region[node] << " to " << r_escape << std::endl;
            } else {
                // highest_escape < region[from]
                // The region is dropping. Don't try to catch falling nodes. This leads to bad cycles.
                // if (trace >= 2) logger << "Dropping attraction " << node << " to " << r_escape << std::endl;
                killQueue.push_back(node);
            }
        }
    }
    lost.clear();
    for(uint idx = 0; idx < killQueue.size(); idx++) {
        int node = killQueue[idx];
        // if (trace >= 3) logger << "Lost support " << node << std::endl;
        strategy[node] = -1;
        region[node] = priority[node];
        const int *_in = ins + ina[node];
        for (int from = *_in; from != -1; from = *++_in) {
            if(disabled[from]) continue;
            // if (trace >= 4) logger << "Lost? " << from << " by " << node << " " << region[from] << std::endl;
            if(strategy[from] >= 0 && strategy[from] != node) continue;             // not in justification graph;
            if(region[from] == priority[from]) {
                continue;                 // don't reset heads; don't reset twice.
            }
            if(region[from] < priority[from]) LOGIC_ERROR;
            killQueue.push_back(from);
        }
    }
    // cannot set a node waiting with a given strategy.
    // recalculate full strategy again instead of fixed edge?
    // The difference is:
    // Knowing the play but not the origin. Iterate over all incoming nodes.
    // Verify all outgoing nodes for a forced attraction as well.
    // Knowing the origin but not the play. Iterate over all outgoing nodes.
    for(auto node : killQueue) {
        escape highest_escape = bestEscape(node);
        int r_escape = highest_escape.region;
        if (r_escape > prio) {
            logger << "Node "<< node << " attracted/forced to " << highest_escape.to << " at " << r_escape << std::endl;
            // At least one non-disabled edge should exist.
            // At least one lower than prio should exist otherwise should have been attracted.
            LOGIC_ERROR;
        }
        // if (trace >= 2) logger << "Best rebuilt support " << node << " at " << r_escape << std::endl;
        setWaiting(node, r_escape);
    }
}

void PPJSolver::setWaiting(int node, int prio) {
    waitingPriority = std::max(waitingPriority, prio);
    if(prio > priority[node]) {
        waiting[prio].push_back(node);
    }
}

bool
PPJSolver::setupRegion(int i, int p, bool mustReset)
{
    // NOTE: a node can be "disabled" and in a region (after dominion attraction)
    if (!mustReset) LOGIC_ERROR;
    // Remove head (to be added again)
    // Remove lost nodes.
    regions[p].erase(std::remove_if(regions[p].begin(), regions[p].end(),
                                    [&](const int x){
            return priority[x] == p || region[x] != p;
        }),
                     regions[p].end());
    std::queue<int> q;
    // Add all escapes to regions vector
    for (int j=i; j>=0 && priority[j] == p; j--) {
        if (disabled[j]) region[j] = -2;
        else if (region[j] == p) {
            if (strategy[j] >= 0 and (disabled[strategy[j]] or region[strategy[j]] != p)) {
                // logger << "Plays outside " << p << " : " << j << " to " << strategy[j] << " at " << region[strategy[j]] << std::endl;
                strategy[j] = -1;
            }
            // if (trace >= 2) logger << "\033[1;37mhead \033[36m" << j << " \033[37mat \033[36m" << p << std::endl;
            regions[p].push_back(j);
            q.push(j);
        }
    }
    for(auto node : waiting[p]) {
        if (disabled[node]) continue;
        if (region[node] >= p) continue;         //has been attracted
        escape best_esc = bestEscape(node);
        if(best_esc.region == p) {
            region[node] = p;
            strategy[node] = best_esc.play;
            regions[p].push_back(node);
            q.push(node);
        }
    }
    waiting[p].clear();

    attract(p, q);
    return regions[p].size() > 0;
}

void
PPJSolver::run()
{
    // obtain highest priority and allocate arrays
    max_prio = priority[n_nodes-1];
    int dom_len = max_prio + 4;
    regions = new std::vector<int>[dom_len];
    region = new int[n_nodes];
    strategy = new int[n_nodes];
    inverse = new int[dom_len];
    waiting = new std::vector<int>[dom_len];
    is_lost = bitset(n_nodes);

    // initialize arrays
    for (int i=0; i<n_nodes; i++) region[i] = disabled[i] ? -2 : priority[i];
    for (int i=0; i<n_nodes; i++) strategy[i] = -1;
    for (int i=0; i<dom_len; i++) inverse[i] = -1;
    for (int i=0; i<n_nodes; i++) if (!disabled[i]) inverse[priority[i]] = i;


    // reset statistics
    promotions = 0;

    /**
     * Two loops: the outer (normal) loop, and the inner (promotion-chain) loop.
     * The outer loop does region setup and attractor on the full region.
     * The inner loop only attracts from the promoted region.
     */

    // start loop at the highest priority
    int p = max_prio;
    while (p >= 0) {
        // get current priority and skip all disabled/attracted nodes
        int i = inverse[p];
        // PPJ: reset region is unused
        if (setupRegion(i, p, true)) {
            // region not empty, maybe promote
            while (true) {
                if (trace >= 2) reportRegion(p);
                int res = getRegionStatus(i, p);

                if (res == -2) {
                    p--;
                    break;
                } else if (res == -1) {
                    // found dominion
                    int dominion = max_prio-(max_prio&1)+2+(p&1);
                    promote(p, dominion);
                    p = max_prio;
                    break;
                } else {
                    // found promotion, promote
                    promote(p, res);
                    p = res;
                    i = inverse[p];
                }
            }
        } else {
            p--;
        }
    }
    for(int i = 0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        int pl = region[i]& 1;
        if ((strategy[i] >= 0) != (owner[i] == pl)) {
            logger << "Strategy of " << i << " is " << strategy[i] << " when owned by " << owner[i] << " and won by " << pl << std::endl;
            LOGIC_ERROR;
        }
        oink->solve(i, pl, owner[i] == pl ? strategy[i] : -1);
    }
    oink->flush();
    delete[] regions;
    delete[] region;
    delete[] strategy;
    delete[] inverse;
    delete[] waiting;

    logger << "solved with " << promotions << " promotions." << std::endl;
}

}
