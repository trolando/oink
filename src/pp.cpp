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
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cassert>

#include "pp.hpp"

namespace pg {

PPSolver::PPSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

PPSolver::~PPSolver()
{
}

void
PPSolver::attract(int prio, std::queue<int> queue)
{
    const int pl = prio & 1;
    auto &rv = regions[prio];

    // if queue is empty, then add all nodes of priority <p> to the queue
    if (queue.empty()) for (int i : rv) queue.push(i);

    // for every node in the queue, try to attract its predecessors
    while (!queue.empty()) {
        // for every node in the queue...
        int cur = queue.front();
        queue.pop();

        // (this iteration is a hot spot...)
        for (auto curedge = ins(cur); *curedge != -1; curedge++) {
            int from = *curedge;
            if (disabled[from] or region[from] > prio) {
                // if not in the subgame of <prio>, skip
                continue;
            } else if (region[from] == prio) {
                // if already in <prio>, check if escape without strategy
                if (owner(from) == pl and strategy[from] == -1) strategy[from] = cur;
            } else if (owner(from) == pl) {
                // if owned by same parity, attract to a-maximal region
                rv.push_back(from);
                region[from] = prio;
                strategy[from] = cur;
                queue.push(from);
                if (trace >= 3) logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << " \033[37mto \033[36m" << prio << "\033[m (via " << label_vertex(cur) << ")" << std::endl;
            } else {
                // if owned by other parity, check all outgoing edges
                bool can_escape = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    // if there is an escape, then it is not attracted
                    if (!disabled[to] and region[to] < prio) {
                        can_escape = true;
                        break;
                    }
                }
                if (can_escape) continue;
                // if all outgoing edges are prio >= p, attract to a-maximal region
                rv.push_back(from);
                region[from] = prio;
                strategy[from] = -1;
                queue.push(from);
                if (trace >= 3) logger << "\033[1;37mforced \033[36m" << label_vertex(from)<< " \033[37mto \033[36m" << prio << "\033[m" << std::endl;
            }
        }
    }
}

void
PPSolver::promote(int from, int to)
{
    assert(from < to);

    if (trace) {
        if (trace >= 2) {
            logger << "\033[1;33mpromoted \033[36m" << from << " \033[37mto \033[36m" << to << "\033[m:";
            for (int n : regions[from]) {
                logger << " \033[37m" << label_vertex(n) << "\033[m";
            }
            logger << std::endl;
        } else {
            logger << "\033[1;33mpromoted \033[36m" << from << " \033[37mto \033[36m" << to << "\033[m" << std::endl;
        }
    }

    // promote all nodes of region <from> to region <to>
    std::queue<int> queue;
    for (int i : regions[from]) {
        region[i] = to;
        queue.push(i);
    }

    regions[to].insert(regions[to].end(), regions[from].begin(), regions[from].end());
    regions[from].clear();

    // attract from the newly promoted nodes
    attract(to, queue);

    promotions++;
}

void
PPSolver::resetRegion(int p)
{
    for (int j : regions[p]) {
        // check if now disabled
        if (disabled[j]) region[j] = -2;
        // reset if currently in region p
        else if (region[j] == p) {
            region[j] = priority(j);
            strategy[j] = -1;
        }
    }
    regions[p].clear();
}

bool
PPSolver::setupRegion(int i, int p, bool mustReset)
{
    // NOTE: a node can be "disabled" and in a region (after dominion attraction)
    
    if (!mustReset) {
        // Check if we must reset anyway (because region is tainted)
        for (int i : regions[p]) {
            if (disabled[i] or region[i] != p) {
                mustReset = true;
                break;
            }
        }
    }

    if (mustReset) {
        if (!regions[p].empty()) resetRegion(p);
    } else {
        // No reset, but remove escapes (to be added)
        regions[p].erase(std::remove_if(regions[p].begin(), regions[p].end(),
                [&](const int x){ return priority(x) == p; }),
                regions[p].end());
    }

    // Add all escapes to regions vector
    for (int j=i; j>=0 && priority(j) == p; j--) {
        if (region[j] == -2) continue;
        else if (disabled[j]) region[j] = -2;
        else if (region[j] == p) {
            if (strategy[j] != -1 and (disabled[strategy[j]] or region[strategy[j]] != p)) strategy[j] = -1;
            regions[p].push_back(j);
        }
    }

    // If the region is empty, return false
    if (regions[p].empty()) return false;

    // Otherwise, attract and return true
    attract(p);
    return true;
}

void
PPSolver::setDominion(int p)
{
    // found a dominion
    const int pl = p&1;
    if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << p << "\033[m";
    for (int i : regions[p]) {
        assert(region[i] == p);
        assert(owner(i) != pl or (strategy[i] != -1 and region[strategy[i]] == p));
        if (trace >= 2) logger << " " << i;
        oink->solve(i, pl, owner(i) == pl ? strategy[i] : -1);
    }
    if (trace) logger << std::endl;
    oink->flush();

    // to recover, all regions containing a disabled node must be reset...
    /*
    for (int i=0; i<nodecount(); i++) {
        if (region[i] == -2) continue;
        if (!disabled[i]) continue;
        int r = region[i];
        logger << "tainted " << r << std::endl;
    }
    */
}

int
PPSolver::getRegionStatus(int i, int p)
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
                int to =*curedge;
                if (disabled[to]) continue;
                if (region[to] < p) return -2; // open
            }
        }
    }
    // closed in the subgame, so find lowest higher region for possible promotion
    int lowest = -1;
    for (int j : regions[p]) {
        if (owner(j) != pl) {
            // losing node, find lowest higher ...
            for (auto curedge = outs(j); *curedge != -1; curedge++) {
                int to = *curedge;
                /** HOT SPOT of the iterator and minor hot spot of obtaining the region **/
                if (disabled[to]) continue;
                int r = region[to];
                if (r > p && (r < lowest || lowest == -1)) lowest = r;
                else if (r < p) LOGIC_ERROR;
            }
        }
    }
    return lowest; // -1 if dominion, <region> otherwise
}

void
PPSolver::reportRegion(int p)
{
    const int pl = p&1;
    logger << "\033[1;33mregion \033[36m" << p << "\033[m";
    for (int n : regions[p]) {
        if (region[n] == p) logger << " \033[37m" << label_vertex(n) << "\033[m";
        if (owner(n) == pl) {
            if (strategy[n] == -1) {
                if (priority(n) != p) logger << "\033[31;1m--\033[m";
            } else {
                if (disabled[strategy[n]] or region[strategy[n]] != p) logger << "->\033[31;1m" << label_vertex(strategy[n]) << "\033[m";
                else logger << "->" << label_vertex(strategy[n]);
            }
        } else {
            bool got = false;
            for (auto curedge = outs(n); *curedge != -1; curedge++) {
                int to = *curedge;
                if (disabled[to]) continue;
                if (region[to] == -2) continue;
                if (region[to] == p) continue;
                if (!got) logger << "(";
                else logger << ",";
                got = true;
                if (region[to] < p) {
                    if (priority(n) != p) logger << "\033[31;1mesc\033[m";
                    else logger << "\033[36m" << region[to] << "\033[m";
                } else logger << "\033[36m" << region[to] << "\033[m";
            }
            if (got) logger << ")";
        }
    }
    logger << std::endl;
}

void
PPSolver::run()
{
    // obtain highest priority and allocate arrays
    max_prio = priority(nodecount()-1);
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

        // PP: always reset region
        if (setupRegion(i, p, true)) {
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
                    break;
                } else {
                    // found promotion, promote
                    promote(p, res);
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

void
PPSolver::printState()
{
    for (int i = nodecount()-1; i>=0; i--) {
        if (region[i] == -2) continue;
        int p = priority(i);
        if (region[i] != p) continue;

        logger << "\033[1;34m|| \033[1;37m" << p << "\033[m (";

        std::vector<int> exits;
        for (int j : regions[p]) {
            logger << "\033[1m" << priority(j) << "\033[m ";
            if (owner(j) == (p&1)) {
                bool escapes = true;
                for (auto curedge = outs(j); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == -2) continue;
                    if (region[to] == p) {
                        escapes = false;
                        break;
                    }
                }
                if (escapes) {
                    for (auto curedge = outs(j); *curedge != -1; curedge++) {
                        int to = *curedge;
                        if (region[to] == -2) continue;
                        logger << "\033[m" << region[to] << "\033[m ";
                        if (region[to] != p) exits.push_back(region[to]);
                    }
                }
            } else {
                for (auto curedge = outs(j); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == -2) continue;
                    logger << "\033[m" << region[to] << "\033[m ";
                    if (region[to] != p) exits.push_back(region[to]);
                }
            }
        }

        if (!regions[p].empty()) logger << "\b";
        logger << ") => {";

        std::sort(exits.begin(), exits.end());
        for (unsigned k=0; k<exits.size(); k++) {
            if (k > 0) {
                if (exits[k-1] == exits[k]) continue;
                logger << ", ";
            }
            logger << exits[k];
        }
        logger << "}";

        // check if the region is closed
        bool fullclosed = true;
        bool subclosed = true;
        bool empty = true;
        for (int j : regions[p]) {
            empty = false;
            if (owner(j) == p%2) {
                // node won by escape owner
                // it is therefore closed if it has an edge to itself
                bool nodeclosed = false;
                for (auto curedge = outs(j); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == p) nodeclosed = true;
                    if (nodeclosed) break;
                }
                if (!nodeclosed) fullclosed = subclosed = false;
            } else {
                // region lost by escape owner
                // it is therefore closed if there are no edges to lower regions
                for (auto curedge = outs(j); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == -2) continue;
                    if (region[to] != p) {
                        fullclosed = false;
                        if (region[to] < p) subclosed = false;
                    }
                    if (!subclosed) break;
                }
            }
            if (!subclosed) break;
        }

        if (empty) {} 
        else if (fullclosed) logger << " (dominion)";
        else if (subclosed) logger << " (closed to " << exits[0] << ")";
        logger << "\033[m";
        logger << std::endl;
    }
}


}
