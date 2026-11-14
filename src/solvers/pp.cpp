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

PPSolver::PPSolver(Oink& oink, Game& game)
    : Solver(oink, game), Z(game.nodecount()), escape_witness(game.nodecount(), -1),
      top_vertex(game.nodecount()-1)
{
    if (game.nodecount() > 0) entries.resize(game.priority(game.nodecount()-1)+1);
}

/**
 * Returns true if no vertex with a priority higher than <p> is left, so there are no regions
 * above <p>. Vertices are only ever disabled, so <top_vertex> only moves down.
 */
bool
PPSolver::isTopmost(int p)
{
    while (top_vertex >= 0 and disabled[top_vertex]) top_vertex--;
    return top_vertex >= 0 and priority(top_vertex) == p;
}

PPSolver::~PPSolver()
{
}

void
PPSolver::attract(int prio, size_t first)
{
    const int pl = prio & 1;
    auto &rv = regions[prio];

    // The region members also form the FIFO. After a promotion, start at
    // the newly appended members because the earlier ones were processed.
    for (size_t i = first; i < rv.size(); ++i) {
        const int cur = rv[i];

        // (this iteration is a hot spot...)
        for (auto curedge = ins(cur); *curedge != -1; curedge++) {
            int from = *curedge;
            if (disabled[from] or region[from] >= prio) {
                // if not in the subgame of <prio>, skip, but remember predecessors in higher
                // regions: after a promotion, only they can be attracted
                if (!disabled[from] and region[from] > prio) entries[prio].emplace_back(from, cur);
                continue;
            } else if (owner(from) == pl) {
                // if owned by same parity, attract to a-maximal region
                rv.push_back(from);
                region[from] = prio;
                strategy[from] = cur;
#ifndef NDEBUG
                if (trace >= 3) logger << "\033[1;37mattracted \033[36m" << label_vertex(from) << " \033[37mto \033[36m" << prio << "\033[m (via " << label_vertex(cur) << ")" << std::endl;
#endif
            } else {
                // Regions and disabled vertices can change between attractions.
                const int witness = escape_witness[from];
                if (witness != -1 && !disabled[witness] && region[witness] < prio) continue;
                // if owned by other parity, check all outgoing edges
                bool can_escape = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    // if there is an escape, then it is not attracted
                    if (!disabled[to] and region[to] < prio) {
                        escape_witness[from] = to;
                        can_escape = true;
                        break;
                    }
                }
                if (can_escape) continue;
                // if all outgoing edges are prio >= p, attract to a-maximal region
                rv.push_back(from);
                region[from] = prio;
                strategy[from] = -1;
#ifndef NDEBUG
                if (trace >= 3) logger << "\033[1;37mforced \033[36m" << label_vertex(from)<< " \033[37mto \033[36m" << prio << "\033[m" << std::endl;
#endif
            }
        }
    }
    // Only original heads can lack a strategy. The priority index avoids
    // revisiting the old target members after a promotion.
    for (int v = inverse[prio]; v >= 0 && priority(v) == prio; --v) {
        if (disabled[v] || region[v] != prio || owner(v) != pl || strategy[v] != -1) continue;
        for (auto edge = outs(v); *edge != -1; ++edge) {
            if (!disabled[*edge] && region[*edge] == prio) {
                strategy[v] = *edge;
                break;
            }
        }
    }
#ifndef NDEBUG
    for (const int v : rv) {
        if (!disabled[v] && region[v] == prio && owner(v) == pl && priority(v) != prio) {
            assert(strategy[v] != -1);
        }
    }
#endif
}

/**
 * Perform promotion of region <from> to region <to>, afterwards maximize the region via attraction
 */
void
PPSolver::promote(int from, int to)
{
    assert(from < to);

#ifndef NDEBUG
    if (trace >= 2) {
        logger << "\033[1;33mpromoted \033[36m" << from << " \033[37mto \033[36m" << to << "\033[m:";
        for (int n : regions[from]) {
            logger << " \033[37m" << label_vertex(n) << "\033[m";
        }
        logger << std::endl;
    } else
#endif
    if (trace) {
        logger << "\033[1;33mpromoted \033[36m" << from << " \033[37mto \033[36m" << to << "\033[m" << std::endl;
    }

    // Promote and append the new attraction seeds to the target region.
    for (int i : regions[from]) region[i] = to;

    regions[to].insert(regions[to].end(), regions[from].begin(), regions[from].end());
    regions[from].clear();

    // The region <from> was closed under attraction, so only its predecessors in the regions
    // between <from> and <to> can be attracted now: those with a lower region still have their
    // escape below <from>. Attract them first, then continue from the newly attracted vertices.
    const int pl = to&1;
    auto &rv = regions[to];
    const size_t first = rv.size();
    for (const auto &entry : entries[from]) {
        const int x = entry.first;
        if (disabled[x]) continue;
        const int r = region[x];
        if (r == to) continue;
        if (r > to) {
            entries[to].push_back(entry); // still in a higher region
            continue;
        }
        if (owner(x) == pl) {
            rv.push_back(x);
            region[x] = to;
            strategy[x] = entry.second;
            continue;
        }
        const int witness = escape_witness[x];
        if (witness != -1 && !disabled[witness] && region[witness] < to) continue;
        bool can_escape = false;
        for (auto curedge = outs(x); *curedge != -1; curedge++) {
            const int y = *curedge;
            if (!disabled[y] and region[y] < to) {
                escape_witness[x] = y;
                can_escape = true;
                break;
            }
        }
        if (can_escape) continue;
        rv.push_back(x);
        region[x] = to;
        strategy[x] = -1;
    }
    entries[from].clear();
    attract(to, first);

    promotions++;
}

/**
 * Reset region <p> by resetting all vertices currently in region <p> to their original priority
 */
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
    entries[p].clear();
}

/**
 * Setup a region, i.e., take the existing region, check if it needs to be reset.
 * Then ensure that all top vertices are in the region...
 */
bool
PPSolver::setupRegion(int i, int p, bool mustReset)
{
    // NOTE: a node can be "disabled" and in a region (after dominion attraction)

    entries[p].clear(); // the attraction below records them again for all members

    if (!mustReset) {
        // Check if we must reset anyway (because region is tainted)
        // When must a region be reset?
        // - some vertices are in a dominion (now disabled)
        // - some vertices are already reset (i.e., attracted by opponent, then reset)
        // - some vertices are attracted to a higher OPPONENT region
        // When does a region not need to be reset?
        // - some vertices are attracted to a higher SAME-PARITY region
        for (int i : regions[p]) {
            if (region[i] == p) {
                Z[i] = true;
            } else if (disabled[i] or region[i] < p or (region[i] > p and (region[i]&1) != (p&1))) {
                mustReset = true;
                break;
            }
        }
    }

    if (mustReset) {
        if (!regions[p].empty()) resetRegion(p);
        Z.reset();
    } else {
        // No reset, but remove escapes (to be added) and remove things in a higher region...
        // The reason is that there may be vertices of priority p that are not yet in the region vector
        regions[p].erase(std::remove_if(regions[p].begin(), regions[p].end(),
                [&](const int x){ return Z[x] != 1; }),
                regions[p].end());
    }

    // Add all escapes to regions vector
    for (int j=i; j>=0 && priority(j) == p; j--) {
        if (region[j] == -2) continue; // already known to be disabled
        else if (disabled[j]) region[j] = -2; // set now to region -2
        else if (region[j] == p) { // is [to be] in region!
            if (Z[j]) {
                // already in the region, but check strategy
                if (strategy[j] != -1 and !Z[strategy[j]]) strategy[j] = -1;
                continue;
            } else {
                regions[p].push_back(j);
                strategy[j] = -1;
            }
        }
    }

    Z.reset();

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
    // Without higher regions, the region was attracted in the whole remaining game, so the
    // flush would attract nothing. The solved vertices stay in the queue of the next flush.
    const bool topmost = isTopmost(p);
    if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << p << "\033[m";
    for (int i : regions[p]) {
        assert(region[i] == p);
        assert(owner(i) != pl or (strategy[i] != -1 and region[strategy[i]] == p));
#ifndef NDEBUG
        if (trace >= 2) logger << " " << i;
#endif
        Solver::solve(i, pl, owner(i) == pl ? strategy[i] : -1);
    }
    if (trace) logger << std::endl;
    if (!topmost) Solver::flush();
}

int
PPSolver::getRegionStatus(int i, int p)
{
    const int pl = p&1;
    int lowest = -1; // the lowest higher region an opponent vertex can escape to

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
                const int r = region[to];
                if (r < p) return -2; // open
                if (r > p && (r < lowest || lowest == -1)) lowest = r;
            }
        }
    }
    // closed in the subgame; without higher regions it is a dominion
    if (isTopmost(p)) return -1;
    // find lowest higher region for possible promotion (the heads were scanned above)
    for (int j : regions[p]) {
        if (owner(j) != pl and priority(j) != p) {
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
#ifndef NDEBUG
                if (trace >= 2) reportRegion(p);
#endif
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
