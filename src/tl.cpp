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

#include "tl.hpp"

namespace pg {

static const int BOT  = 0x80000000; // bottom state for a vertex
static const int DIS  = 0x80000001; // permanently disable vertex/tangle
static const int LEAK = 0x80000002; // when a vertex is not part of a tangle
static const int TANG = 0x80000003; // when a vertex is part of a tangle

#ifndef NDEBUG
#define LOG2(s) if (trace >= 2) { logger << s << std::endl; }
#else
#define LOG2(s) {}
#endif

TLSolver::TLSolver(Oink *oink, Game *game) : Solver(oink, game)
{
    alternating = false; // default
    onthefly = false; // default
}

TLSolver::~TLSolver()
{
}

inline void
TLSolver::attractTo(const int pr, const int pl, int cur)
{
    auto &R = regions[pr];
    // first attract normal vertices
    for (auto curedge = ins(cur); *curedge != -1; curedge++) {
        int from = *curedge;
        int r = region[from];
        if (r == DIS or r > pr) {
            // disabled or in a higher region
            continue; 
        } else if (r == pr) {
            // if already in <pr>, check if escape without strategy
            if (owner(from) == pl and strategy[from] == -1) {
#ifndef NDEBUG
                LOG2("\033[1;37mattracted \033[36m" << label_vertex(from) << "\033[m to \033[1;36m" << pr << "\033[m (via " << label_vertex(cur) << ")")
#endif
                strategy[from] = cur;
            }
        } else if (owner(from) == pl) {
            // it is same parity, lower priority, attract it
#ifndef NDEBUG
            LOG2("\033[1;37mattracted \033[36m" << label_vertex(from) << " \033[mto \033[1;36m" << pr << "\033[m (via " << label_vertex(cur) << ")");
#endif
            region[from] = pr;
            strategy[from] = cur;
            R.push_back(from);
            Q.push(from);
        } else {
            // it is other parity, lower priority, try to force it
            // interpret r as counter (using the negative numbers)
            if (r == BOT or r >= 0) {
                // if lower region without outcount, recompute outcount
                // count number of escapes, including our own priority
                r = 0;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] != DIS and region[to] <= pr) r--; // include our prio!!
                }
            }

            r += 1; // decrease counter (counter uses negative numbers)

            if (r != 0) {
                region[from] = r;
                continue;
            }

#ifndef NDEBUG
            LOG2("\033[1;37mforced \033[36m" << label_vertex(from) << " \033[mto \033[1;36m" << pr << "\033[m");
#endif
            region[from] = pr;
            strategy[from] = -1;
            R.push_back(from);
            Q.push(from);
        }
    }

    // then attract tangles
    auto &vin_cur = vin[cur];
    for (int from : vin_cur) {
        int r = vr[from];
        if (r == DIS or r >= pr) continue; // already attracted higher or same
        const int tangle_prio = vp[from];
        if (tangle_prio >= pr) {
            // higher priority tangle
            continue;
        } else if ((tangle_prio&1) != pl) {
            // i.e., owner(from) == pl
            // we "attract" it (so opponent doesn't attract it by mistake)
            vr[from] = pr;
            vregions[pr].push_back(from);
        } else {
            // i.e., owner(from) != pl, try to force it
            // interpret r as counter (using the negative numbers)
            if (r == BOT or r >= 0) {
                r = 0;
                int* ptr = vout[from], to;
                while ((to=*ptr++) != -1) {
                    if (region[to] != DIS and region[to] <= pr) r--; // include our prio!!
                }
            }

            r += 1; // decrease counter (counter uses negative numbers)

            if (r != 0) {
                vr[from] = r;
                continue;
            }

#ifndef NDEBUG
            int cnt = 0;
#endif
            vr[from] = pr;
            vregions[pr].push_back(from);
            int* ptr = vv[from], v, s;
            while ((v=*ptr++) != -1) {
                s = *ptr++;
                // do NOT overwrite... so only add if not yet in prio or higher
                if (region[v] == DIS or region[v] >= pr) continue; // already solved or in same/higher region
                region[v] = pr;
                strategy[v] = s;
                R.push_back(v);
                Q.push(v);
#ifndef NDEBUG
                cnt++;
#endif
            }

#ifndef NDEBUG
            if (cnt > 0 and trace >= 2) {
                logger << "\033[1;37mforced tangle \033[1;36m" << vp[from] << "\033[m to \033[1;36m" << pr << "\033[m" << std::endl;
            }
#endif
        }
    }
}

void
TLSolver::attract(const int pr)
{
    const int pl = pr&1;

    auto &R = regions[pr];
    const int count = R.size();
    for (int i=0; i<count; i++) attractTo(pr, pl, R[i]);

    while (Q.nonempty()) attractTo(pr, pl, Q.pop());
}

bool
TLSolver::computeRegion(int i)
{
    const int pr = priority(i);
    std::vector<int> &R = regions[pr];
    bool empty_region = true;

    for (int j=i; j>=0 and priority(j)==pr; j--) {
        if (region[j] != DIS and region[j] < pr) {
            region[j] = pr;
            strategy[j] = -1;
            empty_region = false;
            R.push_back(j);
        }
    }

    /**
     * If the region is empty, return False.
     * Otherwise, continue attractor computation.
     */

    if (empty_region) return false;
    else attract(pr);

    /**
     * Report result.
     */

#ifndef NDEBUG
    if (trace >= 3) {
        logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
        for (int n : R) {
            assert(region[n] == pr);
            logger << " " << label_vertex(n);
            if (strategy[n] != -1) logger << "->" << label_vertex(strategy[n]);
            if (strategy[n] != -1) assert(region[strategy[n]] == pr);
        }
        logger << std::endl;
    }
#endif

    return true;
}

int
TLSolver::extractTangles(int i, bool isHighest)
{
    const int pr = priority(i);
    const int pl = pr&1;

    /**
     * Analyze current region to extract tangles.
     * 1) Reduce region by removing escaping vertices (<pl> may not change strategy).
     * 2) Add each bottom SCC of result as new tangle.
     */

    /**
     * Mark escaping vertices with a special value LEAK.
     * We only need to look at the heads, then search backwards from there.
     */

    bool all_heads_leak = true;
    bool a_head_leaks = false;

    auto &R = regions[pr];
    for (auto j : R) {
        if (priority(j) != pr) {
            // only look at the heads
            break;
        } else if (owner(j) == pl) {
            if (strategy[j] == -1) goto current_vertex_leaks;
        } else {
            for (auto curedge = outs(j); *curedge != -1; curedge++) {
                int to = *curedge;
                if (region[to] != DIS and region[to] < pr) goto current_vertex_leaks;
            }
        }
        all_heads_leak = false;
        continue;

current_vertex_leaks:
        region[j] = LEAK;
        Q.push(j);
        a_head_leaks = true;
    }

    /**
     * For the variation where we do not extract tangles when ANY head leaks,
     * uncomment the check for a_head_leaks.
     */

    if (all_heads_leak/* or a_head_leaks*/) {
        // All heads LEAK, reset <region> and clear <Q>.
        while (Q.nonempty()) region[Q.pop()] = pr;
        return -2;
        (void)a_head_leaks; // suppress compiler warning
    }

    // Compute the greatest fixed point with a backward search
    while (Q.nonempty()) {
        unsigned int n = Q.pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (region[from] != pr) continue;
            if (strategy[from] != -1 and (uint) strategy[from] != n) continue;
            region[from] = LEAK;
            Q.push(from);
        }
    }

    if (isHighest) {
        // everything is a dominion
        bool non_empty = false;
        for (int v : R) {
            if (non_empty) {
                if (region[v] == pr) oink->solve(v, pl, strategy[v]);
            } else {
                if (priority(v) != pr) break;
                if (region[v] != pr) continue;
                oink->solve(v, pl, strategy[v]);
                non_empty = true;
            }
        }
        if (non_empty) {
            oink->flush();
            dominions++;
            if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m (highest region)" << std::endl;
            return -1;
        } else {
            // turns out there was no tangle
            for (int i : R) region[i] = pr;
            return -2;
        }
    }

    /**
     * Now all vertices with region[v]==pr are in the reduced region.
     * From every top vertex in the reduced region, we compute the bottom SCC with Tarjan.
     * We start at the top vertices because we know every bottom SCC contains a top vertex.
     */

    int highest_attracting = -1;  // the highest region that attracts one of the tangles
    bool has_tangle = false;      // whether we found a tangle
    bool has_dominion = false;    // whether we found a dominion

    std::vector<int> tangle;      // stores the tangle

    int pre = TANG;               // tarjan counting variable

    for (int j=i; j>=0 and priority(j)==pr; j--) {
        if (region[j] != pr) continue; // either DIS or LEAK or already visited in the search

        // start the tarjan search
        Q.push(j);
        while (Q.nonempty()) {
tarjan_again:
            const unsigned int n = Q.back();
            int min; // lowest successor vertex measure
            if (region[n] == pr) {
                // first time we see it
                // assign next <pre> to it and put it in <tarres>
                min = (region[n] = ++pre);
                tarres.push(n);
            } else {
                min = region[n];
            }

            /**
             * Perform the search step of Tarjan.
             */

            if (owner(n) != pl) {
                for (auto curedge = outs(n); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (region[to] == DIS or region[to] > pr) continue; // not in the tangle
                    // by definition, region[to] cannot be BOT or DIS
                    if (region[to] == pr) {
                        // not visited, add to search stack and go again
                        // push
                        Q.push(to);
                        goto tarjan_again;
                    } else {
                        // visited, update min
                        if (region[to] < min) min = region[to];
                    }
                }
            } else {
                const int s = strategy[n];
                assert(s != -1);
                if (region[s] == pr) {
                    // not visited, add to search stack and go again
                    Q.push(s);
                    goto tarjan_again;
                } else {
                    // visited, update min
                    if (region[s] < min) min = region[s];
                }
            }

            Q.pop();

            /**
             * If we're here, then we pushed no new edge.
             * Check if <n> is the root of an SCC.
             */

            if (min < region[n]) {
                // Not the root! Update measure and go up.
                region[n] = min;
                continue;
            }

            /**
             * Now <n> is the root of an SCC.
             * Every vertex after <n> in <res> is also in this SCC.
             * Also set every node to min.
             */

            for (;;) {
                const unsigned int m = tarres.pop();
                tangle.push_back(m);
                tangle.push_back(strategy[m]);
                region[m] = min; // to compute tangleto
                if (m == n) break;
            }

            /**
             * Now <tangle> contains the current SCC (+ strategy).
             * Every vertex in <tangle> has measure <min>.
             * It is only a real tangle if it contains a cycle, i.e., either 2+ vertices or a self-loop...
             */

            bool is_tangle = (tangle.size() > 2) or
                             ((unsigned int)strategy[n] == n) or
                             (strategy[n] == -1 and game->has_edge(n, n));
            if (!is_tangle) {
                tangle.clear();
                continue;
            }

            /**
             * We have a tangle. Compute the outgoing edges (into <tangleto>) and the next highest region.
             */

            int lowest_escape = 0x7fffffff;
            int esc = -1; // escape node (-1 for none, -2 for more)
            bool bottom_scc = true;
            const auto tangle_end = tangle.end();
            for (auto titer = tangle.begin(); titer != tangle_end;) {
                int v = *titer++;
                int s = *titer++;
                if (s != -1) continue; // not losing
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    int to = *curedge;
                    const int rto = region[to];
                    if (rto == DIS) {
                        // disabled
                        continue; 
                    } else if (bs_exits[to]) {
                        // marked
                        continue;
                    } else if (rto != min) { // not the current SCC, but might be in the same region
                        bs_exits[to] = true;
                        tangleto.push(to);
                        if (rto > pr and rto < lowest_escape) {
                            lowest_escape = rto;
                        } else if (rto < 0) {
                            bottom_scc = false;
                            break;
                        }
                        if (esc == -1) esc = v;
                        else esc = -2;
                    }
                }
                if (bottom_scc == false) break;
            }

            // Unmark exits
            for (unsigned int x=0; x<tangleto.size(); x++) bs_exits[tangleto[x]] = false;

            /**
             * If the tangle is not a bottom SCC, then we skip it.
             * It may have been found already and we don't want duplicates.
             * If it would be a unique tangle, then we'll find it next round!
             */

            if (!bottom_scc) {
                tangle.clear();
                tangleto.clear();
                continue;
            }

            /**
             * If there are no outgoing edges, then we have found a dominion.
             */

            if (tangleto.empty()) {
                // dominion
                if (trace) logger << "\033[1;38;5;201mdominion \033[36m" << pr << "\033[m";
                for (auto titer = tangle.begin(); titer != tangle_end;) {
                    int v = *titer++;
                    int s = *titer++;
                    // only if not already solved
                    if (disabled[v] == 0) oink->solve(v, pl, s);
                }
                if (trace) logger << std::endl;
                dominions++;
                has_dominion = true;
                tangle.clear();
                tangleto.clear();
                continue;
            }

            /**
             * We're not a dominion.
             */

#ifndef NDEBUG
            if (trace >= 2) {
                logger << "\033[1;38;5;198mnew tangle " << pr << "\033[m";
                if (trace >= 3) {
                    for (auto titer = tangle.begin(); titer != tangle_end;) {
                        int v = *titer++;
                        int s = *titer++;
                        logger << " \033[1;36m" << label_vertex(v) << "\033[m";
                        if (s != -1) logger << "->" << label_vertex(s);
                    }
                }

                logger << " with " << tangleto.size() << " exits to";
                for (unsigned int x = 0; x < tangleto.size(); x++) {
                    int r_to = region[tangleto[x]];
                    if (!bs_exits[r_to]) {
                        bs_exits[r_to] = true;
                        logger << " " << r_to;
                    }
                }
                bs_exits.reset();

                logger << std::endl;
            }
#endif
            esc = -1; // just in case
            // add back links to all normal vertices in our [out]
            const int tidx = vp.size();
            for (unsigned int x = 0; x < tangleto.size(); x++) vin[tangleto[x]].push_back(tidx);
            // move tangleto into vout
            int* _vout = new int[tangleto.size()+1];
            std::copy(&tangleto[0], &tangleto[tangleto.size()], _vout);
            _vout[tangleto.size()] = -1;
            vout.push_back(_vout);
            // move tangle into vv
            int* _vv = new int[tangle.size()+1];
            std::copy(tangle.begin(), tangle_end, _vv);
            _vv[tangle.size()] = -1;
            vv.push_back(_vv);
            // and set p to pr and current region to BOT
            vp.push_back(pr);
            vr.push_back(BOT);

            tangles++;
            has_tangle = true;

            /**
             * Finally, update <highest_attracting> and set <has_tangle>.
             */

            if (lowest_escape != 0x7fffffff and highest_attracting < lowest_escape) highest_attracting = lowest_escape;

            if (highest_attracting == lowest_escape) {
                // If we are the highest attracting, then just add to that region then...
                for (auto titer = tangle.begin(); titer != tangle_end;) {
                    int v = *titer++;
                    titer++;
                    region[v] = highest_attracting;
                    regions[highest_attracting].push_back(v);
                }
                if (esc < 0) vr.back() = pr;
            }

            tangleto.clear();
            tangle.clear();
        }
    }

    tarres.clear();

    if (has_dominion) {
        // just flush and return, no need to reset region[]
        oink->flush();
        return -1;
    } else if (has_tangle) {
        // reset region[] of all unattracted vertices
        std::vector<int> &R = regions[pr];
        int reset = onthefly ? BOT : pr;
        for (int i : R) if (region[i] != highest_attracting) region[i] = reset;
        return highest_attracting;
    } else {
        // reset region[] of all vertices to p
        std::vector<int> &R = regions[pr];
        for (int i : R) region[i] = pr;
        return has_tangle ? highest_attracting : -2;
    }
}

void
TLSolver::run()
{
    // get number of nodes and create and initialize inverse array
    int max_prio = game->priority(nodecount()-1);
    inverse = new int[max_prio+1];
    for (int i=0; i<=max_prio; i++) inverse[i] = DIS; // for skipped priorities
    for (int i=0; i<nodecount(); i++) if (!disabled[i]) inverse[priority(i)] = i;

    // allocate arrays
    region = new int[nodecount()];
    strategy = new int[nodecount()];
    vin = new std::vector<int>[nodecount()];
    regions = new std::vector<int>[max_prio+1];
    vregions = new std::vector<int>[max_prio+1];
    Q.resize(nodecount());
    tarres.resize(nodecount()+1);
    tangleto.resize(nodecount());
    bs_exits.resize(nodecount() > (max_prio+1) ? nodecount() : max_prio+1);

    // initialize arrays
    for (int i=0; i<nodecount(); i++) region[i] = disabled[i] ? DIS : BOT;
    for (int i=0; i<nodecount(); i++) strategy[i] = -1;


    tangles = 0;
    iterations = 0;
    turns = alternating ? 1 : 0;
    dominions = 0;

    /**
     * Find highest vertex...
     */

    int top = nodecount()-1;
    while (top >= 0 and region[top] == DIS) top--;
    if (top == -1) return; // empty game?

    bool foundNewTangles = false;        // found new tangles in this iteration (unused with on-the-fly)
    bool skipComputeRegion = false;      // do not recompute next region (after on-the-fly attracting)
    bool skipUntilChange = false;        // after changing player, do not recompute partition until next update
    bool reportedFirstDominion = false;  // did we see the first dominion yet (for trace reporting)
    int player = 0;                      // current player (for alternating strategy)
    int highestAttracting = -1;          // unused with on-the-fly
    int highestEven = -1;                // highest even region
    int highestOdd = -1;                 // highest odd region

    if (trace >= 1) {
        if (alternating) logger << "\033[1;38;5;196mround\033[m \033[1;36m" << turns-1 << "\033[m" << std::endl;
        else logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations << "\033[m" << std::endl;
    }

    int pr = max_prio;

    if (trace >= 2 and alternating) player = 2; // Hah!

    while (true) {
        /**
         * Check if we reached the bottom of the partitioning.
         * - not on-the-fly, new tangles: reset partition, recompute with the new tangles
         * - alternating: keep partition and now extract for other player
         */
        if (pr < 0) {
            if (foundNewTangles) {
                /**
                 * Not on-the-fly, if we reach the bottom and found tangles, reset the partition.
                 */
#ifndef NDEBUG
                /**
                 * Report partition.
                 */
                if (trace >= 3) {
                    for (int pr=max_prio; pr>=0; pr--) {
                        if (regions[pr].empty()) continue;
                        logger << "\033[1;33mregion\033[m \033[1;36m" << pr << "\033[m";
                        std::vector<int> &R = regions[pr];
                        for (int n : R) {
                            logger << " " << label_vertex(n);
                            if (strategy[n] != -1) logger << "->" << label_vertex(strategy[n]);
                        }
                        logger << std::endl;
                    }
                }
#endif
                pr = highestAttracting;
                // reset search below <pr>
                for (int j=inverse[pr]; j>=0; j--) {
                    if (region[j] < pr and region[j] != BOT and region[j] != DIS) region[j] = BOT;
                }
                for (unsigned int j=0; j<vr.size(); j++) {
                    if (vr[j] < pr and vr[j] != BOT and vr[j] != DIS) vr[j] = BOT;
                }
                for (int p=pr-1; p>=0; p--) {
                    regions[p].clear(); // clear vector containing vertices
                    vregions[p].clear(); // clear vector containing tangles
                }
                // continue at <pr>
                attract(pr);
                skipComputeRegion = true; // do not recompute region
                iterations++;
                if (trace >= 1) {
                    logger << "\033[1;38;5;196miteration\033[m \033[1;36m" << iterations<< "\033[m" << std::endl;
                }
                foundNewTangles = false; // reset flag
                skipUntilChange = false; // reset flag
                highestAttracting = -1; // reset value
                if (highestEven < pr) highestEven = -1;
                if (highestOdd < pr) highestOdd = -1;
            } else if (alternating) {
                /**
                 * For alternating tangle learning, change player and go again.
                 */
                if (player == 2) {
                    player = 0;
                } else {
                    player = 1 - player;
                    turns++;
                }
                if (trace >= 1) {
                    logger << "\033[1;38;5;196mround\033[m \033[1;36m" << turns-1 << "\033[m" << std::endl;
                }
                pr = max_prio; // start again from <max_prio>
                skipUntilChange = true; // no need to recompute regions
            } else {
                LOGIC_ERROR;
            }
        }

        /**
         * Check if there actually are vertices with this priority.
         */

        if (inverse[pr] == DIS) {
            pr--;
            continue;
        }

        /**
         * Compute the next region in the partition. (Unless skipping for reasons.)
         */

        if (skipComputeRegion) {
            /* do nothing, we just attracted a tangle */
            skipComputeRegion = false;
        } else if (skipUntilChange) {
            /* do nothing, we just changed the player */
        } else {
            computeRegion(inverse[pr]);
        }

        /**
         * After computing region <pr>, this region may be empty if all <pr> vertices are already in higher regions.
         */

        if (regions[pr].empty()) {
            pr--;
            continue;
        }

        bool isHighest = false;
        if (pr&1) {
            if (highestOdd < pr) {
                highestOdd = pr;
                isHighest = true;
            }
        } else {
            if (highestEven < pr) {
                highestEven = pr;
                isHighest = true;
            }
        }

        /**
         * If we are alternating and this is not our turn, then do not extract tangles.
         */

        if (alternating and (pr&1) != player) {
            pr--;
            continue;
        }

        /**
         * Proceed to extract tangles!
         */

        int res = extractTangles(inverse[pr], isHighest);

        /**
         * No new tangles found, go to next priority.
         */

        if (res == -2) {
            pr--;
            continue;
        }

        /**
         * Check if we have no dominion.
         * In that case, we perform on-the-fly attraction of the new tangle and reset the search at the attracting region.
         */

        if (res != -1) {
            if (trace) {
                logger << "\033[1;38;5;198mtangle\033[m \033[1;36m" << pr << "\033[m => \033[1;36m" << res << "\033[m" << std::endl;
            }

            /**
             * If we are on-the-fly attracting, simply reset the search below <res> and continue at <res>.
             */
            if (onthefly) {
                pr = res;
                // reset search below <pr>
                for (int j=inverse[pr]; j>=0; j--) {
                    if (region[j] < pr and region[j] != BOT and region[j] != DIS) region[j] = BOT;
                }
                for (unsigned int j=0; j<vr.size(); j++) {
                    if (vr[j] < pr and vr[j] != BOT and vr[j] != DIS) vr[j] = BOT;
                }
                for (int p=pr-1; p>=0; p--) {
                    regions[p].clear(); // clear vector containing vertices
                    vregions[p].clear(); // clear vector containing tangles
                }
                // continue at <pr>
                attract(pr);
                skipComputeRegion = true; // do not recompute region
                skipUntilChange = false;
                iterations++;
                continue;
            } else {
                if (res > highestAttracting) highestAttracting = res;
                foundNewTangles = true;
                pr--;
                continue;
            }
        }

        /**
         * If we're here, a dominion was found.
         * Reset the partition and prune tangles.
         */

        if (trace and !reportedFirstDominion) {
            // if this is the first dominion we find
            reportedFirstDominion = true;
            if (alternating) {
                logger << "first dominion in " << tangles << " tangles and " << turns << " turns." << std::endl;
            } else {
                logger << "first dominion in " << tangles << " tangles and " << iterations << " iterations." << std::endl;
            }
        }

        // find new top and reset regions
        int new_top = top;
        while (new_top >= 0 and disabled[new_top]) {
            region[new_top] = DIS;
            new_top--;
        }
        for (int n=new_top; n>=0; n--) {
            if (disabled[n]) {
                region[n] = DIS;
            } else {
                region[n] = BOT;
            }
        }

        if (new_top == -1) {
            // all vertices are now solved, end algorithm
            delete[] region;
            delete[] strategy;
            for (auto &x : vv) delete[] x;
            for (auto &x : vout) delete[] x;
            delete[] vin;
            delete[] regions;
            delete[] vregions;
            delete[] inverse;
            if (alternating) {
                logger << "found " << dominions << " dominions and " << tangles << " tangles." << std::endl;
                logger << "solved in " << turns << " turns." << std::endl;
            } else {
                logger << "found " << dominions << " dominions and " << tangles << " tangles." << std::endl;
                logger << "solved in " << iterations+1 << " iterations." << std::endl;
            }
            return;
        }

        top = new_top;

        // fix inverse
        for (int p=0; p<=max_prio; p++) {
            int n = inverse[p];
            if (n == DIS) continue; // already DIS
            while (true) {
                if (region[n] != DIS) {
                    break; // good!
                } else if (n == 0 or priority(n-1) != p) {
                    n = DIS; // bad!
                    break;
                } else {
                    n--; // try next!
                }
            }
            inverse[p] = n;
        }

        // deactivate tangles that are now permanently bad
        for (unsigned int i=0; i<vp.size(); i++) {
            if (vr[i] == DIS) continue; // already bad
            vr[i] = BOT; // reset to BOT but check if maybe now DIS
            const int vpl = vp[i]&1;
            int* ptr = vout[i], to;
            while ((to=*ptr++) != -1) {
                if (/*region[to] == DIS and*/ game->solved[to] and game->winner[to] != vpl) {
                    vr[i] = DIS;
                    break;
                }
            }
        }

        // clear region vectors
        for (int p=max_prio; p>=0; p--) {
            regions[p].clear(); // clear vector containing vertices
            vregions[p].clear(); // clear vector containing tangles
        }

        // restart partitioning at <max_prio>
        pr = max_prio;
        highestAttracting = -1;
        skipUntilChange = false;
        foundNewTangles = false;
        highestOdd = -1;
        highestEven = -1;
    }
}

}
