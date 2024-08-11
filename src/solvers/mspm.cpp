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

#include <iomanip>
#include "mspm.hpp"

static const bool use_queue = true;

namespace pg {

MSPMSolver::MSPMSolver(Oink& oink, Game& game) : Solver(oink, game)
{
}

MSPMSolver::~MSPMSolver()
{
}

/**
 * Returns true if a progress measure "a" is less than "b"
 * up to and including priority <d>, for player <pl>.
 */
bool
MSPMSolver::pm_less(int* a, int* b, int d, int pl)
{
    // cases where a or b is Top
    if (b[pl] == -1) return a[pl] != -1;
    else if (a[pl] == -1) return false;
    // normal comparison, start with highest priority
    const int start = ((k&1) == pl) ? k-2 : k-1;
    for (int i=start; i>=d; i-=2) {
        if (a[i] == b[i]) continue;
        return a[i] < b[i];
    }
    return false;
}

/**
 * Copy for player <pl>.
 */
void
MSPMSolver::pm_copy(int *dst, int *src, int pl)
{
    for (int i=pl; i<k; i+=2) dst[i] = src[i];
}

/**
 * Write pm to ostream.
 */
void
MSPMSolver::pm_stream(std::ostream &out, int *pm)
{
    bool tope = pm[0] == -1;
    bool topo = pm[1] == -1;
    out << " {";
    if (tope) out << " \033[1;33mTe\033[m";
    else out << " " << pm[0];
    if (topo) out << " \033[1;33mTo\033[m";
    else out << " " << pm[1];
    for (int i=2; i<k; i++) {
        if (i&1) out << " " << (topo ? 0 : pm[i]);
        else     out << " " << (tope ? 0 : pm[i]);
    }
    out << " } ";
}

/**
 * Perform update for player <pl>, node with priority <d>.
 */
void
MSPMSolver::Prog(int *dst, int *src, int d, int pl)
{
    // check if top
    if (src[pl] == -1) {
        dst[pl] = -1;
        return;
    }

    // set every value lower than <d> to 0.
    int i = pl;
    for (; i<d; i+=2) dst[i] = 0;

    int carry = (d == i) ? 1 : 0; // only increase if <d> has parity <pl>

    for (; i<k; i+=2) {
        // increase or copy for same parity
        int v = src[i] + carry;
        if (v > counts[i]) {
            dst[i] = 0;
            carry = 1;
        } else {
            dst[i] = v;
            carry = 0;
        }
    }

    // check if top
    if (carry) dst[pl] = -1;
}

bool
MSPMSolver::lift(int node, int target)
{
    // obtain ptr to current progress measure
    int *pm = pms + k*node;

    // check if already Top for both players
    if (pm[0] == -1 and pm[1] == -1) return false;

    lift_attempt++;

    // initialize stuff
    const int pl_max = owner(node);
    const int pl_min = 1 - pl_max;
    const int d = priority(node);

    if (trace >= 2) {
        logger << "\033[1mupdating node " << node << "/" << d << (owner(node)?" (odd)":" (even)") << "\033[m with current progress measure";
        pm_stream(logger, pm);
        logger << std::endl;
    }

    int best_ch0 = -1, best_ch1 = -1;

    // do max for player <pl_max>
    if (pm[pl_max] != -1) {
        if (trace >= 2) logger << "computing max" << std::endl;
        if (trace >= 2) pm_copy(tmp, pm, 1-pl_max);

        if (target != -1) {
            // just look at target
            Prog(tmp, pms + k*target, d, pl_max);

            if (trace >= 2) {
                logger << "successor node " << target << "/" << priority(target) << " results in";
                pm_stream(logger, tmp);
                logger << std::endl;
            }

            if (pm_less(pm, tmp, d, pl_max)) {
                pm_copy(pm, tmp, pl_max);
                if (pl_max) best_ch1 = target;
                else best_ch0 = target;
            }
        } else {
            for (auto curedge = outs(node); *curedge != -1; curedge++) {
                int to = *curedge;
                if (cover[to] == -2) continue;
                Prog(tmp, pms + k*to, d, pl_max);

                if (trace >= 2) {
                    logger << "successor node " << to << "/" << priority(to) << " results in";
                    pm_stream(logger, tmp);
                    logger << std::endl;
                }

                if (pm_less(pm, tmp, d, pl_max)) {
                    pm_copy(pm, tmp, pl_max);
                    if (pl_max) best_ch1 = to;
                    else best_ch0 = to;
                }
            }
        }
    }

    // do min for player <pl_min>
    if (pm[pl_min] != -1 and (target == -1 or target == strategy[node])) {
        if (trace >= 2) logger << "computing min" << std::endl;
        if (trace >= 2) pm_copy(tmp, pm, 1-pl_min);
        int best_to = -1;
        for (auto curedge = outs(node); *curedge != -1; curedge++) {
            int to = *curedge;
            if (cover[to] == -2) continue;
            Prog(tmp, pms + k*to, d, pl_min);

            if (trace >= 2) {
                logger << "successor node " << to << "/" << priority(to) << " results in";
                pm_stream(logger, tmp);
                logger << std::endl;
            }

            if (best_to == -1 or pm_less(tmp, best, d, pl_min)) {
                for (int i=0; i<k; i++) best[i] = tmp[i];
                best_to = to;
            }
        }

        strategy[node] = best_to;
        // note: sometimes only the strategy changes, but the lowest pm stays the same
        // now "best" contains the smallest Prog, which may be higher than the current min
        if (pm_less(pm, best, d, pl_min)) {
            pm_copy(pm, best, pl_min);
            if (pl_min) best_ch1 = best_to;
            else best_ch0 = best_to;
        }
    }

    bool ch0 = best_ch0 != -1;
    bool ch1 = best_ch1 != -1;

    if (ch0 or ch1) {
        if (trace) {
            logger << "\033[1;32mupdated node " << node << "/" << d << (owner(node)?" (odd)":" (even)") << "\033[m to";
            pm_stream(logger, pm);
            logger << std::endl;
        }
        // check if a changed pm is now Top
        if (ch0 and pm[0] == -1) solve(node, best_ch0);
        if (ch1 and pm[1] == -1) solve(node, best_ch1);
        // increase count and return true
        lift_count++;
        return true;
    } else {
        return false;
    }
}

/**
 * Upon lifting a node to Top, solve is called.
 */
void
MSPMSolver::solve(int node, int str)
{
    int *pm = pms + k*node;
    const int pl = pm[0] == -1 ? 0 : 1;
    if (pm[pl] != -1) LOGIC_ERROR;

    if (trace) logger << "Detected \033[1;31mTop\033[m from " << node << "/" << priority(node) << " to " << str << "/" << priority(str) << std::endl;

    // initialize
    std::queue<int> q;

    // solve <node>
    Solver::solve(node, pl, owner(node) == pl ? str : -1);
    // if (pl == (priority(node)&1)) counts[priority(node)]--;
    cover[node] = -1;
    q.push(node);
    todo_push(node);

    // attract within uncovered area to <node> for <pl>
    while (!q.empty()) {
        int n = q.front();
        q.pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            // logger << "trying edge " << from << "/" << priority(from) << " to " << n << std::endl;
            if (cover[from]) continue;
            if (priority(from) > priority(node)) continue;
            if (owner(from) != pl) {
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (cover[to] < 0) continue; // disabled or already solved
                    // if (cover[to] != 0) LOGIC_ERROR; // should have been covered
                    // ^--- it can escape, happens when alternating!
                    escapes = true;
                    break;
                }
                if (escapes) continue;
            }
            Solver::solve(from, pl, owner(from) == pl ? n : -1);
            // if (winner == (priority(from)&1)) counts[priority(from)]--;
            cover[from] = -1;
            pms[k*from+pl] = -1;
            q.push(from);
            todo_push(from);
        }
    }

    // after lifting, cover everything that is lower
    coverlower(node, ++coverdepth);

    if (trace >= 1) {
        logger << "Cover status:" << std::endl;
        for (int i=0; i<nodecount(); i++) {
            if (cover[i]) logger << i << "/" << priority(i) << ": " << cover[i] << std::endl;
        }
    }
}

void
MSPMSolver::coverlower(int node, int k)
{
    // initialize
    const int pr = priority(node);
    const int pl = 1-(pr&1); // attract for other player

    bool banner = false;

    std::queue<int> q;

    for (int n=node; n<nodecount(); n++) {
        if (cover[n]) continue; // also for "disabled"
        if (priority(n) <= pr) continue; // skip nodes of same priority
        cover[n] = k;
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (cover[from]) continue;
            if (owner(from) != pl) {
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (cover[to]) continue;
                    escapes = true;
                    break;
                }
                if (escapes) continue;
            }
            cover[from] = k;
            if (trace >= 2) {
                if (!banner) logger << "\033[7;31;1mcovering\033[m ";
                banner = true;
                logger << " " << from;
            }
            q.push(from);
        }
    }

    while (!q.empty()) {
        int n = q.front();
        q.pop();
        for (auto curedge = ins(n); *curedge != -1; curedge++) {
            int from = *curedge;
            if (cover[from]) continue;
            if (owner(from) != pl) {
                bool escapes = false;
                for (auto curedge = outs(from); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (cover[to]) continue;
                    escapes = true;
                    break;
                }
                if (escapes) continue;
            }
            cover[from] = k;
            if (trace >= 2) logger << " " << from;
            q.push(from);
        }
    }

    if (banner) logger << std::endl;
}

void
MSPMSolver::run()
{
    // determine k = highest priority + 1
    k = priority(nodecount()-1)+1;

    // now create the data structure, for each node
    pms = new int[k*nodecount()];
    strategy = new int[nodecount()];
    counts = new int[k];
    cover = new int[nodecount()];
    tmp = new int[k];
    best = new int[k];

    // initialize all measures to 0
    for (int i=0; i<k*nodecount(); i++) pms[i] = 0;

    // initialize strategy to -1
    for (int i=0; i<nodecount(); i++) strategy[i] = -1;

    // initialize cover
    for (int i=0; i<nodecount(); i++) cover[i] = disabled[i] ? -2 : 0;

    // initialize counts for each priority
    for (int i=0; i<k; i++) counts[i] = 0;
    for (int i=0; i<nodecount(); i++) if (cover[i] == 0) counts[priority(i)]++;

    // set number of lifts and lift attempts to 0
    lift_count = lift_attempt = 0;

    /**
     * Strategy that updates predecessors then marks updated predecessors for processing.
     * Uses a queue/stack to store the dirty vertices.
     */

    // initialize all nodes as not dirty
    dirty = new int[nodecount()];
    for (int n=0; n<nodecount(); n++) dirty[n] = 0;

    // initialize cover depth to 0
    coverdepth = 0;

    /**
     * First loop, starting with standard "lift" for nodes and "liftR" for predecessors.
     * Can already find "Top" here, so add covered nodes to <deferred>.
     */
    for (int n=nodecount()-1; n>=0; n--) {
        bool lifted = cover[n] == 0 and lift(n, -1);
        if (cover[n] == -1 or lifted) {
            for (auto curedge = ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                if (cover[from] == 0 and lift(from, n)) todo_push(from);
            }
        }
    }

    /**
     * The main loop
     */

    while (!todo.empty()) {
        int n = todo_pop();
        if (cover[n] == -1 or cover[n] == 0) {
            for (auto curedge = ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                if (cover[from] == 0 and lift(from, n)) todo_push(from);
            }
        }

        while (todo.empty() and coverdepth != 0) {
            int m = coverdepth--;
            std::queue<int> q;

            for (int n=0; n<nodecount(); n++) {
                if (cover[n] != m) continue;
                cover[n] = 0;

                if (trace) logger << "\033[7;31;1muncovering\033[m " << n << std::endl;

                // try to see if it is attracted
                const int pl = owner(n);
                bool escapes = false;
                for (auto curedge = outs(n); *curedge != -1; curedge++) {
                    int to = *curedge;
                    if (cover[to] == -2) continue;
                    if (game.isSolved(to)) {
                        if (game.getWinner(to) == pl) {
                            Solver::solve(n, pl, to);
                            cover[n] = -1;
                            pms[k*n+pl] = -1;
                            q.push(n);
                            todo_push(n);
                            break;
                        }
                    } else {
                        escapes = true;
                    }
                }
                if (cover[n] == -1) continue;
                if (!escapes) {
                    Solver::solve(n, 1-pl, -1);
                    cover[n] = -1;
                    pms[k*n+(1-pl)] = -1;
                    q.push(n);
                    todo_push(n);
                    continue;
                }

                // not attracted, lift
                lift(n, -1);

                // TODO: always or only if lifted?
                todo_push(n);
            }

            // attract within uncovered area to <node> for <pl>
            while (!q.empty()) {
                int n = q.front();
                q.pop();
                if (!game.isSolved(n)) LOGIC_ERROR;
                const bool pl = game.getWinner(n);
                for (auto curedge = ins(n); *curedge != -1; curedge++) {
                    int from = *curedge;
                    if (cover[from]) continue;
                    if (owner(from) != pl) {
                        bool escapes = false;
                        for (auto curedge = outs(from); *curedge != -1; curedge++) {
                            int to = *curedge;
                            if (cover[to] < 0) continue;
                            if (cover[to] != 0) LOGIC_ERROR;
                            escapes = true;
                            break;
                        }
                        if (escapes) continue;
                    }
                    Solver::solve(from, pl, owner(from) == pl ? n : -1);
                    // counts[priority(node)]--;
                    cover[from] = -1;
                    pms[k*from+pl] = -1;
                    q.push(from);
                    todo_push(from);
                }
            }
        }
    }

    if (trace) {
        for (int n=0; n<nodecount(); n++) {
            logger << "\033[1;31mnode " << n << "/" << priority(n) << (owner(n)?" (odd)":" (even)") << "\033[m is";
            pm_stream(logger, pms + k*n);
            if (cover[n] >= 0) logger << " cover: " << cover[n];
            logger << std::endl;
        }
    }

    // Now set dominions and derive strategy for even.
    for (int i=0; i<nodecount(); i++) {
        if (disabled[i]) continue;
        int *pm = pms + k*i;
        if ((pm[0] == -1) == (pm[1] == -1)) LOGIC_ERROR;
        const int winner = pm[0] == -1 ? 0 : 1;
        Solver::solve(i, winner, owner(i) == winner ? strategy[i] : -1);
    }

    delete[] pms;
    delete[] strategy;
    delete[] counts;
    delete[] cover;
    delete[] dirty;
    delete[] tmp;
    delete[] best;

    logger << "solved with " << lift_count << " lifts, " << lift_attempt << " lift attempts." << std::endl;
}

}
