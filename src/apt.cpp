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

#include <cassert>
#include <cstring>
#include <unistd.h>

#include "apt.hpp"

namespace pg {

APTSolver::APTSolver(Oink *oink, Game *game) : Solver(oink, game)
{
}

APTSolver::~APTSolver()
{
}

void
APTSolver::run()
{
    unsigned int iterations = 0;

    // allocate and initialize data structures
    int8_t *distraction = new int8_t[n_nodes];
    memset(distraction, 0, sizeof(int8_t[n_nodes]));

    for (;;) {
        iterations++;

        /**
         * We update distractions from low to high.
         * We process per block of priorities.
         * When a block changes, all lower priority vertices are reset.
         *
         * The algorithm refines a partition into two winning regions.
         * APT terminates when the winning regions are inductive, i.e.,
         * no player escapes in one step.
         */
        bool changed = false;
        int cur_pr = priority[0]; // priority of current block
        int first_idx = 0;        // first idx of current block

        for (int i=0; i<n_nodes;i++) {
            if (disabled[i]) continue;
            // on-the-fly compression like Verver's Zielonka trick
            if ((priority[i]&1) != (cur_pr&1)) {
                // new block!
                if (changed) break;
                // the previous block did not change, continue with new block
                cur_pr = priority[i];
                first_idx = i;
            }
            if (distraction[i]) continue; // already a distraction (increases monotonically)
            // now compute if the vertex is attracted to the opponent
            int y;
            if (owner[i] == 0) {
                // see if player Even can go to a vertex currently good for Even
                y = 1;
                const int *_out = outs + outa[i];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    // check if good for Odd (odd priority XOR distraction)
                    // a) pr[to] is Even and it's a distraction
                    // b) pr[to] is Odd and it's not a distraction
                    if ((priority[to] & 1) != distraction[to]) continue;
                    // good for player Even
                    y = 0;
                    break;
                }
            } else {
                // see if player Odd can go to a vertex currently good for Odd
                y = 0;
                const int *_out = outs + outa[i];
                for (int to = *_out; to != -1; to = *++_out) {
                    if (disabled[to]) continue;
                    // check if good for Even (odd priority EQUIV distraction)
                    // a) pr[to] is Even and it's not distraction
                    // b) pr[to] is Odd and it's a distraction
                    if ((priority[to] & 1) == distraction[to]) continue;
                    // good for player Odd
                    y = 1;
                    break;
                }
            }
            // if good-for-player equals vertex-parity, then it's not a distraction
            if ((cur_pr&1) == y) continue;
            changed = true;
            distraction[i] = true;
#ifndef NDEBUG
            if (trace >= 2) logger << "vertex " << label_vertex(i) << " is now a distraction" << std::endl;
#endif
        }

        if (!changed) break; // nothing changed, we're done

        // the previous block changed, reset all distraction from 0 to first_idx
        memset(distraction, 0, sizeof(int8_t[first_idx]));
#ifndef NDEBUG
        if (trace) logger << "restarting after finding distractions of prio " << cur_pr << std::endl;
#endif
    }

    // done -- unfortunately, APT does not compute a strategy!
    for (int i=0; i<n_nodes; i++) {
        if (disabled[i]) continue;
        oink->solve(i, distraction[i] ? 1-(priority[i]&1) : (priority[i]&1), -1);
    }

    // free allocated data structures
    delete[] distraction;

    logger << "solved with " << iterations << " iterations." << std::endl;
}

}
