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
#include <deque>
#include <stack>
#include <cstring> // for memset
#include "oink/oink.hpp"

namespace pg {

/**
 * Tarjan's SCC algorithm, modified to only compute the bottom SCC and avoid disabled nodes.
 */
void
Oink::tarjan(int n, std::vector<int> &res, bool nonempty)
{
    // initialize
    unsigned n_nodes = game->nodecount();
    int *low = new int[n_nodes];
    memset(low, 0, sizeof(int[n_nodes])); // TODO: move this array outside the method
    int pre = 0;

    // search stack "st"
    std::stack<int> st;
    st.push(n);

    while (!st.empty()) { // st is never empty
        int idx = st.top();
        if (low[idx] == 0) {
            // first time we see it
            low[idx] = ++pre;
            res.push_back(idx);
        }
        int min = low[idx];
        bool pushed = false;
        auto curedge = game->outedges() + game->firstout(idx);
        for (int to_idx = *curedge; to_idx != -1; to_idx = *++curedge) {
            if (disabled[to_idx]) continue;
            if (low[to_idx] == 0) {
                // not visited
                st.push(to_idx);
                pushed = true;
                break;
            } else {
                // visited it, update min
                if (low[to_idx] < min) min = low[to_idx];
            }
        }
        if (pushed) continue; // we pushed...

        if (min < low[idx]) {
            low[idx] = min;
            st.pop();
            continue;
        }

        /**
         * At this point, we found a bottom SCC. Now check if it is empty.
         * A SCC is "empty" if it contains no edges, i.e., consists of 1 node without self-loops.
         */

        if (nonempty) {
            if (res.back() == idx and !game->has_edge(idx, idx)) {
                // it has no edges!
                res.pop_back();
                st.pop();
                if (st.empty()) break;
                continue;
            }
        }

        // found bottom SCC; remove nodes not in the bottom SCC...
        if (res.front() != idx) res.erase(res.begin(), std::find(res.begin(), res.end(), idx));
        break;
    }

    delete[] low;
}

/**
 * Find a bottom SCC starting from the first unsolved node.
 */
void
Oink::getBottomSCC(std::vector<int> &scc, bool nonempty)
{
    scc.clear();
    for (int i=0; i<game->nodecount(); i++) {
        if (!disabled[i]) {
            tarjan(i, scc, nonempty);
            return;
        }
    }
}

/**
 * Find a bottom SCC starting from the given node.
 */
void
Oink::getBottomSCC(int start, std::vector<int> &scc, bool nonempty)
{
    scc.clear();
    tarjan(start, scc, nonempty);
}

}
