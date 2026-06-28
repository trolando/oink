/*
 * Copyright 2024 Tom van Dijk
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

#include "oink/solution.hpp"
#include "oink/game.hpp"

namespace pg {

int
Solution::first_strategy_edge(int v) const
{
    if (!has_multi_) return -1;
    const int f = game_->firstout(v);
    const int c = game_->outcount(v);
    const int* o = game_->outedges();
    for (int k=0; k<c; k++) if (edges_.test((std::size_t)f + k)) return o[f + k];
    return -1;
}

bool
Solution::has_strategy_edge_to(int v, int to) const
{
    if (!has_multi_) return false;
    const int f = game_->firstout(v);
    const int c = game_->outcount(v);
    const int* o = game_->outedges();
    for (int k=0; k<c; k++) if (o[f + k] == to and edges_.test((std::size_t)f + k)) return true;
    return false;
}

void
Solution::strategy_targets(int v, std::vector<int>& out) const
{
    const int s = strategy_[v];
    if (s != -1) { out.push_back(s); return; }
    if (has_multi_ and is_solved(v) and winner(v) == game_->owner(v)) {
        const int f = game_->firstout(v);
        const int c = game_->outcount(v);
        const int* o = game_->outedges();
        for (int k=0; k<c; k++) if (edges_.test((std::size_t)f + k)) out.push_back(o[f + k]);
    }
}

}
