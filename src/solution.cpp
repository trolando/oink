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
first_strategy_edge(const Game& game, const Solution& solution, int v)
{
    if (!solution.has_multi()) return -1;
    const int f = game.firstout(v);
    const int c = game.outcount(v);
    const int* o = game.outedges();
    for (int k=0; k<c; k++) if (solution.has_edge_index((std::size_t)f + k)) return o[f + k];
    return -1;
}

bool
has_strategy_edge_to(const Game& game, const Solution& solution, int v, int to)
{
    if (!solution.has_multi()) return false;
    const int f = game.firstout(v);
    const int c = game.outcount(v);
    const int* o = game.outedges();
    for (int k=0; k<c; k++) if (o[f + k] == to and solution.has_edge_index((std::size_t)f + k)) return true;
    return false;
}

void
strategy_targets(const Game& game, const Solution& solution, int v, std::vector<int>& out)
{
    const int s = solution.strategy(v);
    if (s != -1) { out.push_back(s); return; }
    if (solution.has_multi() and solution.is_solved(v) and solution.winner(v) == game.owner(v)) {
        const int f = game.firstout(v);
        const int c = game.outcount(v);
        const int* o = game.outedges();
        for (int k=0; k<c; k++) if (solution.has_edge_index((std::size_t)f + k)) out.push_back(o[f + k]);
    }
}

void
permute(Game& game, Solution& solution, int* mapping)
{
    solution.permute(mapping); // non-destructive; must run before game.permute
    game.permute(mapping);     // destroys mapping in place
}

}
