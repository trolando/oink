/*
 * Copyright 2026 Tom van Dijk
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

/**
 * Unit tests for Game::build_in_to_out, which can reuse the incoming counts.
 */

#include <iostream>
#include <stdexcept>
#include <vector>

#include "oink/game.hpp"

using namespace pg;

/**
 * Check the incoming edges and the in-to-out map of <game>, and return their layout.
 */
static std::vector<int>
signature(const Game& game)
{
    const int n = game.nodecount();
    std::vector<int> result, counts(n, 0);
    std::vector<bool> seen(game.edgeArraySize(), false);
    for (int v = 0; v < n; ++v) {
        for (auto edge = game.outs(v); *edge != -1; ++edge) ++counts[*edge];
    }
    for (int v = 0; v < n; ++v) {
        if (game.incount(v) != counts[v]) throw std::runtime_error("stale incoming count");
        result.push_back(game.firstin(v));
        result.push_back(game.incount(v));
        for (int i = 0; i < game.incount(v); ++i) {
            const int slot = game.firstin(v)+i;
            const int from = game.inedges()[slot];
            const int out = game.inToOut()[slot];
            if (out < game.firstout(from) || out >= game.firstout(from)+game.outcount(from) ||
                game.outedges()[out] != v || seen[out]) {
                throw std::runtime_error("invalid incoming-to-outgoing edge correspondence");
            }
            seen[out] = true;
            result.push_back(from);
            result.push_back(out);
        }
        const int sentinel = game.firstin(v)+game.incount(v);
        if (game.inedges()[sentinel] != -1 || game.inToOut()[sentinel] != -1) {
            throw std::runtime_error("incoming terminator/map was not reset");
        }
    }
    for (int v = 0; v < n; ++v) {
        for (int i = 0; i < game.outcount(v); ++i) {
            if (!seen[game.firstout(v)+i]) throw std::runtime_error("outgoing edge is missing from map");
        }
    }
    return result;
}

/**
 * Building the map with reused counts must give the layout of a full rebuild.
 */
static void
compare(Game& game)
{
    Game reference(game);
    reference.build_in_to_out();
    const auto expected = signature(reference);
    game.build_in_to_out(true);
    if (signature(game) != expected) throw std::runtime_error("reused counts changed edge layout");
    game.build_in_to_out(true);
    if (signature(game) != expected) throw std::runtime_error("repeated map build changed edge layout");
}

static Game
make(int n)
{
    Game game(n, 1000);
    for (int v = 0; v < n; ++v) game.init_vertex(v, v%3, v%2);
    for (int v = 0; v < n; ++v) {
        game.e_start(v);
        // The last vertex has zero indegree; duplicate and self edges remain.
        game.e_add(v, v%(n-1));
        game.e_add(v, (v+1)%(n-1));
        game.e_add(v, (v+1)%(n-1));
        game.e_finish();
    }
    return game;
}

int
test_in_edge_map()
{
    Game a = make(6);
    Game b = make(9);
    // cold fallback, then repeated reuse
    compare(a);
    compare(b);
    // a copy has shorter incoming storage
    Game copy(a);
    compare(copy);
    std::vector<int> permutation(a.nodecount());
    for (int v = 0; v < a.nodecount(); ++v) permutation[v] = a.nodecount()-1-v;
    a.permute(permutation.data());
    compare(a);
    // the maps are scratch, so rebuild after a swap
    a.swap(b);
    compare(a);
    compare(b);
    b.vec_init();
    int target = 0;
    while (target < b.nodecount() && b.vec_has_edge(0, target)) ++target;
    if (target == b.nodecount() || !b.vec_add_edge(0, target)) {
        throw std::runtime_error("edge edit fixture did not change the graph");
    }
    b.vec_finish();
    // the default must still rebuild after structural changes
    b.build_in_to_out();
    signature(b);
    compare(b);
    std::cout << "Incoming edge maps preserve counts, duplicates and exact layout across reuse, copy, permutation, swap and full rebuild." << std::endl;
    return 0;
}
