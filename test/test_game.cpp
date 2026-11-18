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

/**
 * Unit tests for Game dynamic growth.
 */

#include <iostream>
#include <stdexcept>
#include <vector>

#include "oink/game.hpp"
#include "oink/game_builder.hpp"

using namespace pg;

static int failures = 0;

static void
check(const char* name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

int
test_game()
{
    // Start from a single-vertex game and add more vertices. This grows the
    // vertex arrays via v_sizeup; with a 1-vertex start the growth step must
    // still advance (a buggy "size += size/2" would stay at 1 and loop forever).
    Game g(1);
    g.init_vertex(0, 5, 0);
    g.init_vertex(1, 6, 1); // triggers growth from capacity 1
    g.init_vertex(2, 7, 0); // triggers growth again

    check("grew to hold all vertices", g.vertexcount() >= 3);
    check("priorities preserved", g.priority(0) == 5 and g.priority(1) == 6 and g.priority(2) == 7);
    check("owners preserved", g.owner(0) == 0 and g.owner(1) == 1 and g.owner(2) == 0);

    // Streaming rejects a duplicate before storing it, including in a long list.
    {
        Game h(12, 1);
        h.e_start(0);
        for (int to = 0; to < 10; to++) h.e_add(0, to);
        bool rejected = false;
        try {
            h.e_add(0, 2);
        } catch (const std::runtime_error& e) {
            rejected = std::string(e.what()) == "repeated edge from vertex 0 to successor 2";
        }
        check("streaming reports repeated successor", rejected);
        h.e_add(0, 10);
        h.e_finish();
        h.e_start(1);
        h.e_add(1, 2); // the same target from a different source is valid
        h.e_finish();
        check("streaming can continue after rejection", h.outcount(0) == 11 and h.edgecount() == 12);
    }

    // The adjacency-list constructor also checks lists supplied by a caller.
    {
        std::vector<int> priorities(2);
        bitset owners(2);
        std::vector<std::vector<int>> edges{{1}, {1, 0, 1}};
        std::vector<std::string*> labels(2);
        bool rejected = false;
        try {
            Game h(2, 4, priorities, owners, edges, labels);
        } catch (const std::runtime_error& e) {
            rejected = std::string(e.what()) == "repeated edge from vertex 1 to successor 1";
        }
        check("adjacency lists report a repeated self-loop", rejected);
    }

    // build_in_to_out: the in->out edge map must align with the in-edge array, so
    // that for every incoming-edge slot of vertex w the mapped out-edge index
    // belongs to the recorded source and points back to w.
    {
        GameBuilder b(4);
        b.set_priority(0, 0); b.set_owner(0, Player::Even);
        b.set_priority(1, 0); b.set_owner(1, Player::Odd);
        b.set_priority(2, 0); b.set_owner(2, Player::Even);
        b.set_priority(3, 0); b.set_owner(3, Player::Odd);
        b.add_edge(0, 1); b.add_edge(0, 2); b.add_edge(0, 0); // incl. self-loop
        b.add_edge(1, 2);
        b.add_edge(2, 0); b.add_edge(2, 3);
        b.add_edge(3, 3);
        Game h = b.build();
        h.build_in_to_out();

        const int* in = h.inedges();
        const int* in2out = h.inToOut();
        const int* out = h.outedges();
        bool aligned = true;
        long checked = 0;
        for (int w = 0; w < h.nodecount(); w++) {
            const int* s = h.ins(w);
            for (int t = 0; t < h.incount(w); t++) {
                int slot = (int)((s + t) - in);
                int from = in[slot];
                int oidx = in2out[slot];
                // the mapped out-edge must point to w and lie in <from>'s out-block
                if (out[oidx] != w) aligned = false;
                if (oidx < h.firstout(from) || oidx >= h.firstout(from) + h.outcount(from)) aligned = false;
                checked++;
            }
        }
        check("in_to_out maps every in-edge", checked == h.edgecount());
        check("in_to_out endpoints aligned", aligned);
    }

    // strategyTargets dispatches per vertex via the -1 sentinel: a single move,
    // the multi-strategy set, or nothing (loser/unsolved).
    {
        GameBuilder b(3);
        b.set_priority(0, 2); b.set_owner(0, Player::Even); b.add_edge(0, 1); b.add_edge(0, 2);
        b.set_priority(1, 2); b.set_owner(1, Player::Even); b.add_edge(1, 1);
        b.set_priority(2, 2); b.set_owner(2, Player::Odd);  b.add_edge(2, 2);
        Game h = b.build();
        Solution sol(h.nodecount());
        sol.init_multi(h.edgeArraySize());
        sol.solve(0, 0, -1);                       // won by owner, sentinel -> read multi
        sol.add_edge_index(h.firstout(0) + 0);     // edge 0->1
        sol.add_edge_index(h.firstout(0) + 1);     // edge 0->2
        sol.solve(1, 0, 1);                        // won by owner, single move
        sol.solve(2, 0, -1);                       // owner Odd, won by Even -> loser, no strategy

        std::vector<int> t0, t1, t2;
        strategy_targets(h, sol, 0, t0);
        strategy_targets(h, sol, 1, t1);
        strategy_targets(h, sol, 2, t2);
        check("strategyTargets single", (t1 == std::vector<int>{1}));
        check("strategyTargets multi", (t0 == std::vector<int>{1, 2}));
        check("strategyTargets loser empty", t2.empty());
    }

    // Sorting by priority keeps the order of vertices with the same priority, so that every
    // platform solves the same game (std::sort alone leaves their order unspecified).
    {
        const int n = 300;
        Game h(n);
        for (int v = 0; v < n; v++) {
            h.init_vertex(v, (v * 7) % 3, v & 1);
            h.e_start(v);
            h.e_add(v, v);
            h.e_finish();
        }
        std::vector<int> mapping(n);
        h.sort(mapping.data());
        bool sorted = true, stable = true;
        for (int i = 1; i < n; i++) {
            if (h.priority(i-1) > h.priority(i)) sorted = false;
            else if (h.priority(i-1) == h.priority(i) and mapping[i-1] > mapping[i]) stable = false;
        }
        check("sort orders by priority", sorted);
        check("sort keeps the order of equal priorities", stable);
    }

    if (failures) {
        std::cerr << failures << " game test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all game tests passed" << std::endl;
    return 0;
}
