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
main()
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
        h.initMultiStrategy();
        h.solve(0, 0, -1);        // won by owner, sentinel -> read multi
        h.addStrategyEdge(0, 0);  // edge 0->1
        h.addStrategyEdge(0, 1);  // edge 0->2
        h.solve(1, 0, 1);         // won by owner, single move
        h.solve(2, 0, -1);        // owner Odd, won by Even -> loser, no strategy

        check("strategyTargets single", (h.strategyTargets(1) == std::vector<int>{1}));
        check("strategyTargets multi", (h.strategyTargets(0) == std::vector<int>{1, 2}));
        check("strategyTargets loser empty", h.strategyTargets(2).empty());
    }

    if (failures) {
        std::cerr << failures << " game test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all game tests passed" << std::endl;
    return 0;
}
