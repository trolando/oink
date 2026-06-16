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
 * Unit tests for GameBuilder.
 */

#include <iostream>
#include <vector>

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

// Collect the successors of vertex v (raw, -1 terminated).
static std::vector<int>
succ(const Game& g, int v)
{
    std::vector<int> r;
    for (const int* e = g.outs(v); *e != -1; e++) r.push_back(*e);
    return r;
}

int
main()
{
    GameBuilder b(3);
    b.set_priority(0, 1); b.set_owner(0, Player::Even); b.set_label(0, "a");
    b.set_priority(1, 2); b.set_owner(1, Player::Odd); // no label
    b.set_priority(2, 3); b.set_owner(2, Player::Even); b.set_label(2, "c");
    b.add_edge(0, 1); b.add_edge(0, 2);
    b.add_edge(1, 0);
    b.add_edge(2, 2); // self-loop

    check("vertex_count()", b.vertex_count() == 3);

    Game g = b.build();

    check("vertexcount", g.vertexcount() == 3);
    check("edgecount", g.edgecount() == 4);

    check("priority", g.priority(0) == 1 and g.priority(1) == 2 and g.priority(2) == 3);
    check("owner", g.owner(0) == 0 and g.owner(1) == 1 and g.owner(2) == 0);

    const std::string* l0 = g.rawlabel(0);
    const std::string* l1 = g.rawlabel(1);
    const std::string* l2 = g.rawlabel(2);
    check("label 0", l0 != nullptr and *l0 == "a");
    check("label 1 (none)", l1 == nullptr);
    check("label 2", l2 != nullptr and *l2 == "c");

    check("outcount", g.outcount(0) == 2 and g.outcount(1) == 1 and g.outcount(2) == 1);
    check("edges 0", (succ(g, 0) == std::vector<int>{1, 2}));
    check("edges 1", (succ(g, 1) == std::vector<int>{0}));
    check("edges 2", (succ(g, 2) == std::vector<int>{2}));

    // span view over outgoing edges matches the raw iteration
    {
        check("out_edges size", g.out_edges(0).size() == 2);
        std::vector<int> via_span;
        for (int to : g.out_edges(0)) via_span.push_back(to);
        check("out_edges view", (via_span == std::vector<int>{1, 2}));
    }

    // Builder remains usable after build(): build a second, equivalent game.
    Game g2 = b.build();
    check("rebuild vertexcount", g2.vertexcount() == 3 and g2.edgecount() == 4);

    if (failures) {
        std::cerr << failures << " builder test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all builder tests passed" << std::endl;
    return 0;
}
