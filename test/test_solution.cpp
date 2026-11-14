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
 * Unit tests for the Solution class.
 */

#include <iostream>
#include <sstream>
#include <vector>

#include "oink/oink.hpp"
#include "oink/solution.hpp"
#include "verifier.hpp"

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
test_solution()
{
    // A fresh Solution has all vertices unsolved with no strategy.
    Solution s(3);
    check("vertex_count", s.vertex_count() == 3);
    check("initially unsolved", !s.is_solved(0) and !s.is_solved(1) and !s.is_solved(2));
    check("initial strategy -1", s.strategy(0) == -1 and s.strategy(2) == -1);
    check("initially empty solved set", s.solved().count() == 0);

    // solve() records solved/winner/strategy.
    s.solve(0, 1, 2);   // vertex 0 won by Odd with strategy to vertex 2
    s.solve(1, 0, -1);  // vertex 1 won by Even with no strategy
    check("solved after solve", s.is_solved(0) and s.is_solved(1));
    check("unsolved untouched", !s.is_solved(2));
    check("winner recorded", s.winner(0) == 1 and s.winner(1) == 0);
    check("strategy recorded", s.strategy(0) == 2 and s.strategy(1) == -1);
    check("solved set count", s.solved().count() == 2);

    // set_winner / set_strategy update individual fields.
    s.set_winner(2, 1);
    s.set_strategy(2, 0);
    check("set_winner", s.winner(2) == 1);
    check("set_strategy", s.strategy(2) == 0);

    // strategy_data exposes the same values.
    check("strategy_data", s.strategy_data()[0] == 2 and s.strategy_data()[2] == 0);

    // Copy assignment yields an independent copy.
    Solution copy = s;
    copy.solve(2, 0, -1);
    check("copy independent (winner)", s.winner(2) == 1 and copy.winner(2) == 0);

    // resize preserves existing entries and clears new ones.
    s.resize(5);
    check("resize keeps old", s.is_solved(0) and s.strategy(0) == 2);
    check("resize clears new", !s.is_solved(4) and s.strategy(4) == -1);
    check("resize vertex_count", s.vertex_count() == 5);

    // reset clears everything.
    s.reset();
    check("reset unsolves", s.solved().count() == 0);
    check("reset clears strategy", s.strategy(0) == -1 and s.strategy(4) == -1);

    // The multi-strategy edge set (index-based ops are game-agnostic). first/
    // strategy_targets need a Game, so they are covered by the Game/Verifier tests.
    {
        Solution m;
        check("no multi by default", !m.has_multi());
        m.init_multi(10); // index ops are game-free
        check("has multi after init", m.has_multi() and m.solved().count() == 0);
        check("init clears edges", !m.has_edge_index(0) and !m.has_edge_index(9));

        m.add_edge_index(2);
        m.add_edge_index(3);
        m.add_edge_index(7);
        check("add edges", m.has_edge_index(2) and m.has_edge_index(3) and m.has_edge_index(7));
        check("others untouched", !m.has_edge_index(1) and !m.has_edge_index(4));
        m.remove_edge_index(3);
        check("remove one edge", !m.has_edge_index(3) and m.has_edge_index(2));

        // clear_edge_range clears a contiguous block (one vertex's edges)
        m.add_edge_index(3); m.add_edge_index(4); m.add_edge_index(5);
        m.clear_edge_range(2, 4); // indices 2..5
        check("clear_edge_range", !m.has_edge_index(2) and !m.has_edge_index(3)
                                  and !m.has_edge_index(4) and !m.has_edge_index(5));
        check("clear_edge_range keeps outside", m.has_edge_index(7));

        // reset clears the multi-strategy too.
        m.reset();
        check("reset clears multi", !m.has_multi() and !m.has_edge_index(7));
    }

    // A packed update can cross a storage word and must preserve neighboring edges.
    {
        constexpr std::size_t size = 193;
        constexpr uint64_t mask = 0xA965A965A965A965ULL;
        for (const std::size_t start : {0, 1, 31, 63, 64, 65, 127, 128, 129, 192, 193}) {
            for (const unsigned count : {0, 1, 2, 31, 63, 64}) {
                if (start + count > size) continue;
                Solution m;
                m.init_multi(size);
                for (std::size_t edge = 0; edge < size; edge += 2) m.add_edge_index(edge);
                m.replace_edge_bits(start, count, mask);
                for (std::size_t edge = 0; edge < size; ++edge) {
                    const bool expected = edge >= start && edge < start + count
                        ? (mask >> (edge - start)) & 1
                        : edge % 2 == 0;
                    check("replace_edge_bits preserves exact edge set", m.has_edge_index(edge) == expected);
                }
            }
        }
    }

    // A partial solution that is installed before Oink::run() must also be excluded
    // from preprocessing. Vertex 0 has a winning self-loop and is already solved. If
    // the self-loop preprocessor solves it again, its Odd predecessor 1 loses its edge
    // to 0 a second time and is attracted to Even, although Odd wins it via 3.
    for (bool scc : {false, true}) {
        Game game(4);
        const int priorities[] = {0, 1, 2, 3};
        const int owners[] = {0, 1, 1, 0};
        const std::vector<std::vector<int>> edges = {{0}, {0, 3}, {3}, {2}};
        for (int v = 0; v < 4; v++) {
            game.init_vertex(v, priorities[v], owners[v]);
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        game.ensure_sorted();
        std::ostringstream log;
        Oink controller(game, log);
        controller.setSolver("zlk");
        controller.setBottomSCC(scc);
        controller.solution().solve(0, 0, 0);
        bool verified = true;
        try {
            controller.run();
            Verifier(game, controller.solution(), log).verify(true, true, true);
        } catch (const std::exception& e) {
            std::cerr << "partial solution: " << e.what() << std::endl;
            verified = false;
        }
        check("partial solution excluded from preprocessing", verified and controller.solution().winner(1) == 1);
    }

    if (failures) {
        std::cerr << failures << " solution test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all solution tests passed" << std::endl;
    return 0;
}
