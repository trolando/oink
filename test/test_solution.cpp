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

#include "oink/solution.hpp"

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

    if (failures) {
        std::cerr << failures << " solution test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all solution tests passed" << std::endl;
    return 0;
}
