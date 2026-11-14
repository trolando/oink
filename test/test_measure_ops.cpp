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
 * Unit tests for the operations of the progress measures.
 */

#include <iostream>
#include <stdexcept>
#include <vector>

#include "solvers/gpm.hpp"

using namespace pg;

static void
require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

/**
 * A game without edges, with one vertex of each priority 0..maxp.
 */
static Game
make_game(int maxp)
{
    Game game(maxp+1);
    game.vec_init();
    for (int p = 0; p <= maxp; p++) game.init_vertex(p, p, p&1);
    game.vec_finish();
    game.sort();
    return game;
}

int
test_measure_ops()
{
    try {
        auto game = make_game(4);

        // The tuples of Small measures for Odd have one entry per odd priority:
        // with one vertex of priority 1 and one of priority 3, there are 4 values.
        Small odd(game, 1);
        odd.copy_bot(-1);
        for (int i = 0; i < 4; i++) {
            require(!odd.is_top(-1), "odd Small overflowed too early");
            odd.inc(-1);
        }
        require(odd.is_top(-1), "odd Small did not overflow after all tuples");
        auto even_only = make_game(0);
        Small empty(even_only, 1);
        empty.copy_bot(-1);
        empty.inc(-1);
        require(empty.is_top(-1), "empty Small tuple did not increment to Top");
        std::cout << "measure operations passed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
