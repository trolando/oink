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

#include <algorithm>
#include <climits>
#include <cstdint>
#include <iostream>
#include <memory>
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
        for (auto kind : {MeasureKind::Small, MeasureKind::Ordered, MeasureKind::Succinct, MeasureKind::Strahler}) {
            for (int player = 0; player < 2; player++) {
                std::unique_ptr<Measures> source(new_measure(kind, game, player));
                std::unique_ptr<Measures> target(new_measure(kind, game, player));
                // Switching player can resize the tuple storage even for the same graph.
                std::unique_ptr<Measures> opposite(new_measure(kind, game, 1-player));
                opposite->copy_bot(-1);
                opposite->see(-1, 1-player);
                target->set(*opposite);
                require(target->store(-1) == opposite->store(-1), "resized set changed a measure");
                target->set(*source);
                source->copy_bot(-1);
                source->see(-1, player);
                const auto before = source->store(-1);
                source->set(*source);
                require(source->store(-1) == before, "self-set changed a measure");
                for (int step = 0; step < 100; step++) {
                    source->see(-1, step%5);
                    target->copy_from(*source, -1, -2);
                    target->load(-3, source->store(-1));
                    require(target->eq(-2, -3), "copy_from changed a measure");
                    target->copy(-2, -4);
                    require(target->eq(-2, -4), "copy changed a measure");
                    target->set(*source);
                    require(target->store(-1) == source->store(-1), "set changed a measure");
                }
                target->copy_bot(-1);
                target->see(-1, player);
                source->copy_top(-2);
                target->copy_from(*source, -2, -1);
                require(target->is_top(-1), "Top copy_from is not Top");
                target->copy(-1, -3);
                require(target->eq(-1, -3), "Top copies compare unequal");
            }
        }

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

        // Table sizes and row offsets of Small and Ordered, without allocating tables:
        // 65,536 vertices with 32,768 components per slot exceed INT_MAX components.
        const size_t wide = (size_t)(65536+4) * 32768;
        require(wide > (size_t)INT_MAX, "test size does not exceed INT_MAX");
        require(measure_table_size(65536+4, 32768, sizeof(int)) == wide, "table size overflowed");
        require(measure_offset(65535, 32768) == wide - 32768, "row offset of the last vertex overflowed");
        require(measure_offset(-4, 32768) == 0 and measure_offset(-1, 3) == 9, "scratch slot offsets");
        require(measure_table_size(0, 5, sizeof(int)) == 0 and measure_table_size(5, 0, sizeof(int)) == 0,
                "empty table size");
        if (sizeof(size_t) >= 8) {
            require(measure_offset(INT_MAX-4, INT_MAX) == (size_t)INT_MAX * (size_t)INT_MAX,
                    "row offset of the largest slot and width overflowed");
        }
        // the largest table that fits, and one slot more
        const long width = 1024;
        const long slots = (long)std::min<size_t>(SIZE_MAX / sizeof(int) / width, LONG_MAX);
        require(measure_table_size(slots, width, sizeof(int)) == (size_t)slots * width, "largest table size");
        bool rejected = false;
        try {
            measure_table_size(slots+1, width, sizeof(int));
        } catch (const std::length_error&) {
            rejected = true;
        }
        require(rejected, "a table larger than the address space was not rejected");
        rejected = false;
        try {
            measure_table_size(-1, width, sizeof(int));
        } catch (const std::length_error&) {
            rejected = true;
        }
        require(rejected, "a negative table size was not rejected");
        std::cout << "measure operations passed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
