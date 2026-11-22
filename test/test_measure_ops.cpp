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
#include <sstream>
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

/**
 * Small::copy_see gives the same measures as copy and see, also when vertices reach Top
 * (which lowers the bounds of the carries), and when the target is a vertex (which lowers
 * the bounds as copy does). Two copies of the measures follow the same pseudo-random steps,
 * one with copy and see, one with copy_see, and must keep the same values.
 */
static void
check_small_copy_see(const Game& game, int player, int steps, unsigned seed)
{
    Small a(game, player);
    Small b(game, player);
    const int n = game.nodecount(), maxp = game.priority(n-1);
    unsigned x = seed;
    auto next = [&](unsigned m) { x = x * 1103515245u + 12345u; return (int)((x >> 8) % m); };
    for (int s = -4; s < n; s++) {
        a.copy_bot(s);
        b.copy_bot(s);
    }
    for (int i = 0; i < steps; i++) {
        const int from = next(n+4) - 4;
        const int op = next(10);
        // (mostly into a scratch slot, as the lift does)
        int to = op < 4 ? -1 - next(2) : next(n+4) - 4;
        if (to == from) to = from == -1 ? -2 : -1;
        if (op < 6) {
            // the operation under test
            const int p = next(maxp+1);
            a.copy(from, to);
            a.see(to, p);
            b.copy_see(from, to, p);
        } else if (op < 8) {
            // grow the measures, to have carries
            const int p = next(3) == 0 ? player : next(maxp+1);
            for (int k = 1 + next(8); k > 0; k--) {
                a.see(from, p);
                b.see(from, p);
            }
        } else if (op == 8 and next(4) == 0) {
            // lower a bound
            a.copy_top(from);
            b.copy_top(from);
        } else {
            a.inc(from);
            b.inc(from);
        }
        require(a.is_top(to) == b.is_top(to), "copy_see differs in Top from copy and see");
        require(a.is_top(to) or a.store(to) == b.store(to), "copy_see differs from copy and see");
        require(a.compare(to, from) == b.compare(to, from), "copy_see compares differently");
    }
    for (int s = -4; s < n; s++) {
        require(a.is_top(s) == b.is_top(s) and (a.is_top(s) or a.store(s) == b.store(s)),
                "copy_see left different measures");
    }
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

        // copy_see is copy and see, on games with few and with many vertices per priority
        for (int player = 0; player < 2; player++) {
            for (int maxp : {0, 1, 2, 3, 6, 9, 40}) {
                check_small_copy_see(make_game(maxp), player, 20000, 7u + maxp + 100u*player);
            }
            for (int n : {12, 40, 90}) {
                Game few(n);
                few.vec_init();
                for (int v = 0; v < n; v++) few.init_vertex(v, v % 7, v&1);
                few.vec_finish();
                few.sort();
                check_small_copy_see(few, player, 50000, 11u + n + player);
            }
        }

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

/**
 * A game without edges, with <n> vertices, of the priorities in <priorities> in turn.
 */
static Game
make_game_of(int n, const std::vector<int>& priorities)
{
    Game game(n);
    game.vec_init();
    for (int v = 0; v < n; v++) {
        const int p = priorities[v % priorities.size()];
        game.init_vertex(v, p, p&1);
    }
    game.vec_finish();
    game.sort();
    return game;
}

/**
 * Apply the same operations to the measures <a> and <b> (of the same kind, game and player, but
 * codes of different sizes) and require that they keep the same values: store() writes a value
 * independently of the size of the codes. First walk with inc from bottom to Top (all measures
 * of a small tree, at most <walk> steps), then do <steps> pseudo-random operations on the slots
 * -2, -1 and the first (at most six) vertices of the game of <n> vertices. (Ordered::see uses -3,
 * and PMTL reserves -4.)
 */
static void
check_same_values(Measures& a, Measures& b, int n, int maxp, int walk, int steps, unsigned seed)
{
    // (only the values: the codes after the first of a Top slot are not used, and can differ
    // when an empty tree of Succinct overflows its codes to Top, with depth -1)
    auto same = [&](int slot) {
        require(a.is_top(slot) == b.is_top(slot), "codes of different sizes differ in Top");
        if (!a.is_top(slot) and a.store(slot) != b.store(slot)) {
            std::ostringstream msg;
            msg << "codes of different sizes give different measures: ";
            a.stream(msg, slot);
            msg << " and ";
            b.stream(msg, slot);
            throw std::runtime_error(msg.str());
        }
    };
    a.copy_bot(-1);
    b.copy_bot(-1);
    for (int i = 0; i < walk and !a.is_top(-1); i++) {
        a.copy(-1, -2);
        b.copy(-1, -2);
        a.inc(-1);
        b.inc(-1);
        same(-1);
        require(a.compare(-2, -1) == -1 and b.compare(-2, -1) == -1, "inc did not increase");
    }
    same(-1);
    // pseudo-random see, inc, copies and comparisons
    unsigned x = seed;
    auto next = [&](unsigned m) { x = x * 1103515245u + 12345u; return (int)((x >> 8) % m); };
    const int slots = 2 + std::min(n, 6);
    for (int s = -2; s < slots-2; s++) {
        a.copy_bot(s);
        b.copy_bot(s);
    }
    for (int i = 0; i < steps; i++) {
        const int slot = next(slots) - 2, other = next(slots) - 2;
        const int op = next(16);
        if (op < 9) {
            // see a priority, mostly the highest and lowest few (the boundary codes)
            const int r = next(6);
            const int p = r == 0 ? next(maxp+1) : r < 3 ? std::max(0, maxp - next(3)) : std::min(maxp, next(3));
            a.see(slot, p);
            b.see(slot, p);
        } else if (op < 11) {
            a.inc(slot);
            b.inc(slot);
        } else if (op < 13) {
            a.copy(other, slot);
            b.copy(other, slot);
        } else if (op == 13) {
            a.copy_bot(slot);
            b.copy_bot(slot);
        } else if (op == 14 and next(8) == 0) {
            a.copy_top(slot);
            b.copy_top(slot);
        } else {
            a.load(slot, b.store(other));
            b.load(slot, a.store(other));
        }
        same(slot);
        require(a.compare(slot, other) == b.compare(slot, other), "codes of different sizes compare differently");
        require(a.eq(slot, other) == b.eq(slot, other), "codes of different sizes are equal differently");
    }
}

template<class Narrow, class Wide>
static void
check_codes(const Game& game, int player, int walk, int steps)
{
    Narrow narrow(game, player);
    Wide wide(game, player);
    const int n = game.nodecount(), maxp = game.priority(n-1);
    check_same_values(narrow, wide, n, maxp, walk, steps, 1u + (unsigned)maxp*2u + (unsigned)player);
    // a copy keeps the values (copy_from and set between equal types)
    Narrow other(game, player);
    other.copy_from(narrow, -1, 0);
    require(other.store(0) == narrow.store(-1), "copy_from changed a measure");
    other.set(narrow);
    require(other.store(-1) == narrow.store(-1), "set changed a measure");
}

template<class M>
static bool
fits(const Game& game)
{
    try {
        M m0(game, 0), m1(game, 1);
        return true;
    } catch (const pg::Error&) {
        return false;
    }
}

int
test_measure_codes()
{
    try {
        // the size of the codes at the boundaries: the largest priority of each size, and one more
        struct Boundary { MeasureKind kind; int d; int size; };
        const Boundary boundaries[] = {
            {MeasureKind::Ordered, 0, 1}, {MeasureKind::Ordered, 127, 1}, {MeasureKind::Ordered, 128, 4},
            {MeasureKind::Ordered, 32768, 4},
            {MeasureKind::Succinct, 0, 1}, {MeasureKind::Succinct, 251, 1}, {MeasureKind::Succinct, 252, 2},
            {MeasureKind::Succinct, 65531, 2}, {MeasureKind::Succinct, 65532, 4},
            {MeasureKind::Strahler, 0, 1}, {MeasureKind::Strahler, 249, 1}, {MeasureKind::Strahler, 250, 2},
            {MeasureKind::Strahler, 65529, 2}, {MeasureKind::Strahler, 65530, 4},
            {MeasureKind::Small, 2, (int)sizeof(int)},
        };
        for (const auto& b : boundaries) {
            for (int low = 0; low < 2; low++) {
                // a vertex of the highest priority, and one of priority 0 or 1 (both parities)
                const auto game = make_game_of(2, {low, b.d});
                require(measure_code_size(b.kind, game) == b.size, "wrong size of codes at a boundary");
                for (int player = 0; player < 2; player++) {
                    std::unique_ptr<Measures> pm(new_measure(b.kind, game, player));
                    require(pm->code_size() == b.size, "new_measure did not select the size");
                }
            }
        }

        // the narrower codes hold every code up to the boundary, and not one priority more
        // (as the measures check in their constructors)
        for (int low = 0; low < 2; low++) {
            require(fits<Ordered<int8_t>>(make_game_of(2, {low, 127})), "Ordered: 127 does not fit 8 bits");
            require(!fits<Ordered<int8_t>>(make_game_of(2, {low, 128})), "Ordered: 128 fits 8 bits");
            require(fits<Succinct<uint8_t>>(make_game_of(2, {low, 251})), "Succinct: 251 does not fit 8 bits");
            require(!fits<Succinct<uint8_t>>(make_game_of(2, {low, 252})), "Succinct: 252 fits 8 bits");
            require(fits<Succinct<uint16_t>>(make_game_of(2, {low, 65531})), "Succinct: 65531 does not fit 16 bits");
            require(!fits<Succinct<uint16_t>>(make_game_of(2, {low, 65532})), "Succinct: 65532 fits 16 bits");
            require(fits<Strahler<uint8_t>>(make_game_of(2, {low, 249})), "Strahler: 249 does not fit 8 bits");
            require(!fits<Strahler<uint8_t>>(make_game_of(2, {low, 250})), "Strahler: 250 fits 8 bits");
            require(fits<Strahler<uint16_t>>(make_game_of(2, {low, 65529})), "Strahler: 65529 does not fit 16 bits");
            require(!fits<Strahler<uint16_t>>(make_game_of(2, {low, 65530})), "Strahler: 65530 fits 16 bits");
        }

        // the same values in every size: games with 1, 2 and 3 priorities, more priorities, and
        // the largest priorities of the narrower sizes
        struct Shape { int n; std::vector<int> priorities; };
        const Shape small[] = {
            {1, {0}}, {1, {1}}, {7, {0}}, {7, {1}}, {9, {0, 1}}, {200, {0, 1}}, {6, {1, 2}},
            {13, {0, 1, 2}}, {300, {0, 1, 2}}, {40, {1, 2, 3}}, {25, {0, 1, 2, 3, 4, 5, 6, 7}},
            {64, {0, 3, 4, 9}}, {1000, {0, 1, 2, 5, 6}},
        };
        for (const auto& shape : small) {
            const auto game = make_game_of(shape.n, shape.priorities);
            for (int player = 0; player < 2; player++) {
                check_codes<Ordered<int8_t>, Ordered<int32_t>>(game, player, 20000, 5000);
                check_codes<Succinct<uint8_t>, Succinct<uint32_t>>(game, player, 20000, 5000);
                check_codes<Succinct<uint16_t>, Succinct<uint32_t>>(game, player, 1000, 1000);
                check_codes<Strahler<uint8_t>, Strahler<uint32_t>>(game, player, 20000, 5000);
                check_codes<Strahler<uint16_t>, Strahler<uint32_t>>(game, player, 1000, 1000);
            }
        }
        // (with two vertices, the walks reach Top: they pass the deepest codes)
        for (int low = 0; low < 2; low++) {
            for (int n : {2, 5}) {
                for (int player = 0; player < 2; player++) {
                    check_codes<Ordered<int8_t>, Ordered<int32_t>>(make_game_of(n, {low, 126, 127}), player, 100000, 20000);
                    check_codes<Succinct<uint8_t>, Succinct<uint32_t>>(make_game_of(n, {low, 250, 251}), player, 100000, 20000);
                    check_codes<Strahler<uint8_t>, Strahler<uint32_t>>(make_game_of(n, {low, 248, 249}), player, 100000, 20000);
                    check_codes<Succinct<uint16_t>, Succinct<uint32_t>>(make_game_of(n, {low, 65530, 65531}), player, 100000, 5000);
                    check_codes<Strahler<uint16_t>, Strahler<uint32_t>>(make_game_of(n, {low, 65528, 65529}), player, 100000, 5000);
                }
            }
        }
        std::cout << "measure codes passed" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    return 0;
}
