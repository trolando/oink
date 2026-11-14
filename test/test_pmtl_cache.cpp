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
 * Unit tests for PMTL: the blocking escapes of the shortcuts and the work counters.
 */

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "oink/oink.hpp"
#include "solvers/pmtl.hpp"
#include "verifier.hpp"

namespace pg {

class PMTLCacheProbe
{
public:
    static void checkCounters(Game& input, MeasureKind kind, bool tangles)
    {
        struct Counts {
            std::uint64_t lifts, iterations;
            std::string by_player;
        };
        auto run = [&](std::uint64_t seed) {
            Game game(input);
            std::ostringstream log;
            Oink controller(game, log, true);
            PMTLSolver solver(controller, game, kind, tangles);
            solver.lifts = solver.iterations = seed;
            solver.run();
            Verifier(game, controller.solution(), log).verify(true, true, true);
            const auto text = log.str();
            const auto start = text.find("number of lifts:");
            if (start == std::string::npos) throw std::runtime_error("missing per-player lift counts");
            return Counts{solver.lifts, solver.iterations, text.substr(start, text.find('\n', start)-start)};
        };
        const auto baseline = run(0);
        if (baseline.lifts == 0) throw std::runtime_error("counter regression needs actual lifts");
        // Exercise real solver work across both signed and unsigned 32-bit limits.
        for (std::uint64_t seed : {std::uint64_t(std::numeric_limits<int>::max()),
                                   std::uint64_t(std::numeric_limits<std::uint32_t>::max())}) {
            const auto actual = run(seed);
            if (actual.lifts != seed+baseline.lifts || actual.iterations != seed+baseline.iterations ||
                actual.by_player != baseline.by_player) {
                throw std::runtime_error("PMTL work counters overflowed or changed solver work");
            }
        }
    }

    static void check(Game& game, MeasureKind kind, bool tangles)
    {
        std::ostringstream log;
        Oink controller(game, log, true);
        PMTLSolver solver(controller, game, kind, tangles);
        const int n = game.nodecount();
        solver.G.resize(n);
        solver.G.flip();
        solver.Q.resize(n);
        solver.Updated.resize(n);
        solver.Changed[0].resize(n);
        solver.Changed[1].resize(n);
        solver.escape_offset.resize(n);
        solver.escape_stamp.resize(n);
        solver.escape_epoch = std::numeric_limits<unsigned int>::max();
        auto incoming = std::make_unique<std::vector<int>[]>(n);
        auto members = std::make_unique<int[]>(5);
        const int pairs[] = {0, 1, 1, -1, -1};
        std::copy(pairs, pairs+5, members.get());
        incoming[2].push_back(0);
        solver.tin = incoming.get();
        solver.tv.push_back(members.get());
        solver.tpr.push_back(0);

        unsigned expected_epoch = 1;
        for (int player : {0, 0, 1}) {
            std::fill(solver.escape_offset.begin(), solver.escape_offset.end(), 999);
            std::fill(solver.escape_stamp.begin(), solver.escape_stamp.end(),
                      expected_epoch == 1 ? 1 : expected_epoch-1);
            solver.Updated.reset();
            for (int v : {2, 3, 6}) solver.Updated[v] = true;
            solver.Changed[0].reset();
            solver.Changed[1].reset();
            std::unique_ptr<Measures> source(new_measure(kind, game, player));
            std::unique_ptr<Measures> target1(new_measure(kind, game, 1-player));
            std::unique_ptr<Measures> target2(new_measure(kind, game, 1-player));
            // A stable opposing-priority value also covers Ordered's Bottom.
            if (player == 0) {
                for (int v = 0; v < n; ++v) source->see(v, 1);
            }
            solver.shortcuts(player, *source, *target1, *target2);
            if (solver.escape_epoch != expected_epoch || solver.Q.nonempty()) {
                throw std::runtime_error("PMTL shortcut epoch or queue was not completed");
            }
            for (int v = 0; v < n; ++v) {
                const bool excluded = v == 5 || (player == 0 && v == 4);
                if (bool(solver.Updated[v]) != excluded ||
                    bool(solver.Changed[1-player][v]) != excluded ||
                    target1->is_top(v) != excluded || target2->is_top(v) != excluded ||
                    source->is_top(v)) {
                    throw std::runtime_error("PMTL escape cache changed the closure or target measures");
                }
            }
            if (player == 0 && (solver.escape_stamp[4] != expected_epoch || solver.escape_offset[4] != 5)) {
                throw std::runtime_error("PMTL blocker did not advance after its target was attracted");
            }
            ++expected_epoch;
        }
    }
};

}

int
test_pmtl_cache()
{
    using namespace pg;
    Game game(7);
    for (int v = 6; v >= 0; --v) game.init_vertex(v, v < 2 ? 0 : 1, v == 1 || v == 4 ? 1 : 0);
    // The Even tangle {0,1} escapes to 2. LIFO seeds 6 then 3 produce a live
    // blocker at 0; processing seed 2 attracts the tangle and invalidates it.
    // Scanning then passes the duplicate edge to 0 and blocks at vertex 5.
    const std::vector<std::vector<int>> edges = {{1}, {0, 2}, {2}, {3}, {2, 3, 6, 0, 0, 5}, {5}, {6}};
    for (int v = 0; v < 7; ++v) {
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    game.ensure_sorted();
    for (auto kind : {MeasureKind::Small, MeasureKind::Succinct, MeasureKind::Ordered, MeasureKind::Strahler}) {
        for (bool tangles : {false, true}) {
            PMTLCacheProbe::check(game, kind, tangles);
            PMTLCacheProbe::checkCounters(game, kind, tangles);
        }
    }
    std::cout << "PMTL escape cursors preserve tangle and vertex closures across calls, player changes and epoch wrap." << std::endl;
    return 0;
}
