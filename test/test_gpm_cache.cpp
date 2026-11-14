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
 * Unit tests for GPM: the blocking escapes of the shortcuts, the active subgame
 * and Top notifications.
 */

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "oink/oink.hpp"
#include "solvers/gpm.hpp"

namespace pg {

class GPMCacheProbe
{
public:
    static void checkKnownTop(MeasureKind kind)
    {
        Game game(3);
        game.init_vertex(0, 0, 0);
        game.init_vertex(1, 1, 1);
        game.init_vertex(2, 2, 0);
        game.e_start(0);
        game.e_add(0, 0);
        game.e_finish();
        game.e_start(1);
        game.e_add(1, 1);
        game.e_finish();
        game.e_start(2);
        game.e_add(2, 0);
        game.e_add(2, 1);
        game.e_finish();
        game.ensure_sorted();
        std::ostringstream log;
        Oink controller(game, log, true);
        GPMSolver solver(controller, game, kind);
        solver.G.resize(3);
        solver.G.flip();
        solver.Q0 = solver.G;
        solver.Q1.resize(3);
        solver.QQ.resize(3);
        solver.T0.resize(3);
        solver.T1.resize(3);
        solver.L0.resize(3);
        solver.L1.resize(3);
        solver.Queue.resize(3);
        solver.escape_offset.resize(3);
        solver.escape_stamp.resize(3);
        std::unique_ptr<Measures> pm0(new_measure(kind, game, 0));
        std::unique_ptr<Measures> pm1(new_measure(kind, game, 1));
        // Normalize the opposing priority in Ordered, then retain the next wave.
        solver.update(*pm0, 0);
        solver.update(*pm0, 0);
        const bitset source_wave = solver.L0;
        solver.shortcuts(0, *pm0, *pm1);
        if (!solver.T1[1] || !pm1->is_top(1) || !solver.Q1[2]) {
            throw std::runtime_error("fresh Top did not notify its non-Top predecessor");
        }
        // Consume that genuine notification. The unchanged Top target cannot
        // require another lift under fixed bounds.
        solver.update(*pm1, 1);
        if (solver.Q1.any() || solver.T1[2]) {
            throw std::runtime_error("known-Top regression needs a finite, settled predecessor");
        }
        solver.L0 = source_wave;
        solver.shortcuts(0, *pm0, *pm1);
        if (solver.Q1.any()) {
            throw std::runtime_error("unchanged Top scheduled redundant predecessor work");
        }
        if (solver.L0.count() != 1 || !solver.L0[1]) {
            throw std::runtime_error("known-Top filtering changed the full stable complement");
        }
    }

    static void checkDisabled(MeasureKind kind, int player, int priority)
    {
        Game game(2);
        for (int v = 0; v < 2; ++v) game.init_vertex(v, priority, 0);
        game.e_start(0);
        game.e_add(0, 0);
        game.e_finish();
        game.e_start(1);
        game.e_add(1, 0);
        game.e_finish();
        game.ensure_sorted();
        std::ostringstream log;
        Oink controller(game, log, true);
        GPMSolver solver(controller, game, kind);
        // Vertex 1 is a disabled predecessor of the active self-loop at 0,
        // as can occur when Oink solves a bottom SCC of the original graph.
        solver.G.resize(2);
        solver.G[0] = true;
        solver.Q0 = solver.Q1 = solver.G;
        solver.QQ.resize(2);
        solver.T0.resize(2);
        solver.T1.resize(2);
        solver.L0.resize(2);
        solver.L1.resize(2);
        solver.Queue.resize(2);
        solver.escape_offset.resize(2);
        solver.escape_stamp.resize(2);
        std::unique_ptr<Measures> pm0(new_measure(kind, game, 0));
        std::unique_ptr<Measures> pm1(new_measure(kind, game, 1));
        auto& measure = player == 0 ? *pm0 : *pm1;
        for (int wave = 0; wave < 2; ++wave) {
            solver.update(measure, player);
            const auto& pending = player == 0 ? solver.Q0 : solver.Q1;
            const auto& lifted = player == 0 ? solver.L0 : solver.L1;
            if (pending[1] || lifted[1]) throw std::runtime_error("GPM scheduled a disabled predecessor");
        }
        solver.shortcuts(player, *pm0, *pm1);
        if (solver.T0[1] || solver.T1[1] || pm0->is_top(1) || pm1->is_top(1)) {
            throw std::runtime_error("GPM shortcut escaped the active subgame");
        }
    }

    static void check(Game& game, MeasureKind kind)
    {
        std::ostringstream log;
        Oink controller(game, log, true);
        GPMSolver solver(controller, game, kind);
        const int n = game.nodecount();
        solver.G.resize(n);
        solver.G.flip();
        solver.Queue.resize(n);
        solver.Q0.resize(n);
        solver.Q1.resize(n);
        solver.T0.resize(n);
        solver.T1.resize(n);
        solver.L0.resize(n);
        solver.L1.resize(n);
        solver.escape_offset.resize(n);
        solver.escape_stamp.resize(n);
        solver.escape_epoch = std::numeric_limits<unsigned int>::max();

        unsigned expected_epoch = 1;
        for (int player : {0, 0, 1}) {
            // Poison stale entries: wrap and later calls must ignore them.
            std::fill(solver.escape_offset.begin(), solver.escape_offset.end(), 999);
            std::fill(solver.escape_stamp.begin(), solver.escape_stamp.end(),
                      expected_epoch == 1 ? 1 : expected_epoch-1);
            solver.Q0.reset();
            solver.Q1.reset();
            solver.T0.reset();
            solver.T1.reset();
            solver.L0.reset();
            solver.L1.reset();
            auto& target = player == 0 ? solver.L0 : solver.L1;
            target[0] = target[1] = true;
            std::unique_ptr<Measures> pm0(new_measure(kind, game, 0));
            std::unique_ptr<Measures> pm1(new_measure(kind, game, 1));
            // Ordered Bottom changes on the first opposing priority. Give
            // Even a finite value that is already stable under priority 1.
            if (player == 0) {
                for (int v = 0; v < n; ++v) pm0->see(v, 1);
            }
            solver.shortcuts(player, *pm0, *pm1);
            if (solver.escape_epoch != expected_epoch) {
                throw std::runtime_error("shortcut cache did not advance its epoch");
            }
            const auto& other_top = player == 0 ? solver.T1 : solver.T0;
            const auto& other_queue = player == 0 ? solver.Q1 : solver.Q0;
            auto& source = player == 0 ? *pm0 : *pm1;
            auto& other = player == 0 ? *pm1 : *pm0;
            for (int v = 0; v < n; ++v) {
                const bool excluded = player == 0 && (v == 2 || v == 4);
                if (other_top[v] && other_queue[v]) {
                    throw std::runtime_error("new Top vertex remained queued after sibling notifications");
                }
                if (bool(target[v]) != excluded || bool(other_top[v]) != excluded ||
                    other.is_top(v) != excluded || source.is_top(v)) {
                    throw std::runtime_error("shortcut cache changed the attraction closure");
                }
            }
            if (player == 0 && (solver.escape_stamp[2] != expected_epoch || solver.escape_offset[2] != 4)) {
                throw std::runtime_error("expired blocker did not advance past the duplicate edges");
            }
            ++expected_epoch;
        }
    }
};

}

int
test_gpm_cache()
{
    using namespace pg;
    Game game(7);
    for (int v = 6; v >= 0; --v) game.init_vertex(v, 1, v == 2 ? 1 : 0);
    // Vertex 2 first blocks on 3, receives another callback with that blocker
    // live, then resumes past both edges to 3 after the 0<-5<-3 chain grows.
    // Vertex 4 remains a blocker for Even. Its edge to 2 also tests clearing
    // a future Top sibling after an earlier new Top queues it. Odd attracts all.
    const std::vector<std::vector<int>> edges = {{0}, {1}, {0, 1, 3, 3, 4}, {5}, {4, 2}, {0}, {5}};
    for (int v = 0; v < 7; ++v) {
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    game.ensure_sorted();
    for (auto kind : {MeasureKind::Small, MeasureKind::Succinct, MeasureKind::Ordered, MeasureKind::Strahler}) {
        GPMCacheProbe::check(game, kind);
    }
    for (auto kind : {MeasureKind::Small, MeasureKind::Succinct, MeasureKind::Ordered, MeasureKind::Strahler}) {
        for (int player : {0, 1}) {
            for (int priority : {0, 1, 2, 3}) GPMCacheProbe::checkDisabled(kind, player, priority);
        }
    }
    for (auto kind : {MeasureKind::Succinct, MeasureKind::Ordered, MeasureKind::Strahler}) {
        GPMCacheProbe::checkKnownTop(kind);
    }
    std::cout << "GPM updates stay inside the active subgame; shortcut cursors preserve closure across repeated callbacks, player changes and epoch wrap." << std::endl;
    return 0;
}
