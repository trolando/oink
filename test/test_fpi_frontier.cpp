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
 * Exercise winner propagation and resets in FPI and FPIM.
 */

#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lace.h"
#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "verifier.hpp"

static void
check(const char* solver, int invert, int workers)
{
    constexpr int n = 256;
    const int reset_count = workers > 0 ? n : 1;
    const int last = n + reset_count + 1;
    std::vector<int> order(n), successor(n+1);
    std::iota(order.begin(), order.end(), 1);
    // Most changes propagate to a later vertex, which the descending scan of a block
    // only reaches in the next pass. Swapped neighbors exercise newly affected lower
    // lanes and lower words.
    for (int i = 8; i+1 < n; i += 16) std::swap(order[i], order[i+1]);
    for (int i = 62; i+1 < n; i += 64) std::swap(order[i], order[i+1]);
    for (int i = 0; i < n; ++i) successor[order[i]] = i > 0 ? order[i-1] : 0;
    std::ostringstream text;
    text << "parity " << last << ";\n";
    text << "0 " << invert << ' ' << invert << " 0;\n";
    for (int v = 1; v <= n; ++v) {
        text << v << ' ' << 1+invert << ' ' << ((v&1)^invert) << ' ' << successor[v] << ";\n";
    }
    // Once the chain has converged, these distractions reset all its winners.
    // The chain must then converge again, notifying predecessors on both flips.
    // In parallel this also exercises a second priority block large enough to split.
    for (int v = n+1; v < last; ++v) {
        text << v << ' ' << 2+invert << ' ' << (1^invert) << ' ' << last << ";\n";
    }
    text << last << ' ' << 3+invert << ' ' << (1^invert) << ' ' << last << ";\n";
    std::istringstream input(text.str());
    auto game = pg::PGParser::parse_pgsolver(input, false);
    game.ensure_sorted();
    std::ostringstream log;
    pg::Oink controller(game, log, false);
    controller.setRemoveLoops(false);
    controller.setRemoveWCWC(false);
    controller.setSolveSingle(false);
    controller.setWorkers(workers);
    controller.setForceParallel(workers > 0);
    controller.setSolver(solver);
    controller.run();
    if (workers > 0 && !controller.parallel()) {
        throw std::runtime_error("FPI did not run in parallel");
    }
    pg::Verifier(game, controller.solution(), log).verify(true, true, true);
    for (int v = 0; v < game.nodecount(); ++v) {
        const int expected = v <= n ? invert : (1^invert);
        if (controller.solution().winner(v) != expected) {
            throw std::runtime_error(std::string(solver) + ": incorrect chain/reset winner");
        }
    }
}

int
test_fpi_frontier()
{
    try {
        for (const char* solver : {"fpi", "fpim"}) {
            for (int invert : {0, 1}) check(solver, invert, -1);
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}

int
test_fpi_parallel()
{
    constexpr int workers = 4;
    lace_start(workers, 100000);
    int result = 0;
    try {
        for (int invert : {0, 1}) check("fpi", invert, workers);
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        result = 1;
    }
    lace_stop();
    return result;
}
