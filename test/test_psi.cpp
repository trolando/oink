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
 * Unit tests for PSI: vertices that are disabled when the solver starts (here from a
 * partial solution) must not take part in its halting scans. Their state used to be left
 * uninitialised, which shows when new memory is not zero: CTest runs this test with the
 * allocator filling new memory (MallocScribble on macOS, MALLOC_PERTURB_ with glibc; the
 * address sanitizer does so by default).
 */

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "oink/oink.hpp"
#include "verifier.hpp"

using namespace pg;

namespace {

/**
 * A two-vertex cycle 0 <-> 1 of priority 1, won by Odd, followed by <solved> vertices of
 * priority 2 that are Even's, each with a self-loop. Even never stops halting in front of
 * the cycle, so PSI solves the cycle in one major iteration, unless it also scans the
 * solved vertices.
 */
Game
make_game(int solved)
{
    Game game(2+solved);
    for (int v = 0; v < 2+solved; v++) {
        game.init_vertex(v, v < 2 ? 1 : 2, v == 0 ? 1 : 0);
        game.e_start(v);
        game.e_add(v, v < 2 ? 1-v : v);
        game.e_finish();
    }
    game.ensure_sorted();
    return game;
}

/**
 * Solve with PSI and return its "solved with N major iterations, M minor iterations."
 */
std::string
solve(Game& game, int solved, int workers)
{
    std::ostringstream log;
    Oink controller(game, log);
    controller.setRemoveLoops(false);
    controller.setRemoveWCWC(false);
    controller.setSolveSingle(false);
    controller.setWorkers(workers);
    controller.setSolver("psi");
    for (int v = 2; v < 2+solved; v++) controller.solution().solve(v, 0, v);
    controller.run();
    Verifier(game, controller.solution(), log).verify(true, true, true);
    const std::string text = log.str();
    const size_t at = text.find("solved with ");
    if (at == std::string::npos) throw std::runtime_error("PSI did not report its iterations");
    return text.substr(at, text.find('\n', at)-at);
}

}

int
test_psi()
{
    int failures = 0;
    for (int workers : {-1, 4}) {
        Game reference = make_game(0);
        const std::string expected = solve(reference, 0, workers);
        Game game = make_game(30);
        const std::string actual = solve(game, 30, workers);
        if (actual != expected) {
            std::cerr << "PSI with " << workers << " workers and a partial solution: " << actual
                      << " (without the solved vertices: " << expected << ")" << std::endl;
            failures++;
        }
    }
    if (failures) {
        std::cerr << failures << " PSI test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all PSI tests passed" << std::endl;
    return 0;
}
