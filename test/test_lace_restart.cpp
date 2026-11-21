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
 * Unit test for parallel solving with Oink::run more than once in the same process: each run
 * starts Lace and stops it again, and Lace must not be seen as running after it stopped.
 */

#include <iostream>
#include <sstream>

#include "lace.h"
#include "oink/oink.hpp"
#include "verifier.hpp"

using namespace pg;

namespace {

/**
 * A ring of <n> vertices with priorities 0..n-1, alternating owners.
 */
Game
make_ring(int n)
{
    Game game(n);
    for (int v = 0; v < n; v++) {
        game.init_vertex(v, v, v & 1);
        game.e_start(v);
        game.e_add(v, (v+1) % n);
        game.e_finish();
    }
    game.ensure_sorted();
    return game;
}

}

int
test_lace_restart()
{
    for (const char* solver : {"zlk", "psi", "tlpar"}) {
        for (int round = 0; round < 2; round++) {
            Game game = make_ring(100);
            std::ostringstream log;
            Oink controller(game, log);
            controller.setWorkers(2);
            controller.setSolver(solver);
            controller.run();
            Verifier(game, controller.solution(), log).verify(true, true, true);
            if (lace_is_running()) {
                std::cerr << solver << ": Lace still running after run " << round << std::endl;
                return 1;
            }
            if (log.str().find("initialized Lace") == std::string::npos) {
                std::cerr << solver << ": run " << round << " did not start Lace" << std::endl;
                return 1;
            }
        }
    }
    std::cout << "Lace restart test passed" << std::endl;
    return 0;
}
