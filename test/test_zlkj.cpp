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
 * Unit tests for ZLKJ: its registry entry, and all small games.
 */

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "verifier.hpp"

using namespace pg;

int
test_zlkj()
{
    if (Solvers::family("zlkj") != Solvers::Family::Zielonka ||
        Solvers::isExperimental("zlkj") || Solvers::hasOptions("zlkj")) {
        throw std::runtime_error("incorrect ZLKJ registry metadata");
    }

    // Exhaust every total graph, ownership and priority assignment up to three
    // vertices. Four priorities exercise both players and same-parity levels.
    // Verify both complete winning strategies without preprocessing shortcuts.
    unsigned cases = 0;
    for (int n = 1; n <= 3; ++n) {
        const int masks = (1<<n)-1;
        int graphs = 1;
        int priorities = 1;
        for (int v = 0; v < n; ++v) {
            graphs *= masks;
            priorities *= 4;
        }
        for (int graph = 0; graph < graphs; ++graph) {
            for (int owners = 0; owners < (1<<n); ++owners) {
                for (int prios = 0; prios < priorities; ++prios) {
                    Game game(n);
                    int g = graph;
                    int p = prios;
                    for (int v = 0; v < n; ++v) {
                        game.init_vertex(v, p%4, (owners>>v)&1);
                        p /= 4;
                        const int successors = g%masks+1;
                        g /= masks;
                        game.e_start(v);
                        for (int w = 0; w < n; ++w) {
                            if (successors & (1<<w)) game.e_add(v, w);
                        }
                        game.e_finish();
                    }
                    game.ensure_sorted();
                    std::ostringstream log;
                    Oink controller(game, log, true);
                    controller.setRemoveLoops(false);
                    controller.setRemoveWCWC(false);
                    controller.setSolveSingle(false);
                    controller.setWorkers(-1);
                    controller.setSolver("zlkj");
                    controller.run();
                    Verifier(game, controller.solution(), log).verify(true, true, true);
                    ++cases;
                }
            }
        }
    }
    std::cout << "Verified ZLKJ strategies on " << cases << " exhaustive small games." << std::endl;
    return 0;
}
