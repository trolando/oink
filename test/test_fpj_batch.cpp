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
 * Unit tests for FPJ: batch resets, restarts with tracked changes, and all small games.
 */

#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "solvers/fpj.hpp"
#include "verifier.hpp"

using namespace pg;

static void
check_batch()
{
    for (int invert : {0, 1}) {
        Game game(3);
        const int priorities[] = {0, 1, 3};
        for (int v = 0; v < 3; ++v) {
            game.init_vertex(v, priorities[v]+invert, invert);
            game.e_start(v);
            // 1 -> 2 -> 0: the descending scan evaluates 2 before 1, as one batch
            game.e_add(v, v == 1 ? 2 : 0);
            game.e_finish();
        }
        game.ensure_sorted();
        std::ostringstream log;
        Oink controller(game, log, true);
        FPJSolver solver(controller, game);
        solver.run();
        Verifier(game, controller.solution(), log).verify(true, true, true);
        if (solver.iterations != 1) throw std::runtime_error("new batch justifications were reset");
    }
}

/**
 * A forwards chain (changes propagate to later vertices, against the descending scan of
 * a block) exercises repeated evaluation and activates the frontier
 * when sparse evaluation is enabled. The initially attractive higher-priority
 * parameter later changes winner; its dependents must be reconsidered,
 * including across word boundaries.
 */
static void
check_restarts(int invert)
{
    constexpr int n = 256;
    std::vector<int> order(n), successor(n+1);
    std::iota(order.begin(), order.end(), 1);
    for (int i = 8; i+1 < n; i += 16) std::swap(order[i], order[i+1]);
    for (int i = 62; i+1 < n; i += 64) std::swap(order[i], order[i+1]);
    for (int i = 0; i < n; ++i) successor[order[i]] = i > 0 ? order[i-1] : 0;
    std::ostringstream text;
    text << "parity " << n+2 << ";\n0 " << invert << ' ' << invert << " 0;\n";
    for (int v = 1; v <= n; ++v) {
        const int player = v == n/2 ? invert : ((v&1)^invert);
        text << v << ' ' << 1+invert << ' ' << player << ' ';
        if (v == n/2) text << n+1 << ',';
        text << successor[v] << ";\n";
    }
    text << n+1 << ' ' << 2+invert << ' ' << (1^invert) << ' ' << n+2 << ";\n";
    text << n+2 << ' ' << 3+invert << ' ' << (1^invert) << ' ' << n+2 << ";\n";
    std::istringstream input(text.str());
    auto game = PGParser::parse_pgsolver(input, false);
    game.ensure_sorted();
    std::ostringstream log;
    Oink controller(game, log, true);
    controller.setRemoveLoops(false);
    controller.setRemoveWCWC(false);
    controller.setSolveSingle(false);
    controller.setWorkers(-1);
    controller.setSolver("fpj");
    controller.run();
    Verifier(game, controller.solution(), log).verify(true, true, true);
    for (int v = 0; v < game.nodecount(); ++v) {
        const int expected = game.priority(v) < 2+invert ? invert : (1^invert);
        if (controller.solution().winner(v) != expected) {
            throw std::runtime_error("FPJ frontier missed a reset or hypothesis change");
        }
    }
}

int
test_fpj_batch()
{
    check_batch();
    for (int invert : {0, 1}) check_restarts(invert);

    // Exhaust all total graphs, ownerships and priority assignments on at most
    // three vertices. Verify the complete strategies without preprocessing.
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
                    controller.setSolver("fpj");
                    controller.run();
                    Verifier(game, controller.solution(), log).verify(true, true, true);
                    ++cases;
                }
            }
        }
    }
    std::cout << "Verified FPJ on " << cases << " exhaustive total games." << std::endl;
    return 0;
}
