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
 * Unit tests for PMA: it must solve as PMTL without tangle learning.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "oink/solvers.hpp"
#include "verifier.hpp"

using namespace pg;

struct Result {
    std::vector<std::pair<int, int>> vertices;
    std::string summary;
};

/**
 * Solve a copy of <input> without preprocessing; return the solution and the statistics.
 */
static Result
solve(const Game& input, const std::string& id, const std::string& options)
{
    Game game(input);
    std::ostringstream log;
    Oink controller(game, log, Solvers::needsInEdges(id));
    controller.setRemoveLoops(false);
    controller.setRemoveWCWC(false);
    controller.setSolveSingle(false);
    controller.setWorkers(-1);
    controller.setSolver(id);
    controller.setSolverOptions(options);
    controller.run();
    Verifier(game, controller.solution(), log).verify(true, true, true);
    Result result;
    for (int v = 0; v < game.nodecount(); ++v) {
        result.vertices.emplace_back(controller.solution().winner(v), controller.solution().strategy(v));
    }
    const auto text = log.str();
    const auto start = text.find("found ");
    if (start == std::string::npos) throw std::runtime_error("missing PMTL statistics");
    const auto last = text.find("number of lifts:", start);
    if (last == std::string::npos) throw std::runtime_error("missing lift statistics");
    result.summary = text.substr(start, text.find('\n', last)-start);
    return result;
}

static void
compare(const Game& game, const std::string& pma, const std::string& pmtl, const std::string& options)
{
    const auto a = solve(game, pma, options);
    const auto b = solve(game, pmtl, "notangles,"+options);
    if (a.vertices != b.vertices || a.summary != b.summary) {
        throw std::runtime_error(pma+" differs from "+pmtl+"=notangles,"+options);
    }
    if (a.summary.find("and 0 tangles.") == std::string::npos) {
        throw std::runtime_error(pma+" learned a tangle");
    }
}

int
test_pma(const char* fixture)
{
    try {
        std::ifstream input(fixture);
        Game game = PGParser::parse_pgsolver(input, false);
        game.ensure_sorted();
        // Two disjoint winning cycles can be certified as closed regions
        // without SCC extraction, including when tangle learning is disabled.
        std::istringstream cycles("parity 3;\n0 0 0 1;\n1 2 1 0;\n2 1 1 3;\n3 3 0 2;\n");
        Game closed = PGParser::parse_pgsolver(cycles, false);
        closed.ensure_sorted();
        bool exercised_learning = false;
        for (const std::string kind : {"small", "succinct", "ordered", "strahler"}) {
            const auto pma = "pma-"+kind;
            const auto pmtl = "pmtl-"+kind;
            if (Solvers::family(pma) != Solvers::Family::ProgressMeasures || !Solvers::hasOptions(pma)) {
                throw std::runtime_error("missing PMA registry metadata");
            }
            for (const std::string options : {"", "up", "no-shortcuts", "up,no-shortcuts"}) {
                std::ostringstream errors;
                if (!Solvers::validateOptions(pma, options, errors)) {
                    throw std::runtime_error(errors.str());
                }
                compare(game, pma, pmtl, options);
                compare(closed, pma, pmtl, options);
                const auto result = solve(closed, pma, options);
                if (result.summary.find("found 0 dominions") != std::string::npos) {
                    throw std::runtime_error(pma+" failed to recognize closed-region dominions");
                }
            }
            const auto learning = solve(game, pmtl, "");
            exercised_learning |= learning.summary.find("and 0 tangles.") == std::string::npos;
            const auto override = solve(game, pma, "tangles");
            if (override.vertices != learning.vertices || override.summary != learning.summary) {
                throw std::runtime_error("explicit tangle configuration differs");
            }
        }
        if (!exercised_learning) throw std::runtime_error("fixture did not exercise tangle learning");
        compare(game, "pma", "pmtl-ordered", "");
        compare(game, "pma-small", "pmtl-ordered", "ordered");
        std::cout << "PMA solvers match PMTL without tangles, including configuration overrides." << std::endl;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}
