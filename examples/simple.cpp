/*
 * Copyright 2024 Tom van Dijk
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

/*
 * Minimal example of using Oink as a library: load a parity game in PGSolver
 * format, solve it, and print the winner (and strategy) of each vertex.
 *
 * Build (in-tree): the CMake build produces the `oink-example-simple` binary.
 * Usage:
 *     oink-example-simple [game.pg] [solver]
 * The game is read from the file argument, or from stdin if none is given.
 * The solver defaults to "tl" (tangle learning); run `oink --solvers` for the
 * full list of solver ids.
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "oink/game.hpp"
#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "oink/version.hpp"

int
main(int argc, char** argv)
{
    // 1. Load a parity game in PGSolver format, from a file argument or stdin.
    pg::Game game;
    if (argc > 1) {
        std::ifstream file(argv[1], std::ios_base::binary);
        if (!file) {
            std::cerr << "could not open " << argv[1] << std::endl;
            return 1;
        }
        game = pg::PGParser::parse_pgsolver(file, true);
    } else {
        game = pg::PGParser::parse_pgsolver(std::cin, true);
    }

    // 2. Solve the game. Oink writes its progress to the stream you pass it; we
    //    send that to a stringstream so this program's own output stays clean.
    //    Pass std::cout instead to watch the solver work.
    const std::string solver = argc > 2 ? argv[2] : "tl";
    std::stringstream solver_log;
    pg::Oink oink(game, solver_log);
    oink.setSolver(solver);
    oink.run();

    // 3. Read the result. Oink owns the solution; after solving, every vertex is
    //    won by player 0 (Even) or player 1 (Odd). strategy(v) gives the move to
    //    play when the owner of v is the winner, or -1 otherwise.
    const pg::Solution& sol = oink.solution();
    long won0 = 0, won1 = 0;
    for (int v = 0; v < game.vertexcount(); v++) {
        if (!sol.is_solved(v)) continue;
        if (sol.winner(v) == 0) won0++;
        else won1++;
    }

    std::cout << "oink " << pg::version() << ": solved " << game.vertexcount()
              << " vertices with solver \"" << solver << "\"" << std::endl;
    std::cout << won0 << " won by Even (player 0), "
              << won1 << " won by Odd (player 1)" << std::endl;

    // For small games, show the per-vertex outcome and winning strategy.
    if (game.vertexcount() <= 20) {
        for (int v = 0; v < game.vertexcount(); v++) {
            if (!sol.is_solved(v)) continue;
            const int strategy = sol.strategy(v);
            std::cout << "  vertex " << v
                      << " (owner " << game.owner(v)
                      << ", priority " << game.priority(v) << ")"
                      << " won by player " << sol.winner(v);
            if (strategy != -1) std::cout << ", strategy -> " << strategy;
            std::cout << std::endl;
        }
    }

    return 0;
}
