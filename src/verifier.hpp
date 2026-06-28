/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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

#ifndef VERIFIER_HPP
#define VERIFIER_HPP

#include <istream>
#include "oink/game.hpp"
#include "oink/solution.hpp"

namespace pg {

class Verifier
{
public:
    // The game must already be sorted by priority (call game.ensure_sorted()
    // before constructing the Verifier); verification reads the game (structure)
    // and the solution (solved/winner/strategy).
    Verifier(const Game& game, const Solution& solution, std::ostream &logger)
        : game(game), solution(solution), logger(logger) { }

    // TODO: make a std::exception class for verification exceptions?

    /**
     * Verify the game strategy.
     */
    void verify(bool fullgame=true, bool even=true, bool odd=true);

    /**
     * Return the number of checked strategies in the game.
     */
    int numberOfStrategies(void) { return n_strategies; }

protected:
    const Game& game;
    const Solution& solution;
    std::ostream &logger;
    int n_strategies = 0;

    // winner of a solved vertex (0 or 1), or -1 if unsolved
    [[nodiscard]] int getWinner(int v) const { return solution.is_solved(v) ? solution.winner(v) : -1; }
};

}

#endif 
