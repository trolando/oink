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

#include "display.hpp"

#include <iostream>

#include "common/terminal.hpp"

namespace bench {

PlainDisplay::PlainDisplay(const Options& options, const std::vector<SolverRun>& solvers) : options(options), solvers(solvers)
{
    std::cerr << "Oink benchmark | " << options.repeats << (options.calibrated() ? "+ calibrated" : "") << " measured, ";
    std::cerr << options.warmups << " warm-up rounds | order seed " << options.seed << '\n';
    std::cerr << "game";
    for (const auto& solver : solvers) std::cerr << '\t' << ui::clean(solver.name);
    std::cerr << std::endl;
}

void
PlainDisplay::completed(const Case& game, size_t)
{
    std::cerr << ui::clean(game.name);
    for (size_t s = 0; s < solvers.size(); s++) {
        const Cell& cell = game.cells[s];
        std::cerr << '\t';
        if (game.fully_preprocessed) {
            std::cerr << "preprocessed";
        } else if (game.input_failed) {
            std::cerr << "input_error";
        } else if (game.prep_failed) {
            std::cerr << "preparation_error";
        } else if (cell.failure) {
            std::cerr << status_name(*cell.failure);
        } else if (cell.complete(options.repeats)) {
            const Stats& st = cell.stats();
            std::cerr << duration(st.median) << " [IQR " << duration(st.q3 - st.q1) << ']';
        } else {
            std::cerr << "incomplete";
        }
    }
    std::cerr << std::endl;
    if (!game.error.empty()) std::cerr << "  " << ui::clean(game.error) << std::endl;
}

}
