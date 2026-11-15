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
 * bench_solvers: benchmark parity game solvers on files or random games.
 * See test/bench/README.md and bench_solvers --help.
 */

#include <exception>
#include <iostream>

#include "bench.hpp"
#include "options.hpp"
#include "report.hpp"
#include "runner.hpp"

int
main(int argc, char** argv)
{
    using namespace bench;
    try {
        Options options;
        std::vector<SolverRun> solvers;
        if (!parse_command_line(argc, argv, options, solvers, std::cout)) return 0;
        validate(options, solvers);
        const std::string executable = executable_path(argv[0]);
        const auto games = plan_games(options, executable);
        return run_benchmark(options, solvers, games, executable);
    } catch (const std::exception& error) {
        std::cerr << "bench_solvers: " << error.what() << std::endl;
        return 3;
    }
}
