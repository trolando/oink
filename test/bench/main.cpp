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

#include <csignal>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>

#include "bench.hpp"
#include "experiment.hpp"
#include "options.hpp"
#include "report.hpp"
#include "runner.hpp"
#include "setup.hpp"

int
main(int argc, char** argv)
{
    using namespace bench;
    try {
        Options options;
        std::vector<SolverRun> solvers;
        if (!parse_command_line(argc, argv, options, solvers, std::cout, &std::cerr)) return 0;
        const std::string executable = executable_path(argv[0]);
        volatile std::sig_atomic_t cancelled = 0;

        // an existing experiment is resumed with its own plan
        std::unique_ptr<Experiment> experiment;
        std::vector<PlannedGame> games;
        if (!options.resume && !options.experiment.empty() && std::filesystem::exists(options.experiment)) {
            if (!confirm_resume(options.experiment, &cancelled)) return cancelled ? 130 : 0;
            options.resume = true;
        }
        if (options.resume) {
            experiment = Experiment::resume(options.experiment, options, solvers, games, executable);
        } else {
            if (options.interactive && !setup_benchmark(options, solvers, &cancelled)) return cancelled ? 130 : 0;
            validate(options, solvers);
            games = plan_games(options, executable);
            if (!options.experiment.empty()) experiment = Experiment::create(options.experiment, options, solvers, games, executable);
        }
        if (experiment) experiment->configure(options);
        validate(options, solvers);
        return run_benchmark(options, solvers, games, executable, experiment.get());
    } catch (const std::exception& error) {
        std::cerr << "bench_solvers: " << error.what() << std::endl;
        return 3;
    }
}
