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
 * The benchmark loop of bench_solvers.
 *
 * The games are read or generated one at a time, in the order of the plan
 * (shuffled with the schedule seed unless --sort). For every game:
 *
 *   1. it is read, parsed and sorted (outside any timing) and prepared once;
 *   2. the rounds run: warm-ups, three pilot trials when calibrating, and the
 *      measured repetitions; every round runs the solvers in a shuffled order,
 *      each trial on a fresh copy of the game, verified;
 *   3. its statistics are reported and shown.
 *
 * A failed trial (solver error, invalid solution, timeout, crash) fails the
 * cell of that solver on that game; its later trials are reported as skipped.
 */

#ifndef OINK_BENCH_RUNNER_HPP
#define OINK_BENCH_RUNNER_HPP

#include <optional>
#include <string>
#include <vector>

#include "bench.hpp"

namespace bench {

/** One game of the plan: a file, or the seed of a random game. */
struct PlannedGame
{
    std::string path;
    std::optional<uint32_t> seed;
};

/**
 * The games of a benchmark: the parity game files among <options.inputs>
 * (canonical, without duplicates, sorted and then shuffled with the schedule
 * seed unless --sort), or the seeds of the random games.
 */
std::vector<PlannedGame> plan_games(const Options& options, const std::string& executable);

class Experiment;

/**
 * Run the benchmark of <solvers> on <games>; returns the exit status. The
 * executable is reported in the metadata. With a named <experiment>, its
 * finished games are not measured again, and every game that finishes is
 * committed to it.
 */
int run_benchmark(Options& options, std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& games, const std::string& executable, Experiment* experiment = nullptr);

}

#endif
