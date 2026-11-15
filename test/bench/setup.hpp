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
 * The setup screen of bench_solvers -i: the shared solver chooser, with the
 * settings of the benchmark (r repetitions, w warm-ups, t timeout, b schedule
 * seed, p sampling target, P preprocessing) and R for the reference solver.
 */

#ifndef OINK_BENCH_SETUP_HPP
#define OINK_BENCH_SETUP_HPP

#include <csignal>
#include <vector>

#include "bench.hpp"
#include "common/chooser.hpp"

namespace bench {

/** The chooser configuration whose settings edit <options> directly. */
ui::ChooserConfig setup_config(Options& options);

/**
 * Show the setup screen, starting from <solvers>; on start, the chosen roster
 * and reference are stored. Returns false if the user cancelled; SIGINT
 * cancels too, and sets <*quit>.
 */
bool setup_benchmark(Options& options, std::vector<SolverRun>& solvers, volatile std::sig_atomic_t* quit);

}

#endif
