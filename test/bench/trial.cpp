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

#include "trial.hpp"

#include <chrono>
#include <ctime>
#include <sstream>
#include <stdexcept>

#include "lace.h"
#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "verifier.hpp"

namespace bench {

using Clock = std::chrono::steady_clock;

static double
seconds_since(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

// The longest captured solver log kept with a failure.
static const size_t MAX_LOG = 262144;

std::unique_ptr<Prepared>
prepare(const pg::Game& original, const Options& options, std::ostream& log)
{
    (void)log;
    if (options.preprocess) throw std::runtime_error("preprocessing is not supported");
    auto prepared = std::make_unique<Prepared>();
    prepared->original = &original;
    prepared->residual = &original;
    prepared->partial = pg::Solution(original.nodecount());
    prepared->info.residual_vertices = original.nodecount();
    prepared->info.residual_edges = original.edgecount();
    return prepared;
}

/** Add the captured solver log to the message of a failed trial. */
static void
add_log(Trial& trial, const std::ostringstream& log)
{
    std::string text = log.str();
    if (text.empty()) return;
    if (text.size() > MAX_LOG) text = "[...]\n" + text.substr(text.size() - MAX_LOG);
    trial.message += "\n" + text;
}

Trial
execute(const Prepared& prepared, const SolverRun& run, const Options& options, uint32_t seed)
{
    Trial result;
    // The solver log is captured only with --trace; otherwise it is discarded.
    std::ostringstream captured;
    std::ostream muted(nullptr);
    std::ostream& log = options.trace >= 0 ? static_cast<std::ostream&>(captured) : muted;

    auto start = Clock::now();
    try {
        // a fresh copy, with the in-edges built by Oink's constructor: not timed
        pg::Game game(*prepared.residual);
        game.set_random_seed(seed);
        pg::Oink solver(game, log, pg::Solvers::needsInEdges(run.id));
        // the trials solve the residual game as it is, without preprocessing
        solver.setSolveSingle(false);
        solver.setRemoveLoops(false);
        solver.setRemoveWCWC(false);
        solver.setBottomSCC(options.scc);
        solver.setWorkers(options.workers);
        solver.setForceParallel(options.force_parallel);
        solver.setSolver(run.id);
        solver.setTrace(options.trace < 0 ? 0 : options.trace);
        if (!run.config.empty()) solver.setSolverOptions(run.config);
        result.setup = seconds_since(start);

        // the measurement: Oink::run, including SCC decomposition and lazy in-edges
        const std::clock_t cpu_start = std::clock();
        start = Clock::now();
        try {
            solver.run();
        } catch (const std::exception& error) {
            result.solve = seconds_since(start);
            result.status = Status::SolverError;
            result.message = error.what();
            add_log(result, captured);
            return result;
        }
        result.solve = seconds_since(start);
        const std::clock_t cpu_end = std::clock();
        if (cpu_start != (std::clock_t)-1 && cpu_end != (std::clock_t)-1) {
            result.cpu = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        }
        result.effective_solver = solver.effectiveSolver().value_or("none");
        result.effective_config = solver.effectiveSolver() ? solver.effectiveSolverOptions() : "";
        result.parallel = solver.parallel();
        result.workers = options.workers < 0 ? -1 : (int)lace_worker_count();

        // verify every strategy of the solution, outside the measurement
        start = Clock::now();
        try {
            std::ostringstream verify_log;
            pg::Verifier verifier(game, solver.solution(), verify_log);
            verifier.verify(true, true, true);
        } catch (const std::exception& error) {
            result.verify = seconds_since(start);
            result.status = Status::InvalidSolution;
            result.message = error.what();
            add_log(result, captured);
            return result;
        }
        result.verify = seconds_since(start);
        result.status = Status::Ok;
    } catch (const std::exception& error) {
        result.setup = seconds_since(start);
        result.status = Status::SetupError;
        result.message = error.what();
        add_log(result, captured);
    }
    return result;
}

}
