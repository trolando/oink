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
#include <limits>
#include <map>
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

/**
 * The maps of the vertices and edges of the residual game to the original
 * game. Edges are matched by their targets, in order, so parallel edges and
 * any reordering of the successors are handled.
 */
static void
map_edges(Prepared& prepared)
{
    const pg::Game& original = *prepared.original;
    const pg::Game& residual = *prepared.residual;
    std::vector<int> inverse(original.nodecount(), -1);
    for (int v = 0; v < residual.nodecount(); v++) inverse[prepared.vertices[v]] = v;
    prepared.edges.assign((size_t)residual.edgeArraySize(), std::numeric_limits<size_t>::max());
    std::map<int, std::vector<size_t>> pending; // residual target -> original edges, in order
    for (int v = 0; v < residual.nodecount(); v++) {
        const int ov = prepared.vertices[v];
        if (original.owner(ov) != residual.owner(v) || (original.priority(ov) & 1) != (residual.priority(v) & 1)) {
            throw std::runtime_error("the residual game differs from the original in an owner or parity");
        }
        pending.clear();
        for (int j = 0; j < original.outcount(ov); j++) {
            const int target = original.outs(ov)[j];
            if (inverse[target] >= 0) pending[inverse[target]].push_back((size_t)original.firstout(ov) + (size_t)j);
        }
        std::map<int, size_t> used;
        for (int k = 0; k < residual.outcount(v); k++) {
            const int target = residual.outs(v)[k];
            auto& edges = pending[target];
            const size_t next = used[target]++;
            if (next >= edges.size()) throw std::runtime_error("the residual game has an edge that the original does not have");
            prepared.edges[(size_t)residual.firstout(v) + (size_t)k] = edges[next];
        }
        for (const auto& entry : pending) {
            if (used[entry.first] != entry.second.size()) throw std::runtime_error("the residual game lost an edge of the original");
        }
        if (residual.outcount(v) == 0) throw std::runtime_error("the residual game has a dead end");
    }
}

std::unique_ptr<Prepared>
prepare(const pg::Game& original, const Options& options, std::ostream& log)
{
    auto prepared = std::make_unique<Prepared>();
    prepared->original = &original;
    const int n = original.nodecount();
    prepared->partial = pg::Solution(n);
    if (!options.preprocess) {
        // the solvers solve the original game itself
        prepared->residual = &original;
        prepared->info.residual_vertices = n;
        prepared->info.residual_edges = original.edgecount();
        return prepared;
    }

    // preprocess a copy: Oink without a solver only preprocesses
    pg::Game working(original);
    const auto start = Clock::now();
    {
        pg::Oink preprocessor(working, log);
        if (options.inflate) {
            preprocessor.setInflate();
        } else if (options.compress) {
            preprocessor.setCompress();
        } else {
            preprocessor.setRenumber();
        }
        preprocessor.setSolveSingle(options.single);
        preprocessor.setRemoveLoops(options.loops);
        preprocessor.setRemoveWCWC(options.wcwc);
        preprocessor.setTrace(options.trace < 0 ? 0 : options.trace);
        preprocessor.run();
        prepared->partial = preprocessor.solution();
    }
    prepared->info.seconds = seconds_since(start);
    const pg::Solution& partial = prepared->partial;
    if (partial.has_multi()) throw std::runtime_error("preprocessing returned a multi-strategy");

    // what preprocessing solved must be a correct partial solution of the original
    pg::Verifier(original, partial, log).verify(false, true, true);

    const int solved = (int)partial.solved().count();
    prepared->info.solved_vertices = solved;
    if (solved == n) return prepared; // no residual: nothing for the solvers

    pg::bitset mask(n);
    for (int v = 0; v < n; v++) {
        if (!partial.is_solved(v)) mask.set(v);
    }
    // the unsolved part must be closed: no unsolved vertex can move to a vertex its owner won
    for (int v = 0; v < n; v++) {
        if (!mask[v]) continue;
        for (const int* e = working.outs(v); *e != -1; e++) {
            if (!mask[*e] && partial.winner(*e) == working.owner(v)) throw std::runtime_error("the preprocessed part is not closed under attraction");
        }
    }

    if (solved == 0) {
        prepared->owned = std::make_unique<pg::Game>(working);
        prepared->vertices.resize(n);
        for (int v = 0; v < n; v++) prepared->vertices[v] = v;
    } else {
        std::vector<int> mapping;
        prepared->owned = working.extract_subgame(mask, mapping);
        // sort the residual game; mapping[v] is then the original of vertex v
        std::vector<int> order(prepared->owned->nodecount());
        prepared->owned->sort(order.data());
        prepared->vertices.resize(order.size());
        for (size_t v = 0; v < order.size(); v++) prepared->vertices[v] = mapping[order[v]];
    }
    prepared->residual = prepared->owned.get();
    map_edges(*prepared);
    prepared->info.residual_vertices = prepared->residual->nodecount();
    prepared->info.residual_edges = prepared->residual->edgecount();
    prepared->info.residual_hash = hex64(fingerprint(*prepared->residual));
    return prepared;
}

pg::Solution
lift(const Prepared& prepared, const pg::Game& trial, const pg::Solution& solution)
{
    const pg::Game& residual = *prepared.residual;
    const int n = residual.nodecount();
    if (solution.vertex_count() != n || (int)solution.solved().count() != n) throw std::runtime_error("the solver did not solve every vertex of the residual game");
    pg::Solution full = prepared.partial;
    if (solution.has_multi()) full.init_multi((size_t)prepared.original->edgeArraySize());
    for (int v = 0; v < n; v++) {
        if (trial.priority(v) != residual.priority(v) || trial.owner(v) != residual.owner(v) || trial.firstout(v) != residual.firstout(v) || trial.outcount(v) != residual.outcount(v)) {
            throw std::runtime_error("the solver changed the residual game");
        }
        const int target = solution.strategy(v);
        if (target < -1 || target >= n) throw std::runtime_error("invalid strategy in the residual game");
        full.solve(prepared.vertices[v], solution.winner(v), target < 0 ? -1 : prepared.vertices[target]);
        if (!solution.has_multi()) continue;
        for (int k = 0; k < trial.outcount(v); k++) {
            if (trial.outs(v)[k] != residual.outs(v)[k]) throw std::runtime_error("the solver changed the residual game");
            const size_t edge = (size_t)trial.firstout(v) + (size_t)k;
            if (solution.has_edge_index(edge)) full.add_edge_index(prepared.edges[edge]);
        }
    }
    return full;
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
execute(const Prepared& prepared, const SolverRun& run, const Options& options, uint32_t seed, const PhaseUpdate& update)
{
    Trial result;
    Phase phase = Phase::Setup;
    // The solver log is captured only with --trace; otherwise it is discarded.
    std::ostringstream captured;
    std::ostream muted(nullptr);
    std::ostream& log = options.trace >= 0 ? static_cast<std::ostream&>(captured) : muted;

    auto start = Clock::now();
    try {
        if (update) update(phase, result);
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
        phase = Phase::Solve;
        if (update) update(phase, result);
        const std::clock_t cpu_start = std::clock();
        start = Clock::now();
        solver.run();
        result.solve = seconds_since(start);
        const std::clock_t cpu_end = std::clock();
        if (cpu_start != (std::clock_t)-1 && cpu_end != (std::clock_t)-1) {
            result.cpu = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        }
        phase = Phase::Verify;
        start = Clock::now();
        if (update) update(phase, result);
        result.effective_solver = solver.effectiveSolver().value_or("none");
        result.effective_config = solver.effectiveSolver() ? solver.effectiveSolverOptions() : "";
        result.parallel = solver.parallel();
        result.workers = options.workers < 0 ? -1 : (int)lace_worker_count();

        // verify every strategy of the solution, outside the measurement
        std::ostringstream verify_log;
        if (prepared.identity()) {
            pg::Verifier verifier(game, solver.solution(), verify_log);
            verifier.verify(true, true, true);
        } else {
            // against the original game, with what preprocessing solved
            const pg::Solution full = lift(prepared, game, solver.solution());
            pg::Verifier verifier(*prepared.original, full, verify_log);
            verifier.verify(true, true, true);
        }
        result.verify = seconds_since(start);
        result.status = Status::Ok;
    } catch (const std::exception& error) {
        const double elapsed = seconds_since(start);
        if (phase == Phase::Setup) result.setup = elapsed;
        else if (phase == Phase::Solve) result.solve = elapsed;
        else if (result.verify < 0) result.verify = elapsed;
        record_failure(result, phase, error);
        add_log(result, captured);
    }
    return result;
}

}
