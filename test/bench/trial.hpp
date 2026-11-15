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
 * Running one trial: one solver on a fresh copy of the game, with the solve
 * time measured around Oink::run, followed by the verification of the
 * solution against the original game.
 *
 * The solvers solve the residual game of the once-per-game preparation. With
 * preprocessing off, the residual is the original game itself. With
 * preprocessing on, it is the part that preprocessing left unsolved, and the
 * solution of a trial is mapped back to the original game for verification.
 */

#ifndef OINK_BENCH_TRIAL_HPP
#define OINK_BENCH_TRIAL_HPP

#include <memory>
#include <ostream>
#include <vector>

#include "bench.hpp"
#include "oink/game.hpp"
#include "oink/solution.hpp"

namespace bench {

/**
 * The input of the trials of one game.
 */
struct Prepared
{
    const pg::Game* original = nullptr; // the parsed, sorted game
    const pg::Game* residual = nullptr; // the game the solvers solve; null if preprocessing solved all
    std::unique_ptr<pg::Game> owned;    // the residual, if it is not the original
    pg::Solution partial;               // what preprocessing solved, for the original
    std::vector<int> vertices;          // residual vertex -> original vertex
    std::vector<size_t> edges;          // residual edge index -> original edge index
    PrepInfo info;

    bool identity() const { return residual == original; }
};

/**
 * Prepare the trials of game <original>, which must be sorted. With
 * preprocessing off, the residual is <original> itself. The preprocessing
 * log goes to <log>. Throws if the preparation fails.
 */
std::unique_ptr<Prepared> prepare(const pg::Game& original, const Options& options, std::ostream& log);

/**
 * Run solver <run> on a fresh copy of the residual game, with solver seed
 * <seed>, and verify its solution. Never throws: failures are in the status.
 */
Trial execute(const Prepared& prepared, const SolverRun& run, const Options& options, uint32_t seed);

}

#endif
