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
 * The final summary of a benchmark: per solver the coverage, the sum of its
 * per-game medians, the geometric mean of its medians over the fastest median
 * of each game (the virtual best solver), the virtual best solver counts and,
 * with a reference, the geometric mean speedup and wins/ties/losses against
 * the reference.
 */

#ifndef OINK_BENCH_SUMMARY_HPP
#define OINK_BENCH_SUMMARY_HPP

#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "bench.hpp"

namespace bench {

// Point estimates within 1% count as ties (W/T/L and VBS).
const double TIE_BAND = 0.01;

struct SolverTotals
{
    size_t complete = 0;            // verified games
    size_t failed = 0;              // failures, not counting interruptions
    double sum = 0;                 // the sum of the medians of the verified games
    size_t paired = 0;              // games verified for both this solver and the reference
    size_t wins = 0;                // faster than the reference by more than 1%
    size_t ties = 0;
    size_t losses = 0;
    size_t reference_only = 0;      // verified for the reference only
    size_t solver_only = 0;         // verified for this solver only
    double log_ratios = 0;          // the sum of log(reference/solver) over the pairs
    size_t vbs_wins = 0;            // games where this solver is uniquely fastest
    size_t vbs_ties = 0;            // games where it shares the fastest median (1% band)
    size_t vbs_games = 0;           // games verified or timed out, with a fastest median
    double log_vbs = 0;             // the sum of log(solver/fastest) over these games

    std::optional<double> geomean() const;
    /** The geometric mean of the median over the fastest median (a timeout counts as the limit). */
    std::optional<double> vbs_ratio() const;
};

struct Totals
{
    size_t planned = 0;             // games planned
    size_t started = 0;             // games with results
    size_t preprocessed = 0;        // games solved by preprocessing
    size_t timed_planned = 0;       // games that need a solver
    size_t prep_failures = 0;       // games whose preparation failed
    bool complete_coverage = true;  // every solver verified every game
    std::optional<size_t> reference;
    double timeout = 0;             // the time limit of a trial; 0: none
    size_t vbs_games = 0;           // games with a fastest median (of the VBS ratios)
    std::vector<SolverTotals> solvers;
    std::vector<size_t> order;      // full coverage first, then by sum
};

Totals compute_totals(const std::vector<Case>& games, size_t planned, size_t solver_count, int repeats, std::optional<size_t> reference, double timeout);

/** What the "Run details" part of the summary reports. */
struct RunDetails
{
    uint64_t schedule_seed = 0;
    std::optional<uint32_t> game_seed; // for random games
    bool preprocess = false;
    std::string jsonl;
    std::string csv;
    bool interrupted = false;
};

/** Print the summary table, the column guide and the run details. */
void print_summary(std::ostream& out, const Totals& totals, const std::vector<SolverRun>& solvers, const RunDetails& details, bool styled);

}

#endif
