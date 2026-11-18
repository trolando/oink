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
 * Core types of bench_solvers: the settings of a benchmark, the outcome of a
 * single solver call (a trial), the observations of one solver on one game (a
 * cell) and of all solvers on one game (a case), and small helpers for
 * statistics, formatting and seeds that have no I/O.
 */

#ifndef OINK_BENCH_HPP
#define OINK_BENCH_HPP

#include <cstdint>
#include <exception>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "common/solver_run.hpp"

namespace pg {
class Game;
}

namespace bench {

/**
 * The settings of a benchmark, from the command line or the setup screen.
 */
struct Options
{
    // Solving
    int workers = -1;               // Lace workers, -1 for sequential, 0 for autodetect
    bool force_parallel = false;    // parallel solvers also on small games
    int trace = -1;                 // trace level of the captured solver log; -1: no log
    bool scc = false;               // solve one bottom SCC at a time (inside the timing)

    // Once-per-game preprocessing (outside the timing)
    bool preprocess = false;
    bool inflate = false;
    bool compress = false;
    bool single = true;
    bool loops = true;
    bool wcwc = true;

    // Sampling
    int repeats = 7;                // minimum measured repetitions per cell
    int warmups = 1;                // warm-up rounds per game
    int max_repeats = 10000;        // cap when sampling has a solve-time target
    double min_solve_time = 0.001;  // accumulated solve time per cell; 0: fixed count
    double timeout = 30;            // wall seconds per solve; 0: no solve limit
    double safety_timeout = 300;    // wall seconds each for worker setup and verification
    bool isolate = true;            // run the trials in a worker process
    uint64_t seed = 0;              // the schedule seed (--bench-seed)
    bool vary_solver_seed = false;  // vary the solver RNG seed between rounds
    bool sort = false;              // keep the sorted file order
    bool recursive = false;         // search input directories recursively

    // Random games, when no inputs are given
    int count = 100;
    int size = 100;
    int maxp = 100;
    int maxe = 400;
    uint32_t game_seed = 0;         // the seed of the series of random games
    bool game_seed_given = false;   // --seed was given (else derived from the schedule seed)

    // Inputs and external solvers ("name:command")
    std::vector<std::string> inputs;
    std::vector<std::string> externals;

    // Display
    bool live = false;              // the live table
    bool review = false;            // keep the table open after the benchmark
    bool interactive = false;       // setup screen, live table and review
    bool heatmap = false;           // colour all timings by their rank in the game
    bool plain = false;             // no escape sequences

    // Reports
    std::string jsonl;
    std::string csv;
    std::string label;
    std::string artifacts = "benchmark-artifacts";
    std::string exports = "benchmark-extracts";
    bool overwrite = false;
    std::optional<std::string> reference; // a solver name, "first", or "none"

    // Named experiments
    std::string experiment;
    bool resume = false;

    std::vector<std::string> argv;  // the command line, for the metadata

    bool adaptive_sampling() const { return min_solve_time > 0; }
};

/** The solvers benchmarked when none are selected. */
std::vector<SolverRun> default_solvers();

/**
 * The index of the reference solver: a selected name, a selected id with one
 * variant, or "first"; nothing for "none" or no reference. Throws for an
 * unselected or ambiguous reference.
 */
std::optional<size_t> reference_index(const std::optional<std::string>& reference, const std::vector<SolverRun>& solvers);

/** Check the settings; throws std::runtime_error with an explanation. */
void validate(const Options& options, const std::vector<SolverRun>& solvers);

/**
 * The outcome of one trial.
 */
enum class Status { Ok, SolverError, InvalidSolution, SetupError, Timeout, SafetyTimeout, MemoryFailure, Crash, Interrupted };
enum class Phase { None, Setup, Solve, Verify };

const char* status_name(Status status);
const char* phase_name(Phase phase);
std::string failure_label(Status status, Phase phase);

/** One solver call: the measured times and what actually ran. */
struct Trial
{
    Status status = Status::SetupError;
    Phase failure_phase = Phase::None;
    // A negative time means the measurement was unavailable.
    double setup = -1;              // copy the game and construct Oink, with the in-edges
    double solve = -1;              // Oink::run, the measured time
    double verify = -1;             // map the solution back and verify it
    double cpu = -1;                // CPU seconds of the process during the solve
    bool parallel = false;          // whether the solver ran in parallel
    int workers = -1;               // Lace workers, -1 when sequential
    std::string effective_solver;   // the solver that ran, after a sequential fallback
    std::string effective_config;
    std::string message;            // the error, and the captured log of a failure
};

using PhaseUpdate = std::function<void(Phase, const Trial&)>;
void record_failure(Trial& trial, Phase phase, const std::exception& error);

/** Mean, median, quartiles and median absolute deviation of a sample. */
struct Stats
{
    size_t samples = 0;
    double mean = 0;
    double median = 0;
    double q1 = 0;
    double q3 = 0;
    double mad = 0;
};

/** Quartiles by linear interpolation at (n-1)p. */
Stats summarize(std::vector<double> values);

/**
 * The observations of one solver on one game.
 */
struct Cell
{
    int warmups = 0;                // successful warm-up trials
    int required = 0;               // the current target if extended beyond --repeat
    std::optional<Status> failure;  // the first failure; later trials are skipped
    Phase failure_phase = Phase::None;

    int target(int repeats) const { return required ? required : repeats; }
    bool complete(int repeats) const { return !failure && (int)values.size() == target(repeats); }
    bool capped(const Options& o) const { return complete(o.repeats) && solve_total < o.min_solve_time; }

    void add_sample(double seconds);
    void extend_sampling(const Options& o);
    const std::vector<double>& samples() const { return values; }
    double accumulated_solve() const { return solve_total; }

    /** The statistics of the measured values (cached; values are only appended). */
    const Stats& stats() const;

private:
    std::vector<double> values;
    double solve_total = 0;
    mutable Stats cache;
};

/** What the once-per-game preparation did to a game. */
struct PrepInfo
{
    long residual_vertices = 0;
    long residual_edges = 0;
    long solved_vertices = 0;       // solved by preprocessing
    double seconds = 0;             // preprocessing time (not part of any trial)
    std::string residual_hash;
};

/**
 * One game of the benchmark, with a cell per solver.
 */
struct Case
{
    size_t index = 0;               // position in the schedule
    std::string name;               // the file name, or "game N seed=S"
    std::string path;               // the canonical path; empty for a random game
    std::optional<uint32_t> seed;   // the seed of a random game
    std::string hash;               // fingerprint of the parsed game
    long vertices = 0;
    long edges = 0;
    double load = 0;                // reading and decompressing
    double parse = 0;               // parsing and sorting
    std::string error;              // why the game could not be read or prepared
    bool input_failed = false;      // the game could not be read
    bool prep_failed = false;       // the preparation failed
    bool fully_preprocessed = false; // preprocessing solved the whole game
    PrepInfo prep;
    std::vector<Cell> cells;

    /** Whether the solvers ran (or should have run) on this game. */
    bool timed() const { return !fully_preprocessed; }

    /** The status of the cell of solver <s> for the reports. */
    std::string status(size_t s, int repeats) const;
};

/** Format seconds as "12.3 us", "4.56 ms" or "1.23 s". */
std::string duration(double seconds);

/** Format seconds briefly for a setting, such as "30 s" or "1 ms". */
std::string seconds_setting(double seconds);

/**
 * Colours of the 256-colour palette for <values> by their rank: green for the
 * fastest 5%, red for the slowest 5%, yellow to orange between. Equal values
 * share a colour, all-equal values are neutral (yellow), and NaN or negative
 * values get -1 (no colour).
 */
std::vector<int> heatmap_colors(const std::vector<double>& values);

/** FNV-1a fingerprint of a game: its size, priorities, owners and edges. */
uint64_t fingerprint(const pg::Game& game);

/** FNV-1a fingerprint of a file, or nothing if it cannot be read. */
std::optional<uint64_t> file_fingerprint(const std::string& path);

/** A 64-bit number as 16 hexadecimal digits. */
std::string hex64(uint64_t value);

/** The splitmix64 finalizer, to derive independent seeds. */
uint64_t mix(uint64_t x);

/** The seed of the random game series, derived from the schedule seed. */
uint32_t series_seed(uint64_t schedule_seed);

/** The RNG seed of the solvers on a game in a round (round 0 unless varied). */
uint32_t trial_seed(uint64_t schedule_seed, uint64_t game, int round);

/** The seed that shuffles the solvers of a round of a game. */
uint64_t round_seed(uint64_t schedule_seed, uint64_t game, int round);

}

#endif
