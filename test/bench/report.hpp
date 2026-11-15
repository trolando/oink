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
 * The reports of a benchmark.
 *
 * JSON Lines: one event per line, written and flushed as it happens: the
 * metadata of the run, every game, every warm-up, calibration and measured
 * trial, the statistics per game and solver, the summary and an end marker.
 * A run that is killed still leaves complete lines.
 *
 * CSV: one row per game and solver with the median, quartiles and median
 * absolute deviation of the measured solve times.
 */

#ifndef OINK_BENCH_REPORT_HPP
#define OINK_BENCH_REPORT_HPP

#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "bench.hpp"

namespace bench {

/** The program, build and machine, for the metadata. */
std::string build_json();
std::string machine_json();
std::string executable_path(const char* argv0);
std::string host_name();
std::string utc_now();

/** Totals of the benchmark per solver, for the summary. */
struct Totals;

class Report
{
public:
    /**
     * Open the reports of <options> (none if their paths are empty). Existing
     * files are only replaced with --overwrite.
     */
    Report(const Options& options, const std::vector<SolverRun>& solvers);

    /**
     * Write the metadata event of a run of <planned> games; <session> is the
     * session of a named experiment, if any.
     */
    void metadata(size_t planned, const std::string& executable, const std::string& session = "");

    /** Write one event, an object built with Json. */
    void event(const std::string& json);

    void game(const Case& game);
    void sample(const Case& game, size_t solver, size_t position, int round, const char* phase, uint32_t seed, int generation, const Trial& trial, bool executed, const std::string& artifact);
    void calibration(const Case& game, size_t solver, double pilot, int target, bool capped);

    /** The result per solver of a finished game (JSON and CSV). */
    void results(const Case& game);

    /** The summary events and the end marker. */
    void summary(const Totals& totals, bool interrupted, size_t started);

    /**
     * Also copy every event to <out> (null to stop). Named experiments keep the
     * events of a game apart this way.
     */
    void copy_to(std::ostream* out) { copy = out; }

    const std::string& jsonl_path() const { return jsonl_path_; }
    const std::string& csv_path() const { return csv_path_; }

private:
    const Options& options;
    const std::vector<SolverRun>& solvers;
    std::ofstream jsonl;
    std::ofstream csv;
    std::ostream* copy = nullptr;
    std::string jsonl_path_;
    std::string csv_path_;
};

/** The CSV header and the rows of a finished game. */
std::string csv_header();
std::string csv_rows(const Case& game, const std::vector<SolverRun>& solvers, int repeats);

/** The result events of a finished game. */
std::vector<std::string> result_events(const Case& game, const std::vector<SolverRun>& solvers, int repeats);

}

#endif
