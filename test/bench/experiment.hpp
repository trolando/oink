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
 * Named experiments: a benchmark in a directory that can be resumed.
 *
 * The directory holds
 *
 *   plan.txt          the frozen plan: settings, solvers, games in their order,
 *                     the program and the host;
 *   done.txt          one line per finished game with all its observations;
 *   games/N.jsonl     the events of finished game N;
 *   sessions/N.jsonl  everything that happened in session N, also the
 *                     trials of a game that was interrupted;
 *   failures/         the games and logs of failed trials;
 *   results.jsonl     metadata, the events of all finished games and the
 *   results.csv       summary, and the statistics, rebuilt after every session;
 *   status.json       how far the experiment is.
 *
 * The unit of progress is a whole game: an interrupted game is measured again
 * from its warm-ups when the experiment resumes; its partial trials stay in
 * the session log only. Failures are final and are not retried. A resume
 * checks that the program, the host and the input files are unchanged.
 */

#ifndef OINK_BENCH_EXPERIMENT_HPP
#define OINK_BENCH_EXPERIMENT_HPP

#include <csignal>
#include <istream>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "bench.hpp"
#include "runner.hpp"

namespace bench {

struct Totals;

/** The settings of <options> as command-line arguments, for the plan. */
std::vector<std::string> settings_arguments(const Options& options);

/** One finished game as a line of done.txt, and back (false if the line is incomplete). */
std::string encode_case(const Case& game);
bool decode_case(const std::string& line, Case& game);

class Experiment
{
public:
    /**
     * Create the experiment directory <directory> with the plan of
     * <options>, <solvers> and <games>. Throws if it exists.
     */
    static std::unique_ptr<Experiment> create(const std::string& directory, const Options& options, const std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& games, const std::string& executable);

    /**
     * Open the experiment in <directory> to resume it: its settings, solvers
     * and games replace <options>, <solvers> and <games>. Throws if the
     * program, the host or an input changed.
     */
    static std::unique_ptr<Experiment> resume(const std::string& directory, Options& options, std::vector<SolverRun>& solvers, std::vector<PlannedGame>& games, const std::string& executable);

    ~Experiment();

    /** Point the reports of <options> to this session; no live table. */
    void configure(Options& options) const;

    /** The finished game <index>, if any. */
    const Case* finished(size_t index) const;

    /** Commit finished game <game> with its events <events>. */
    void commit(const Case& game, const std::string& events);

    /** Rebuild results.jsonl, results.csv and status.json from the finished games. */
    void publish(const Options& options, const std::vector<SolverRun>& solvers, const std::string& executable, const std::vector<PlannedGame>& games);

    const std::string& directory() const { return root; }
    int session() const { return session_; }
    size_t finished_count() const { return done.size(); }
    int lock_fd() const { return lock; }

private:
    std::string root;
    int session_ = 0;
    int lock = -1;
    std::map<size_t, Case> done;

    Experiment(const std::string& directory);
    void acquire_lock();
    void start_session();
    void load_done(size_t solvers);
};

/**
 * Whether to resume the existing experiment <directory>: asks on a terminal
 * (Enter resumes); without one, an existing experiment needs --resume.
 * Returns false if the user declines.
 */
bool confirm_resume(const std::string& directory, volatile std::sig_atomic_t* quit);

}

#endif
