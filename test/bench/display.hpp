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
 * What bench_solvers shows while it runs. The runner calls the display only
 * between trials, never during one, so drawing and input handling do not
 * overlap with a measurement.
 */

#ifndef OINK_BENCH_DISPLAY_HPP
#define OINK_BENCH_DISPLAY_HPP

#include <optional>
#include <string>
#include <vector>

#include "bench.hpp"

namespace bench {

/** The state of the benchmark for the display. */
struct Progress
{
    const std::vector<Case>& games; // the games so far; the last one may be running
    size_t planned;                 // games planned
    size_t done;                    // finished slots (cells, or trials without calibration)
    size_t slots;                   // planned slots
    std::string activity;           // what runs next
};

class Display
{
public:
    virtual ~Display() = default;

    /** Show the state between two trials and handle queued input. */
    virtual void update(const Progress& progress, bool force) = 0;

    /** A game has finished. */
    virtual void completed(const Case& game, size_t planned) = 0;

    /** Whether the user paused the measurements (only the live table can). */
    virtual bool paused() const { return false; }

    /** Wait a moment for input while paused. */
    virtual void wait() { }

    /** The game the user asked to extract, if any. */
    virtual std::optional<size_t> take_extract() { return std::nullopt; }

    /** Show a message, such as where an extracted game was saved. */
    virtual void notice(const std::string&) { }

    /** Stop the run timer: the measurements are over. */
    virtual void stop_timer() { }

    /** Leave the screen, before the summary is printed. */
    virtual void close() { }
};

/**
 * Plain output on stderr: a header, and a line per finished game with the
 * median and interquartile range of every solver.
 */
class PlainDisplay : public Display
{
public:
    PlainDisplay(const Options& options, const std::vector<SolverRun>& solvers);
    void update(const Progress&, bool) override { }
    void completed(const Case& game, size_t planned) override;

private:
    const Options& options;
    const std::vector<SolverRun>& solvers;
};

}

#endif
