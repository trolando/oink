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
 * The live table on the terminal (see table.hpp), with keyboard input between
 * trials and, after the benchmark, the review with mouse support.
 */

#ifndef OINK_BENCH_LIVE_HPP
#define OINK_BENCH_LIVE_HPP

#include <chrono>
#include <csignal>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "common/terminal.hpp"
#include "display.hpp"
#include "table.hpp"

namespace bench {

class LiveDisplay : public Display
{
public:
    /** Take over the terminal; q and Ctrl-C set <*quit>. */
    LiveDisplay(const Options& options, const std::vector<SolverRun>& solvers, volatile std::sig_atomic_t* quit);

    void update(const Progress& progress, bool force) override;
    void completed(const Case&, size_t) override { }
    bool paused() const override { return view.paused && !stopped; }
    void wait() override;
    std::optional<size_t> take_extract() override;
    void notice(const std::string& text) override;
    void stop_timer() override;
    void close() override;

    /**
     * The review after the benchmark: browse and sort the table until q.
     * Extraction requests go to <extract_game>.
     */
    void review(const Progress& progress, const std::function<void(size_t)>& extract_game);

private:
    const Options& options;
    const std::vector<SolverRun>& solvers;
    volatile std::sig_atomic_t* quit;
    ui::Terminal terminal;
    TableView view;
    std::chrono::steady_clock::time_point started;
    std::optional<std::chrono::steady_clock::time_point> stopped;
    std::chrono::steady_clock::time_point last_draw;
    int last_width = 0;
    int last_height = 0;
    std::optional<size_t> extract;  // a requested extraction
    std::string message;            // a notice instead of the activity

    double elapsed() const;
    bool input(const Progress& progress);
    void draw(const Progress& progress);
};

}

#endif
