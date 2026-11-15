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

#include "live.hpp"

#include <chrono>

namespace bench {


using Clock = std::chrono::steady_clock;

LiveDisplay::LiveDisplay(const Options& options, const std::vector<SolverRun>& solvers, volatile std::sig_atomic_t* quit)
    : options(options), solvers(solvers), quit(quit), started(Clock::now()), last_draw(Clock::now())
{
    view.heatmap = options.heatmap;
}

double
LiveDisplay::elapsed() const
{
    const auto end = stopped ? *stopped : Clock::now();
    return std::chrono::duration<double>(end - started).count();
}

/** Read the queued keys; true if any changed the view. */
bool
LiveDisplay::input(const Progress& progress)
{
    const TableData data{options, solvers, progress.games, progress.planned, progress.done, progress.slots, elapsed(), progress.activity};
    bool changed = false;
    while (terminal.pending()) {
        const ui::Key key = terminal.read(0);
        if (key.type == ui::Key::None) break;
        const TableAction action = table_key(view, data, key, !stopped);
        message.clear(); // a notice lasts until the next key
        if (action == TableAction::Quit && quit != nullptr) *quit = 1;
        if (action == TableAction::Extract) {
            const auto game = selected_game(view);
            if (game) extract = game;
        }
        if (action != TableAction::None) changed = true;
    }
    return changed;
}

void
LiveDisplay::draw(const Progress& progress)
{
    int width = 0;
    int height = 0;
    terminal.size(width, height);
    const std::string activity = message.empty() ? progress.activity : message;
    const TableData data{options, solvers, progress.games, progress.planned, progress.done, progress.slots, elapsed(), activity};
    terminal.draw(render_table(view, data, width, height));
    last_width = width;
    last_height = height;
    last_draw = Clock::now();
}

void
LiveDisplay::update(const Progress& progress, bool force)
{
    const bool changed = input(progress);
    int width = 0;
    int height = 0;
    terminal.size(width, height);
    const bool resized = width != last_width || height != last_height;
    // without news, redraw at most ten times a second
    const double since = std::chrono::duration<double>(Clock::now() - last_draw).count();
    if (!force && !changed && !resized && since < 0.1) return;
    draw(progress);
}

void
LiveDisplay::wait()
{
    terminal.pending(100);
}

std::optional<size_t>
LiveDisplay::take_extract()
{
    auto result = extract;
    extract.reset();
    return result;
}

void
LiveDisplay::notice(const std::string& text)
{
    message = text;
}

void
LiveDisplay::stop_timer()
{
    if (!stopped) stopped = Clock::now();
    view.paused = false;
}

void
LiveDisplay::review(const Progress& progress, const std::function<void(size_t)>& extract_game)
{
    view.review = true;
    view.follow = false;
    terminal.set_mouse(true);
    draw(progress);
    while (quit == nullptr || !*quit) {
        const ui::Key key = terminal.read(100);
        bool changed = false;
        if (key.type != ui::Key::None) {
            const TableData data{options, solvers, progress.games, progress.planned, progress.done, progress.slots, elapsed(), progress.activity};
            message.clear();
            const TableAction action = table_key(view, data, key, false);
            if (action == TableAction::Quit) break;
            if (action == TableAction::Extract) {
                const auto game = selected_game(view);
                if (game) extract_game(*game);
            }
            changed = action != TableAction::None;
        }
        int width = 0;
        int height = 0;
        terminal.size(width, height);
        if (changed || width != last_width || height != last_height) draw(progress);
    }
    terminal.set_mouse(false);
}

void
LiveDisplay::close()
{
    terminal.close();
}

}
