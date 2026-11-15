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
 * The live table of bench_solvers: a row per game and a column per solver,
 * with the totals and the progress below. The selected game is marked with a
 * cyan >, the fastest complete median of a row is green and failures are red;
 * in heatmap mode all timings are coloured by their rank in the row.
 *
 * After the benchmark, the review can sort the rows by a solver column or by
 * name, with the keys or by clicking the headers.
 *
 * The view state, the rendering and the key handling have no I/O, so the unit
 * tests can drive them; live.cpp connects them to the terminal.
 */

#ifndef OINK_BENCH_TABLE_HPP
#define OINK_BENCH_TABLE_HPP

#include <optional>
#include <string>
#include <vector>

#include "bench.hpp"
#include "common/terminal.hpp"

namespace bench {

/** What the table shows. */
struct TableData
{
    const Options& options;
    const std::vector<SolverRun>& solvers;
    const std::vector<Case>& games;
    size_t planned;                 // games planned
    size_t done;                    // finished slots
    size_t slots;                   // planned slots
    double elapsed;                 // seconds since the start (frozen at the end)
    std::string activity;           // the status line
};

/** The state of the view: selection, scrolling, columns and sorting. */
struct TableView
{
    size_t page = 0;                // the page of solver columns
    size_t selected = 0;            // the selected row, a position in the display order
    size_t first = 0;               // the first visible row
    bool follow = true;             // the selection follows the newest game
    bool heatmap = false;
    bool paused = false;            // the user paused the measurements
    bool review = false;            // the benchmark is over: sorting and mouse
    bool sorted = false;
    bool descending = false;
    std::optional<size_t> sort_column; // the sorted solver; nothing: sorted by name
    std::vector<size_t> order;      // display position -> game

    // the layout of the last frame, for mouse clicks and the 1-9 keys
    int name_width = 0;
    int cell_width = 0;
    size_t column_first = 0;
    size_t column_last = 0;
    size_t columns = 1;             // solver columns per page
    size_t body = 1;                // rows of games
};

/** The rows of the screen before the games, and the row of the headers. */
const int TABLE_HEADER_ROW = 4;

/** Keep the display order in step with the games, and follow the newest game. */
void sync_view(TableView& view, size_t games);

/** The game of the selected row, if there is one. */
std::optional<size_t> selected_game(const TableView& view);

/**
 * Sort the rows by solver <column> (by name if nothing), or restore the
 * original order with <original>. Sorting the same column again reverses the
 * direction. Failed and incomplete cells stay last in either direction, and
 * equal values keep the benchmark order.
 */
void sort_rows(TableView& view, const std::vector<Case>& games, std::optional<size_t> column, bool original, int repeats);

/** Render the table to <height> lines of <width> columns. */
std::vector<std::string> render_table(TableView& view, const TableData& data, int width, int height);

enum class TableAction { None, Changed, Quit, Extract };

/**
 * Apply a key or mouse event. While <measuring>, browsing pauses the
 * measurements and Space resumes them; the mouse and sorting only work in
 * the review.
 */
TableAction table_key(TableView& view, const TableData& data, const ui::Key& key, bool measuring);

}

#endif
