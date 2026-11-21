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
 * Unit tests for the live table of bench_solvers: its layout and colours, the
 * keys and the mouse, and sorting in the review. The table has no I/O, so the
 * tests render it to lines and press keys.
 */

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "bench/table.hpp"
#include "common/terminal.hpp"

using namespace bench;

static int failures = 0;

static void
check(const char* name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

/** A finished game with one measured value per solver (NaN: a timeout). */
static Case
game(const std::string& name, const std::vector<double>& medians)
{
    Case result;
    result.name = name;
    result.cells.resize(medians.size());
    for (size_t s = 0; s < medians.size(); s++) {
        if (std::isnan(medians[s])) {
            result.cells[s].failure = Status::Timeout;
        } else {
            result.cells[s].add_sample(medians[s]);
        }
    }
    return result;
}

static std::vector<SolverRun>
roster(size_t count)
{
    std::vector<SolverRun> result;
    for (size_t i = 0; i < count; i++) result.emplace_back("s" + std::to_string(i), "");
    return result;
}

static size_t
find_row(const std::vector<std::string>& lines, const std::string& text)
{
    for (size_t i = 0; i < lines.size(); i++) {
        if (ui::strip(lines[i]).find(text) != std::string::npos) return i;
    }
    return lines.size();
}

static void
test_render()
{
    Options o;
    o.repeats = 1;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto solvers = roster(3);
    std::vector<Case> games = {game("alpha", {0.002, 0.001, nan}), game("beta", {0.003, 0.004, 0.005})};
    games.back().worker_peak_rss = 3 * 1024 * 1024;
    TableView view;
    const TableData data{o, solvers, games, 4, 6, 12, 1.5, "beta | s1 | repeat 1/1"};
    const auto lines = render_table(view, data, 100, 30);
    check("one line per row but the last", lines.size() == 29);
    check("title", ui::strip(lines[0]).find("OINK  /  SOLVER BENCHMARK") == 0);
    check("headers", ui::strip(lines[TABLE_HEADER_ROW - 1]).find("GAME") == 0 && ui::strip(lines[TABLE_HEADER_ROW - 1]).find("s2") != std::string::npos);
    const size_t alpha = find_row(lines, "alpha");
    const size_t beta = find_row(lines, "beta");
    check("rows in order", alpha == (size_t)TABLE_HEADER_ROW && beta == alpha + 1);
    check("the newest game is selected", ui::strip(lines[beta]).find("> beta") == 0);
    check("best median green", lines[alpha].find(std::string(ui::style::good) + "1.00 ms") != std::string::npos);
    check("failure red", lines[alpha].find(std::string(ui::style::bad) + "timeout") != std::string::npos);
    const size_t sums = find_row(lines, "SUM MEDIANS");
    check("sums of the complete cells", sums < lines.size() && ui::strip(lines[sums]).find("5.00 ms") != std::string::npos);
    check("coverage", find_row(lines, "verified games") == sums + 1 && ui::strip(lines[sums + 1]).find("1/4") != std::string::npos);
    check("selected row mean", ui::strip(lines[sums + 2]).find("selected mean") == 0 &&
          ui::strip(lines[sums + 2]).find("|3.00 ms") != std::string::npos);
    check("selected row sample count", ui::strip(lines[sums + 3]).find("selected solves") == 0 &&
          ui::strip(lines[sums + 3]).find("|1") != std::string::npos);
    check("reuse RSS", ui::strip(lines[sums + 4]).find("Worker process peak RSS  3.0 MiB") == 0);
    check("progress", ui::strip(lines[sums + 5]).find("50%  6/12") != std::string::npos);
    check("activity", ui::strip(lines[sums + 6]).find("beta | s1 | repeat 1/1") == 0);
    check("keys", ui::strip(lines[sums + 7]).find("Space pause/resume") == 0);

    o.fresh_process = true;
    games.back().cells[0].max_trial_rss = 5 * 1024 * 1024;
    games.back().cells[0].max_solve_rss = 4 * 1024 * 1024;
    view.review = true;
    const auto review = render_table(view, data, 100, 30);
    check("fresh RSS in review", find_row(review, "|5.0 MiB") == find_row(review, "trial peak RSS") &&
          find_row(review, "|4.0 MiB") == find_row(review, "solve-end RSS") &&
          find_row(review, "solve-end RSS") < review.size());
    check("fresh has no reuse RSS", find_row(review, "Worker process peak RSS") == review.size());
    view.review = false;
    o.fresh_process = false;

    view.heatmap = true;
    const auto heat = render_table(view, data, 100, 30);
    check("heatmap: fastest of beta green", heat[beta].find(ui::style::color256(82) + "3.00 ms") != std::string::npos);
    check("heatmap: slowest of beta red", heat[beta].find(ui::style::color256(196) + "5.00 ms") != std::string::npos);
    check("heatmap mode shown", ui::strip(heat[1]).find("heatmap") != std::string::npos);

    // a running game shows its phase
    Case running = game("gamma", {0.001, 0.001, 0.001});
    running.cells[0] = Cell();
    running.cells[0].warmups = 1;
    running.cells[1] = Cell();
    running.cells[2].required = 5;
    games.push_back(running);
    TableView live;
    const TableData more{o, solvers, games, 4, 7, 12, 1.5, ""};
    const auto partial = render_table(live, more, 100, 30);
    const size_t gamma = find_row(partial, "gamma");
    check("warm-up cell", ui::strip(partial[gamma]).find("warm 1") != std::string::npos);
    check("incomplete cell count", ui::strip(partial[gamma]).find("1.00 ms (1)") != std::string::npos);
}

static void
test_keys()
{
    Options o;
    o.repeats = 1;
    const auto solvers = roster(2);
    std::vector<Case> games;
    for (int i = 0; i < 40; i++) games.push_back(game("g" + std::to_string(i), {0.001 * (40 - i), 0.001}));
    TableView view;
    const TableData data{o, solvers, games, 40, 0, 80, 0, ""};
    render_table(view, data, 100, 30);
    check("follows the newest game", view.selected == 39 && view.follow);
    check("up selects and pauses", table_key(view, data, ui::Key(ui::Key::Up), true) == TableAction::Changed && view.selected == 38 && view.paused && !view.follow);
    check("space resumes", table_key(view, data, ui::Key::character(' '), true) == TableAction::Changed && !view.paused);
    table_key(view, data, ui::Key::character('j'), true);
    check("down to the last game follows again", view.selected == 39 && view.follow);
    check("h toggles the heatmap and pauses", table_key(view, data, ui::Key::character('h'), true) == TableAction::Changed && view.heatmap && view.paused);
    check("e asks for an extraction", table_key(view, data, ui::Key::character('e'), true) == TableAction::Extract && selected_game(view) == size_t(39));
    check("q quits", table_key(view, data, ui::Key::character('q'), true) == TableAction::Quit);
    check("no sorting while measuring", table_key(view, data, ui::Key::character('1'), true) == TableAction::None && !view.sorted);
    ui::Key click(ui::Key::Mouse);
    click.x = 40;
    click.y = TABLE_HEADER_ROW;
    check("no mouse while measuring", table_key(view, data, click, true) == TableAction::None);

    // the viewport moves only at its edges
    for (int i = 0; i < 5; i++) table_key(view, data, ui::Key(ui::Key::Up), true);
    render_table(view, data, 100, 30);
    const size_t first = view.first;
    table_key(view, data, ui::Key(ui::Key::Down), true);
    render_table(view, data, 100, 30);
    check("viewport stays while the selection moves inside", view.first == first);
}

static void
test_review()
{
    Options o;
    o.repeats = 1;
    const auto solvers = roster(12);
    std::vector<Case> games;
    // solver 5 is fastest on c, then a, b; it failed on d
    games.push_back(game("b", std::vector<double>(12, 0.003)));
    games.push_back(game("d", std::vector<double>(12, 0.001)));
    games.push_back(game("a", std::vector<double>(12, 0.002)));
    games.push_back(game("c", std::vector<double>(12, 0.004)));
    games[1].cells[5] = Cell();
    games[1].cells[5].failure = Status::Crash;
    games[3].cells[5] = Cell();
    games[3].cells[5].add_sample(0.0005);
    TableView view;
    view.review = true;
    const TableData data{o, solvers, games, 4, 12, 12, 1, ""};
    render_table(view, data, 80, 24);
    check("columns per page", view.columns == 4 && view.column_first == 0 && view.column_last == 4);

    // page two and a click on its second header sorts solver 5
    table_key(view, data, ui::Key(ui::Key::Right), false);
    render_table(view, data, 80, 24);
    check("second page", view.column_first == 4 && view.column_last == 8);
    ui::Key click(ui::Key::Mouse);
    click.button = 0;
    click.y = TABLE_HEADER_ROW;
    click.x = view.name_width + 3 + view.cell_width + 1;
    check("click on a header sorts", table_key(view, data, click, false) == TableAction::Changed && view.sorted && view.sort_column == size_t(5));
    check("ascending, failures last", view.order == std::vector<size_t>({3, 2, 0, 1}));
    auto lines = render_table(view, data, 80, 24);
    check("sort marker", ui::strip(lines[TABLE_HEADER_ROW - 1]).find("s5") != std::string::npos && ui::strip(lines[TABLE_HEADER_ROW - 1]).find(" ^") != std::string::npos);
    table_key(view, data, click, false);
    check("again: descending, failures still last", view.descending && view.order == std::vector<size_t>({0, 2, 3, 1}));
    click.release = true;
    check("a release does not sort again", table_key(view, data, click, false) == TableAction::None && view.descending);

    // keys: 1-9 for the visible columns, g by name, o the original order
    table_key(view, data, ui::Key::character('2'), false);
    check("2 sorts the second visible column", view.sort_column == size_t(5) && !view.descending);
    table_key(view, data, ui::Key::character('g'), false);
    check("g sorts by name", !view.sort_column && view.order == std::vector<size_t>({2, 0, 3, 1}));
    lines = render_table(view, data, 80, 24);
    check("name sort marker", ui::strip(lines[TABLE_HEADER_ROW - 1]).find("GAME ^") == 0);
    table_key(view, data, ui::Key::character('o'), false);
    check("o restores the order", !view.sorted && view.order == std::vector<size_t>({0, 1, 2, 3}));

    // a click on a row selects it; the wheel moves the selection
    ui::Key row(ui::Key::Mouse);
    row.button = 0;
    row.y = TABLE_HEADER_ROW + 3;
    row.x = 5;
    table_key(view, data, row, false);
    check("click selects a row", view.selected == 2 && selected_game(view) == size_t(2));
    ui::Key wheel(ui::Key::Mouse);
    wheel.button = 64;
    table_key(view, data, wheel, false);
    check("wheel up", view.selected == 0);
    wheel.button = 65;
    table_key(view, data, wheel, false);
    check("wheel down", view.selected == 3);
    check("review keys shown", ui::strip(render_table(view, data, 100, 24).back()).find("Click headers/1-9 sort") != std::string::npos);
}

int
test_table()
{
    test_render();
    test_keys();
    test_review();
    if (failures) std::cerr << failures << " table test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
