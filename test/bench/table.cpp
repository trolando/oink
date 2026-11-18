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

#include "table.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace bench {

using ui::fit;
namespace style = ui::style;

void
sync_view(TableView& view, size_t games)
{
    while (view.order.size() < games) view.order.push_back(view.order.size());
    if (view.follow && games > 0) view.selected = games - 1;
    if (games > 0 && view.selected >= games) view.selected = games - 1;
}

std::optional<size_t>
selected_game(const TableView& view)
{
    if (view.selected >= view.order.size()) return std::nullopt;
    return view.order[view.selected];
}

void
sort_rows(TableView& view, const std::vector<Case>& games, std::optional<size_t> column, bool original, int repeats)
{
    view.descending = !original && view.sorted && view.sort_column == column ? !view.descending : false;
    view.sorted = !original;
    view.sort_column = column;
    view.order.resize(games.size());
    std::iota(view.order.begin(), view.order.end(), 0);
    if (view.sorted) {
        const bool descending = view.descending;
        std::stable_sort(view.order.begin(), view.order.end(), [&](size_t a, size_t b) {
            if (!column) return descending ? games[a].name > games[b].name : games[a].name < games[b].name;
            const Cell& x = games[a].cells[*column];
            const Cell& y = games[b].cells[*column];
            const bool left = games[a].timed() && x.complete(repeats);
            const bool right = games[b].timed() && y.complete(repeats);
            // failed and incomplete cells stay last
            if (left != right) return left;
            if (!left) return false;
            return descending ? x.stats().median > y.stats().median : x.stats().median < y.stats().median;
        });
    }
    view.selected = 0;
    view.first = 0;
    view.follow = false;
}

/** Move the selection by <delta> rows; moving up stops following the newest game. */
static void
move(TableView& view, size_t count, long delta)
{
    if (count == 0) return;
    long target = (long)view.selected + delta;
    if (target < 0) target = 0;
    if (target >= (long)count) target = (long)count - 1;
    view.selected = (size_t)target;
    view.follow = !view.review && view.selected + 1 == count;
}

TableAction
table_key(TableView& view, const TableData& data, const ui::Key& key, bool measuring)
{
    const size_t count = data.games.size();
    // browsing during the measurements pauses them
    auto browse = [&] {
        if (measuring) view.paused = true;
        return TableAction::Changed;
    };
    switch (key.type) {
        case ui::Key::Up:
            move(view, count, -1);
            return browse();
        case ui::Key::Down:
            move(view, count, 1);
            return browse();
        case ui::Key::PageUp:
            move(view, count, -(long)view.body);
            return browse();
        case ui::Key::PageDown:
            move(view, count, (long)view.body);
            return browse();
        case ui::Key::Home:
            move(view, count, -(long)count);
            return browse();
        case ui::Key::End:
            move(view, count, (long)count);
            return browse();
        case ui::Key::Left:
            if (view.page > 0) view.page--;
            return browse();
        case ui::Key::Right:
            view.page++;
            return browse();
        case ui::Key::Mouse:
            break;
        case ui::Key::Char:
            break;
        default:
            return TableAction::None;
    }

    if (key.type == ui::Key::Mouse) {
        if (!view.review || key.release) return TableAction::None;
        if (key.button == 64 || key.button == 65) {
            // the wheel moves the selection three rows
            move(view, count, key.button == 65 ? 3 : -3);
            return TableAction::Changed;
        }
        if (key.button != 0) return TableAction::None;
        if (key.y == TABLE_HEADER_ROW) {
            if (key.x <= view.name_width) {
                sort_rows(view, data.games, std::nullopt, false, data.options.repeats);
                return TableAction::Changed;
            }
            if (view.cell_width > 0 && key.x >= view.name_width + 3) {
                const size_t column = view.column_first + (size_t)((key.x - view.name_width - 3) / view.cell_width);
                if (column < view.column_last) {
                    sort_rows(view, data.games, column, false, data.options.repeats);
                    return TableAction::Changed;
                }
            }
            return TableAction::None;
        }
        const int top = TABLE_HEADER_ROW + 1;
        if (key.y >= top && key.y < top + (int)view.body) {
            const size_t row = view.first + (size_t)(key.y - top);
            if (row < count) {
                view.selected = row;
                view.follow = false;
                return TableAction::Changed;
            }
        }
        return TableAction::None;
    }

    const int c = key.ch;
    if (c == 'q' || c == 'Q') return TableAction::Quit;
    if (c == 'k') {
        move(view, count, -1);
        return browse();
    }
    if (c == 'j') {
        move(view, count, 1);
        return browse();
    }
    if (c == '[') {
        if (view.page > 0) view.page--;
        return browse();
    }
    if (c == ']') {
        view.page++;
        return browse();
    }
    if (c == 'h') {
        view.heatmap = !view.heatmap;
        return browse();
    }
    if (c == 'e') {
        if (count == 0) return TableAction::None;
        if (measuring) view.paused = true;
        return TableAction::Extract;
    }
    if (c == ' ') {
        if (measuring) view.paused = !view.paused;
        return TableAction::Changed;
    }
    if (!view.review) return TableAction::None;
    if (c == 'g') {
        sort_rows(view, data.games, std::nullopt, false, data.options.repeats);
        return TableAction::Changed;
    }
    if (c == 'o') {
        sort_rows(view, data.games, std::nullopt, true, data.options.repeats);
        return TableAction::Changed;
    }
    if (c >= '1' && c <= '9') {
        // the visible columns of the current page
        const size_t column = view.column_first + (size_t)(c - '1');
        if (column < view.column_last) {
            sort_rows(view, data.games, column, false, data.options.repeats);
            return TableAction::Changed;
        }
    }
    return TableAction::None;
}

/** The text and style of the cell of solver <s> in game <game>. */
static void
cell(const Case& game, size_t s, int repeats, bool heatmap, int color, double best, std::string& text, std::string& code)
{
    const Cell& c = game.cells[s];
    if (game.fully_preprocessed) {
        text = "preprocessed";
        code = style::dim;
    } else if (game.input_failed) {
        text = "input error";
        code = style::bad;
    } else if (game.prep_failed) {
        text = "prep error";
        code = style::bad;
    } else if (c.failure) {
        if (c.failure_phase == Phase::Setup) text = "setup ";
        if (c.failure_phase == Phase::Solve) text = "solve ";
        if (c.failure_phase == Phase::Verify) text = "verify ";
        if (*c.failure == Status::SafetyTimeout) text += "safety";
        else if (*c.failure == Status::MemoryFailure) text += "memory";
        else text += status_name(*c.failure);
        code = style::bad;
    } else if (c.values.empty()) {
        if (!c.pilots.empty()) text = "pilot " + std::to_string(c.pilots.size());
        else if (c.warmups) text = "warm " + std::to_string(c.warmups);
        else text = "-";
        code = style::dim;
    } else {
        const Stats& st = c.stats();
        text = duration(st.median);
        if (!c.complete(repeats)) {
            text += " " + std::to_string(c.values.size()) + "/" + std::to_string(c.target(repeats));
        } else if (heatmap) {
            code = style::color256(color);
        } else if (st.median == best) {
            code = style::good;
        }
    }
}

std::vector<std::string>
render_table(TableView& view, const TableData& data, int width, int height)
{
    const Options& o = data.options;
    const auto& solvers = data.solvers;
    const auto& games = data.games;
    sync_view(view, games.size());

    // the layout: the game names, then pages of solver columns
    const int name_width = std::min(30, std::max(16, width / 4));
    const size_t columns = std::max<size_t>(1, (size_t)std::max(0, width - name_width - 2) / 12);
    const size_t pages = (solvers.size() + columns - 1) / columns;
    if (view.page >= pages) view.page = pages ? pages - 1 : 0;
    const size_t first_column = view.page * columns;
    const size_t last_column = std::min(first_column + columns, solvers.size());
    const int cell_width = std::min(16, std::max(1, (width - name_width - 2) / std::max(1, (int)(last_column - first_column))));
    view.name_width = name_width;
    view.cell_width = cell_width;
    view.column_first = first_column;
    view.column_last = last_column;
    view.columns = columns;

    std::vector<std::string> lines;
    auto line = [&](const std::string& text) { lines.push_back(text); };

    line(std::string(style::title) + "OINK  /  SOLVER BENCHMARK" + style::reset + "  " + fit(o.label, std::max(0, width - 27)));
    std::string mode = "median Oink::run | " + std::to_string(o.repeats) + (o.calibrated() ? "+ calibrated" : "") + " repeats | " + std::to_string(o.warmups) + " warm-ups | ";
    mode += view.heatmap ? "heatmap: per-game ranks (green/red 5% tails)" : "verified";
    line(fit(mode, width));
    const std::string rule((size_t)width, '-');
    line(style::dim + rule + style::reset);

    // the headers, with the sort direction
    const std::string direction = view.descending ? "v" : "^";
    std::string header = fit("GAME" + std::string(view.sorted && !view.sort_column ? " " + direction : ""), name_width) + " |";
    for (size_t s = first_column; s < last_column; s++) {
        if (view.sorted && view.sort_column == s) {
            header += std::string(style::cyan) + fit(solvers[s].name, cell_width - 3) + " " + direction + style::reset + style::bold + " ";
        } else {
            header += fit(solvers[s].name, cell_width - 1) + " ";
        }
    }
    line(style::bold + header + style::reset);

    // the games
    const size_t body = (size_t)std::max(1, height - 11);
    view.body = body;
    view.first = ui::keep_visible(view.first, view.selected, games.size(), body);
    for (size_t i = 0; i < body; i++) {
        const size_t row = view.first + i;
        if (row >= games.size()) {
            line("");
            continue;
        }
        const Case& game = games[view.order[row]];
        const bool selected = row == view.selected;
        std::string text = (selected ? style::title : "") + fit((selected ? "> " : "  ") + game.name, name_width) + style::reset + " |";
        std::vector<double> medians(solvers.size(), std::numeric_limits<double>::quiet_NaN());
        double best = std::numeric_limits<double>::infinity();
        for (size_t s = 0; s < solvers.size(); s++) {
            const Cell& c = game.cells[s];
            if (!game.timed() || !c.complete(o.repeats)) continue;
            medians[s] = c.stats().median;
            best = std::min(best, medians[s]);
        }
        const auto colors = view.heatmap ? heatmap_colors(medians) : std::vector<int>(solvers.size(), -1);
        for (size_t s = first_column; s < last_column; s++) {
            std::string value;
            std::string code;
            cell(game, s, o.repeats, view.heatmap, colors[s], best, value, code);
            text += code + fit(value, cell_width - 1) + " " + style::reset;
        }
        line(text);
    }
    line(style::dim + rule + style::reset);

    // the totals of the complete cells, coloured by rank in heatmap mode
    std::vector<double> totals(solvers.size(), 0);
    std::vector<size_t> coverage(solvers.size(), 0);
    size_t timed = 0;
    for (const auto& game : games) {
        if (!game.timed()) continue;
        timed++;
        for (size_t s = 0; s < solvers.size(); s++) {
            if (!game.cells[s].complete(o.repeats)) continue;
            totals[s] += game.cells[s].stats().median;
            coverage[s]++;
        }
    }
    size_t preprocessed = games.size() - timed;
    std::vector<double> ranked(solvers.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t s = 0; s < solvers.size(); s++) {
        if (timed && coverage[s] == timed) ranked[s] = totals[s];
    }
    const auto total_colors = view.heatmap ? heatmap_colors(ranked) : std::vector<int>(solvers.size(), -1);
    std::string sums = fit("SUM MEDIANS (complete)", name_width) + " |";
    std::string counts = fit("verified games", name_width) + " |";
    const size_t timed_planned = data.planned > preprocessed ? data.planned - preprocessed : 0;
    for (size_t s = first_column; s < last_column; s++) {
        sums += style::color256(total_colors[s]) + fit(coverage[s] ? duration(totals[s]) : "-", cell_width - 1) + " " + style::reset + style::bold;
        counts += fit(std::to_string(coverage[s]) + "/" + std::to_string(timed_planned), cell_width - 1) + " ";
    }
    line(style::bold + sums + style::reset);
    line(counts);

    // progress, activity and keys
    const double fraction = data.slots ? std::min(1.0, (double)data.done / (double)data.slots) : 0;
    const int bar = std::max(8, width - 44);
    const int filled = (int)(fraction * bar);
    line(std::string(style::cyan) + "[" + std::string((size_t)filled, '=') + std::string((size_t)(bar - filled), ' ') + "]" + style::reset + " " + std::to_string((int)(fraction * 100)) + "%  " + std::to_string(data.done) + "/" + std::to_string(data.slots) + "  " + duration(data.elapsed));
    line(fit((view.paused ? "PAUSED | " : "") + data.activity, width));
    const std::string page = "[" + std::to_string(view.page + 1) + "/" + std::to_string(std::max<size_t>(pages, 1)) + "]";
    std::string keys;
    if (view.review) {
        keys = page + (width < 80 ? " 1-9 sort | o reset | q quit" : " Click headers/1-9 sort | g game | o reset | arrows | h colors | e extract | q quit");
    } else {
        keys = "Space pause/resume | arrows browse " + page + " | h colors | e extract | q quit";
    }
    line(fit(keys, width));
    return lines;
}

}
