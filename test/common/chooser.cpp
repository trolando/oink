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

#include "chooser.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

#include "oink/solvers.hpp"

using pg::Solvers;

namespace ui {

/** <text> in the style <code>, when colours are used. */
static std::string
paint(const std::string& text, const char* code, bool color)
{
    if (!color || code[0] == '\0') return text;
    return code + text + style::reset;
}

static std::string
lower(std::string text)
{
    for (char& c : text) c = (char)std::tolower((unsigned char)c);
    return text;
}

/** Split <text> into lines of at most <width> columns, at spaces where possible. */
static std::vector<std::string>
wrap(const std::string& text, int width)
{
    std::vector<std::string> lines;
    std::istringstream words(text);
    std::string word;
    std::string line;
    while (words >> word) {
        if (!line.empty() && (int)(line.size() + 1 + word.size()) > width) {
            lines.push_back(line);
            line.clear();
        }
        while ((int)word.size() > width) {
            if (!line.empty()) {
                lines.push_back(line);
                line.clear();
            }
            lines.push_back(word.substr(0, (size_t)width));
            word.erase(0, (size_t)width);
        }
        if (!line.empty()) line += ' ';
        line += word;
    }
    if (!line.empty()) lines.push_back(line);
    return lines;
}

bool
valid_configuration(const std::string& id, const std::string& config, std::string& error)
{
    std::ostringstream errors;
    try {
        if (Solvers::validateOptions(id, config, errors)) return true;
        error = errors.str();
    } catch (const std::exception& e) {
        error = e.what();
    }
    while (!error.empty() && std::isspace((unsigned char)error.back())) error.pop_back();
    for (char& c : error) {
        if (c == '\n') c = ' ';
    }
    if (error.empty()) error = "invalid configuration for " + id;
    return false;
}

/*------------------------------------------------------------------------*/

LineEditor::LineEditor(std::string title, std::string hint, std::string text, Accept accept)
    : title(std::move(title)), hint(std::move(hint)), text_(std::move(text)), accept(std::move(accept))
{
    cursor = text_.size();
}

LineEditor::Result
LineEditor::handle(const Key& key)
{
    switch (key.type) {
        case Key::Escape:
            return Result::Cancelled;
        case Key::Enter:
            message_.clear();
            if (!accept || accept(text_, message_)) return Result::Saved;
            break;
        case Key::Left:
            if (cursor > 0) cursor--;
            break;
        case Key::Right:
            if (cursor < text_.size()) cursor++;
            break;
        case Key::Home:
            cursor = 0;
            break;
        case Key::End:
            cursor = text_.size();
            break;
        case Key::Backspace:
            if (cursor > 0) {
                cursor--;
                text_.erase(cursor, 1);
            }
            break;
        case Key::Delete:
            if (cursor < text_.size()) text_.erase(cursor, 1);
            break;
        case Key::ClearLine:
            text_.clear();
            cursor = 0;
            break;
        case Key::Char:
            if (key.ch >= 32 && key.ch < 127) {
                text_.insert(cursor, 1, (char)key.ch);
                cursor++;
                message_.clear();
            }
            break;
        default:
            break;
    }
    return Result::Editing;
}

std::vector<std::string>
LineEditor::render(int width, bool color, int& cursor_row, int& cursor_column) const
{
    std::vector<std::string> rows;
    rows.push_back(paint(fit("OINK  /  " + title, width), style::title, color));
    rows.push_back("");
    for (const auto& line : wrap(hint, width)) rows.push_back(fit(line, width));
    rows.push_back("");
    // scroll the text horizontally to keep the cursor visible
    const size_t room = width > 5 ? (size_t)width - 5 : 1;
    const size_t start = cursor > room ? cursor - room : 0;
    cursor_row = (int)rows.size() + 1;
    cursor_column = (int)(cursor - start) + 3;
    rows.push_back("> " + fit(text_.substr(start), width - 2));
    rows.push_back("");
    rows.push_back(paint(fit(message_, width), style::error, color));
    rows.push_back(fit("Enter: save   Esc: cancel   Ctrl-U: clear", width));
    return rows;
}

/*------------------------------------------------------------------------*/

Chooser::Chooser(std::vector<SolverRun> selected, ChooserConfig config)
    : Chooser(std::move(selected), std::move(config), Solvers::getSolverIDsByFamily())
{
}

Chooser::Chooser(std::vector<SolverRun> selected, ChooserConfig config, std::vector<std::string> ids)
    : config(std::move(config)), ids(std::move(ids)), roster(std::move(selected))
{
    refresh();
}

/**
 * Recompute the solvers shown in AVAILABLE after a change of the filter or of
 * the experimental flag, keeping the cursor on the same solver if it is shown.
 */
void
Chooser::refresh()
{
    const std::string keep = choice < shown.size() ? shown[choice] : "";
    const std::string needle = lower(filter);
    shown.clear();
    for (const auto& id : ids) {
        if (!experimental && Solvers::isExperimental(id)) continue;
        const std::string text = lower(id + " " + Solvers::desc(id) + " " + Solvers::familyName(Solvers::family(id)));
        if (text.find(needle) != std::string::npos) shown.push_back(id);
    }
    const auto found = std::find(shown.begin(), shown.end(), keep);
    if (found != shown.end()) {
        choice = (size_t)(found - shown.begin());
    } else if (shown.empty()) {
        choice = 0;
    } else {
        choice = std::min(choice, shown.size() - 1);
    }
}

/** The solver under the cursor of the active pane, or "" if there is none. */
std::string
Chooser::current() const
{
    if (right) return roster.empty() ? "" : roster[row].id;
    return shown.empty() ? "" : shown[choice];
}

bool
Chooser::duplicate(const std::string& name, size_t except) const
{
    for (size_t i = 0; i < roster.size(); i++) {
        if (i != except && roster[i].name == name) return true;
    }
    return false;
}

void
Chooser::remove(size_t index)
{
    if (reference_ && *reference_ == roster[index].name) reference_.reset();
    notice_ = "Removed " + roster[index].name;
    roster.erase(roster.begin() + (long)index);
    if (row >= roster.size() && row > 0) row = roster.empty() ? 0 : roster.size() - 1;
}

void
Chooser::edit(const std::string& title, const std::string& hint, const std::string& text, LineEditor::Accept accept)
{
    editor.emplace(title, hint, text, std::move(accept));
}

/**
 * Edit the configuration of <run>: of roster row <target>, or of a new variant
 * if <target> is past the end of the roster.
 */
void
Chooser::configure(const SolverRun& run, size_t target)
{
    if (!Solvers::hasOptions(run.id)) {
        notice_ = run.id + " has no solver-specific options.";
        return;
    }
    const bool existing = target < roster.size();
    const std::string title = (existing ? "CONFIGURE " : "NEW VARIANT / ") + run.id;
    edit(title, Solvers::optionsHelp(run.id), run.config, [this, run, target](const std::string& text, std::string& error) {
        const SolverRun updated(run.id, text);
        if (duplicate(updated.name, target)) {
            error = "That solver/configuration is already selected.";
            return false;
        }
        if (!valid_configuration(run.id, text, error)) return false;
        if (target < roster.size()) {
            if (reference_ && *reference_ == roster[target].name) reference_ = updated.name;
            roster[target] = updated;
        } else {
            roster.push_back(updated);
            row = roster.size() - 1;
            right = true;
        }
        notice_ = "Saved " + updated.name;
        return true;
    });
}

/** Check the roster before starting. */
Chooser::Result
Chooser::start()
{
    if (roster.empty()) {
        notice_ = "Choose at least one solver first.";
        return Result::Continue;
    }
    for (size_t i = 0; i < roster.size(); i++) {
        std::string error;
        if (!valid_configuration(roster[i].id, roster[i].config, error)) {
            row = i;
            right = true;
            notice_ = error;
            return Result::Continue;
        }
    }
    return Result::Start;
}

Chooser::Result
Chooser::handle(const Key& key)
{
    if (editor) {
        if (editor->handle(key) != LineEditor::Result::Editing) editor.reset();
        return Result::Continue;
    }

    size_t& active = right ? row : choice;
    const size_t count = right ? roster.size() : shown.size();
    switch (key.type) {
        case Key::Up:
            if (active > 0) active--;
            return Result::Continue;
        case Key::Down:
            if (active + 1 < count) active++;
            return Result::Continue;
        case Key::Home:
            active = 0;
            return Result::Continue;
        case Key::End:
            if (count) active = count - 1;
            return Result::Continue;
        case Key::PageUp:
            active = active > body ? active - body : 0;
            return Result::Continue;
        case Key::PageDown:
            if (count) active = std::min(count - 1, active + body);
            return Result::Continue;
        case Key::Tab:
        case Key::Left:
        case Key::Right:
            right = !right;
            return Result::Continue;
        case Key::Escape:
            return Result::Cancel;
        case Key::Delete:
        case Key::Backspace:
            if (right && !roster.empty()) remove(row);
            return Result::Continue;
        case Key::Enter:
            if (right) {
                if (!roster.empty()) configure(roster[row], row);
                return Result::Continue;
            }
            break; // add, as with Space
        case Key::Char:
            break;
        default:
            return Result::Continue;
    }

    const int c = key.type == Key::Enter ? '\n' : key.ch;
    if (c == 'q') return Result::Cancel;
    if (c == 's') return start();
    if (c == '\n' || c == ' ') {
        if (right) {
            if (!roster.empty()) remove(row);
            return Result::Continue;
        }
        if (shown.empty()) return Result::Continue;
        const std::string id = shown[choice];
        const auto found = std::find_if(roster.begin(), roster.end(), [&](const SolverRun& r) {
            return r.id == id && r.config.empty();
        });
        if (found == roster.end()) {
            roster.emplace_back(id, "");
            notice_ = "Added " + id;
        } else if (c == ' ') {
            remove((size_t)(found - roster.begin()));
        } else {
            row = (size_t)(found - roster.begin());
            right = true;
            notice_ = "Already selected";
        }
        return Result::Continue;
    }
    if (c == '-') {
        if (right && !roster.empty()) remove(row);
        return Result::Continue;
    }
    if (c == 'x') {
        experimental = !experimental;
        refresh();
        notice_ = experimental ? "Experimental solvers shown." : "Experimental solvers hidden; selected solvers are kept.";
        return Result::Continue;
    }
    if (c == 'c' || c == 'd') {
        if (right) {
            if (!roster.empty()) configure(roster[row], c == 'd' ? roster.size() : row);
        } else if (!shown.empty()) {
            configure(SolverRun(shown[choice], ""), roster.size());
        }
        return Result::Continue;
    }
    if (c == '/') {
        edit("FILTER SOLVERS", "Match the solver ID, description or family; empty shows all.", filter, [this](const std::string& text, std::string&) {
            filter = text;
            choice = 0;
            right = false;
            refresh();
            return true;
        });
        return Result::Continue;
    }
    if (c == 'R' && config.reference) {
        if (!right || roster.empty()) {
            notice_ = "Highlight a solver in SELECTED, then press R to set or clear the reference.";
        } else if (reference_ && *reference_ == roster[row].name) {
            reference_.reset();
            notice_ = "Reference comparisons disabled.";
        } else {
            reference_ = roster[row].name;
            notice_ = "Reference: " + *reference_;
        }
        return Result::Continue;
    }
    for (const auto& setting : config.settings) {
        if (c != setting.key) continue;
        if (setting.toggle) {
            notice_ = setting.toggle();
        } else {
            edit(setting.title, setting.hint, setting.text(), setting.apply);
        }
        return Result::Continue;
    }
    return Result::Continue;
}

std::vector<std::string>
Chooser::render(int width, int height, bool color, int& cursor_row, int& cursor_column)
{
    if (editor) return editor->render(width, color, cursor_row, cursor_column);
    cursor_row = 0;
    cursor_column = 0;

    // The settings line shows the edited settings; the status line the toggles.
    std::string settings_line;
    std::string status = std::string("x Experimental: ") + (experimental ? "shown" : "hidden");
    for (const auto& setting : config.settings) {
        if (setting.toggle) {
            status += " | " + setting.label + ": " + setting.show();
        } else {
            if (!settings_line.empty()) settings_line += "   ";
            settings_line += setting.label + " " + setting.show();
        }
    }
    if (config.reference) status += " | Reference: " + reference_.value_or("none");

    const int fixed = settings_line.empty() ? 10 : 11;
    body = (size_t)std::max(1, height - fixed - 1);

    // Width of the panes: SELECTED is compact, AVAILABLE gets the rest, but no
    // more than its unfiltered entries need, so the divider stays in place when
    // filtering.
    const std::string selected_title = "SELECTED (" + std::to_string(roster.size()) + ")";
    int selected_needed = (int)selected_title.size();
    for (size_t i = 0; i < roster.size(); i++) {
        const bool ref = reference_ && *reference_ == roster[i].name;
        const int length = 3 + (int)std::to_string(i + 1).size() + (int)roster[i].name.size() + (ref ? 6 : 0);
        selected_needed = std::max(selected_needed, length);
    }
    int available_needed = (int)std::string("AVAILABLE").size();
    for (const auto& id : ids) {
        if (!experimental && Solvers::isExperimental(id)) continue;
        available_needed = std::max(available_needed, 21 + (int)Solvers::desc(id).size());
        available_needed = std::max(available_needed, 12 + (int)Solvers::familyName(Solvers::family(id)).size());
    }
    const int third = (width - 3) / 3;
    const int selected_width = std::clamp(selected_needed, std::min(14, third), third);
    const int left_width = std::min(available_needed, width - 3 - selected_width);
    const int right_width = width - 3 - left_width;

    // AVAILABLE with a header row per family; the cursor skips the headers.
    struct Row {
        std::string text;
        bool header;
        size_t index;
    };
    std::vector<Row> rows;
    std::vector<size_t> position(shown.size(), 0);
    for (size_t i = 0; i < shown.size(); i++) {
        const auto family = Solvers::family(shown[i]);
        if (i == 0 || Solvers::family(shown[i-1]) != family) rows.push_back({Solvers::familyName(family), true, 0});
        position[i] = rows.size();
        rows.push_back({shown[i], false, i});
    }
    const size_t visual = shown.empty() ? 0 : position[choice];
    left_first = choice ? keep_visible(left_first, visual, rows.size(), body) : 0;
    right_first = keep_visible(right_first, row, roster.size(), body);

    std::vector<std::string> lines;
    lines.push_back(paint(fit("OINK  /  " + config.title, width), style::title, color));
    if (!settings_line.empty()) lines.push_back(fit(settings_line, width));
    lines.push_back(fit(status, width));
    lines.push_back(std::string((size_t)width, '-'));
    const std::string family = shown.empty() ? "" : Solvers::familyName(Solvers::family(shown[choice]));
    const std::string heading = "AVAILABLE / " + (filter.empty() ? family : filter);
    lines.push_back(paint(fit(heading, left_width), right ? "" : style::title, color) + " | " + paint(fit(selected_title, right_width), right ? style::title : "", color));
    for (size_t i = 0; i < body; i++) {
        std::string lhs;
        const char* left_style = "";
        const size_t l = left_first + i;
        if (l < rows.size()) {
            const Row& entry = rows[l];
            if (entry.header) {
                lhs = "-- " + entry.text + " --";
                left_style = style::title;
            } else {
                const bool checked = std::any_of(roster.begin(), roster.end(), [&](const SolverRun& r) { return r.id == entry.text; });
                const bool active = entry.index == choice && !right;
                lhs = std::string(active ? "> " : "  ") + (checked ? "[x] " : "[ ] ") + fit(entry.text, 14) + " " + Solvers::desc(entry.text);
                if (active) left_style = style::reverse;
            }
        }
        std::string rhs;
        const char* right_style = "";
        const size_t r = right_first + i;
        if (r < roster.size()) {
            const bool active = r == row && right;
            const bool ref = reference_ && *reference_ == roster[r].name;
            rhs = std::string(active ? "> " : "  ") + std::to_string(r + 1) + " " + (ref ? "[ref] " : "") + roster[r].name;
            if (active) right_style = style::reverse;
        }
        lines.push_back(paint(fit(lhs, left_width), left_style, color) + " | " + paint(fit(rhs, right_width), right_style, color));
    }
    lines.push_back(std::string((size_t)width, '-'));

    const std::string id = current();
    if (id.empty()) {
        lines.push_back(fit("No matching solvers; x shows experimental entries", width));
        lines.push_back(fit("", width));
    } else {
        lines.push_back(fit(id + ": " + Solvers::desc(id) + (Solvers::isExperimental(id) ? " (experimental)" : ""), width));
        if (Solvers::hasOptions(id)) {
            lines.push_back(fit("c: configure / d: variant | Options: " + Solvers::optionsHelp(id), width));
        } else {
            lines.push_back(fit("No solver-specific configuration", width));
        }
    }
    lines.push_back(paint(fit(notice_, width), style::notice, color));
    if (right) {
        lines.push_back(fit("Up/Down select   Tab available   Enter configure   Space/Del remove", width));
    } else {
        lines.push_back(fit("Up/Down select   Tab selected   Enter add   Space toggle   / filter", width));
    }
    std::string keys = "s START";
    for (const auto& setting : config.settings) keys += std::string("   ") + setting.key + " " + setting.label;
    if (config.reference) keys += "   R Reference";
    keys += "   q Cancel";
    lines.push_back(paint(fit(keys, width), style::bold, color));
    return lines;
}

/*------------------------------------------------------------------------*/

bool
choose(Chooser& chooser, volatile std::sig_atomic_t* quit)
{
    volatile std::sig_atomic_t local = 0;
    volatile std::sig_atomic_t* flag = quit != nullptr ? quit : &local;
    InterruptGuard guard(flag);
    Terminal terminal;
    const bool color = color_allowed(2, false);
    bool redraw = true;
    int last_width = 0;
    int last_height = 0;
    while (!*flag) {
        int width = 0;
        int height = 0;
        terminal.size(width, height);
        if (redraw || width != last_width || height != last_height) {
            int cursor_row = 0;
            int cursor_column = 0;
            const auto lines = chooser.render(width, height, color, cursor_row, cursor_column);
            terminal.draw(lines, cursor_row, cursor_column);
            last_width = width;
            last_height = height;
        }
        const Key key = terminal.read(100);
        redraw = key.type != Key::None;
        if (!redraw) continue;
        const auto result = chooser.handle(key);
        if (result == Chooser::Result::Start) return true;
        if (result == Chooser::Result::Cancel) return false;
    }
    return false;
}

}
