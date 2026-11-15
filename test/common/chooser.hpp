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
 * The interactive solver chooser of test_solvers -i and bench_solvers -i.
 *
 * The screen has two panes. AVAILABLE lists the solvers of the registry by
 * family, with their descriptions; experimental solvers are hidden until x is
 * pressed, and / filters the list. SELECTED is the roster in the chosen order,
 * including configured variants such as tlpar[parmin=1,chunk=8]. A
 * configuration is checked by the parser of the solver itself.
 *
 * The caller can add settings of its own (bench_solvers adds repetitions,
 * warm-ups, the timeout and so on) and a reference solver.
 *
 * The Chooser class holds the state and renders it to lines of text; it does no
 * I/O, so the unit tests drive it with keys. choose() runs it on the terminal.
 */

#ifndef OINK_TEST_CHOOSER_HPP
#define OINK_TEST_CHOOSER_HPP

#include <csignal>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "solver_run.hpp"
#include "terminal.hpp"

namespace ui {

/**
 * A line editor for one value: the text, the cursor and keys to edit them.
 */
class LineEditor
{
public:
    /**
     * Validates the text when Enter is pressed: return false and set the error
     * message to keep the editor open.
     */
    using Accept = std::function<bool(const std::string& text, std::string& error)>;

    LineEditor(std::string title, std::string hint, std::string text, Accept accept);

    enum class Result { Editing, Saved, Cancelled };
    Result handle(const Key& key);

    /** Render the editor screen; sets the cursor position. */
    std::vector<std::string> render(int width, bool color, int& cursor_row, int& cursor_column) const;

    const std::string& text() const { return text_; }
    const std::string& message() const { return message_; }

private:
    std::string title;
    std::string hint;
    std::string text_;
    std::string message_;
    size_t cursor;
    Accept accept;
};

/**
 * A setting that the chooser shows and edits besides the solvers, such as the
 * number of repetitions of bench_solvers. Its key opens a line editor; for a
 * toggle, the key flips the setting at once.
 */
struct Setting
{
    char key = 0;                       // the key that edits the setting
    std::string label;                  // its name in the settings line and the key help
    std::string title;                  // the title of the editor
    std::string hint;                   // the explanation shown in the editor
    std::function<std::string()> show;  // the current value, for the settings line
    std::function<std::string()> text;  // the initial text in the editor
    LineEditor::Accept apply;           // validate and apply the edited text
    std::function<std::string()> toggle; // for a toggle: flip it, and return a notice
};

struct ChooserConfig
{
    std::string title = "SOLVER SELECTION";
    std::vector<Setting> settings;
    bool reference = false;             // whether R marks a reference solver
};

class Chooser
{
public:
    enum class Result { Continue, Start, Cancel };

    /**
     * Start with the roster <selected>. The available solvers are <ids>, by
     * default all registered solvers ordered by family.
     */
    Chooser(std::vector<SolverRun> selected, ChooserConfig config);
    Chooser(std::vector<SolverRun> selected, ChooserConfig config, std::vector<std::string> ids);

    // The editors refer to the chooser, so it stays in place.
    Chooser(const Chooser&) = delete;
    Chooser& operator=(const Chooser&) = delete;

    /** Handle one key; Start or Cancel ends the chooser. */
    Result handle(const Key& key);

    /**
     * Render the screen of <width> by <height>, with or without colours. Sets
     * the cursor position (row 0 hides it). This also scrolls the panes.
     */
    std::vector<std::string> render(int width, int height, bool color, int& cursor_row, int& cursor_column);

    const std::vector<SolverRun>& selected() const { return roster; }
    const std::optional<std::string>& reference() const { return reference_; }
    void set_reference(const std::optional<std::string>& name) { reference_ = name; }

    /** The visible entries of AVAILABLE (filtered, without the family headers). */
    const std::vector<std::string>& available() const { return shown; }
    size_t available_cursor() const { return choice; }
    size_t selected_cursor() const { return row; }
    bool selected_pane() const { return right; }
    bool experimental_shown() const { return experimental; }
    bool editing() const { return editor.has_value(); }
    const std::string& notice() const { return notice_; }

private:
    ChooserConfig config;
    std::vector<std::string> ids;       // all solvers, by family
    std::vector<std::string> shown;     // the solvers shown in AVAILABLE
    std::vector<SolverRun> roster;      // SELECTED
    std::optional<std::string> reference_;
    std::string filter;
    std::string notice_;
    bool experimental = false;          // show experimental solvers
    bool right = false;                 // the SELECTED pane has the focus
    size_t choice = 0;                  // cursor in AVAILABLE
    size_t row = 0;                     // cursor in SELECTED
    size_t left_first = 0;              // first visible row of AVAILABLE
    size_t right_first = 0;             // first visible row of SELECTED
    size_t body = 10;                   // number of rows of the panes, from the last render
    std::optional<LineEditor> editor;   // the open editor, if any

    void refresh();
    std::string current() const;
    bool duplicate(const std::string& name, size_t except) const;
    void remove(size_t index);
    void configure(const SolverRun& run, size_t target);
    void edit(const std::string& title, const std::string& hint, const std::string& text, LineEditor::Accept accept);
    Result start();
};

/**
 * Show the chooser on the terminal until the user starts or cancels; true to
 * start. SIGINT and SIGTERM cancel it and set <*quit> (if given).
 */
bool choose(Chooser& chooser, volatile std::sig_atomic_t* quit = nullptr);

/**
 * Check configuration <config> of solver <id> with the solver's own parser;
 * false with the parser's message in <error> if it is rejected.
 */
bool valid_configuration(const std::string& id, const std::string& config, std::string& error);

}

#endif
