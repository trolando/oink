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
 * Unit tests for the interactive solver chooser of test_solvers -i and
 * bench_solvers -i. The chooser has no I/O, so the tests press keys and
 * inspect its state and the rendered screen.
 */

#include <algorithm>
#include <iostream>
#include <string>

#include "common/chooser.hpp"
#include "oink/solvers.hpp"

using namespace ui;
using pg::Solvers;

static int failures = 0;

static void
check(const char* name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

static void
press(Chooser& chooser, const Key& key)
{
    chooser.handle(key);
}

static void
press(Chooser& chooser, char c)
{
    chooser.handle(Key::character(c));
}

static void
type(Chooser& chooser, const std::string& text)
{
    for (char c : text) press(chooser, c);
}

static std::vector<std::string>
screen(Chooser& chooser, int width, int height, bool color = false)
{
    int row = 0;
    int column = 0;
    return chooser.render(width, height, color, row, column);
}

static bool
contains(const std::vector<std::string>& lines, const std::string& text)
{
    for (const auto& line : lines) {
        if (strip(line).find(text) != std::string::npos) return true;
    }
    return false;
}

/** Move the cursor of AVAILABLE to solver <id>. */
static bool
go_to(Chooser& chooser, const std::string& id)
{
    if (chooser.selected_pane()) press(chooser, Key(Key::Tab));
    press(chooser, Key(Key::Home));
    for (size_t i = 0; i < chooser.available().size(); i++) {
        if (chooser.available()[chooser.available_cursor()] == id) return true;
        press(chooser, Key(Key::Down));
    }
    return false;
}

static void
test_available()
{
    Chooser chooser({}, ChooserConfig());
    bool experimental = false;
    for (const auto& id : chooser.available()) experimental = experimental || Solvers::isExperimental(id);
    check("experimental solvers hidden by default", !experimental && !chooser.available().empty());
#ifdef OINK_HAVE_EP
    const auto& available = chooser.available();
    for (const char* id : {"ep-clipped", "ep-sign", "gpm-energy", "pmtl-energy"}) {
        check(id, std::find(available.begin(), available.end(), id) != available.end());
    }
    check("bare EP is experimental", Solvers::isExperimental("ep"));
#endif
    const size_t standard = chooser.available().size();
    press(chooser, 'x');
    check("x shows experimental solvers", chooser.experimental_shown() && chooser.available().size() > standard);
#ifdef OINK_HAVE_EP
    check("x shows experimental gpa-energy", Solvers::isExperimental("gpa-energy") &&
          std::find(available.begin(), available.end(), "gpa-energy") != available.end());
#endif
    press(chooser, 'x');
    check("x hides them again", chooser.available().size() == standard);

    // the list is ordered by family, as the registry orders it
    const auto ids = Solvers::getSolverIDsByFamily();
    size_t last = 0;
    bool ordered = true;
    for (const auto& id : chooser.available()) {
        size_t at = 0;
        while (at < ids.size() && ids[at] != id) at++;
        if (at < last) ordered = false;
        last = at;
    }
    check("available ordered by family", ordered);

    // the filter matches the id, description and family name
    press(chooser, '/');
    check("slash opens the filter editor", chooser.editing());
    type(chooser, "recursion");
    press(chooser, Key(Key::Enter));
    bool matches = !chooser.available().empty();
    for (const auto& id : chooser.available()) {
        if (Solvers::family(id) != Solvers::Family::Zielonka) matches = false;
    }
    check("filter by family name", !chooser.editing() && matches);
    press(chooser, '/');
    press(chooser, Key(Key::ClearLine));
    press(chooser, Key(Key::Enter));
    check("empty filter shows all", chooser.available().size() == standard);
    press(chooser, '/');
    type(chooser, "no such solver");
    press(chooser, Key(Key::Enter));
    check("filter without matches", chooser.available().empty() && contains(screen(chooser, 100, 30), "No matching solvers"));
    press(chooser, Key(Key::Enter)); // nothing to add
    check("enter with an empty list adds nothing", chooser.selected().empty());
}

static void
test_roster()
{
    Chooser chooser({SolverRun("tl", "")}, ChooserConfig());
    check("roster from the command line", chooser.selected().size() == 1 && chooser.selected()[0].name == "tl");

    check("cursor finds zlk", go_to(chooser, "zlk"));
    press(chooser, Key(Key::Enter));
    check("enter adds a solver", chooser.selected().size() == 2 && chooser.selected()[1].name == "zlk");
    press(chooser, Key(Key::Enter));
    check("enter on a selected solver jumps to it", chooser.selected().size() == 2 && chooser.selected_pane() && chooser.selected_cursor() == 1);
    go_to(chooser, "zlk");
    press(chooser, ' ');
    check("space removes a selected default", chooser.selected().size() == 1);
    press(chooser, ' ');
    check("space adds it again", chooser.selected().size() == 2);

    // in SELECTED, Delete, Backspace, - and Space remove the row
    press(chooser, Key(Key::Tab));
    press(chooser, Key(Key::Home));
    press(chooser, Key(Key::Delete));
    check("delete removes the row", chooser.selected().size() == 1 && chooser.selected()[0].name == "zlk");
    press(chooser, '-');
    check("minus removes the row", chooser.selected().empty());

    // solvers without options cannot be configured
    go_to(chooser, "tl");
    press(chooser, 'c');
    check("no options for tl", !chooser.editing() && chooser.notice().find("no solver-specific options") != std::string::npos);
}

static void
test_configure()
{
    ChooserConfig config;
    config.reference = true;
    Chooser chooser({SolverRun("tl", "")}, config);

    // a new variant of tlpar, checked by the parser of tlpar itself
    check("cursor finds tlpar", go_to(chooser, "tlpar"));
    press(chooser, 'c');
    check("c opens the configuration editor", chooser.editing());
    type(chooser, "bogus=1");
    press(chooser, Key(Key::Enter));
    check("invalid configuration rejected", chooser.editing() && chooser.selected().size() == 1);
    press(chooser, Key(Key::ClearLine));
    type(chooser, "parmin=1,chunk=8");
    press(chooser, Key(Key::Enter));
    check("valid configuration saved", !chooser.editing() && chooser.selected().size() == 2);
    check("variant name", chooser.selected()[1].name == "tlpar[parmin=1,chunk=8]" && chooser.selected()[1].config == "parmin=1,chunk=8");
    check("focus moves to the new variant", chooser.selected_pane() && chooser.selected_cursor() == 1);

    // an identical second variant is rejected; Escape cancels
    press(chooser, 'd');
    check("d opens a new variant editor", chooser.editing());
    press(chooser, Key(Key::Enter));
    check("duplicate variant rejected", chooser.editing() && chooser.selected().size() == 2);
    press(chooser, Key(Key::Escape));
    check("escape cancels the editor", !chooser.editing() && chooser.selected().size() == 2);

    // the reference follows configuration edits and is cleared on removal
    press(chooser, 'R');
    check("R sets the reference", chooser.reference() && *chooser.reference() == "tlpar[parmin=1,chunk=8]");
    check("reference marked", contains(screen(chooser, 120, 30), "[ref] tlpar[parmin=1,chunk=8]"));
    press(chooser, Key(Key::Enter));
    press(chooser, Key(Key::ClearLine));
    type(chooser, "parmin=2,chunk=8");
    press(chooser, Key(Key::Enter));
    check("enter reconfigures the row", chooser.selected()[1].name == "tlpar[parmin=2,chunk=8]");
    check("reference follows the edit", chooser.reference() && *chooser.reference() == "tlpar[parmin=2,chunk=8]");
    press(chooser, 'R');
    check("R clears the reference", !chooser.reference());
    press(chooser, 'R');
    press(chooser, Key(Key::Delete));
    check("removal clears the reference", !chooser.reference() && chooser.selected().size() == 1);

    // R needs a row in SELECTED
    press(chooser, Key(Key::Tab));
    press(chooser, 'R');
    check("R in AVAILABLE only explains", !chooser.reference() && chooser.notice().find("Highlight") != std::string::npos);
}

static void
test_settings()
{
    int reps = 7;
    bool prep = false;
    ChooserConfig config;
    config.title = "BENCHMARK SETUP";
    Setting repetitions;
    repetitions.key = 'r';
    repetitions.label = "Reps";
    repetitions.title = "MEASURED REPETITIONS";
    repetitions.hint = "Enter a positive integer";
    repetitions.show = [&] { return std::to_string(reps); };
    repetitions.text = [&] { return std::to_string(reps); };
    repetitions.apply = [&](const std::string& text, std::string& error) {
        if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos || std::stoi(text) < 1) {
            error = "enter a positive integer";
            return false;
        }
        reps = std::stoi(text);
        return true;
    };
    config.settings.push_back(repetitions);
    Setting preprocess;
    preprocess.key = 'P';
    preprocess.label = "Preprocessing";
    preprocess.show = [&] { return std::string(prep ? "ON" : "OFF"); };
    preprocess.toggle = [&] {
        prep = !prep;
        return std::string(prep ? "Preprocessing ON" : "Preprocessing OFF");
    };
    config.settings.push_back(preprocess);

    Chooser chooser({SolverRun("tl", "")}, config);
    auto lines = screen(chooser, 100, 30);
    check("settings line", contains(lines, "Reps 7"));
    check("toggle in the status line", contains(lines, "Preprocessing: OFF"));
    check("key help lists the settings", contains(lines, "s START   r Reps   P Preprocessing   q Cancel"));

    press(chooser, 'r');
    check("r opens the editor", chooser.editing());
    lines = screen(chooser, 100, 30);
    check("editor title", contains(lines, "OINK  /  MEASURED REPETITIONS"));
    press(chooser, Key(Key::ClearLine));
    type(chooser, "x");
    press(chooser, Key(Key::Enter));
    check("invalid setting rejected", chooser.editing() && reps == 7 && contains(screen(chooser, 100, 30), "enter a positive integer"));
    press(chooser, Key(Key::Backspace));
    type(chooser, "12");
    press(chooser, Key(Key::Enter));
    check("setting applied", !chooser.editing() && reps == 12 && contains(screen(chooser, 100, 30), "Reps 12"));

    press(chooser, 'P');
    check("toggle flips at once", prep && !chooser.editing() && chooser.notice() == "Preprocessing ON");
}

static void
test_start()
{
    Chooser empty({}, ChooserConfig());
    check("start needs a solver", empty.handle(Key::character('s')) == Chooser::Result::Continue && !empty.notice().empty());
    check("q cancels", empty.handle(Key::character('q')) == Chooser::Result::Cancel);
    check("escape cancels", empty.handle(Key(Key::Escape)) == Chooser::Result::Cancel);

    // a bad configuration from the command line is caught before starting
    Chooser bad({SolverRun("tl", ""), SolverRun("tlpar", "bogus")}, ChooserConfig());
    check("bad configuration blocks the start", bad.handle(Key::character('s')) == Chooser::Result::Continue);
    check("bad row selected", bad.selected_pane() && bad.selected_cursor() == 1);

    Chooser good({SolverRun("tl", ""), SolverRun("tlpar", "parmin=1,chunk=4")}, ChooserConfig());
    check("start", good.handle(Key::character('s')) == Chooser::Result::Start);
}

static void
test_render()
{
    Chooser chooser({SolverRun("tl", ""), SolverRun("zlk", "")}, ChooserConfig());
    for (int width : {60, 80, 100, 200}) {
        const auto lines = screen(chooser, width, 30, true);
        check("screen fits the height", (int)lines.size() < 30);
        bool fits = true;
        for (const auto& line : lines) {
            if ((int)strip(line).size() != width) fits = false;
        }
        check("every line has the width", fits);
        // the divider between the panes
        const std::string heading = strip(lines[3]);
        const size_t divider = heading.find(" | SELECTED");
        check("pane headings", heading.find("AVAILABLE") == 0 && divider != std::string::npos);
        if (width <= 80) {
            check("narrow: AVAILABLE gets most of the width", divider > (size_t)(width - 3) / 2);
        } else {
            check("wide: AVAILABLE no wider than its entries", divider < 110);
        }
    }
    const auto lines = screen(chooser, 100, 30, false);
    check("title", strip(lines[0]).find("OINK  /  SOLVER SELECTION") == 0);
    check("family headers", contains(lines, "-- Tangle learning --"));
    check("selected count", contains(lines, "SELECTED (2)"));
    check("checked entries", contains(lines, "[x] tl "));
    bool plain = true;
    for (const auto& line : lines) {
        if (line.find('\033') != std::string::npos) plain = false;
    }
    check("no escape sequences without colour", plain);

    // the cursor stays visible in a short window
    press(chooser, Key(Key::End));
    const auto small = screen(chooser, 100, 14);
    const std::string last = chooser.available().back();
    check("end scrolls to the last solver", contains(small, "> [ ] " + last) || contains(small, "> [x] " + last));
}

int
test_chooser()
{
    test_available();
    test_roster();
    test_configure();
    test_settings();
    test_start();
    test_render();
    if (failures) std::cerr << failures << " chooser test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
