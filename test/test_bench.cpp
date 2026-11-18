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
 * Unit tests for bench_solvers: the statistics, the colours, the summary with
 * its VBS and reference comparisons, the command line, the schedule, single
 * trials and the worker process that isolates them.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "bench/bench.hpp"
#include "bench/experiment.hpp"
#include "bench/json.hpp"
#include "bench/options.hpp"
#include "bench/runner.hpp"
#include "bench/setup.hpp"
#include "bench/summary.hpp"
#include "bench/trial.hpp"
#include "bench/worker.hpp"
#include "common/chooser.hpp"
#include "common/terminal.hpp"
#include "oink/game.hpp"
#include "oink/game_builder.hpp"

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

static bool
near(double a, double b)
{
    return std::abs(a - b) <= 1e-12 * std::max(1.0, std::abs(b));
}

static void
test_statistics()
{
    const Stats odd = summarize({5, 1, 3, 2, 4});
    check("odd median", odd.samples == 5 && near(odd.median, 3));
    check("odd quartiles", near(odd.q1, 2) && near(odd.q3, 4));
    check("odd mad", near(odd.mad, 1));
    const Stats even = summarize({4, 1, 3, 2});
    check("even median interpolates", near(even.median, 2.5));
    check("even quartiles interpolate", near(even.q1, 1.75) && near(even.q3, 3.25));
    check("even mad", near(even.mad, 1));
    const Stats empty = summarize({});
    check("empty sample", empty.samples == 0 && empty.median == 0);

    Cell cell;
    cell.values = {3, 1, 2};
    check("cell statistics", near(cell.stats().median, 2));
    cell.values.push_back(10);
    check("cell statistics follow new values", near(cell.stats().median, 2.5));
    check("cell target", cell.target(7) == 7 && !cell.complete(7) && cell.complete(4));
    cell.required = 4;
    check("calibrated target", cell.target(7) == 4 && cell.complete(7));
    cell.failure = Status::Timeout;
    check("failed cell is not complete", !cell.complete(7));

    check("calibration rounds up", calibrated_repetitions(0.001, 0.0003, 7, 100) == 7);
    check("calibration target", calibrated_repetitions(0.001, 0.00001, 7, 1000) == 100);
    check("calibration cap", calibrated_repetitions(0.001, 1e-9, 7, 50) == 50);
    check("calibration of a slow game", calibrated_repetitions(0.001, 2.0, 3, 50) == 3);

    check("duration us", duration(0.0000123) == "12.3 us");
    check("duration small us", duration(0.00000123) == "1.23 us");
    check("duration ms", duration(0.00456) == "4.56 ms");
    check("duration s", duration(1.5) == "1.50 s");
    check("setting seconds", seconds_setting(30) == "30 s" && seconds_setting(0.001) == "1 ms" && seconds_setting(0) == "off");
}

static void
test_colors()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto colors = heatmap_colors({3, 1, 2, nan, 4});
    check("fastest green", colors[1] == 82);
    check("slowest red", colors[4] == 196);
    check("middle between", colors[0] != 82 && colors[0] != 196 && colors[2] != 82 && colors[2] != 196);
    check("NaN uncoloured", colors[3] == -1);
    const auto ties = heatmap_colors({1, 1, 2});
    check("ties share a colour", ties[0] == ties[1] && ties[0] == 82 && ties[2] == 196);
    const auto equal = heatmap_colors({2, 2, 2});
    check("all equal is neutral", equal[0] == 226 && equal[1] == 226 && equal[2] == 226);
    const auto none = heatmap_colors({nan, -1});
    check("nothing to rank", none[0] == -1 && none[1] == -1);
}

static void
test_reference()
{
    const std::vector<SolverRun> roster = {SolverRun("tl", ""), SolverRun("tlpar", "parmin=1"), SolverRun("tlpar", "parmin=2")};
    check("no reference", !reference_index(std::nullopt, roster));
    check("reference none", !reference_index(std::string("none"), roster));
    check("reference first", reference_index(std::string("first"), roster) == size_t(0));
    check("reference by id", reference_index(std::string("tl"), roster) == size_t(0));
    check("reference by name", reference_index(std::string("tlpar[parmin=2]"), roster) == size_t(2));
    bool ambiguous = false;
    try {
        reference_index(std::string("tlpar"), roster);
    } catch (const std::runtime_error&) {
        ambiguous = true;
    }
    check("ambiguous reference rejected", ambiguous);
    bool unselected = false;
    try {
        reference_index(std::string("zlk"), roster);
    } catch (const std::runtime_error&) {
        unselected = true;
    }
    check("unselected reference rejected", unselected);
}

/** A game with cells whose medians are <medians> (NaN: failed). */
static Case
game_with(size_t index, const std::vector<double>& medians)
{
    Case game;
    game.index = index;
    game.name = "g" + std::to_string(index);
    game.cells.resize(medians.size());
    for (size_t s = 0; s < medians.size(); s++) {
        if (std::isnan(medians[s])) {
            game.cells[s].failure = Status::Timeout;
        } else {
            game.cells[s].values = {medians[s]};
        }
    }
    return game;
}

static void
test_totals()
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<Case> games;
    games.push_back(game_with(0, {1.0, 2.0, 1.005})); // solver 0 and 2 tie for fastest
    games.push_back(game_with(1, {4.0, 2.0, nan}));   // solver 1 uniquely fastest
    games.push_back(game_with(2, {3.0, 3.0, 9.0}));   // 0 and 1 tie
    const auto totals = compute_totals(games, 3, 3, 1, size_t(0), 10);
    check("coverage", totals.solvers[0].complete == 3 && totals.solvers[2].complete == 2 && !totals.complete_coverage);
    check("failures", totals.solvers[2].failed == 1 && totals.solvers[0].failed == 0);
    check("sums", near(totals.solvers[0].sum, 8) && near(totals.solvers[1].sum, 7));
    check("vbs", totals.solvers[0].vbs_wins == 0 && totals.solvers[0].vbs_ties == 2);
    check("vbs unique", totals.solvers[1].vbs_wins == 1 && totals.solvers[1].vbs_ties == 1);
    check("vbs tie band", totals.solvers[2].vbs_ties == 1);
    // against reference solver 0: ratios 0.5, 2, 1
    const SolverTotals& t = totals.solvers[1];
    check("pairs", t.paired == 3 && t.wins == 1 && t.losses == 1 && t.ties == 1);
    check("geometric mean", t.geomean() && near(*t.geomean(), 1.0));
    check("one-sided pairs", totals.solvers[2].paired == 2 && totals.solvers[2].reference_only == 1);
    check("order: full coverage by sum, then the rest", totals.order == std::vector<size_t>({1, 0, 2}));
    // the fastest medians are 1, 2 and 3; the timeout of solver 2 counts as the limit 10
    check("vbs ratio", totals.solvers[0].vbs_games == 3 && near(*totals.solvers[0].vbs_ratio(), std::cbrt(2.0)));
    check("vbs ratio with a timeout", totals.solvers[2].vbs_games == 3 && near(*totals.solvers[2].vbs_ratio(), std::cbrt(1.005 * 5 * 3)));

    const auto none = compute_totals(games, 3, 3, 1, std::nullopt, 0);
    check("no reference, no pairs", none.solvers[1].paired == 0 && !none.solvers[1].geomean());
    check("vbs ratio without a time limit", none.solvers[2].vbs_games == 2 && near(*none.solvers[2].vbs_ratio(), std::sqrt(1.005 * 3)));

    // the printed summary, plain and styled
    const std::vector<SolverRun> roster = {SolverRun("tl", ""), SolverRun("zlk", ""), SolverRun("npp", "")};
    RunDetails details;
    details.schedule_seed = 42;
    std::ostringstream plain;
    print_summary(plain, totals, roster, details, false);
    check("plain summary has no escapes", plain.str().find('\033') == std::string::npos);
    check("summary names the reference", plain.str().find("Reference: tl.") != std::string::npos);
    check("summary has the comparison columns", plain.str().find("geom. speedup") != std::string::npos && plain.str().find("W/T/L") != std::string::npos);
    check("summary has the replay seed", plain.str().find("--bench-seed=42") != std::string::npos);
    check("summary has the vbs ratio", plain.str().find("geom. VBS") != std::string::npos && plain.str().find("1.260x") != std::string::npos);
    check("summary explains the timeout", plain.str().find("counts as the time limit") != std::string::npos);
    std::ostringstream plain_none;
    print_summary(plain_none, none, roster, details, false);
    check("no comparison columns without a reference", plain_none.str().find("W/T/L") == std::string::npos && plain_none.str().find("VBS W/T") != std::string::npos);

    std::vector<Case> full;
    full.push_back(game_with(0, {1.0, 2.0, 3.0}));
    full.push_back(game_with(1, {1.0, 2.0, 3.0}));
    std::ostringstream styled;
    print_summary(styled, compute_totals(full, 2, 3, 1, std::nullopt, 0), roster, details, true);
    check("styled sums are coloured", styled.str().find("\033[38;5;82m") != std::string::npos && styled.str().find("\033[38;5;196m") != std::string::npos);
    check("styled vbs ratios are coloured", styled.str().find("\033[38;5;82m      1.000x") != std::string::npos && styled.str().find("\033[38;5;196m      3.000x") != std::string::npos);
}

static void
test_json()
{
    check("json escapes", json_string("a\"b\\c\n\x01") == "\"a\\\"b\\\\c\\n\\u0001\"");
    Json json;
    json.add("s", "x").add("i", 3).add("b", true).add("d", 0.5).null("n").add("o", std::optional<int>());
    check("json object", json.str() == "{\"s\":\"x\",\"i\":3,\"b\":true,\"d\":0.5,\"n\":null,\"o\":null}");
    check("json non-finite", json_number(std::numeric_limits<double>::infinity()) == "null");
}

/** Parse a command line given as a list of words. */
static bool
parse(const std::vector<std::string>& words, Options& options, std::vector<SolverRun>& solvers)
{
    std::vector<const char*> argv = {"bench_solvers"};
    for (const auto& word : words) argv.push_back(word.c_str());
    std::ostringstream out;
    return parse_command_line((int)argv.size(), argv.data(), options, solvers, out);
}

static void
test_command_line()
{
    {
        Options o;
        std::vector<SolverRun> solvers;
        parse({"some/dir"}, o, solvers);
        check("default roster", solvers.size() == 7 && solvers[0].name == "tl" && solvers[6].name == "zlk");
        check("default sampling", o.repeats == 7 && o.warmups == 1 && near(o.min_solve_time, 0.001) && o.max_repeats == 10000);
        check("inputs", o.inputs == std::vector<std::string>({"some/dir"}));
        check("isolation by default", o.isolate == Worker::supported() && (o.isolate ? near(o.timeout, 30) : o.timeout == 0));
        check("setup and verification safety limit by default", o.safety_timeout == (o.isolate ? 300 : 0));
        check("series seed from the schedule seed", o.game_seed == series_seed(o.seed));
    }
    {
        Options o;
        std::vector<SolverRun> solvers;
        parse({"--zlk", "--tlpar=parmin=1,chunk=8", "--tl", "--quick", "--bench-seed=9", "--seed=4", "--reference=tl"}, o, solvers);
        check("solvers in the given order", solvers.size() == 3 && solvers[0].name == "zlk" && solvers[1].name == "tlpar[parmin=1,chunk=8]" && solvers[2].name == "tl");
        check("quick", o.repeats == 1 && o.warmups == 0 && o.min_solve_time == 0);
        check("seeds", o.seed == 9 && o.game_seed == 4 && o.game_seed_given);
        check("reference option", o.reference && *o.reference == "tl");
    }
    {
        Options o;
        std::vector<SolverRun> solvers;
        parse({"--tlpar", "-c", "parmin=2,chunk=4", "--in-process", "--quick", "--repeat=3"}, o, solvers);
        check("-c configures the single solver", solvers.size() == 1 && solvers[0].name == "tlpar[parmin=2,chunk=4]");
        check("in-process has no timeout", !o.isolate && o.timeout == 0 && o.safety_timeout == 0);
        check("explicit setting overrides quick", o.repeats == 3);
    }
    {
        Options o;
        std::vector<SolverRun> solvers;
        check("help returns false", !parse({"--help"}, o, solvers));
    }
    bool unknown = false;
    try {
        Options o;
        std::vector<SolverRun> solvers;
        parse({"--no-such-option"}, o, solvers);
    } catch (const std::runtime_error&) {
        unknown = true;
    }
    check("unknown option rejected", unknown);

    bool bad_config = false;
    try {
        Options o;
        std::vector<SolverRun> solvers;
        parse({"--tlpar=bogus"}, o, solvers);
        validate(o, solvers);
    } catch (const std::runtime_error&) {
        bad_config = true;
    }
    check("invalid configuration rejected", bad_config);
}

static void
test_schedule()
{
    Options o;
    o.count = 5;
    o.game_seed = 77;
    const auto a = plan_games(o, "");
    const auto b = plan_games(o, "");
    check("random plan", a.size() == 5 && a[0].seed && a[0].path.empty());
    bool same = true;
    for (size_t i = 0; i < a.size(); i++) same = same && a[i].seed == b[i].seed;
    check("random plan is reproducible", same);
    o.count = 1;
    check("a single random game uses the series seed", plan_games(o, "")[0].seed == uint32_t(77));

    check("trial seeds depend on the round", trial_seed(1, 2, 0) != trial_seed(1, 2, 1));
    check("trial seeds depend on the game", trial_seed(1, 2, 0) != trial_seed(1, 3, 0));
    check("round seeds", round_seed(1, 2, 0) != round_seed(1, 2, 1));

    Options files;
    files.inputs = {OINK_TEST_DATA_DIR};
    files.seed = 5;
    const auto first = plan_games(files, "");
    const auto again = plan_games(files, "");
    files.seed = 6;
    const auto other = plan_games(files, "");
    files.sort = true;
    const auto sorted = plan_games(files, "");
    bool reproducible = first.size() == again.size();
    for (size_t i = 0; reproducible && i < first.size(); i++) reproducible = first[i].path == again[i].path;
    bool differs = false;
    for (size_t i = 0; i < first.size() && i < other.size(); i++) differs = differs || first[i].path != other[i].path;
    bool ordered = true;
    for (size_t i = 1; i < sorted.size(); i++) ordered = ordered && sorted[i-1].path < sorted[i].path;
    check("directory inputs found", first.size() > 100);
    check("only parity games", std::none_of(first.begin(), first.end(), [](const PlannedGame& g) { return g.path.find("dot.sh") != std::string::npos; }));
    check("shuffle is reproducible", reproducible);
    check("shuffle depends on the seed", differs);
    check("--sort keeps the sorted order", ordered);
}

/** A small game: an even and an odd cycle connected by a choice. */
static pg::Game
small_game()
{
    pg::GameBuilder b(4);
    b.set_priority(0, 0);
    b.set_owner(0, pg::Player::Even);
    b.add_edge(0, 1);
    b.add_edge(0, 2);
    b.set_priority(1, 2);
    b.set_owner(1, pg::Player::Even);
    b.add_edge(1, 0);
    b.set_priority(2, 3);
    b.set_owner(2, pg::Player::Odd);
    b.add_edge(2, 3);
    b.set_priority(3, 4);
    b.set_owner(3, pg::Player::Even);
    b.add_edge(3, 2);
    b.add_edge(3, 3);
    pg::Game game = b.build();
    game.ensure_sorted();
    return game;
}

static void
test_trials()
{
    const pg::Game game = small_game();
    Options o;
    o.preprocess = false;
    std::ostringstream log;
    const auto prepared = prepare(game, o, log);
    check("prepared without preprocessing", prepared->identity() && prepared->info.residual_vertices == 4);
    for (const char* id : {"tl", "zlk", "fpi", "npp", "psi"}) {
        const Trial trial = execute(*prepared, SolverRun(id, ""), o, 1);
        check("trial verified", trial.status == Status::Ok && trial.solve >= 0 && trial.effective_solver == id);
    }
    const Trial bad = execute(*prepared, SolverRun("mspm", "bogus"), o, 1);
    check("bad options fail the trial", bad.status == Status::SolverError && !bad.message.empty());

#ifndef _WIN32
    volatile std::sig_atomic_t quit = 0;
    ui::InterruptGuard guard(&quit);
    int calls = 0;
    Worker worker([&](size_t solver, uint32_t seed, const PhaseUpdate& update) {
        calls++; // counted in the child only
        if (solver == 4) std::this_thread::sleep_for(std::chrono::seconds(10));
        Trial trial = execute(*prepared, SolverRun("tl", ""), o, seed, [&](Phase phase, const Trial& progress) {
            update(phase, progress);
            if ((solver == 2 || solver == 3) && phase == Phase::Setup) std::this_thread::sleep_for(std::chrono::milliseconds(300));
            if (solver == 1 && phase == Phase::Solve) std::raise(SIGKILL);
            if (solver == 2 && phase == Phase::Solve) std::this_thread::sleep_for(std::chrono::seconds(10));
            if (solver == 3 && phase == Phase::Verify) std::this_thread::sleep_for(std::chrono::milliseconds(300));
            if (solver == 5 && phase == Phase::Verify) std::this_thread::sleep_for(std::chrono::seconds(10));
            if (solver == 6 && phase == Phase::Verify) throw std::bad_alloc();
        });
        if (trial.status == Status::Ok) trial.message = "seed " + std::to_string(seed);
        return trial;
    });
    const Trial ok = worker.run(0, 12, 5, 5, &quit);
    check("worker trial", ok.status == Status::Ok && ok.message == "seed 12" && ok.effective_solver == "tl");
    const Trial again = worker.run(0, 13, 5, 5, &quit);
    check("worker serves several trials", again.status == Status::Ok && worker.generation() == 1);
    check("trials run in the child", calls == 0);
    const Trial crash = worker.run(1, 0, 5, 5, &quit);
    check("worker crash", crash.status == Status::Crash && crash.failure_phase == Phase::Solve && crash.message.find("signal") != std::string::npos && !worker.running());
    const auto start = std::chrono::steady_clock::now();
    const Trial slow = worker.run(2, 0, 0.2, 3, &quit);
    const double waited = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    check("solve timeout after slow setup", slow.status == Status::Timeout && slow.failure_phase == Phase::Solve &&
          slow.setup >= 0.3 && slow.solve < 0 && waited < 5 && worker.generation() == 2);
    const Trial after = worker.run(3, 14, 0.2, 3, &quit);
    check("fresh worker accepts slow setup and verification", after.status == Status::Ok && after.setup >= 0.3 && after.verify >= 0.3 &&
          after.solve < 0.2 && worker.generation() == 3);
    const Trial setup = worker.run(4, 0, 0, 0.2, &quit);
    check("setup safety limit", setup.status == Status::SafetyTimeout && setup.failure_phase == Phase::Setup && !worker.running());
    const Trial verify = worker.run(5, 0, 0, 0.2, &quit);
    check("verification safety limit", verify.status == Status::SafetyTimeout && verify.failure_phase == Phase::Verify &&
          verify.setup >= 0 && verify.solve >= 0 && verify.verify < 0 && !worker.running());
    const Trial memory = worker.run(6, 0, 0, 3, &quit);
    check("verification allocation failure", memory.status == Status::MemoryFailure && memory.failure_phase == Phase::Verify);
    quit = 1;
    const Trial stopped = worker.run(0, 15, 5, 5, &quit);
    check("interrupted worker", stopped.status == Status::Interrupted);
    worker.stop();
#endif
}

/**
 * A game for preprocessing: two vertices with self-loops that preprocessing
 * solves, and a residual of mixed parity (12 edges, 7 in the residual).
 */
static pg::Game
preprocessing_game()
{
    pg::GameBuilder b(6);
    const int priority[] = {1, 2, 3, 4, 5, 6};
    const pg::Player owner[] = {pg::Player::Odd, pg::Player::Even, pg::Player::Odd, pg::Player::Even, pg::Player::Odd, pg::Player::Even};
    for (int v = 0; v < 6; v++) {
        b.set_priority(v, priority[v]);
        b.set_owner(v, owner[v]);
    }
    b.add_edge(0, 0); // odd self-loop, won by odd
    b.add_edge(0, 1);
    b.add_edge(1, 1); // even self-loop, won by even
    b.add_edge(2, 1); // odd avoids the even self-loop
    b.add_edge(2, 3);
    b.add_edge(3, 4);
    b.add_edge(3, 5);
    b.add_edge(4, 3);
    b.add_edge(4, 5);
    b.add_edge(5, 2);
    b.add_edge(5, 3);
    b.add_edge(5, 0); // even avoids the odd self-loop
    pg::Game game = b.build();
    game.ensure_sorted();
    return game;
}

static void
test_preprocessing()
{
    const pg::Game game = preprocessing_game();
    Options o;
    o.preprocess = true;
    o.single = true;
    o.loops = true;
    o.wcwc = false; // keep a residual: the cycles 3-5 and 4-5 are winner-controlled
    std::ostringstream log;
    const auto prepared = prepare(game, o, log);
    check("residual after preprocessing", prepared->residual != nullptr && !prepared->identity());
    check("self-loops solved", prepared->info.solved_vertices == 2 && prepared->info.residual_vertices == 4);
    check("residual edge count", prepared->info.residual_edges == 7);
    bool distinct = true;
    for (size_t i = 0; i < prepared->edges.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            if (prepared->edges[i] != std::numeric_limits<size_t>::max() && prepared->edges[i] == prepared->edges[j]) distinct = false;
        }
    }
    check("edges map to distinct original edges", distinct);
    for (const char* id : {"tl", "zlk", "fpjm", "fpim", "npp"}) {
        const Trial trial = execute(*prepared, SolverRun(id, ""), o, 1);
        if (trial.status != Status::Ok) std::cerr << id << ": " << trial.message << std::endl;
        check("residual solution verified on the original", trial.status == Status::Ok);
    }

    // with every reduction, preprocessing solves this game completely
    o.wcwc = true;
    const auto solved = prepare(game, o, log);
    check("fully preprocessed", solved->residual == nullptr && solved->info.solved_vertices == 6);

    // without preprocessing the residual is the original
    o.preprocess = false;
    const auto off = prepare(game, o, log);
    check("preprocessing off", off->identity() && off->info.solved_vertices == 0);
}

/** Type <text> into the chooser and press Enter, after clearing the editor. */
static void
enter(ui::Chooser& chooser, char key, const std::string& text)
{
    chooser.handle(ui::Key::character(key));
    chooser.handle(ui::Key(ui::Key::ClearLine));
    for (char c : text) chooser.handle(ui::Key::character(c));
    chooser.handle(ui::Key(ui::Key::Enter));
}

static void
test_setup()
{
    Options o;
    o.seed = 5;
    o.game_seed = series_seed(5);
    ui::Chooser chooser({SolverRun("tl", ""), SolverRun("zlk", "")}, setup_config(o));
    int row = 0;
    int column = 0;
    std::string screen;
    for (const auto& line : chooser.render(120, 30, false, row, column)) screen += line + "\n";
    check("setup title", screen.find("OINK  /  BENCHMARK SETUP") != std::string::npos);
    check("setup settings", screen.find("Reps 7   Warm-ups 1   Timeout 30 s   Seed 5   Sampling 1 ms") != std::string::npos);
    check("setup toggles", screen.find("Preprocessing: OFF | Reference: none") != std::string::npos);
    check("setup keys", screen.find("s START   r Reps   w Warm-ups   t Timeout   b Seed   p Sampling   P Preprocessing   R Reference   q Cancel") != std::string::npos);

    enter(chooser, 'r', "3");
    check("r sets the repetitions", o.repeats == 3 && !chooser.editing());
    enter(chooser, 'w', "x");
    check("bad warm-ups rejected", o.warmups == 1 && chooser.editing());
    chooser.handle(ui::Key(ui::Key::Escape));
    enter(chooser, 'w', "0");
    check("w sets the warm-ups", o.warmups == 0);
    enter(chooser, 't', "2.5");
    check("t sets the timeout", near(o.timeout, 2.5));
    enter(chooser, 't', "-1");
    check("negative timeout rejected", near(o.timeout, 2.5) && chooser.editing());
    chooser.handle(ui::Key(ui::Key::Escape));
    enter(chooser, 'b', "99");
    check("b sets the schedule seed and the derived series seed", o.seed == 99 && o.game_seed == series_seed(99));
    enter(chooser, 'p', "0");
    check("p sets the sampling target", o.min_solve_time == 0);
    chooser.handle(ui::Key::character('P'));
    check("P toggles preprocessing", o.preprocess);
    chooser.handle(ui::Key(ui::Key::Tab));
    chooser.handle(ui::Key::character('R'));
    check("R sets the reference", chooser.reference() && *chooser.reference() == "tl");
    screen.clear();
    for (const auto& line : chooser.render(120, 30, false, row, column)) screen += line + "\n";
    check("settings shown", screen.find("Reps 3   Warm-ups 0   Timeout 2.5 s   Seed 99   Sampling off") != std::string::npos);
    check("toggles shown", screen.find("Preprocessing: ON | Reference: tl") != std::string::npos);

    Options local;
    local.isolate = false;
    local.timeout = 0;
    ui::Chooser in_process({SolverRun("tl", "")}, setup_config(local));
    enter(in_process, 't', "5");
    check("no timeout without workers", local.timeout == 0 && in_process.editing());
}

static void
test_experiment_records()
{
    // a finished game survives the round trip through done.txt exactly
    Case game;
    game.index = 7;
    game.name = "a \"quoted\" name";
    game.path = "/some/path with spaces/game.pg";
    game.seed = 1234;
    game.hash = "0123456789abcdef";
    game.vertices = 10;
    game.edges = 22;
    game.load = 0.000125;
    game.parse = 1.0 / 3.0;
    game.error = "";
    game.prep.residual_vertices = 8;
    game.prep.residual_edges = 17;
    game.prep.solved_vertices = 2;
    game.prep.seconds = 0.1;
    game.prep.residual_hash = "fedcba9876543210";
    game.cells.resize(3);
    game.cells[0].values = {1e-6, 2.5e-6, 1.0 / 7.0};
    game.cells[0].warmups = 1;
    game.cells[1].failure = Status::Timeout;
    game.cells[1].failure_phase = Phase::Solve;
    game.cells[2].required = 40;
    game.cells[2].pilots = {0.1, 0.2, 0.3};
    game.cells[2].values.assign(40, 0.25);
    const std::string line = encode_case(game);
    Case back;
    check("decode", decode_case(line, back));
    check("decoded game", back.index == 7 && back.name == game.name && back.path == game.path && back.seed == game.seed && back.hash == game.hash);
    check("decoded sizes and times", back.vertices == 10 && back.edges == 22 && back.load == game.load && back.parse == game.parse);
    check("decoded preparation", back.prep.residual_vertices == 8 && back.prep.solved_vertices == 2 && back.prep.residual_hash == game.prep.residual_hash);
    check("decoded cells", back.cells.size() == 3 && back.cells[0].values == game.cells[0].values && back.cells[0].warmups == 1);
    check("decoded failure", back.cells[1].failure == Status::Timeout && back.cells[1].failure_phase == Phase::Solve);
    check("decoded calibration", back.cells[2].required == 40 && back.cells[2].pilots == game.cells[2].pilots && back.cells[2].values.size() == 40);
    check("a cut line is rejected", !decode_case(line.substr(0, line.size() - 10), back));
    game.seed.reset();
    check("no seed", decode_case(encode_case(game), back) && !back.seed);

    // the saved settings give the same options when they are parsed again
    Options o;
    o.repeats = 5;
    o.warmups = 2;
    o.min_solve_time = 0.0123;
    o.max_repeats = 77;
    o.timeout = 12.5;
    o.safety_timeout = 91.25;
    o.seed = 18446744073709551000ULL;
    o.vary_solver_seed = true;
    o.preprocess = true;
    o.single = true;
    o.loops = false;
    o.wcwc = true;
    o.scc = true;
    o.workers = 3;
    o.force_parallel = true;
    o.trace = 1;
    o.reference = std::string("tl");
    o.label = "label with spaces";
    o.count = 9;
    o.size = 30;
    o.maxp = 12;
    o.maxe = 90;
    o.game_seed = 4242;
    const auto args = settings_arguments(o);
    std::vector<std::string> words(args.begin(), args.end());
    Options parsed;
    std::vector<SolverRun> solvers;
    parse(words, parsed, solvers);
    check("saved sampling", parsed.repeats == 5 && parsed.warmups == 2 && parsed.min_solve_time == o.min_solve_time && parsed.max_repeats == 77);
    check("saved timeout and seeds", parsed.timeout == 12.5 && parsed.seed == o.seed && parsed.vary_solver_seed && parsed.game_seed == 4242 && parsed.game_seed_given);
    check("saved safety timeout", parsed.safety_timeout == o.safety_timeout);
    check("saved recipe", parsed.preprocess && parsed.single && !parsed.loops && parsed.wcwc && !parsed.inflate && !parsed.compress);
    check("saved solving", parsed.scc && parsed.workers == 3 && parsed.force_parallel && parsed.trace == 1);
    check("saved reference and label", parsed.reference && *parsed.reference == "tl" && parsed.label == o.label);
    check("saved random games", parsed.count == 9 && parsed.size == 30 && parsed.maxp == 12 && parsed.maxe == 90);
    o.preprocess = false;
    Options off;
    const auto off_args = settings_arguments(o);
    parse(std::vector<std::string>(off_args.begin(), off_args.end()), off, solvers);
    check("saved: no preprocessing", !off.preprocess);
}

int
test_bench()
{
    test_statistics();
    test_colors();
    test_reference();
    test_totals();
    test_json();
    test_command_line();
    test_schedule();
    test_trials();
    test_preprocessing();
    test_setup();
    test_experiment_records();
    if (failures) std::cerr << failures << " bench test(s) failed" << std::endl;
    return failures ? 1 : 0;
}
