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

#include "options.hpp"

#include <algorithm>
#include <limits>
#include <random>
#include <stdexcept>

#include "tools/cxxopts.hpp"
#include "common/external_solver.hpp"
#include "oink/solvers.hpp"
#include "worker.hpp"

using pg::Solvers;

namespace bench {

/** The short user guide at the end of --help; see test/bench/README.md. */
static const char* const GUIDE = R"(
Benchmark the selected solvers (default: tl, tlj, tlq, npp, fpi, fpj, zlk) on
the given parity game files and directories, or on random games without them.
Every trial solves a fresh copy of the game and its solution is verified; the
metric is the wall time of Oink::run. The trials of a game are interleaved in
a shuffled order per round: warm-up rounds, three pilot trials that choose the
sample count for --min-solve-time, then the measured rounds. Preprocessing is
off unless --preprocess: then it runs once per game, outside the timing, the
solvers solve what it leaves, and their solutions are mapped back and verified
against the whole game. On POSIX systems the trials run in a worker process,
so crashes are caught and --timeout applies.

The per-game medians are printed to stderr as each game completes, and the
summary (sum of medians, the geometric mean over the virtual best solver (VBS),
VBS wins/ties and, with --reference, the geometric mean speedup and W/T/L) to
stdout. --jsonl writes every trial, --csv the
statistics per game and solver.

With -i, a setup screen comes first: AVAILABLE lists the solvers by family
(x shows the experimental ones, / filters), SELECTED is the roster. Enter
adds a solver or configures a selected one, d adds a configured variant,
Space or Delete removes, R marks the reference; r, w, t, b and p edit the
repetitions, warm-ups, timeout, schedule seed and sampling target, P toggles
preprocessing, and s starts. The live table and the review follow.

With --live, the table of medians updates as the trials finish: the selected
game is marked with a cyan >, the fastest complete median of a row is green,
failures are red, and h (or --heatmap) colours every timing by its rank in the
row. The keys work between trials: up/down (k/j) select a game, left/right
([ ]) page the solver columns, e extracts the selected game; browsing pauses
the measurements and Space resumes them; q or Ctrl-C stops the benchmark and
keeps the results so far. With --review the table stays open at the end: click
a solver header (or press 1-9) to sort the games by it, again to reverse; g
sorts by name and o restores the order. Without a terminal, or with --plain,
TERM=dumb or NO_COLOR, the output is plain.

A named experiment (--experiment=DIR) saves its plan, every finished game and
the reports in DIR, and shows only a progress bar. If it is interrupted,
--resume=DIR (or --experiment=DIR again, on a terminal) continues it with the
saved games, solvers and settings: finished games are kept, an interrupted
game is measured again from its warm-ups. The program, the host and the input
files must be unchanged. -i with --experiment shows the setup screen first.

Exit status: 0 on success, 1 if any trial failed, 130 when interrupted, 3 for
an invalid command line. The user guide is test/bench/README.md.
)";

/** Whether option <key> was given on the command line. */
static bool
given(const cxxopts::ParseResult& result, const char* key)
{
    return result.count(key) != 0;
}

bool
parse_command_line(int argc, const char* const* argv, Options& o, std::vector<SolverRun>& solvers, std::ostream& out)
{
    cxxopts::Options opts("bench_solvers", "Benchmark parity game solvers");
    opts.custom_help("[OPTIONS...] [FILES/DIRECTORIES...]");
    opts.add_options()
        ("help", "Print the help")
        ;
    opts.add_options("Benchmark")
        ("repeat", "Measured repetitions per game and solver (the minimum, with calibration)", cxxopts::value<int>()->default_value("7"))
        ("warmup", "Warm-up rounds per game, not part of the results", cxxopts::value<int>()->default_value("1"))
        ("min-solve-time", "Calibrate the sample count to about this many seconds of solving per game and solver; 0 disables", cxxopts::value<double>()->default_value("0.001"))
        ("max-repeat", "The maximum calibrated sample count", cxxopts::value<int>()->default_value("10000"))
        ("quick", "One repetition, no warm-up and no calibration (for testing)")
        ("timeout", "Wall seconds per solve in the worker process; 0 disables (default 30)", cxxopts::value<double>())
        ("safety-timeout", "Wall seconds each for worker setup and verification (default 300)", cxxopts::value<double>())
        ("in-process", "Run the trials in this process: no crash isolation and no timeout")
        ("bench-seed", "Seed of the schedule: file order, solver order per round and solver seeds (default: random, reported)", cxxopts::value<uint64_t>())
        ("vary-solver-seed", "Vary the solver RNG seed between rounds (default: fixed per game)")
        ("reference", "Compare the solvers with this one: a selected name, its id, or first (default: none)", cxxopts::value<std::string>())
        ("sort", "Benchmark the files in sorted order instead of shuffled")
        ("recursive", "Search input directories recursively")
        ;
    opts.add_options("Display")
        ("i,interactive", "Choose the solvers and settings in a setup screen, then show the live table and the review")
        ("live", "Show the live table on the terminal (plain output without one)")
        ("review", "Keep the live table open after the benchmark, to sort, browse and extract")
        ("heatmap", "Colour all timings of the live table by their rank in the game (h toggles)")
        ("plain", "Plain output without escape sequences")
        ("exports", "Directory for the games extracted from the live table (e)", cxxopts::value<std::string>()->default_value("benchmark-extracts"))
        ;
    opts.add_options("Experiments")
        ("experiment", "Run a named experiment in this directory, with simple progress; offers to resume it if it exists", cxxopts::value<std::string>())
        ("resume", "Resume the named experiment in this directory, with its saved games, solvers and settings", cxxopts::value<std::string>())
        ;
    opts.add_options("Reports")
        ("jsonl", "Write the metadata and every trial as JSON Lines", cxxopts::value<std::string>())
        ("csv", "Write the statistics per game and solver as CSV", cxxopts::value<std::string>())
        ("overwrite", "Replace existing report files")
        ("label", "A label for the run, shown and reported", cxxopts::value<std::string>())
        ("artifacts", "Directory for the games and logs of failed trials", cxxopts::value<std::string>()->default_value("benchmark-artifacts"))
        ;
    opts.add_options("Preprocessing")
        ("preprocess", "Preprocess every game once before its trials, outside the timing: single-parity games, self-loops and winner-controlled winning cycles (default: off)", cxxopts::value<bool>()->implicit_value("true"))
        ("no-preprocess", "Do not preprocess (the default)")
        ("single", "Preprocessing recipe: solve single-parity games (implies --preprocess)")
        ("loops", "Preprocessing recipe: solve self-loops (implies --preprocess)")
        ("wcwc", "Preprocessing recipe: solve winner-controlled winning cycles (implies --preprocess)")
        ("inflate", "Preprocessing recipe: inflate the priorities instead of renumbering them")
        ("compress", "Preprocessing recipe: compress the priorities instead of renumbering them")
        ;
    opts.add_options("Random games")
        ("count", "Number of random games", cxxopts::value<int>()->default_value("100"))
        ("seed", "Seed of the series of random games (default: derived from the schedule seed)", cxxopts::value<uint32_t>())
        ("size", "Vertices of each random game", cxxopts::value<int>()->default_value("100"))
        ("maxp", "Maximum priority of each random game (default: <size>)", cxxopts::value<int>())
        ("maxe", "Maximum number of edges of each random game (default: 4*<size>)", cxxopts::value<int>())
        ;
    opts.add_options("Solving")
        ("w,workers", "Workers for parallel solvers: -1 sequential, 0 autodetect", cxxopts::value<int>()->default_value("-1"))
        ("force-parallel", "Run parallel solvers in parallel also on games below their size threshold")
        ("scc", "Solve one bottom SCC at a time (part of the measured time)")
        ("t,trace", "Capture the solver log at this trace level, kept with failures", cxxopts::value<int>())
        ("c,configure", "Configuration of the single selected solver (otherwise use --<solver>=<options>)", cxxopts::value<std::string>())
        ;
    opts.add_options("Solvers")
        ("all", "Benchmark all solvers")
        ("e,external", "External solver as name:command, e.g. 'mysolv:mysolver %I %O'", cxxopts::value<std::vector<std::string>>())
        ;
    // only solvers with options take a configuration: --tlpar=parmin=1,chunk=8
    for (const auto& id : Solvers::getSolverIDs()) {
        if (Solvers::hasOptions(id)) {
            opts.add_options("Solvers")(id, Solvers::desc(id) + "; options: " + Solvers::optionsHelp(id), cxxopts::value<std::string>()->implicit_value(""));
        } else {
            opts.add_options("Solvers")(id, Solvers::desc(id));
        }
    }
    opts.allow_unrecognised_options();

    const auto result = opts.parse(argc, argv);
    if (given(result, "help")) {
        out << opts.help({"", "Benchmark", "Display", "Experiments", "Reports", "Preprocessing", "Random games", "Solving", "Solvers"}) << GUIDE << std::endl;
        return false;
    }

    o.argv.assign(argv, argv + argc);

    // named experiments
    if (given(result, "experiment") && given(result, "resume")) throw std::runtime_error("choose --experiment or --resume");
    if (given(result, "resume")) {
        // the saved plan decides everything
        for (const auto& arg : result.arguments()) {
            if (arg.key() != "resume" && arg.key() != "plain") throw std::runtime_error("--resume uses the saved games, solvers and settings; only --plain may be added");
        }
        if (!result.unmatched().empty()) throw std::runtime_error("--resume uses the saved games; give no inputs");
        o.experiment = result["resume"].as<std::string>();
        o.resume = true;
        o.plain = given(result, "plain");
        return true;
    }
    if (given(result, "experiment")) {
        o.experiment = result["experiment"].as<std::string>();
        for (const char* key : {"jsonl", "csv", "overwrite"}) {
            if (given(result, key)) throw std::runtime_error(std::string("a named experiment writes its own reports; --") + key + " cannot be used with --experiment");
        }
        if (given(result, "in-process") && Worker::supported()) throw std::runtime_error("a named experiment runs its trials in worker processes; --in-process cannot be used with --experiment");
    }
    if ((given(result, "experiment") || given(result, "resume")) && o.experiment.empty()) throw std::runtime_error("an experiment needs a directory name");
    for (const auto& item : result.unmatched()) {
        if (item.size() > 1 && item[0] == '-') throw std::runtime_error("unknown option: " + item);
        o.inputs.push_back(item);
    }

    // sampling
    o.repeats = result["repeat"].as<int>();
    o.warmups = result["warmup"].as<int>();
    o.min_solve_time = result["min-solve-time"].as<double>();
    o.max_repeats = result["max-repeat"].as<int>();
    if (given(result, "quick")) {
        if (!given(result, "repeat")) o.repeats = 1;
        if (!given(result, "warmup")) o.warmups = 0;
        if (!given(result, "min-solve-time")) o.min_solve_time = 0;
    }

    // isolation: worker processes where the platform has them
    o.isolate = Worker::supported() && !given(result, "in-process");
    if (!o.isolate) o.timeout = o.safety_timeout = 0;
    if (given(result, "timeout")) o.timeout = result["timeout"].as<double>();
    if (given(result, "safety-timeout")) o.safety_timeout = result["safety-timeout"].as<double>();

    if (given(result, "bench-seed")) {
        o.seed = result["bench-seed"].as<uint64_t>();
    } else {
        std::random_device random;
        o.seed = ((uint64_t)random() << 32) ^ (uint64_t)random();
    }
    o.vary_solver_seed = given(result, "vary-solver-seed");
    if (given(result, "reference")) o.reference = result["reference"].as<std::string>();
    o.sort = given(result, "sort");
    o.recursive = given(result, "recursive");

    // display
    o.live = given(result, "live");
    o.review = given(result, "review");
    o.heatmap = given(result, "heatmap");
    o.plain = given(result, "plain");
    o.exports = result["exports"].as<std::string>();
    o.interactive = given(result, "interactive");
    if (o.interactive) {
        o.live = true;
        o.review = true;
    }
    if (o.review && !o.live) throw std::runtime_error("--review needs the live table (--live)");

    // reports
    if (given(result, "jsonl")) o.jsonl = result["jsonl"].as<std::string>();
    if (given(result, "csv")) o.csv = result["csv"].as<std::string>();
    o.overwrite = given(result, "overwrite");
    if (given(result, "label")) o.label = result["label"].as<std::string>();
    o.artifacts = result["artifacts"].as<std::string>();

    // preprocessing: --preprocess gives the default recipe; recipe flags choose their own
    const bool recipe = given(result, "single") || given(result, "loops") || given(result, "wcwc") || given(result, "inflate") || given(result, "compress");
    o.preprocess = given(result, "preprocess") ? result["preprocess"].as<bool>() : recipe;
    if (given(result, "no-preprocess")) {
        if (given(result, "preprocess") && o.preprocess) throw std::runtime_error("choose --preprocess or --no-preprocess");
        o.preprocess = false;
    }
    if (recipe) {
        if (!o.preprocess) throw std::runtime_error("preprocessing is off, but preprocessing steps were given");
        o.single = given(result, "single");
        o.loops = given(result, "loops");
        o.wcwc = given(result, "wcwc");
        o.inflate = given(result, "inflate");
        o.compress = given(result, "compress");
    }

    // random games
    o.count = result["count"].as<int>();
    o.size = result["size"].as<int>();
    o.maxp = given(result, "maxp") ? result["maxp"].as<int>() : o.size;
    if (given(result, "maxe")) {
        o.maxe = result["maxe"].as<int>();
    } else {
        if (o.size > std::numeric_limits<int>::max() / 4) throw std::runtime_error("--size is too large for the default --maxe");
        o.maxe = 4 * o.size;
    }
    if (o.size > 0 && (int64_t)o.maxe > (int64_t)o.size * o.size) o.maxe = (int)std::min<int64_t>((int64_t)o.size * o.size, std::numeric_limits<int>::max());
    o.game_seed_given = given(result, "seed");
    o.game_seed = o.game_seed_given ? result["seed"].as<uint32_t>() : series_seed(o.seed);

    // solving
    o.workers = result["workers"].as<int>();
    o.force_parallel = given(result, "force-parallel");
    o.scc = given(result, "scc");
    if (given(result, "trace")) o.trace = result["trace"].as<int>();

    // the solvers, in the order given (with --all, all solvers first)
    solvers.clear();
    auto add = [&](const std::string& id, const std::string& config) {
        const SolverRun run(id, config);
        for (const auto& r : solvers) {
            if (r.name == run.name) return; // given twice
        }
        solvers.push_back(run);
    };
    const auto ids = Solvers::getSolverIDs();
    if (given(result, "all")) {
        for (const auto& id : ids) add(id, "");
    }
    for (const auto& arg : result.arguments()) {
        if (ids.count(arg.key())) add(arg.key(), Solvers::hasOptions(arg.key()) ? arg.value() : "");
    }
    if (given(result, "external")) {
        for (const auto& spec : result["external"].as<std::vector<std::string>>()) {
            o.externals.push_back(spec);
            add(register_external_solver(spec), "");
        }
    }
    if (given(result, "configure")) {
        if (solvers.size() != 1) throw std::runtime_error("-c/--configure needs exactly one selected solver; use --<solver>=<options>");
        if (!solvers[0].config.empty()) throw std::runtime_error("configuration given both with --" + solvers[0].id + "= and with -c");
        solvers[0] = SolverRun(solvers[0].id, result["configure"].as<std::string>());
    }
    if (solvers.empty()) solvers = default_solvers();
    return true;
}

}
