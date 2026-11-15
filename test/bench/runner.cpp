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

#include "runner.hpp"

#include <algorithm>
#include <chrono>
#include <csignal>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>

#include "common/terminal.hpp"
#include "display.hpp"
#include "json.hpp"
#include "lace.h"
#include "live.hpp"
#include "oink/game.hpp"
#include "oink/io.hpp"
#include "oink/pgparser.hpp"
#include "report.hpp"
#include "summary.hpp"
#include "trial.hpp"
#include "worker.hpp"

namespace fs = std::filesystem;

namespace bench {

using Clock = std::chrono::steady_clock;

static double
seconds_since(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

/**
 * Whether a file found in an input directory is a parity game: a known
 * extension, or no extension and a first word "parity".
 */
static bool
parity_file(const fs::path& path)
{
    const std::string ext = path.extension().string();
    if (ext.empty()) {
        std::ifstream in(path);
        std::string keyword;
        in >> keyword;
        return keyword == "parity";
    }
    return ext == ".pg" || ext == ".gm" || ext == ".pgsolver" || ext == ".gz" || ext == ".bz2" || ext == ".xz";
}

std::vector<PlannedGame>
plan_games(const Options& o, const std::string& executable)
{
    std::vector<PlannedGame> plan;
    if (o.inputs.empty()) {
        std::mt19937 series(o.game_seed);
        for (int i = 0; i < o.count; i++) {
            const uint32_t seed = o.count == 1 ? o.game_seed : (uint32_t)series();
            plan.push_back({"", seed});
        }
        return plan;
    }

    std::vector<fs::path> files;
    for (const auto& input : o.inputs) {
        const fs::path path(input);
        if (!fs::exists(path)) throw std::runtime_error("input path not found: " + input);
        if (fs::is_directory(path)) {
            if (o.recursive) {
                for (const auto& entry : fs::recursive_directory_iterator(path)) {
                    if (entry.is_regular_file() && parity_file(entry.path())) files.push_back(fs::canonical(entry.path()));
                }
            } else {
                for (const auto& entry : fs::directory_iterator(path)) {
                    if (entry.is_regular_file() && parity_file(entry.path())) files.push_back(fs::canonical(entry.path()));
                }
            }
        } else if (fs::is_regular_file(path)) {
            // an explicit file can have any name
            files.push_back(fs::canonical(path));
        } else {
            throw std::runtime_error("input is not a regular file: " + input);
        }
    }
    std::sort(files.begin(), files.end());
    files.erase(std::unique(files.begin(), files.end()), files.end());
    if (files.empty()) throw std::runtime_error("no parity game files found");

    // the reports must not replace an input or the program
    for (const auto& report : {o.jsonl, o.csv}) {
        if (report.empty()) continue;
        const auto absolute = fs::weakly_canonical(fs::absolute(report));
        if (std::find(files.begin(), files.end(), absolute) != files.end()) throw std::runtime_error("the report would replace an input game: " + report);
        if (!executable.empty() && absolute == fs::weakly_canonical(fs::absolute(executable))) throw std::runtime_error("the report would replace the program");
    }

    if (!o.sort) {
        std::mt19937_64 order(o.seed);
        std::shuffle(files.begin(), files.end(), order);
    }
    for (const auto& file : files) plan.push_back({file.string(), std::nullopt});
    return plan;
}

namespace {

/** Starts Lace for the trials in this process, and stops it again. */
struct LaceRunning
{
    bool running = false;

    explicit LaceRunning(int workers)
    {
        if (workers < 0) return;
        lace_start((unsigned)workers, 10000000UL, 0);
        running = true;
    }

    ~LaceRunning()
    {
        if (running) lace_stop();
    }
};

/** A name for the directories of this run: the UTC time and a random tag. */
std::string
make_run_id()
{
    const std::time_t now = std::time(nullptr);
    std::tm value{};
#ifdef _WIN32
    gmtime_s(&value, &now);
#else
    gmtime_r(&now, &value);
#endif
    std::ostringstream out;
    out << std::put_time(&value, "%Y%m%d-%H%M%S") << '-' << std::hex << std::setw(4) << std::setfill('0') << (std::random_device()() & 0xffff);
    return out.str();
}

class Runner
{
public:
    Runner(Options& options, std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& plan, const std::string& executable)
        : o(options), solvers(solvers), plan(plan), executable(executable)
    {
    }

    int run();

private:
    Options& o;
    std::vector<SolverRun>& solvers;
    const std::vector<PlannedGame>& plan;
    std::string executable;
    volatile std::sig_atomic_t quit = 0;
    std::vector<Case> games;
    std::unique_ptr<Report> report;
    std::unique_ptr<Display> display;
    std::string run_id;
    size_t slots = 0;           // planned slots
    size_t done = 0;            // finished slots
    size_t position = 0;        // schedule position of the next sample
    bool failures = false;      // any failure, reflected in the exit status
    size_t current = 0;         // the game being measured
    pg::Game* current_game = nullptr; // and its graph, while it is in memory
    size_t extractions = 0;
    size_t pauses = 0;
    double paused_seconds = 0;

    /** Slots of a cell: one with calibration, else one per repetition. */
    size_t per_cell() const { return o.calibrated() ? 1 : (size_t)o.repeats; }

    Progress progress(const std::string& activity) const { return Progress{games, plan.size(), done, slots, activity}; }
    std::unique_ptr<pg::Game> load(Case& game, const PlannedGame& planned);
    void measure(Case& game, pg::Game& original, const Prepared& prepared);
    bool service(const std::string& activity, bool force);
    void handle_extract();
    void extract(size_t index);
    std::string extract_game(size_t index);
    std::string save_failure(const Case& game, size_t solver, pg::Game& original, const Trial& trial);
    void skip_game(Case& game);
    void report_preparation(const Case& game);
};

/** Read or generate a game, parse and sort it, and fill in what we know of it. */
std::unique_ptr<pg::Game>
Runner::load(Case& game, const PlannedGame& planned)
{
    auto start = Clock::now();
    std::unique_ptr<pg::Game> result;
    if (planned.seed) {
        game.seed = planned.seed;
        game.name = "game " + std::to_string(game.index + 1) + " seed=" + std::to_string(*planned.seed);
        result = std::make_unique<pg::Game>();
        result->set_random_seed(*planned.seed);
        result->init_random_game(o.size, o.maxp, o.maxe - o.size);
        result->ensure_sorted();
        result->renumber();
        game.parse = seconds_since(start);
    } else {
        game.path = planned.path;
        game.name = fs::path(planned.path).filename().string();
        auto in = pg::open_input(planned.path);
        game.load = seconds_since(start);
        start = Clock::now();
        // priorities above INT_MAX are renumbered; self-loops are kept
        result = std::make_unique<pg::Game>(pg::PGParser::parse_pgsolver_renumber(*in, false));
        result->ensure_sorted();
        game.parse = seconds_since(start);
    }
    game.vertices = result->nodecount();
    game.edges = result->edgecount();
    game.hash = hex64(fingerprint(*result));
    return result;
}

/**
 * Show the progress between two trials, handle the input, and wait while the
 * user paused the measurements. Returns false when the benchmark must stop.
 */
bool
Runner::service(const std::string& activity, bool force)
{
    display->update(progress(activity), force);
    handle_extract();
    if (!display->paused() || quit) return !quit;

    // paused by browsing: no trials until Space resumes; the pause is reported
    pauses++;
    Json json;
    json.add("event", "ui_pause");
    json.add("at_utc", utc_now());
    json.add("after_schedule_position", position);
    json.add("game_index", current);
    report->event(json.str());
    const auto start = Clock::now();
    while (display->paused() && !quit) {
        display->update(progress("Measurements paused. Browse or extract, then Space resumes."), true);
        handle_extract();
        if (display->paused() && !quit) display->wait();
    }
    const double seconds = seconds_since(start);
    paused_seconds += seconds;
    Json end;
    end.add("event", "ui_pause_end");
    end.add("at_utc", utc_now());
    end.add("duration_seconds", seconds);
    end.add("resumed", quit == 0);
    end.add("after_schedule_position", position);
    report->event(end.str());
    if (!quit) display->update(progress(activity), true);
    return !quit;
}

void
Runner::handle_extract()
{
    const auto index = display->take_extract();
    if (index) extract(*index);
}

/** Extract game <index>, and tell the user and the report where it went. */
void
Runner::extract(size_t index)
{
    try {
        const std::string path = extract_game(index);
        display->notice("Saved " + path);
        Json json;
        json.add("event", "extract");
        json.add("game_index", index);
        json.add("path", path);
        report->event(json.str());
    } catch (const std::exception& error) {
        display->notice(std::string("Extraction failed: ") + error.what());
    }
}

/**
 * Write game <index> to the extraction directory, with a JSON Lines sidecar
 * that has the metadata, the game, its results and the measured times so far.
 * With preprocessing, the residual game, the partial solution and the vertex
 * map are saved too. A game that is no longer in memory is read (or
 * generated) again and must still have its fingerprint. Returns the path.
 */
std::string
Runner::extract_game(size_t index)
{
    const Case& game = games.at(index);
    if (game.input_failed) throw std::runtime_error("this input could not be read; there is no game to extract");
    std::unique_ptr<pg::Game> loaded;
    pg::Game* graph = index == current ? current_game : nullptr;
    if (graph == nullptr) {
        Case copy;
        copy.index = game.index;
        loaded = load(copy, plan[index]);
        if (copy.hash != game.hash) throw std::runtime_error("the game changed since it was benchmarked: " + game.path);
        graph = loaded.get();
    }

    const fs::path dir = fs::path(o.exports) / run_id;
    fs::create_directories(dir);
    extractions++;
    const std::string stem = (dir / ("game-" + std::to_string(index + 1) + "-" + game.hash + "-" + std::to_string(extractions))).string();
    {
        std::ofstream out(stem + ".pg");
        graph->write_pgsolver(out);
        if (!out) throw std::runtime_error("cannot write " + stem + ".pg");
    }

    Options sidecar = o;
    sidecar.jsonl = stem + ".jsonl";
    sidecar.csv.clear();
    sidecar.overwrite = false;
    Report record(sidecar, solvers);
    record.metadata(1, executable);
    record.game(game);
    if (o.preprocess) {
        // the measured input: what preprocessing leaves of the game
        std::ostringstream log;
        auto prepared = prepare(*graph, o, log);
        {
            std::ofstream out(stem + "-partial.sol");
            pg::write_solution(*graph, prepared->partial, out);
        }
        if (prepared->owned) {
            std::ofstream out(stem + "-residual.pg");
            prepared->owned->write_pgsolver(out);
            std::ofstream map(stem + "-residual-vertices.txt");
            for (int v : prepared->vertices) map << v << '\n';
        }
    }
    record.results(game);
    std::vector<std::string> times;
    for (size_t s = 0; s < solvers.size(); s++) {
        std::vector<std::string> values;
        for (double v : game.cells[s].values) values.push_back(json_number(v));
        Json entry;
        entry.add("solver", solvers[s].name);
        entry.raw("seconds", json_array(values));
        times.push_back(entry.str());
    }
    Json json;
    json.add("event", "extraction");
    json.add("extracted_utc", utc_now());
    json.add("graph_file", stem + ".pg");
    json.add("source_report", o.jsonl);
    json.raw("recorded_solve_seconds", json_array(times));
    record.event(json.str());
    return stem + ".pg";
}

/** Save the game and the log of a failed trial; returns the path of the game. */
std::string
Runner::save_failure(const Case& game, size_t solver, pg::Game& original, const Trial& trial)
{
    if (o.artifacts.empty()) return "";
    try {
        const fs::path dir = fs::path(o.artifacts) / run_id;
        fs::create_directories(dir);
        const std::string stem = (dir / ("case-" + std::to_string(game.index) + "-" + game.hash + "-s" + std::to_string(solver) + "-r" + std::to_string(position))).string();
        std::ofstream graph(stem + ".pg");
        original.write_pgsolver(graph);
        std::ofstream log(stem + ".log");
        log << solvers[solver].name << ": " << status_name(trial.status) << '\n' << trial.message << '\n';
        if (!graph || !log) return "";
        return stem + ".pg";
    } catch (const std::exception&) {
        return "";
    }
}

/** The once-per-game preparation: what it solved and what is left. */
void
Runner::report_preparation(const Case& game)
{
    Json json;
    json.add("event", "preparation");
    json.add("game_index", game.index);
    json.add("enabled", o.preprocess);
    json.add("status", game.prep_failed ? "preparation_error" : "ok");
    if (game.prep_failed) {
        json.null("residual_vertices");
        json.null("residual_edges");
        json.null("solved_vertices");
        json.null("preprocessing_seconds");
        json.null("residual_hash");
    } else {
        json.add("residual_vertices", game.prep.residual_vertices);
        json.add("residual_edges", game.prep.residual_edges);
        json.add("solved_vertices", game.prep.solved_vertices);
        json.add("preprocessing_seconds", game.prep.seconds);
        json.add("residual_hash", game.prep.residual_hash);
    }
    json.add("original_hash", game.hash);
    json.add("message", game.error);
    report->event(json.str());
}

/** A game without trials (unreadable, failed preparation, or fully preprocessed). */
void
Runner::skip_game(Case& game)
{
    if (!game.fully_preprocessed) {
        for (auto& cell : game.cells) cell.failure = Status::SetupError;
        failures = true;
    }
    done += solvers.size() * per_cell();
    report->results(game);
    display->completed(game, plan.size());
}

void
Runner::measure(Case& game, pg::Game& original, const Prepared& prepared)
{
    const uint64_t hash = std::stoull(game.hash, nullptr, 16);
    // the worker is forked with the prepared game; Lace runs in the worker
    Worker worker([&](size_t s, uint32_t seed) { return execute(prepared, solvers[s], o, seed); },
                  [&] {
                      if (o.workers >= 0) lace_start((unsigned)o.workers, 10000000UL, 0);
                  },
                  [&] {
                      if (o.workers >= 0) lace_stop();
                  });

    const int pilots = o.pilots();
    int rounds = o.warmups + pilots + o.repeats; // grows when a cell is calibrated
    bool first = true;
    for (int round = 0; round < rounds && !quit; round++) {
        const bool warmup = round < o.warmups;
        const bool pilot = !warmup && round < o.warmups + pilots;
        const bool measured = !warmup && !pilot;
        const char* phase = warmup ? "warmup" : pilot ? "calibration" : "measured";
        int repetition = round - o.warmups - pilots + 1;
        if (warmup) repetition = round + 1;
        else if (pilot) repetition = round - o.warmups + 1;

        // the solvers of a round run in a shuffled order
        std::vector<size_t> order(solvers.size());
        std::iota(order.begin(), order.end(), 0);
        std::mt19937_64 shuffle(round_seed(o.seed, hash, round));
        std::shuffle(order.begin(), order.end(), shuffle);

        for (size_t s : order) {
            if (quit) break;
            Cell& cell = game.cells[s];
            if (measured && repetition > cell.target(o.repeats)) continue;
            const uint32_t seed = trial_seed(o.seed, hash, o.vary_solver_seed ? round : 0);
            const bool slot = measured && (!o.calibrated() || repetition == cell.target(o.repeats));

            if (cell.failure) {
                // a failed cell is not run again
                Trial skipped;
                skipped.status = *cell.failure;
                report->sample(game, s, position++, repetition, phase, seed, worker.generation(), skipped, false, "");
                if (slot) done++;
                continue;
            }

            const int total = warmup ? o.warmups : pilot ? pilots : cell.target(o.repeats);
            const std::string activity = game.name + " | " + solvers[s].name + " | " + (warmup ? "warm-up " : pilot ? "pilot " : "repeat ") + std::to_string(repetition) + "/" + std::to_string(total);
            if (!service(activity, first)) break;
            first = false;

            Trial trial;
            if (o.isolate) {
                trial = worker.run(s, seed, o.timeout, &quit);
            } else {
                trial = execute(prepared, solvers[s], o, seed);
            }

            std::string artifact;
            if (trial.status == Status::Ok) {
                if (pilot) {
                    cell.pilots.push_back(trial.solve);
                    if ((int)cell.pilots.size() == pilots) {
                        // the pilots fix the sample count of this cell
                        const double estimate = summarize(cell.pilots).median;
                        cell.required = calibrated_repetitions(o.min_solve_time, estimate, o.repeats, o.max_repeats);
                        rounds = std::max(rounds, o.warmups + pilots + cell.required);
                        const bool capped = o.min_solve_time / std::max(estimate, 1e-9) > o.max_repeats;
                        report->calibration(game, s, estimate, cell.required, capped);
                    }
                } else if (warmup) {
                    cell.warmups++;
                } else {
                    cell.values.push_back(trial.solve);
                }
            } else {
                cell.failure = trial.status;
                if (trial.status != Status::Interrupted) {
                    failures = true;
                    artifact = save_failure(game, s, original, trial);
                }
                worker.stop(); // the next trial gets a fresh worker
            }
            report->sample(game, s, position++, repetition, phase, seed, worker.generation(), trial, true, artifact);
            if (slot) done++;
        }
    }
}

int
Runner::run()
{
    ui::InterruptGuard guard(&quit);
    run_id = make_run_id();
    slots = plan.size() * solvers.size() * per_cell();
    report = std::make_unique<Report>(o, solvers);
    report->metadata(plan.size(), executable);
    LiveDisplay* live = nullptr;
    if (o.live && !o.plain && ui::Terminal::available()) {
        auto table = std::make_unique<LiveDisplay>(o, solvers, &quit);
        live = table.get();
        display = std::move(table);
    } else {
        display = std::make_unique<PlainDisplay>(o, solvers);
    }
    // with isolation, Lace runs in the workers; the benchmark itself never forks with threads
    std::unique_ptr<LaceRunning> lace;
    if (!o.isolate) lace = std::make_unique<LaceRunning>(o.workers);

    for (size_t index = 0; index < plan.size() && !quit; index++) {
        games.emplace_back();
        Case& game = games.back();
        game.index = index;
        game.cells.resize(solvers.size());

        std::unique_ptr<pg::Game> original;
        try {
            original = load(game, plan[index]);
        } catch (const std::exception& error) {
            game.input_failed = true;
            game.error = error.what();
            report->game(game);
            skip_game(game);
            continue;
        }
        report->game(game);

        service(game.name + " | preparing", true);
        if (quit) {
            for (auto& cell : game.cells) cell.failure = Status::Interrupted;
            report->results(game);
            break;
        }
        std::unique_ptr<Prepared> prepared;
        try {
            std::ostringstream log;
            prepared = prepare(*original, o, log);
            game.prep = prepared->info;
            if (prepared->identity()) game.prep.residual_hash = game.hash;
        } catch (const std::exception& error) {
            game.prep_failed = true;
            game.error = error.what();
        }
        report_preparation(game);
        if (game.prep_failed) {
            skip_game(game);
            continue;
        }
        if (prepared->residual == nullptr) {
            game.fully_preprocessed = true;
            skip_game(game);
            continue;
        }

        current = index;
        current_game = original.get();
        measure(game, *original, *prepared);
        current_game = nullptr;
        report->results(game);
        if (!quit) display->completed(game, plan.size());
        service(game.name + " complete", true);
    }

    const bool interrupted = quit != 0;
    display->stop_timer();
    handle_extract();
    if (live != nullptr && o.review && !interrupted) {
        live->review(progress("Finished. Browse with the arrows; click the headers to sort. e extracts, q exits."), [&](size_t index) {
            extract(index);
        });
        quit = 0; // leaving the review is a normal end
    }
    display->close();
    Json interaction;
    interaction.add("event", "interaction_summary");
    interaction.add("pause_count", pauses);
    interaction.add("paused_seconds", paused_seconds);
    interaction.add("extractions", extractions);
    report->event(interaction.str());

    const auto totals = compute_totals(games, plan.size(), solvers.size(), o.repeats, reference_index(o.reference, solvers));
    report->summary(totals, interrupted, games.size());
    RunDetails details;
    details.schedule_seed = o.seed;
    if (o.inputs.empty()) details.game_seed = o.game_seed;
    details.preprocess = o.preprocess;
    details.jsonl = o.jsonl;
    details.csv = o.csv;
    details.interrupted = interrupted;
    print_summary(std::cout, totals, solvers, details, ui::color_allowed(1, o.plain));
    if (failures && !o.artifacts.empty() && fs::exists(fs::path(o.artifacts) / run_id)) {
        std::cout << "Failed trials: " << (fs::path(o.artifacts) / run_id).string() << '\n';
    }
    if (pauses) std::cout << "Browsing paused the measurements " << pauses << " times (" << duration(paused_seconds) << "). Use a headless run for measurements to publish.\n";
    if (extractions) std::cout << "Extracted games: " << (fs::path(o.exports) / run_id).string() << '\n';
    std::cout << std::flush;
    if (interrupted) return 130;
    return failures ? 1 : 0;
}

}

int
run_benchmark(Options& options, std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& games, const std::string& executable)
{
    Runner runner(options, solvers, games, executable);
    return runner.run();
}

}
