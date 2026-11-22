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

#include "common/external_solver.hpp"
#include "common/terminal.hpp"
#include "display.hpp"
#include "experiment.hpp"
#include "json.hpp"
#include "lace.h"
#include "live.hpp"
#include "oink/game.hpp"
#include "oink/io.hpp"
#include "oink/pgparser.hpp"
#include "options.hpp"
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

/**
 * The size of the task deque of each Lace worker for a game of <vertices> vertices: a task for
 * every vertex, as Oink::run gives it (ZLK's parallel attractor can have that many spawned at
 * once), but at least the 10M tasks that the benchmark always used. Lace reserves the deques
 * lazily, so a large deque costs address space but no memory up front.
 */
size_t
lace_deque_size(long vertices)
{
    return std::max<size_t>(10000000UL, (size_t)std::max(vertices, 0L) + 1000);
}

/** Starts Lace for the trials in this process, and stops it again. */
struct LaceRunning
{
    int workers = -1;
    size_t deque = 0; // the size of the task deques, 0 if Lace does not run

    explicit LaceRunning(int workers, long vertices = 0) : workers(workers)
    {
        if (workers < 0) return;
        deque = lace_deque_size(vertices);
        lace_start((unsigned)workers, deque, 0);
    }

    /** Restart Lace with larger deques if a game of <vertices> vertices needs them. */
    void fit(long vertices)
    {
        if (deque == 0 or lace_deque_size(vertices) <= deque) return;
        lace_stop();
        deque = lace_deque_size(vertices);
        lace_start((unsigned)workers, deque, 0);
    }

    ~LaceRunning()
    {
        if (deque != 0) lace_stop();
    }
};

/** A name for the directories of this run: the UTC time and a random tag. */
std::string
make_run_id()
{
    const std::time_t now = std::time(nullptr);
    std::tm value{};
#ifdef _WIN32
    value = *std::gmtime(&now); // the benchmark is single-threaded here
#else
    gmtime_r(&now, &value);
#endif
    std::ostringstream out;
    out << std::put_time(&value, "%Y%m%d-%H%M%S") << '-' << std::hex << std::setw(4) << std::setfill('0') << (std::random_device()() & 0xffff);
    return out.str();
}

void
identify(Case& game, const PlannedGame& planned)
{
    if (planned.seed) {
        game.seed = planned.seed;
        game.name = "game " + std::to_string(game.index + 1) + " seed=" + std::to_string(*planned.seed);
    } else {
        game.path = planned.path;
        game.name = fs::path(planned.path).filename().string();
    }
}

/** Read or generate a game, parse and sort it, and fill in what we know of it. */
std::unique_ptr<pg::Game>
load(Case& game, const Options& o)
{
    auto start = Clock::now();
    std::unique_ptr<pg::Game> result;
    if (game.seed) {
        result = std::make_unique<pg::Game>();
        result->set_random_seed(*game.seed);
        result->init_random_game(o.size, o.maxp, o.maxe - o.size);
        result->ensure_sorted();
        result->renumber();
        game.parse = seconds_since(start);
    } else {
        auto in = pg::open_input(game.path);
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

std::unique_ptr<Prepared>
prepare_case(Case& game, const pg::Game& original, const Options& o)
{
    std::ostringstream log;
    auto prepared = prepare(original, o, log);
    game.prep = prepared->info;
    if (prepared->identity()) game.prep.residual_hash = game.hash;
    game.fully_preprocessed = prepared->residual == nullptr;
    return prepared;
}

class Runner
{
public:
    Runner(Options& options, std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& plan, const std::string& executable, Experiment* experiment)
        : o(options), solvers(solvers), plan(plan), executable(executable), experiment(experiment)
    {
    }

    int run();

private:
    Options& o;
    std::vector<SolverRun>& solvers;
    const std::vector<PlannedGame>& plan;
    std::string executable;
    Experiment* experiment;     // the named experiment, if any
    volatile std::sig_atomic_t quit = 0;
    std::vector<Case> games;
    std::unique_ptr<Report> report;
    std::unique_ptr<Display> display;
    std::ostringstream round_events;
    std::string game_events;
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
    std::unique_ptr<LaceRunning> lace; // Lace for all trials, without isolation

    /** One progress slot per cell with a solve-time target, else one per trial. */
    size_t per_cell() const { return o.adaptive_sampling() ? 1 : (size_t)o.repeats; }

    Progress progress(const std::string& activity) const { return Progress{games, plan.size(), done, slots, activity}; }
    void measure(Case& game, pg::Game* original, const Prepared* prepared);
    bool service(const std::string& activity, bool force);
    void handle_extract();
    void extract(size_t index);
    std::string extract_game(size_t index);
    std::string save_failure(const Case& game, size_t solver, pg::Game* original, const Trial& trial);
    void skip_game(Case& game, Status status = Status::SetupError, Phase phase = Phase::None);
    void report_preparation(const Case& game);
    Trial run_fresh(Worker& worker, const Case& game, const SolverRun* solver, uint32_t seed);
    bool bench_fresh(Case& game);
    bool bench_game(size_t index);
};

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
        identify(copy, plan[index]);
        loaded = load(copy, o);
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
        for (double v : game.cells[s].samples()) values.push_back(json_number(v));
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
Runner::save_failure(const Case& game, size_t solver, pg::Game* original, const Trial& trial)
{
    if (o.artifacts.empty()) return "";
    try {
        const fs::path dir = fs::path(o.artifacts) / run_id;
        fs::create_directories(dir);
        const std::string stem = (dir / ("case-" + std::to_string(game.index) + "-" + game.hash + "-s" + std::to_string(solver) + "-r" + std::to_string(position))).string();
        std::ofstream log(stem + ".log");
        log << solvers[solver].name << ": " << failure_label(trial.status, trial.failure_phase) << '\n' << trial.message << '\n';
        if (!log) return "";
        if (trial.failure_phase == Phase::Load) return stem + ".log";
        if (!game.path.empty() && original == nullptr) {
            const auto ext = fs::path(game.path).extension().string();
            const std::string path = stem + ".pg" + (ext == ".gz" || ext == ".xz" || ext == ".bz2" ? ext : "");
            fs::copy_file(game.path, path);
            return path;
        }
        std::unique_ptr<pg::Game> generated;
        if (original == nullptr) {
            Case copy;
            copy.seed = game.seed;
            generated = load(copy, o);
            original = generated.get();
        }
        std::ofstream graph(stem + ".pg");
        original->write_pgsolver(graph);
        if (!graph) return "";
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
Runner::skip_game(Case& game, Status status, Phase phase)
{
    if (!game.fully_preprocessed) {
        for (auto& cell : game.cells) {
            cell.failure = status;
            cell.failure_phase = phase;
        }
        if (status != Status::Interrupted) failures = true;
    }
    done += solvers.size() * per_cell();
    report->results(game);
    display->completed(game, plan.size());
}

void
Runner::measure(Case& game, pg::Game* original, const Prepared* prepared)
{
    const uint64_t hash = std::stoull(game.hash, nullptr, 16);
    // without isolation, the trials run here: give Lace deques that fit this game
    if (lace) lace->fit(game.vertices);
    // the worker is forked with the prepared game; Lace runs in the worker
    Worker worker([&](size_t s, uint32_t seed, const PhaseUpdate& update) { return execute(*prepared, solvers[s], o, seed, update); },
                  [&] {
                      if (o.workers >= 0) lace_start((unsigned)o.workers, lace_deque_size(game.vertices), 0);
                  },
                  [&] {
                      if (o.workers >= 0) lace_stop();
                  });

    int rounds = o.warmups + o.repeats;
    const auto observe = [&](int round, const RoundSample& sample) {
        Cell& cell = game.cells[sample.solver];
        const bool warmup = round < o.warmups;
        if (sample.executed && o.isolate && !o.fresh_process) {
            game.worker_peak_rss = std::max(game.worker_peak_rss, sample.peak_rss);
        }
        if (!sample.executed) {
            if (sample.status == Status::Timeout) cell.skipped_timeouts++;
        } else if (sample.status == Status::Ok) {
            if (warmup) {
                cell.warmups++;
            } else {
                cell.add_sample(sample.solve);
                if (o.fresh_process) cell.max_trial_rss = std::max(cell.max_trial_rss, sample.peak_rss);
                cell.max_solve_rss = std::max(cell.max_solve_rss, sample.solve_end_rss);
                cell.extend_sampling(o);
                rounds = std::max(rounds, o.warmups + cell.target(o.repeats));
            }
        } else if (sample.status != Status::Timeout || !warmup) {
            cell.failure = sample.status;
            cell.failure_phase = sample.phase;
            if (sample.status == Status::Timeout) cell.solve_timeouts++;
            if (sample.status != Status::Interrupted) failures = true;
        }
        position++;
        if (!warmup && (!o.adaptive_sampling() || round - o.warmups + 1 == cell.target(o.repeats))) done++;
    };

    int start = 0;
    if (experiment != nullptr) {
        game_events = round_events.str(); // game and preparation from this session
        round_events.str("");
        start = experiment->restore_rounds(game.index, solvers.size(), game_events, observe);
    }
    bool first = true;
    for (int round = start; round < rounds && !quit; round++) {
        const bool warmup = round < o.warmups;
        const char* phase = warmup ? "warmup" : "measured";
        const int repetition = (warmup ? round : round - o.warmups) + 1;
        std::vector<RoundSample> samples;
        bool completed = true;

        // the solvers of a round run in a shuffled order
        std::vector<size_t> order(solvers.size());
        std::iota(order.begin(), order.end(), 0);
        std::mt19937_64 shuffle(round_seed(o.seed, hash, round));
        std::shuffle(order.begin(), order.end(), shuffle);

        for (size_t s : order) {
            Cell& cell = game.cells[s];
            if (!warmup && repetition > cell.target(o.repeats)) continue;
            const uint32_t seed = trial_seed(o.seed, hash, o.vary_solver_seed ? round : 0);

            const bool timeout_skip = cell.failure == Status::Timeout && o.skip_after_timeouts > 0 &&
                                      cell.solve_timeouts >= o.skip_after_timeouts;
            const bool executed = !timeout_skip && (!cell.failure || cell.failure == Status::Timeout);
            const int total = warmup ? o.warmups : o.repeats;
            std::string activity = game.name + " | " + solvers[s].name + " | " + (warmup ? "warm-up " : "repeat ") + std::to_string(repetition);
            if (repetition <= total) activity += "/" + std::to_string(total);
            if (quit || !service(activity, first)) {
                completed = false;
                break;
            }
            first = false;

            Trial trial;
            std::string artifact;
            if (!executed) {
                trial.status = *cell.failure;
                trial.failure_phase = cell.failure_phase;
                if (timeout_skip) {
                    const std::string limit = o.skip_after_timeouts == 1 ? "a timeout" : std::to_string(o.skip_after_timeouts) + " timeouts";
                    trial.message = "not repeated after " + limit;
                }
            } else if (o.fresh_process) {
                trial = run_fresh(worker, game, &solvers[s], seed);
            } else if (o.isolate) {
                trial = worker.run(s, seed, o.timeout, o.safety_timeout, &quit);
            } else {
                trial = execute(*prepared, solvers[s], o, seed);
            }

            if (executed && trial.status != Status::Ok) {
                if (trial.status != Status::Interrupted) artifact = save_failure(game, s, original, trial);
                worker.stop(); // the next trial gets a fresh worker
            }
            report->sample(game, s, position, repetition, phase, seed, worker.generation(), trial, executed, artifact);
            const RoundSample sample{s, executed, trial.status, trial.failure_phase, trial.solve, trial.peak_rss, trial.solve_end_rss};
            observe(round, sample);
            if (experiment != nullptr) samples.push_back(sample);
            if (trial.status == Status::Interrupted) {
                completed = false;
                break;
            }
        }
        if (completed && experiment != nullptr) {
            const auto events = round_events.str();
            experiment->record_round(round, samples, events);
            game_events += events;
            round_events.str("");
        }
    }
}

/** A fresh worker either inspects the game or runs one selected solver. */
Trial
Runner::run_fresh(Worker& worker, const Case& game, const SolverRun* solver, uint32_t seed)
{
    std::string external;
    if (solver != nullptr) {
        for (const auto& spec : o.externals) {
            if (spec.substr(0, spec.find(':')) == solver->id) external = spec;
        }
    }
    std::vector<std::string> command = {
        executable, game.path, game.seed ? std::to_string(*game.seed) : "", game.hash, game.prep.residual_hash,
        solver ? solver->id : "", solver ? solver->config : "", external,
    };
    const auto settings = settings_arguments(o);
    command.insert(command.end(), settings.begin(), settings.end());
    return worker.run(0, seed, o.timeout, o.safety_timeout, &quit, command, experiment ? experiment->lock_fd() : -1);
}

bool
Runner::bench_fresh(Case& game)
{
    service(game.name + " | inspecting input", true);
    Worker inspector({});
    Trial inspection = run_fresh(inspector, game, nullptr, 0);
    if (inspection.status == Status::Ok) {
        Case inspected;
        if (!decode_case(inspection.message, inspected)) throw std::runtime_error("invalid fresh worker inspection result");
        game.hash = std::move(inspected.hash);
        game.vertices = inspected.vertices;
        game.edges = inspected.edges;
        game.prep = std::move(inspected.prep);
        game.fully_preprocessed = inspected.fully_preprocessed;
    } else {
        game.input_failed = inspection.failure_phase == Phase::Load;
        game.prep_failed = !game.input_failed;
        game.error = inspection.message;
    }
    game.load = inspection.load;
    game.parse = inspection.parse;
    report->game(game);
    if (!game.input_failed) report_preparation(game);
    if (inspection.status != Status::Ok || game.fully_preprocessed) {
        skip_game(game, inspection.status, inspection.failure_phase);
        return !quit;
    }

    current = game.index;
    measure(game, nullptr, nullptr);
    report->results(game);
    if (!quit) display->completed(game, plan.size());
    service(game.name + " complete", true);
    return !quit;
}

/** Benchmark game <index>; false if it was interrupted before it finished. */
bool
Runner::bench_game(size_t index)
{
    games.emplace_back();
    Case& game = games.back();
    game.index = index;
    game.cells.resize(solvers.size());
    identify(game, plan[index]);
    if (o.fresh_process) return bench_fresh(game);

    std::unique_ptr<pg::Game> original;
    try {
        original = load(game, o);
    } catch (const std::exception& error) {
        game.input_failed = true;
        game.error = error.what();
        report->game(game);
        skip_game(game);
        return true;
    }
    report->game(game);

    service(game.name + " | preparing", true);
    if (quit) {
        for (auto& cell : game.cells) cell.failure = Status::Interrupted;
        report->results(game);
        return false;
    }
    std::unique_ptr<Prepared> prepared;
    try {
        prepared = prepare_case(game, *original, o);
    } catch (const std::exception& error) {
        game.prep_failed = true;
        game.error = error.what();
    }
    report_preparation(game);
    if (game.prep_failed) {
        skip_game(game);
        return true;
    }
    if (game.fully_preprocessed) {
        skip_game(game);
        return true;
    }

    current = index;
    current_game = original.get();
    measure(game, original.get(), prepared.get());
    current_game = nullptr;
    const bool finished = !quit;
    report->results(game);
    if (finished) display->completed(game, plan.size());
    service(game.name + " complete", true);
    return finished;
}

int
Runner::run()
{
    ui::InterruptGuard guard(&quit);
    run_id = make_run_id();
    slots = plan.size() * solvers.size() * per_cell();
    report = std::make_unique<Report>(o, solvers);
    report->metadata(plan.size(), executable, experiment != nullptr ? std::to_string(experiment->session()) : "");
    LiveDisplay* live = nullptr;
    if (experiment != nullptr) {
        display = std::make_unique<ExperimentDisplay>(o, experiment->directory(), experiment->session(), executable);
    } else if (o.live && !o.plain && ui::Terminal::available()) {
        auto table = std::make_unique<LiveDisplay>(o, solvers, &quit);
        live = table.get();
        display = std::move(table);
    } else {
        display = std::make_unique<PlainDisplay>(o, solvers);
    }
    // with isolation, Lace runs in the workers; the benchmark itself never forks with threads
    if (!o.isolate) lace = std::make_unique<LaceRunning>(o.workers);

    for (size_t index = 0; index < plan.size() && !quit; index++) {
        const Case* saved = experiment != nullptr ? experiment->finished(index) : nullptr;
        if (saved != nullptr) {
            // finished in an earlier session of the experiment
            games.push_back(*saved);
            done += solvers.size() * per_cell();
            for (const auto& cell : saved->cells) {
                if (cell.failure) failures = true;
                if (saved->timed() && !saved->input_failed && !saved->prep_failed) {
                    position += o.warmups + cell.target(o.repeats);
                }
            }
            display->completed(games.back(), plan.size());
            continue;
        }
        // an experiment keeps the events of every game apart
        if (experiment != nullptr) {
            game_events.clear();
            round_events.str("");
            report->copy_to(&round_events);
        }
        const bool finished = bench_game(index);
        report->copy_to(nullptr);
        if (experiment != nullptr) {
            if (finished) experiment->commit(games.back(), game_events + round_events.str());
            else games.pop_back(); // only finished games enter the published results
        }
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

    const auto totals = compute_totals(games, plan.size(), solvers.size(), o.repeats, reference_index(o.reference, solvers), o.timeout);
    report->summary(totals, interrupted, games.size());
    RunDetails details;
    details.schedule_seed = o.seed;
    if (o.inputs.empty()) details.game_seed = o.game_seed;
    details.preprocess = o.preprocess;
    if (experiment == nullptr) {
        details.jsonl = o.jsonl;
        details.csv = o.csv;
    }
    details.interrupted = interrupted;
    print_summary(std::cout, totals, solvers, details, ui::color_allowed(1, o.plain));
    if (experiment != nullptr) {
        experiment->publish(o, solvers, executable, plan);
        const fs::path root(experiment->directory());
        std::cout << "Experiment results: " << (root / "results.jsonl").string() << " and results.csv (" << experiment->finished_count() << "/" << plan.size() << " games)\n";
        if (experiment->finished_count() < plan.size()) std::cout << "Continue with: " << executable << " --resume=" << root.string() << '\n';
    }
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
run_trial_worker(int argc, char** argv)
{
    if (argc < 12) throw std::runtime_error("invalid fresh worker arguments");
    return Worker::serve_once(std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]),
                              [&](size_t, uint32_t seed, const PhaseUpdate& update) {
        Trial result;
        Case game;
        game.load = game.parse = -1;
        Phase phase = Phase::Load;
        auto setup_start = Clock::now();
        try {
            Options options;
            std::vector<SolverRun> unused_solvers;
            std::ostringstream unused_help;
            std::vector<const char*> args = {argv[0]};
            args.insert(args.end(), argv + 12, argv + argc);
            parse_command_line((int)args.size(), args.data(), options, unused_solvers, unused_help);
            game.path = argv[5];
            if (*argv[6]) game.seed = (uint32_t)std::stoul(argv[6]);
            const std::string expected_hash = argv[7];
            const std::string expected_residual = argv[8];
            const auto original = load(game, options);
            if (!expected_hash.empty() && game.hash != expected_hash) {
                throw std::runtime_error("the game changed since it was inspected: " + game.path);
            }
            result.load = game.load;
            result.parse = game.parse;
            phase = Phase::Setup;
            update(phase, result);
            setup_start = Clock::now();
            const bool inspecting = !*argv[9];
            std::unique_ptr<Prepared> prepared;
            if (inspecting && !options.preprocess) {
                // The inspection only needs the graph metadata in this case.
                game.prep.residual_vertices = game.vertices;
                game.prep.residual_edges = game.edges;
                game.prep.residual_hash = game.hash;
            } else {
                prepared = prepare_case(game, *original, options);
            }

            if (inspecting) {
                result.status = Status::Ok;
                result.message = encode_case(game);
                return result;
            }
            if (prepared->residual == nullptr || expected_residual != game.prep.residual_hash) {
                throw std::runtime_error("the prepared game changed since it was inspected");
            }
            if (*argv[11]) register_external_solver(argv[11]);
            const SolverRun solver(argv[9], argv[10]);
            LaceRunning lace(options.workers, game.vertices);
            result.setup = seconds_since(setup_start);
            return execute(*prepared, solver, options, seed, update, result);
        } catch (const std::exception& error) {
            result.load = game.load;
            result.parse = game.parse;
            if (phase == Phase::Setup) result.setup = seconds_since(setup_start);
            record_failure(result, phase, error);
            return result;
        }
    });
}

int
run_benchmark(Options& options, std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& games, const std::string& executable, Experiment* experiment)
{
    Runner runner(options, solvers, games, executable, experiment);
    return runner.run();
}

}
