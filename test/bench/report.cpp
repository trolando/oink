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

#include "report.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "json.hpp"
#include "summary.hpp"
#include "oink/version.hpp"

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <sys/utsname.h>
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <mach-o/dyld.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#ifdef __FreeBSD__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

#ifndef OINK_BENCH_BUILD_TYPE
#define OINK_BENCH_BUILD_TYPE ""
#endif

namespace fs = std::filesystem;

namespace bench {

std::string
utc_now()
{
    const std::time_t now = std::time(nullptr);
    std::tm value{};
#ifdef _WIN32
    value = *std::gmtime(&now); // the benchmark is single-threaded here
#else
    gmtime_r(&now, &value);
#endif
    std::ostringstream out;
    out << std::put_time(&value, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

std::string
build_json()
{
    Json json;
    json.add("oink_version", pg::version());
#if defined(__VERSION__)
    json.add("compiler", __VERSION__);
#elif defined(_MSC_FULL_VER)
    json.add("compiler", "MSVC " + std::to_string(_MSC_FULL_VER));
#else
    json.add("compiler", "unknown");
#endif
    json.add("cpp", (long)__cplusplus);
    json.add("build_type", OINK_BENCH_BUILD_TYPE);
#ifdef NDEBUG
    json.add("assertions", false);
#else
    json.add("assertions", true);
#endif
    return json.str();
}

std::string
host_name()
{
#ifdef _WIN32
    const char* name = std::getenv("COMPUTERNAME");
    return name != nullptr ? name : "";
#else
    char name[256] = {0};
    if (gethostname(name, sizeof(name) - 1) != 0) return "";
    return name;
#endif
}

/** The model name of the CPU, where the system tells it. */
static std::string
cpu_model()
{
#if defined(__linux__)
    std::ifstream cpuinfo("/proc/cpuinfo");
    std::string line;
    while (std::getline(cpuinfo, line)) {
        if (line.rfind("model name", 0) != 0) continue;
        const auto colon = line.find(':');
        return colon == std::string::npos ? line : line.substr(colon + 2);
    }
    return "";
#elif defined(__APPLE__) || defined(__FreeBSD__)
    char model[256] = {0};
    size_t size = sizeof(model) - 1;
#ifdef __APPLE__
    const char* key = "machdep.cpu.brand_string";
#else
    const char* key = "hw.model";
#endif
    if (sysctlbyname(key, model, &size, nullptr, 0) != 0) return "";
    return model;
#else
    return "";
#endif
}

std::string
machine_json()
{
    Json json;
    json.add("host", host_name());
#ifndef _WIN32
    utsname host{};
    if (uname(&host) == 0) {
        json.add("os", std::string(host.sysname) + " " + host.release);
        json.add("architecture", host.machine);
    }
#else
    json.add("os", "Windows");
#endif
    json.add("cpu", cpu_model());
    json.add("hardware_threads", (long)std::thread::hardware_concurrency());
#ifdef __linux__
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0) {
        std::vector<std::string> cpus;
        for (int i = 0; i < CPU_SETSIZE; i++) {
            if (CPU_ISSET(i, &mask)) cpus.push_back(std::to_string(i));
        }
        json.raw("affinity", json_array(cpus));
    }
#endif
    Json environment;
    for (const char* key : {"LACE_PIN", "LACE_WORKERS", "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES"}) {
        const char* value = std::getenv(key);
        if (value != nullptr) environment.add(key, value);
    }
    json.raw("runtime_environment", environment.str());
    return json.str();
}

std::string
executable_path(const char* argv0)
{
    std::error_code error;
#if defined(__linux__)
    const auto path = fs::read_symlink("/proc/self/exe", error);
    if (!error) return path.string();
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> path(size + 1, 0);
    if (_NSGetExecutablePath(path.data(), &size) == 0) {
        const auto canonical = fs::canonical(path.data(), error);
        if (!error) return canonical.string();
    }
#elif defined(__FreeBSD__)
    int mib[4] = {CTL_KERN, KERN_PROC, KERN_PROC_PATHNAME, -1};
    char path[4096] = {0};
    size_t size = sizeof(path) - 1;
    if (sysctl(mib, 4, path, &size, nullptr, 0) == 0) return path;
#elif defined(_WIN32)
    char path[4096] = {0};
    const DWORD n = GetModuleFileNameA(nullptr, path, sizeof(path) - 1);
    if (n > 0 && n < sizeof(path) - 1) return path;
#endif
    if (argv0 == nullptr) return "";
    const auto canonical = fs::canonical(argv0, error);
    return error ? std::string(argv0) : canonical.string();
}

/** Empty steady_clock intervals: a diagnostic of the clock, never subtracted. */
static std::string
clock_probe_json()
{
    using Clock = std::chrono::steady_clock;
    std::vector<double> observations;
    observations.reserve(10000);
    for (int i = 0; i < 10000; i++) {
        const auto a = Clock::now();
        const auto b = Clock::now();
        observations.push_back(std::chrono::duration<double>(b - a).count());
    }
    const Stats stats = summarize(observations);
    Json json;
    json.add("kind", "empty steady_clock interval; diagnostic only, never subtracted");
    json.add("samples", (long)observations.size());
    json.add("min_seconds", *std::min_element(observations.begin(), observations.end()));
    json.add("median_seconds", stats.median);
    json.add("q1_seconds", stats.q1);
    json.add("q3_seconds", stats.q3);
    return json.str();
}

static std::string
csv_string(const std::string& text)
{
    std::string out = "\"";
    for (char c : text) {
        if (c == '"') out += '"';
        out += c;
    }
    return out + '"';
}

/** Open <path> for a report, refusing to replace a file unless <overwrite>. */
static void
open_report(std::ofstream& out, const std::string& path, bool overwrite)
{
    if (fs::exists(path) && !overwrite) throw std::runtime_error("the report already exists (use --overwrite): " + path);
    const auto parent = fs::path(path).parent_path();
    if (!parent.empty()) fs::create_directories(parent);
    out.open(path, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
}

Report::Report(const Options& options, const std::vector<SolverRun>& solvers) : options(options), solvers(solvers)
{
    jsonl_path_ = options.jsonl;
    csv_path_ = options.csv;
    if (!jsonl_path_.empty()) open_report(jsonl, jsonl_path_, options.overwrite);
    if (!csv_path_.empty()) {
        open_report(csv, csv_path_, options.overwrite);
        csv << csv_header();
        csv.flush();
    }
}

void
Report::event(const std::string& json)
{
    if (jsonl.is_open()) {
        jsonl << json << '\n';
        jsonl.flush();
        if (!jsonl) throw std::runtime_error("cannot write the JSON Lines report " + jsonl_path_);
    }
    if (copy != nullptr) {
        *copy << json << '\n';
        copy->flush();
    }
}

void
Report::metadata(size_t planned, const std::string& executable, const std::string& session)
{
    const Options& o = options;
    Json json;
    json.add("event", "metadata");
    json.add("schema_version", 7);
    json.add("tool", "bench_solvers");
    json.add("started_utc", utc_now());
    json.add("oink_version", pg::version());
    json.add("label", o.label);
    json.add("experiment", o.experiment);
    if (!session.empty()) json.raw("experiment_session", session);
    const auto ref = reference_index(o.reference, solvers);
    if (ref) json.add("reference", solvers[*ref].name);
    else json.null("reference");
    json.raw("build", build_json());
    json.raw("clock_probe", clock_probe_json());
    json.add("binary", executable);
    const auto binary = file_fingerprint(executable);
    json.add("binary_fnv1a64", binary ? hex64(*binary) : std::string("unavailable"));
    json.add("metric", "wall seconds of Oink::run on a fresh copy of the residual game; SCC solving and lazy in-edges included; reading, preparation, copying, Oink construction with in-edges and verification excluded");
    json.add("clock", "steady_clock");
    json.add("clock_period_seconds", (double)std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den);
    Json sampling;
    sampling.add("target_solve_seconds", o.min_solve_time);
    sampling.add("max_repetitions", o.adaptive_sampling() ? o.max_repeats : o.repeats);
    json.raw("sampling", sampling.str());
    json.add("repetitions", o.repeats);
    json.add("warmup_rounds", o.warmups);
    json.add("order_seed", o.seed);
    json.add("solver_rng_policy", o.vary_solver_seed ? "vary by round" : "fixed per game");
    json.add("cpu_metric", "process CPU time of the trial process during the solve; excludes external subprocesses");
    json.add("random_series_seed", o.game_seed);
    json.add("games_planned", planned);
    if (o.timeout > 0) json.add("solve_wall_limit_seconds", o.timeout);
    else json.null("solve_wall_limit_seconds");
    if (o.isolate) json.add("safety_wall_limit_seconds", o.safety_timeout);
    else json.null("safety_wall_limit_seconds");
    json.add("isolation", o.isolate ? "worker process per game, restarted after a failure" : "in_process");
    json.add("solver_logging", o.trace >= 0 ? "captured at trace level " + std::to_string(o.trace) : std::string("disabled"));
    json.add("workers_requested", o.workers);
    json.add("force_parallel", o.force_parallel);
    json.raw("machine", machine_json());
    Json prep;
    prep.add("enabled", o.preprocess);
    prep.add("once_per_game", true);
    prep.add("timed", false);
    prep.add("renumber", o.preprocess && !o.inflate && !o.compress);
    prep.add("inflate", o.preprocess && o.inflate);
    prep.add("compress", o.preprocess && o.compress);
    prep.add("single", o.preprocess && o.single);
    prep.add("loops", o.preprocess && o.loops);
    prep.add("wcwc", o.preprocess && o.wcwc);
    json.raw("preprocessing", prep.str());
    json.add("bottom_scc_solving", o.scc);
    Json random;
    random.add("count", o.count);
    random.add("size", o.size);
    random.add("maxp", o.maxp);
    random.add("maxe", o.maxe);
    json.raw("random_parameters", random.str());
    std::vector<std::string> argv;
    for (const auto& arg : o.argv) argv.push_back(json_string(arg));
    json.raw("argv", json_array(argv));
    std::vector<std::string> list;
    for (const auto& solver : solvers) {
        Json entry;
        entry.add("id", solver.id);
        entry.add("configuration", solver.config);
        entry.add("name", solver.name);
        list.push_back(entry.str());
    }
    json.raw("solvers", json_array(list));
    event(json.str());
}

void
Report::game(const Case& game)
{
    Json json;
    json.add("event", "game");
    json.add("index", game.index);
    json.add("name", game.name);
    json.add("path", game.path);
    json.add("graph_fnv1a64", game.hash);
    json.add("vertices", game.vertices);
    json.add("edges", game.edges);
    json.add("seed", game.seed);
    json.add("load_seconds", game.load);
    json.add("parse_sort_seconds", game.parse);
    json.add("error", game.error);
    event(json.str());
}

void
Report::sample(const Case& game, size_t solver, size_t position, int round, const char* phase, uint32_t seed, int generation, const Trial& trial, bool executed, const std::string& artifact)
{
    if (!jsonl.is_open() && copy == nullptr) return;
    Json json;
    json.add("event", "sample");
    json.add("game_index", game.index);
    json.add("solver", solvers[solver].name);
    json.add("phase", phase);
    json.add("schedule_position", position);
    json.add("round", round);
    json.add("warmup", std::string(phase) != "measured");
    json.add("game_rng_seed", seed);
    json.add("worker_generation", generation);
    json.add("status", executed ? status_name(trial.status) : "skipped");
    if (trial.failure_phase != Phase::None) json.add("failure_phase", phase_name(trial.failure_phase));
    else json.null("failure_phase");
    json.add("executed", executed);
    auto timing = [&](const char* name, double seconds) {
        if (executed && seconds >= 0) json.add(name, seconds);
        else json.null(name);
    };
    timing("setup_seconds", trial.setup);
    timing("solve_seconds", trial.solve);
    timing("verification_seconds", trial.verify);
    timing("cpu_seconds", trial.cpu);
    json.add("parallel", trial.parallel);
    json.add("effective_solver", trial.effective_solver);
    json.add("effective_configuration", trial.effective_config);
    json.add("actual_workers", trial.workers);
    json.add("message", trial.message);
    json.add("artifact", artifact);
    event(json.str());
}

std::string
csv_header()
{
    return "game,graph_fnv1a64,vertices,edges,game_seed,solver,configuration,status,samples,median_seconds,q1_seconds,q3_seconds,mad_seconds,target_samples,residual_vertices,residual_edges,residual_fnv1a64,preprocessing_seconds,mean_seconds,accumulated_solve_seconds,sampling_capped,failure_phase\n";
}

std::string
csv_rows(const Case& game, const std::vector<SolverRun>& solvers, const Options& options)
{
    const int repeats = options.repeats;
    std::ostringstream out;
    out << std::setprecision(17);
    for (size_t s = 0; s < solvers.size(); s++) {
        const Cell& cell = game.cells[s];
        const Stats& st = cell.stats();
        const bool complete = game.timed() && cell.complete(repeats);
        out << csv_string(game.name) << ',' << game.hash << ',' << game.vertices << ',' << game.edges << ',';
        if (game.seed) out << *game.seed;
        out << ',' << csv_string(solvers[s].id) << ',' << csv_string(solvers[s].config) << ',' << game.status(s, repeats) << ',' << st.samples << ',';
        if (complete) out << st.median << ',' << st.q1 << ',' << st.q3 << ',' << st.mad;
        else out << ",,,";
        out << ',' << (game.timed() ? cell.target(repeats) : 0) << ',';
        if (!game.input_failed && !game.prep_failed) {
            out << game.prep.residual_vertices << ',' << game.prep.residual_edges << ',' << game.prep.residual_hash << ',' << game.prep.seconds;
        } else {
            out << ",,,";
        }
        out << ',';
        if (complete) out << st.mean;
        out << ',' << cell.accumulated_solve() << ',' << (complete && cell.capped(options) ? "true" : "false");
        out << ',' << phase_name(cell.failure_phase) << '\n';
    }
    return out.str();
}

std::vector<std::string>
result_events(const Case& game, const std::vector<SolverRun>& solvers, const Options& options)
{
    const int repeats = options.repeats;
    std::vector<std::string> events;
    for (size_t s = 0; s < solvers.size(); s++) {
        const Cell& cell = game.cells[s];
        const Stats& st = cell.stats();
        const bool complete = game.timed() && cell.complete(repeats);
        Json json;
        json.add("event", "result");
        json.add("game_index", game.index);
        json.add("solver", solvers[s].name);
        json.add("status", game.status(s, repeats));
        if (cell.failure_phase != Phase::None) json.add("failure_phase", phase_name(cell.failure_phase));
        else json.null("failure_phase");
        json.add("samples", st.samples);
        json.add("target_samples", game.timed() ? cell.target(repeats) : 0);
        json.add("accumulated_solve_seconds", cell.accumulated_solve());
        json.add("sampling_capped", complete && cell.capped(options));
        if (complete) {
            json.add("mean_seconds", st.mean);
            json.add("median_seconds", st.median);
            json.add("q1_seconds", st.q1);
            json.add("q3_seconds", st.q3);
            json.add("mad_seconds", st.mad);
        } else {
            json.null("mean_seconds");
            json.null("median_seconds");
            json.null("q1_seconds");
            json.null("q3_seconds");
            json.null("mad_seconds");
        }
        events.push_back(json.str());
    }
    return events;
}

void
Report::results(const Case& game)
{
    for (const auto& line : result_events(game, solvers, options)) event(line);
    if (csv.is_open()) {
        csv << csv_rows(game, solvers, options);
        csv.flush();
        if (!csv) throw std::runtime_error("cannot write the CSV report " + csv_path_);
    }
}

void
Report::summary(const Totals& totals, bool interrupted, size_t started)
{
    for (size_t s = 0; s < solvers.size(); s++) {
        const SolverTotals& t = totals.solvers[s];
        const bool ref = totals.reference.has_value();
        Json json;
        json.add("event", "summary");
        json.add("solver", solvers[s].name);
        json.add("games_planned", totals.planned);
        json.add("fully_preprocessed_games", totals.preprocessed);
        json.add("games_timed_planned", totals.timed_planned);
        json.add("complete", t.complete);
        json.add("failed", t.failed);
        if (t.complete) json.add("sum_game_medians_seconds", t.sum);
        else json.null("sum_game_medians_seconds");
        if (ref) json.add("reference", solvers[*totals.reference].name);
        else json.null("reference");
        json.add("paired_games", ref ? std::optional<size_t>(t.paired) : std::nullopt);
        json.add("reference_only_complete", ref ? std::optional<size_t>(t.reference_only) : std::nullopt);
        json.add("solver_only_complete", ref ? std::optional<size_t>(t.solver_only) : std::nullopt);
        json.add("geomean_speedup", t.geomean());
        json.add("wins", ref ? std::optional<size_t>(t.wins) : std::nullopt);
        json.add("ties", ref ? std::optional<size_t>(t.ties) : std::nullopt);
        json.add("losses", ref ? std::optional<size_t>(t.losses) : std::nullopt);
        json.add("geomean_vbs_ratio", t.vbs_ratio());
        json.add("vbs_ratio_games", t.vbs_games);
        json.add("vbs_wins", t.vbs_wins);
        json.add("vbs_ties", t.vbs_ties);
        json.add("point_estimate_tie_band", TIE_BAND);
        json.add("interrupted", interrupted);
        event(json.str());
    }
    Json end;
    end.add("event", "end");
    end.add("ended_utc", utc_now());
    end.add("interrupted", interrupted);
    end.add("games_started", started);
    event(end.str());
}

}
