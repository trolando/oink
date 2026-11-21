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

#include "experiment.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "common/terminal.hpp"
#include "json.hpp"
#include "options.hpp"
#include "report.hpp"
#include "summary.hpp"
#include "worker.hpp"

#ifndef _WIN32
#include <fcntl.h>
#include <poll.h>
#include <sys/file.h>
#include <unistd.h>
#endif

namespace fs = std::filesystem;

namespace bench {

/** A number with all its digits. */
static std::string
exact(double value)
{
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

std::vector<std::string>
settings_arguments(const Options& o)
{
    std::vector<std::string> args;
    args.push_back("--repeat=" + std::to_string(o.repeats));
    args.push_back("--warmup=" + std::to_string(o.warmups));
    args.push_back("--min-solve-time=" + exact(o.min_solve_time));
    args.push_back("--max-repeat=" + std::to_string(o.max_repeats));
    args.push_back("--skip-after-timeouts=" + std::to_string(o.skip_after_timeouts));
    if (o.isolate) {
        args.push_back("--timeout=" + exact(o.timeout));
        args.push_back("--safety-timeout=" + exact(o.safety_timeout));
        if (o.fresh_process) args.push_back("--fresh-process");
    } else {
        args.push_back("--in-process");
    }
    args.push_back("--bench-seed=" + std::to_string(o.seed));
    if (o.vary_solver_seed) args.push_back("--vary-solver-seed");
    if (!o.preprocess) {
        args.push_back("--no-preprocess");
    } else if (o.single && o.loops && o.wcwc && !o.inflate && !o.compress) {
        args.push_back("--preprocess");
    } else {
        // a custom recipe
        args.push_back("--preprocess");
        if (o.single) args.push_back("--single");
        if (o.loops) args.push_back("--loops");
        if (o.wcwc) args.push_back("--wcwc");
        if (o.inflate) args.push_back("--inflate");
        if (o.compress) args.push_back("--compress");
    }
    if (o.scc) args.push_back("--scc");
    args.push_back("--workers=" + std::to_string(o.workers));
    if (o.force_parallel) args.push_back("--force-parallel");
    if (o.trace >= 0) args.push_back("--trace=" + std::to_string(o.trace));
    if (o.reference) args.push_back("--reference=" + *o.reference);
    if (!o.label.empty()) args.push_back("--label=" + o.label);
    args.push_back("--count=" + std::to_string(o.count));
    args.push_back("--size=" + std::to_string(o.size));
    args.push_back("--maxp=" + std::to_string(o.maxp));
    args.push_back("--maxe=" + std::to_string(o.maxe));
    args.push_back("--seed=" + std::to_string(o.game_seed));
    return args;
}

std::string
encode_case(const Case& g)
{
    std::ostringstream out;
    out << std::setprecision(17);
    out << "game " << g.index << ' ' << std::quoted(g.name) << ' ' << std::quoted(g.path) << ' ';
    if (g.seed) out << *g.seed;
    else out << '-';
    out << ' ' << std::quoted(g.hash) << ' ' << g.vertices << ' ' << g.edges << ' ' << g.load << ' ' << g.parse;
    out << ' ' << std::quoted(g.error) << ' ' << g.input_failed << ' ' << g.prep_failed << ' ' << g.fully_preprocessed;
    out << ' ' << g.prep.residual_vertices << ' ' << g.prep.residual_edges << ' ' << g.prep.solved_vertices << ' ' << g.prep.seconds << ' ' << std::quoted(g.prep.residual_hash);
    out << ' ' << g.cells.size();
    for (const auto& cell : g.cells) {
        out << ' ' << cell.required << ' ' << cell.warmups << ' ' << (cell.failure ? (int)*cell.failure : -1) << ' ' << (int)cell.failure_phase;
        out << ' ' << cell.samples().size();
        for (double v : cell.samples()) out << ' ' << v;
    }
    out << " rss " << g.worker_peak_rss;
    for (const auto& cell : g.cells) out << ' ' << cell.max_trial_rss << ' ' << cell.max_solve_rss;
    out << " timeouts";
    for (const auto& cell : g.cells) out << ' ' << cell.solve_timeouts << ' ' << cell.skipped_timeouts;
    out << " end";
    return out.str();
}

bool
decode_case(const std::string& line, Case& g)
{
    std::istringstream in(line);
    std::string word;
    std::string seed;
    in >> word;
    if (word != "game") return false;
    in >> g.index >> std::quoted(g.name) >> std::quoted(g.path) >> seed >> std::quoted(g.hash) >> g.vertices >> g.edges >> g.load >> g.parse;
    in >> std::quoted(g.error) >> g.input_failed >> g.prep_failed >> g.fully_preprocessed;
    in >> g.prep.residual_vertices >> g.prep.residual_edges >> g.prep.solved_vertices >> g.prep.seconds >> std::quoted(g.prep.residual_hash);
    if (!in) return false;
    g.seed.reset();
    if (seed != "-") g.seed = (uint32_t)std::stoul(seed);
    size_t cells = 0;
    in >> cells;
    if (!in || cells > 100000) return false;
    g.cells.assign(cells, Cell());
    for (auto& cell : g.cells) {
        int failure = -1;
        int phase = 0;
        size_t count = 0;
        in >> cell.required >> cell.warmups >> failure >> phase >> count;
        if (!in || failure < -1 || failure > (int)Status::Interrupted || phase < 0 || phase > (int)Phase::Load || count > 100000000) return false;
        if (failure >= 0) cell.failure = (Status)failure;
        cell.failure_phase = (Phase)phase;
        for (size_t i = 0; i < count; i++) {
            double v;
            if (!(in >> v)) return false;
            cell.add_sample(v);
        }
    }
    in >> word;
    if (word == "rss") {
        in >> g.worker_peak_rss;
        for (auto& cell : g.cells) in >> cell.max_trial_rss >> cell.max_solve_rss;
        in >> word;
    }
    if (word == "timeouts") {
        for (auto& cell : g.cells) in >> cell.solve_timeouts >> cell.skipped_timeouts;
        in >> word;
    }
    return (bool)in && word == "end";
}

/** Replace <to> by <from>; atomically where the system can. */
static void
replace_file(const fs::path& from, const fs::path& to)
{
#ifdef _WIN32
    // rename may refuse to replace a file here
    std::error_code ignored;
    fs::remove(to, ignored);
#endif
    fs::rename(from, to);
}

/** Write <text> to <path> through a temporary file, so the file is never half written. */
static void
write_file(const fs::path& path, const std::string& text)
{
    const fs::path temporary = path.string() + ".tmp";
    {
        std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
        out << text;
        out.flush();
        if (!out) throw std::runtime_error("cannot write " + temporary.string());
    }
    replace_file(temporary, path);
}

static std::string
read_file(const fs::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path.string());
    std::ostringstream text;
    text << in.rdbuf();
    return text.str();
}

Experiment::Experiment(const std::string& directory) : root(fs::absolute(directory).lexically_normal().string())
{
}

Experiment::~Experiment()
{
#ifndef _WIN32
    if (lock >= 0) ::close(lock);
#endif
}

/**
 * Only one benchmark at a time works on an experiment. The lock is inherited
 * by the worker processes, so a worker that outlives a killed benchmark still
 * holds it.
 */
void
Experiment::acquire_lock()
{
#ifndef _WIN32
    lock = ::open((fs::path(root) / "lock").c_str(), O_CREAT | O_RDWR, 0644);
    if (lock < 0) throw std::runtime_error("cannot open the lock of the experiment");
    fcntl(lock, F_SETFD, FD_CLOEXEC);
    if (flock(lock, LOCK_EX | LOCK_NB) != 0) throw std::runtime_error("the experiment is busy: another benchmark (or a worker of one) still runs on it");
#endif
}

void
Experiment::start_session()
{
    const fs::path sessions = fs::path(root) / "sessions";
    fs::create_directories(sessions);
    session_ = 1;
    while (fs::exists(sessions / (std::to_string(session_) + ".jsonl"))) session_++;
}

void
Experiment::configure(Options& o) const
{
    o.experiment = root;
    o.jsonl = (fs::path(root) / "sessions" / (std::to_string(session_) + ".jsonl")).string();
    o.csv.clear();
    o.overwrite = false;
    o.artifacts = (fs::path(root) / "failures").string();
    // simple progress only: no live table
    o.live = false;
    o.review = false;
    o.heatmap = false;
}

std::unique_ptr<Experiment>
Experiment::create(const std::string& directory, const Options& o, const std::vector<SolverRun>& solvers, const std::vector<PlannedGame>& games, const std::string& executable)
{
    if (fs::exists(directory)) throw std::runtime_error("the experiment already exists: " + directory);
    std::unique_ptr<Experiment> experiment(new Experiment(directory));
    const fs::path root(experiment->root);
    fs::create_directories(root / "games");
    experiment->acquire_lock();

    const auto binary = file_fingerprint(executable);
    std::ostringstream plan;
    plan << "bench_solvers experiment 3\n";
    plan << "created " << std::quoted(utc_now()) << '\n';
    plan << "binary " << std::quoted(executable) << ' ' << (binary ? hex64(*binary) : std::string("unavailable")) << '\n';
    plan << "host " << std::quoted(host_name()) << '\n';
    plan << "directory " << std::quoted(fs::current_path().string()) << '\n';
    for (const auto& arg : settings_arguments(o)) plan << "setting " << std::quoted(arg) << '\n';
    for (const auto& spec : o.externals) plan << "external " << std::quoted(spec) << '\n';
    for (const auto& solver : solvers) plan << "solver " << std::quoted(solver.id) << ' ' << std::quoted(solver.config) << '\n';
    for (const auto& game : games) {
        plan << "game " << std::quoted(game.path) << ' ';
        if (game.seed) {
            plan << "- " << *game.seed;
        } else {
            const auto hash = file_fingerprint(game.path);
            if (!hash) throw std::runtime_error("cannot read " + game.path);
            plan << hex64(*hash) << " -";
        }
        plan << '\n';
    }
    plan << "end\n";
    write_file(root / "plan.txt", plan.str());
    write_file(root / "done.txt", "");
    experiment->start_session();
    return experiment;
}

std::unique_ptr<Experiment>
Experiment::resume(const std::string& directory, Options& options, std::vector<SolverRun>& solvers, std::vector<PlannedGame>& games, const std::string& executable)
{
    if (!fs::is_directory(directory) || !fs::exists(fs::path(directory) / "plan.txt")) throw std::runtime_error("not an experiment: " + directory);
    std::unique_ptr<Experiment> experiment(new Experiment(directory));
    const fs::path root(experiment->root);
    experiment->acquire_lock();

    std::istringstream plan(read_file(root / "plan.txt"));
    std::string line;
    std::getline(plan, line);
    if (line != "bench_solvers experiment 3") throw std::runtime_error("unknown experiment format: " + (root / "plan.txt").string());
    std::vector<std::string> args = {"bench_solvers"};
    std::vector<SolverRun> roster;
    std::vector<PlannedGame> planned;
    std::vector<std::string> hashes;
    std::string binary_hash;
    std::string host;
    std::string working_directory;
    bool complete = false;
    while (std::getline(plan, line)) {
        std::istringstream in(line);
        std::string key;
        in >> key;
        if (key == "end") {
            complete = true;
            break;
        }
        if (key == "binary") {
            std::string path;
            in >> std::quoted(path) >> binary_hash;
        } else if (key == "host") {
            in >> std::quoted(host);
        } else if (key == "directory") {
            in >> std::quoted(working_directory);
        } else if (key == "setting") {
            std::string arg;
            in >> std::quoted(arg);
            args.push_back(arg);
        } else if (key == "external") {
            std::string spec;
            in >> std::quoted(spec);
            args.push_back("--external=" + spec);
        } else if (key == "solver") {
            std::string id;
            std::string config;
            in >> std::quoted(id) >> std::quoted(config);
            roster.emplace_back(id, config);
        } else if (key == "game") {
            PlannedGame game;
            std::string hash;
            std::string seed;
            in >> std::quoted(game.path) >> hash >> seed;
            if (seed != "-") game.seed = (uint32_t)std::stoul(seed);
            planned.push_back(game);
            hashes.push_back(hash);
        }
        if (!in && key != "created") throw std::runtime_error("damaged experiment plan: " + line);
    }
    if (!complete) throw std::runtime_error("the experiment plan is incomplete: " + (root / "plan.txt").string());

    // the same program on the same host, or the measurements are not comparable
    const auto binary = file_fingerprint(executable);
    if (!binary || hex64(*binary) != binary_hash) throw std::runtime_error("the program differs from the one that started the experiment (rebuilt?); start a new experiment");
    if (host != host_name()) throw std::runtime_error("the experiment started on host " + host + ", not on " + host_name());

    // the saved settings; the options of this command line are ignored, except --plain
    const bool plain = options.plain;
    const auto command_line = options.argv;
    std::vector<const char*> argv;
    for (const auto& arg : args) argv.push_back(arg.c_str());
    std::vector<SolverRun> ignored;
    std::ostringstream help;
    Options saved;
    parse_command_line((int)argv.size(), argv.data(), saved, ignored, help);
    saved.plain = plain;
    saved.argv = command_line;
    saved.experiment = experiment->root;
    saved.resume = true;
    options = saved;
    solvers = roster;
    games = planned;
    // relative paths of external solvers mean what they meant at the start
    if (!working_directory.empty() && fs::is_directory(working_directory)) fs::current_path(working_directory);

    experiment->load_done(solvers.size());
    for (size_t i = 0; i < games.size(); i++) {
        if (games[i].seed || experiment->finished(i) != nullptr) continue;
        const auto hash = file_fingerprint(games[i].path);
        if (!hash || hex64(*hash) != hashes[i]) throw std::runtime_error("an input changed or disappeared since the experiment started: " + games[i].path);
    }
    experiment->start_session();
    return experiment;
}

void
Experiment::load_done(size_t solvers)
{
    const auto path = fs::path(root) / "done.txt";
    auto contents = read_file(path);
    if (!contents.empty() && contents.back() != '\n') {
        const auto last = contents.find_last_of('\n');
        contents.resize(last == std::string::npos ? 0 : last + 1);
        fs::resize_file(path, contents.size());
    }
    std::istringstream in(contents);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        Case game;
        // a line cut short by a crash during the commit is ignored: that game runs again
        if (!decode_case(line, game)) continue;
        if (game.cells.size() != solvers) throw std::runtime_error("done.txt does not match the plan of the experiment");
        done[game.index] = game;
    }
}

const Case*
Experiment::finished(size_t index) const
{
    const auto found = done.find(index);
    return found == done.end() ? nullptr : &found->second;
}

int
Experiment::restore_rounds(size_t index, size_t solvers, std::string& events,
                           const std::function<void(int, const RoundSample&)>& observe)
{
    const auto path = fs::path(root) / "games" / (std::to_string(index) + ".rounds");
    int round = 0;
    if (fs::exists(path)) {
        const auto size = fs::file_size(path);
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("cannot read " + path.string());
        std::streamoff end = 0;
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream header(line);
            size_t bytes = 0;
            size_t count = 0;
            int number = 0;
            const std::streamoff begin = in.tellg();
            if (!(header >> number >> bytes >> count) || number != round || count > solvers ||
                begin < 0 || bytes == 0 || bytes > size - (uintmax_t)begin) break;
            std::vector<RoundSample> samples;
            for (size_t i = 0; i < count; i++) {
                RoundSample sample;
                int status = 0;
                int phase = 0;
                if (!(header >> sample.solver >> sample.executed >> status >> phase >> sample.solve >> sample.peak_rss >> sample.solve_end_rss) ||
                    sample.solver >= solvers || status < 0 || status > (int)Status::Interrupted ||
                    phase < 0 || phase > (int)Phase::Load) break;
                sample.status = (Status)status;
                sample.phase = (Phase)phase;
                samples.push_back(sample);
            }
            if (samples.size() != count) break;
            std::string chunk(bytes, '\0');
            if (!in.read(chunk.data(), bytes) || chunk.back() != '\n') break;
            end = in.tellg();
            for (const auto& sample : samples) observe(round, sample);
            events += chunk;
            round++;
        }
        if (in.bad()) throw std::runtime_error("cannot read " + path.string());
        in.close();
        // A process killed during the last write may have left a partial round.
        if ((uintmax_t)end < size) fs::resize_file(path, end);
    }
    rounds.open(path, std::ios::binary | std::ios::app);
    if (!rounds) throw std::runtime_error("cannot open " + path.string());
    return round;
}

void
Experiment::record_round(int round, const std::vector<RoundSample>& samples, const std::string& events)
{
    rounds << std::setprecision(17) << round << ' ' << events.size() << ' ' << samples.size();
    for (const auto& s : samples) {
        rounds << ' ' << s.solver << ' ' << s.executed << ' ' << (int)s.status << ' ' << (int)s.phase
               << ' ' << s.solve << ' ' << s.peak_rss << ' ' << s.solve_end_rss;
    }
    rounds << '\n' << events << std::flush;
    if (!rounds) throw std::runtime_error("cannot write the completed rounds of the experiment");
}

void
Experiment::commit(const Case& game, const std::string& events)
{
    const fs::path root_path(root);
    write_file(root_path / "games" / (std::to_string(game.index) + ".jsonl"), events);
    {
        std::ofstream out(root_path / "done.txt", std::ios::app | std::ios::binary);
        out << encode_case(game) << '\n';
        out.flush();
        if (!out) throw std::runtime_error("cannot write " + (root_path / "done.txt").string());
    }
#ifndef _WIN32
    // the commit survives a crash of the machine
    const int fd = ::open((root_path / "done.txt").c_str(), O_RDONLY);
    if (fd >= 0) {
        fsync(fd);
        ::close(fd);
    }
#endif
    done[game.index] = game;
    if (rounds.is_open()) rounds.close();
    fs::remove(root_path / "games" / (std::to_string(game.index) + ".rounds"));
}

void
Experiment::publish(const Options& options, const std::vector<SolverRun>& solvers, const std::string& executable, const std::vector<PlannedGame>& games)
{
    const fs::path root_path(root);
    std::vector<Case> finished_games;
    for (const auto& entry : done) finished_games.push_back(entry.second);
    const bool complete = finished_games.size() == games.size();

    // the JSON Lines results: this session's metadata, every finished game, the summary
    Options aggregate = options;
    aggregate.jsonl = (root_path / "results.jsonl.tmp").string();
    aggregate.csv.clear();
    aggregate.overwrite = true;
    {
        Report report(aggregate, solvers);
        report.metadata(games.size(), executable, std::to_string(session_));
        for (const auto& game : finished_games) {
            std::istringstream events(read_file(root_path / "games" / (std::to_string(game.index) + ".jsonl")));
            std::string line;
            while (std::getline(events, line)) {
                if (!line.empty()) report.event(line);
            }
        }
        const auto totals = compute_totals(finished_games, games.size(), solvers.size(), options.repeats, reference_index(options.reference, solvers), options.timeout);
        report.summary(totals, !complete, finished_games.size());
    }
    replace_file(aggregate.jsonl, root_path / "results.jsonl");

    std::string csv = csv_header();
    for (const auto& game : finished_games) csv += csv_rows(game, solvers, options);
    write_file(root_path / "results.csv", csv);

    Json status;
    status.add("updated_utc", utc_now());
    status.add("session", session_);
    status.add("games_planned", games.size());
    status.add("games_finished", finished_games.size());
    status.add("complete", complete);
    status.add("results", "results.jsonl");
    status.add("statistics", "results.csv");
    write_file(root_path / "status.json", status.str() + "\n");
}

bool
confirm_resume(const std::string& directory, volatile std::sig_atomic_t* quit)
{
    if (!ui::is_tty(0) || !ui::is_tty(2)) throw std::runtime_error("the experiment already exists; use --resume=" + directory + " to continue it");
    ui::InterruptGuard guard(quit);
    std::cerr << "Experiment already exists: " << ui::clean(directory) << '\n';
    std::cerr << "Resume uses its saved games, solvers and settings; the other options of this command are ignored.\n";
    while (!*quit) {
        std::cerr << "Resume the experiment? [Y/n] " << std::flush;
#ifndef _WIN32
        // wait for a line, but let Ctrl-C through
        pollfd input{0, POLLIN, 0};
        while (!*quit && poll(&input, 1, 100) <= 0) {
        }
        if (*quit) {
            std::cerr << '\n';
            return false;
        }
#endif
        std::string answer;
        if (!std::getline(std::cin, answer)) return false;
        answer.erase(0, answer.find_first_not_of(" \t\r"));
        answer.erase(answer.find_last_not_of(" \t\r") + 1);
        for (char& c : answer) c = (char)std::tolower((unsigned char)c);
        if (answer.empty() || answer == "y" || answer == "yes") return true;
        if (answer == "n" || answer == "no" || answer == "q") return false;
    }
    return false;
}

}
