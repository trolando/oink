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

#include "bench.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <system_error>

#include "common/chooser.hpp"
#include "oink/game.hpp"

namespace bench {

std::vector<SolverRun>
default_solvers()
{
    std::vector<SolverRun> result;
    for (const char* id : {"tl", "tlj", "tlq", "npp", "fpi", "fpj", "zlk"}) result.emplace_back(id, "");
    return result;
}

std::optional<size_t>
reference_index(const std::optional<std::string>& reference, const std::vector<SolverRun>& solvers)
{
    if (!reference || *reference == "none") return std::nullopt;
    const std::string& name = *reference;
    if (name == "first") {
        if (solvers.empty()) return std::nullopt;
        return 0;
    }
    for (size_t i = 0; i < solvers.size(); i++) {
        if (solvers[i].name == name) return i;
    }
    // an id with exactly one selected variant
    std::optional<size_t> match;
    for (size_t i = 0; i < solvers.size(); i++) {
        if (solvers[i].id != name) continue;
        if (match) throw std::runtime_error("the reference matches several configurations; use the full solver[configuration] name: " + name);
        match = i;
    }
    if (!match) throw std::runtime_error("the reference solver is not selected: " + name);
    return match;
}

void
validate(const Options& o, const std::vector<SolverRun>& solvers)
{
    if (solvers.empty()) throw std::runtime_error("select at least one solver");
    for (size_t i = 0; i < solvers.size(); i++) {
        std::string error;
        if (!ui::valid_configuration(solvers[i].id, solvers[i].config, error)) {
            throw std::runtime_error("invalid configuration for " + solvers[i].id + ": " + error);
        }
        for (size_t j = 0; j < i; j++) {
            if (solvers[i].name == solvers[j].name) throw std::runtime_error("solver selected twice: " + solvers[i].name);
        }
    }
    reference_index(o.reference, solvers);
    if (o.repeats < 1) throw std::runtime_error("the number of repetitions must be positive");
    if (o.warmups < 0) throw std::runtime_error("the number of warm-up rounds cannot be negative");
    if (!std::isfinite(o.timeout) || o.timeout < 0) throw std::runtime_error("the timeout must be zero (no limit) or positive");
    if (!std::isfinite(o.safety_timeout) || (o.isolate ? o.safety_timeout <= 0 : o.safety_timeout != 0)) {
        throw std::runtime_error("isolated trials need a positive --safety-timeout; in-process trials have no limits");
    }
    if (!std::isfinite(o.min_solve_time) || o.min_solve_time < 0) throw std::runtime_error("the sampling target must be zero (off) or positive");
    if (o.max_repeats < 1 || (o.calibrated() && o.max_repeats < o.repeats)) throw std::runtime_error("the repetition cap must be at least the number of repetitions");
    const int rounds = o.calibrated() ? o.max_repeats : o.repeats;
    if (o.warmups > std::numeric_limits<int>::max() - rounds - o.pilots()) throw std::runtime_error("too many sampling rounds");
    if (o.timeout > 0 && !o.isolate) throw std::runtime_error("a timeout needs isolated trials; use --timeout=0 with --in-process");
    if (o.count < 1 || o.size < 1 || o.maxp < 0 || o.maxe < o.size) throw std::runtime_error("invalid random game parameters (count, size > 0, maxp >= 0, maxe >= size)");
    if (!o.preprocess && (o.inflate || o.compress)) throw std::runtime_error("--inflate and --compress are preprocessing steps; they need --preprocess");
    if (o.inflate && o.compress) throw std::runtime_error("choose --inflate or --compress");
    if (!o.jsonl.empty() && o.jsonl == o.csv) throw std::runtime_error("the JSON Lines and CSV reports need different paths");
}

const char*
status_name(Status status)
{
    switch (status) {
        case Status::Ok: return "ok";
        case Status::SolverError: return "solver_error";
        case Status::InvalidSolution: return "invalid_solution";
        case Status::SetupError: return "setup_error";
        case Status::Timeout: return "timeout";
        case Status::SafetyTimeout: return "safety_timeout";
        case Status::MemoryFailure: return "memory_failure";
        case Status::Crash: return "crash";
        case Status::Interrupted: return "interrupted";
    }
    return "unknown";
}

const char*
phase_name(Phase phase)
{
    switch (phase) {
        case Phase::None: return "";
        case Phase::Setup: return "setup";
        case Phase::Solve: return "solve";
        case Phase::Verify: return "verification";
    }
    return "";
}

std::string
failure_label(Status status, Phase phase)
{
    std::string label = status_name(status);
    if (phase != Phase::None) label += std::string(" (") + phase_name(phase) + ')';
    return label;
}

void
record_failure(Trial& trial, Phase phase, const std::exception& error)
{
    trial.failure_phase = phase;
    trial.status = phase == Phase::Solve ? Status::SolverError :
                   phase == Phase::Verify ? Status::InvalidSolution : Status::SetupError;
    const auto* system = dynamic_cast<const std::system_error*>(&error);
    if (dynamic_cast<const std::bad_alloc*>(&error) ||
        (system && system->code() == std::errc::not_enough_memory)) {
        trial.status = Status::MemoryFailure;
    }
    trial.message = error.what();
}

Stats
summarize(std::vector<double> values)
{
    Stats result;
    result.samples = values.size();
    if (values.empty()) return result;
    std::sort(values.begin(), values.end());
    auto quantile = [&](double p) {
        const double position = p * (double)(values.size() - 1);
        const size_t i = (size_t)position;
        const size_t j = std::min(i + 1, values.size() - 1);
        return values[i] + (values[j] - values[i]) * (position - (double)i);
    };
    result.median = quantile(0.5);
    result.q1 = quantile(0.25);
    result.q3 = quantile(0.75);
    for (double& v : values) v = std::abs(v - result.median);
    std::sort(values.begin(), values.end());
    result.mad = quantile(0.5);
    return result;
}

const Stats&
Cell::stats() const
{
    if (cache.samples != values.size()) cache = summarize(values);
    return cache;
}

std::string
Case::status(size_t s, int repeats) const
{
    if (fully_preprocessed) return "preprocessed";
    if (input_failed) return "input_error";
    if (prep_failed) return "preparation_error";
    const Cell& cell = cells[s];
    if (cell.failure) return status_name(*cell.failure);
    return cell.complete(repeats) ? "ok" : "incomplete";
}

std::string
duration(double seconds)
{
    std::ostringstream out;
    out << std::fixed;
    if (seconds < 0.001) {
        out << std::setprecision(seconds < 0.00001 ? 2 : 1) << seconds * 1e6 << " us";
    } else if (seconds < 1) {
        out << std::setprecision(seconds < 0.01 ? 2 : 1) << seconds * 1e3 << " ms";
    } else {
        out << std::setprecision(2) << seconds << " s";
    }
    return out.str();
}

std::string
seconds_setting(double seconds)
{
    if (seconds <= 0) return "off";
    std::ostringstream out;
    if (seconds < 1) {
        out << seconds * 1e3 << " ms";
    } else {
        out << seconds << " s";
    }
    return out.str();
}

std::vector<int>
heatmap_colors(const std::vector<double>& values)
{
    // green, yellow-green, yellow, orange, red
    static const int palette[] = {82, 118, 154, 190, 226, 220, 214, 208, 202, 196};
    std::vector<int> colors(values.size(), -1);
    std::vector<size_t> order;
    for (size_t i = 0; i < values.size(); i++) {
        if (std::isfinite(values[i]) && values[i] >= 0) order.push_back(i);
    }
    if (order.empty()) return colors;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return values[a] < values[b]; });
    if (values[order.front()] == values[order.back()]) {
        for (size_t i : order) colors[i] = 226;
        return colors;
    }
    for (size_t first = 0; first < order.size();) {
        // a run of equal values shares the rank of its middle
        size_t last = first + 1;
        while (last < order.size() && values[order[last]] == values[order[first]]) last++;
        double rank = ((double)first + (double)(last - 1)) / 2 / (double)(order.size() - 1);
        if (first == 0) rank = 0;
        else if (last == order.size()) rank = 1;
        const double scaled = std::max(0.0, std::min(1.0, (rank - 0.05) / 0.90));
        const int color = palette[(size_t)std::round(scaled * 9)];
        for (size_t i = first; i < last; i++) colors[order[i]] = color;
        first = last;
    }
    return colors;
}

int
calibrated_repetitions(double target, double pilot, int minimum, int maximum)
{
    const double count = std::ceil(target / std::max(pilot, 1e-9));
    return (int)std::max((double)minimum, std::min((double)maximum, count));
}

namespace {

const uint64_t FNV_OFFSET = 14695981039346656037ULL;
const uint64_t FNV_PRIME = 1099511628211ULL;

struct Fnv
{
    uint64_t h = FNV_OFFSET;

    void byte(unsigned char c)
    {
        h ^= c;
        h *= FNV_PRIME;
    }

    void word(uint64_t x)
    {
        for (int i = 0; i < 8; i++) {
            byte((unsigned char)(x & 255));
            x >>= 8;
        }
    }
};

}

uint64_t
fingerprint(const pg::Game& game)
{
    Fnv fnv;
    fnv.word((uint64_t)game.nodecount());
    fnv.word((uint64_t)game.edgecount());
    for (int v = 0; v < game.nodecount(); v++) {
        fnv.word((uint64_t)(int64_t)game.priority(v));
        fnv.word((uint64_t)game.owner(v));
        for (const int* e = game.outs(v); *e != -1; e++) fnv.word((uint64_t)(int64_t)*e);
        fnv.word(UINT64_MAX);
    }
    return fnv.h;
}

std::optional<uint64_t>
file_fingerprint(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) return std::nullopt;
    Fnv fnv;
    char data[65536];
    while (in) {
        in.read(data, sizeof(data));
        for (std::streamsize i = 0; i < in.gcount(); i++) fnv.byte((unsigned char)data[i]);
    }
    if (in.bad()) return std::nullopt;
    return fnv.h;
}

std::string
hex64(uint64_t value)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

uint64_t
mix(uint64_t x)
{
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

uint32_t
series_seed(uint64_t schedule_seed)
{
    return (uint32_t)mix(schedule_seed ^ 0xa0761d6478bd642fULL);
}

uint32_t
trial_seed(uint64_t schedule_seed, uint64_t game, int round)
{
    return (uint32_t)mix(schedule_seed ^ game ^ ((uint64_t)round + 0x9e3779b97f4a7c15ULL));
}

uint64_t
round_seed(uint64_t schedule_seed, uint64_t game, int round)
{
    return mix(schedule_seed ^ game ^ ((uint64_t)round + 0xd1b54a32d192ed03ULL));
}

}
