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

#include "setup.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace bench {

// the largest repetition or warm-up count
static const uint64_t INT_LIMIT = (uint64_t)std::numeric_limits<int>::max();

/** <text> as an integer from <minimum> to <maximum>; false with a message otherwise. */
static bool
parse_integer(const std::string& text, uint64_t minimum, uint64_t maximum, uint64_t& value, std::string& error)
{
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos) {
        error = "enter a nonnegative integer";
        return false;
    }
    try {
        value = std::stoull(text);
    } catch (const std::exception&) {
        error = "the number is too large";
        return false;
    }
    if (value < minimum || value > maximum) {
        error = "enter a number from " + std::to_string(minimum) + " to " + std::to_string(maximum);
        return false;
    }
    return true;
}

/** <text> as nonnegative seconds; false with a message otherwise. */
static bool
parse_seconds(const std::string& text, double& value, std::string& error)
{
    try {
        size_t used = 0;
        value = std::stod(text, &used);
        if (used == text.size() && std::isfinite(value) && value >= 0) return true;
    } catch (const std::exception&) {
    }
    error = "enter nonnegative seconds, such as 0.5";
    return false;
}

static std::string
number(double value)
{
    std::ostringstream out;
    out << value;
    return out.str();
}

ui::ChooserConfig
setup_config(Options& o)
{
    ui::ChooserConfig config;
    config.title = "BENCHMARK SETUP";
    config.reference = true;

    ui::Setting repeats;
    repeats.key = 'r';
    repeats.label = "Reps";
    repeats.title = "MEASURED REPETITIONS";
    repeats.hint = "Measured repetitions per game and solver; with a sampling target, the minimum.";
    repeats.show = [&o] { return std::to_string(o.repeats); };
    repeats.text = repeats.show;
    repeats.apply = [&o](const std::string& text, std::string& error) {
        uint64_t value = 0;
        if (!parse_integer(text, 1, INT_LIMIT, value, error)) return false;
        o.repeats = (int)value;
        if (o.max_repeats < o.repeats) o.max_repeats = o.repeats;
        return true;
    };
    config.settings.push_back(repeats);

    ui::Setting warmups;
    warmups.key = 'w';
    warmups.label = "Warm-ups";
    warmups.title = "WARM-UP ROUNDS";
    warmups.hint = "Warm-up rounds per game; their trials are verified but not part of the results.";
    warmups.show = [&o] { return std::to_string(o.warmups); };
    warmups.text = warmups.show;
    warmups.apply = [&o](const std::string& text, std::string& error) {
        uint64_t value = 0;
        if (!parse_integer(text, 0, INT_LIMIT, value, error)) return false;
        o.warmups = (int)value;
        return true;
    };
    config.settings.push_back(warmups);

    ui::Setting timeout;
    timeout.key = 't';
    timeout.label = "Timeout";
    timeout.title = "SOLVE TIMEOUT";
    timeout.hint = "Wall seconds per Oink::run; 0 disables the solve limit.";
    timeout.show = [&o] { return seconds_setting(o.timeout); };
    timeout.text = [&o] { return number(o.timeout); };
    timeout.apply = [&o](const std::string& text, std::string& error) {
        double value = 0;
        if (!parse_seconds(text, value, error)) return false;
        if (value > 0 && !o.isolate) {
            error = "a timeout needs worker processes; this run has none (--in-process, or no POSIX)";
            return false;
        }
        o.timeout = value;
        return true;
    };
    config.settings.push_back(timeout);

    ui::Setting seed;
    seed.key = 'b';
    seed.label = "Seed";
    seed.title = "SCHEDULE SEED";
    seed.hint = "The seed of the file order, the solver order of every round and the solver seeds.";
    seed.show = [&o] { return std::to_string(o.seed); };
    seed.text = seed.show;
    seed.apply = [&o](const std::string& text, std::string& error) {
        uint64_t value = 0;
        if (!parse_integer(text, 0, std::numeric_limits<uint64_t>::max(), value, error)) return false;
        o.seed = value;
        if (!o.game_seed_given) o.game_seed = series_seed(o.seed);
        return true;
    };
    config.settings.push_back(seed);

    ui::Setting sampling;
    sampling.key = 'p';
    sampling.label = "Sampling";
    sampling.title = "SOLVE-TIME SAMPLING TARGET";
    sampling.hint = "Accumulate this many measured solve seconds per game and solver, up to --max-repeat; 0 uses a fixed count.";
    sampling.show = [&o] { return seconds_setting(o.min_solve_time); };
    sampling.text = [&o] { return number(o.min_solve_time); };
    sampling.apply = [&o](const std::string& text, std::string& error) {
        double value = 0;
        if (!parse_seconds(text, value, error)) return false;
        o.min_solve_time = value;
        return true;
    };
    config.settings.push_back(sampling);

    ui::Setting preprocess;
    preprocess.key = 'P';
    preprocess.label = "Preprocessing";
    preprocess.show = [&o] { return std::string(o.preprocess ? "ON" : "OFF"); };
    preprocess.toggle = [&o] {
        o.preprocess = !o.preprocess;
        return std::string(o.preprocess ? "Preprocessing ON: once per game, outside the solve time." : "Preprocessing OFF.");
    };
    config.settings.push_back(preprocess);
    return config;
}

bool
setup_benchmark(Options& o, std::vector<SolverRun>& solvers, volatile std::sig_atomic_t* quit)
{
    ui::Chooser chooser(solvers, setup_config(o));
    // start with the reference of the command line, if it is selected
    bool resolved = true;
    try {
        const auto index = reference_index(o.reference, solvers);
        if (index) chooser.set_reference(solvers[*index].name);
    } catch (const std::exception&) {
        resolved = false;
    }
    if (!ui::choose(chooser, quit)) return false;
    solvers = chooser.selected();
    // an unselected reference of the command line stays, to be checked at the start
    if (resolved || chooser.reference()) o.reference = chooser.reference();
    return true;
}

}
