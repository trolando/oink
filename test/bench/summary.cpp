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

#include "summary.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

#include "common/terminal.hpp"

namespace bench {

std::optional<double>
SolverTotals::geomean() const
{
    if (paired == 0) return std::nullopt;
    return std::exp(log_ratios / (double)paired);
}

Totals
compute_totals(const std::vector<Case>& games, size_t planned, size_t solver_count, int repeats, std::optional<size_t> reference)
{
    Totals totals;
    totals.planned = planned;
    totals.started = games.size();
    totals.reference = reference;
    totals.solvers.resize(solver_count);
    for (const auto& game : games) {
        if (game.fully_preprocessed) totals.preprocessed++;
        if (game.prep_failed) totals.prep_failures++;
    }
    totals.timed_planned = planned - totals.preprocessed;
    totals.complete_coverage = games.size() == planned;

    for (const auto& game : games) {
        if (!game.timed()) continue;
        // the virtual best solver: the fastest median, with a 1% band for ties
        double best = std::numeric_limits<double>::infinity();
        for (const auto& cell : game.cells) {
            if (cell.complete(repeats)) best = std::min(best, cell.stats().median);
        }
        std::vector<size_t> leaders;
        for (size_t s = 0; s < solver_count; s++) {
            const Cell& cell = game.cells[s];
            if (cell.complete(repeats) && cell.stats().median <= best * (1 + TIE_BAND)) leaders.push_back(s);
        }
        for (size_t s : leaders) {
            if (leaders.size() == 1) totals.solvers[s].vbs_wins++;
            else totals.solvers[s].vbs_ties++;
        }

        const bool ref_ok = reference && game.cells[*reference].complete(repeats);
        for (size_t s = 0; s < solver_count; s++) {
            SolverTotals& t = totals.solvers[s];
            const Cell& cell = game.cells[s];
            const bool ok = cell.complete(repeats);
            if (!ok || cell.stats().median <= 0) totals.complete_coverage = false;
            if (ok) {
                t.complete++;
                t.sum += cell.stats().median;
            }
            if (cell.failure && *cell.failure != Status::Interrupted) t.failed++;
            if (!reference) continue;
            if (ref_ok && !ok) t.reference_only++;
            if (ok && !ref_ok) t.solver_only++;
            if (ok && ref_ok) {
                const double a = game.cells[*reference].stats().median;
                const double b = cell.stats().median;
                if (a > 0 && b > 0) {
                    const double ratio = a / b;
                    t.log_ratios += std::log(ratio);
                    t.paired++;
                    if (ratio > 1 + TIE_BAND) t.wins++;
                    else if (ratio < 1 / (1 + TIE_BAND)) t.losses++;
                    else t.ties++;
                }
            }
        }
    }

    totals.order.resize(solver_count);
    std::iota(totals.order.begin(), totals.order.end(), 0);
    std::stable_sort(totals.order.begin(), totals.order.end(), [&](size_t a, size_t b) {
        const bool full_a = totals.solvers[a].complete == totals.timed_planned;
        const bool full_b = totals.solvers[b].complete == totals.timed_planned;
        if (full_a != full_b) return full_a;
        return totals.solvers[a].sum < totals.solvers[b].sum;
    });
    return totals;
}

void
print_summary(std::ostream& out, const Totals& totals, const std::vector<SolverRun>& solvers, const RunDetails& details, bool styled)
{
    const char* heading = styled ? ui::style::title : "";
    const char* bold = styled ? ui::style::bold : "";
    const char* dim = styled ? ui::style::dim : "";
    const char* reset = styled ? ui::style::reset : "";
    const bool reference = totals.reference.has_value();

    size_t name_width = 22;
    for (const auto& solver : solvers) name_width = std::max(name_width, solver.name.size() + 2);
    const int digits = (int)std::to_string(totals.planned).size();
    const int coverage_width = std::max(10, 3 * digits + 4);
    const int wtl_width = std::max(9, 3 * digits + 2);
    size_t rule_width = name_width + (size_t)coverage_width + 14 + 2 + (size_t)std::max(7, 2 * digits + 1);
    if (reference) rule_width += 14 + (totals.complete_coverage ? 0 : 10) + 2 + (size_t)wtl_width;
    const std::string rule(rule_width, '-');

    out << '\n' << heading << "Benchmark summary - sum of per-game medians (seconds)" << reset << "\n\n";
    out << bold << std::left << std::setw((int)name_width) << "solver" << std::right << std::setw(coverage_width) << "complete" << std::setw(14) << "sum medians";
    if (reference) out << std::setw(14) << "geom. speedup";
    if (reference && !totals.complete_coverage) out << std::setw(10) << "pairs";
    if (reference) out << "  " << std::left << std::setw(wtl_width) << "W/T/L";
    out << "  VBS W/T" << reset << '\n';
    out << dim << rule << reset << '\n';

    // colour the sums of the solvers with full coverage by their rank
    std::vector<double> ranked(solvers.size(), std::numeric_limits<double>::quiet_NaN());
    for (size_t s = 0; s < solvers.size(); s++) {
        const SolverTotals& t = totals.solvers[s];
        if (totals.timed_planned > 0 && t.complete == totals.timed_planned) ranked[s] = t.sum;
    }
    const auto colors = styled ? heatmap_colors(ranked) : std::vector<int>(solvers.size(), -1);

    for (size_t s : totals.order) {
        const SolverTotals& t = totals.solvers[s];
        std::string coverage = std::to_string(t.complete) + "/" + std::to_string(totals.timed_planned);
        if (t.failed) coverage += " !" + std::to_string(t.failed);
        std::ostringstream sum;
        sum << std::right << std::setw(14);
        if (t.complete) sum << std::fixed << std::setprecision(6) << t.sum;
        else sum << "-";
        out << std::left << std::setw((int)name_width) << ui::clean(solvers[s].name) << std::right << std::setw(coverage_width) << coverage;
        out << ui::style::color256(colors[s]) << sum.str() << (colors[s] >= 0 ? reset : "");
        if (reference) {
            std::ostringstream speedup;
            if (t.geomean()) {
                speedup << std::fixed << std::setprecision(3) << *t.geomean() << "x";
            } else {
                speedup << "-";
            }
            out << std::setw(14) << speedup.str();
            if (!totals.complete_coverage) out << std::setw(10) << t.paired;
            const std::string pairwise = s == *totals.reference ? "reference" : std::to_string(t.wins) + "/" + std::to_string(t.ties) + "/" + std::to_string(t.losses);
            out << "  " << std::left << std::setw(wtl_width) << pairwise << std::right;
        }
        out << "  " << t.vbs_wins << '/' << t.vbs_ties << '\n';
    }
    out << dim << rule << reset << "\n\n";

    out << bold << "Column guide" << reset << '\n' << dim;
    out << "  Complete       Verified games / games needing a solver; !N counts failures.\n";
    out << "  Sum medians    Sum of the per-game medians of the solve time, in seconds.\n";
    if (reference) {
        out << "  Geom. speedup  Geometric mean of the paired ratios against the reference (>1x is faster).\n";
        out << "  W/T/L          Wins/ties/losses against the reference, with a 1% tie band.\n";
    }
    out << "  VBS W/T        Uniquely fastest / tied within 1% of the best median of the game.\n";
    if (!totals.complete_coverage) {
        if (reference) out << "  Pairs          Games verified for both the solver and the reference.\n";
        out << "  Coverage       Solvers with incomplete coverage follow the others.\n";
    }
    out << reset;

    out << '\n' << bold << "Run details" << reset << '\n';
    out << "  Schedule seed: " << details.schedule_seed << " (replay with --bench-seed=" << details.schedule_seed << ")\n";
    if (details.game_seed) out << "  Game-series seed: " << *details.game_seed << '\n';
    out << "  Preprocessing: " << (details.preprocess ? "ON (once per game, excluded from the solve time)" : "OFF") << '\n';
    if (totals.prep_failures) out << "  Shared preparation failures: " << totals.prep_failures << '\n';
    if (totals.preprocessed) out << "  Fully solved by preprocessing: " << totals.preprocessed << '/' << totals.planned << "; no solver trials for these games.\n";
    if (reference) out << "  Reference: " << ui::clean(solvers[*totals.reference].name) << ". Speedup and W/T/L compare each solver with this reference.\n";
    if (!details.jsonl.empty()) out << "  Raw samples: " << details.jsonl << '\n';
    if (!details.csv.empty()) out << "  Per-game statistics: " << details.csv << '\n';
    if (details.interrupted) out << "  Interrupted; the completed results are kept.\n";
}

}
