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

#include "ep.hpp"
#include "ep_arena.hpp"

namespace pg::ep {

// A potential supports Even on an edge v->u when x[v] <= w[v]+x[u].
// Keep the greatest subgame where Even has a supported move and every Odd
// move is supported. The inequalities telescope, so all its cycles are
// positive. Reverse the inequalities for Odd. Each edge is visited O(1) times.
std::vector<int> Arena::supported_dominions(const Prices& weight, const Prices& low, const Prices& high,
                                            const Integer& unit, std::vector<int>& moves) const
{
    std::vector<int> winners(size(), -1);
    bool found = false;
    for (int player : {0, 1}) {
        const auto& price = player ? high : low;
        Prices threshold(size());
        for (int v = 0; v < size(); ++v) threshold[v] = price[v] - unit * weight[v];
        auto supports = [&](int v, int to) {
            return player ? price[to] <= threshold[v] : price[to] >= threshold[v];
        };

        std::vector<std::vector<int>> incoming(size());
        std::vector<std::size_t> remaining(size());
        std::vector<bool> bad(size());
        std::vector<int> removed;
        for (int v = 0; v < size(); ++v) {
            for (int to : successors(v)) {
                if (supports(v, to)) {
                    incoming[to].push_back(v);
                    ++remaining[v];
                } else if (owners[v] != player) bad[v] = true;
            }
            if (remaining[v] == 0) bad[v] = true;
            if (bad[v]) removed.push_back(v);
        }
        for (std::size_t i = 0; i < removed.size(); ++i) {
            for (int v : incoming[removed[i]]) {
                if (bad[v]) continue;
                if (owners[v] != player || --remaining[v] == 0) {
                    bad[v] = true;
                    removed.push_back(v);
                }
            }
        }
        for (int v = 0; v < size(); ++v) {
            if (bad[v]) continue;
            found = true;
            winners[v] = player;
            if (owners[v] == player) {
                for (int to : successors(v)) {
                    if (!bad[to] && supports(v, to)) {
                        moves[v] = to;
                        break;
                    }
                }
            }
        }
    }
    if (!found) return {};
    return winners;
}

} // namespace pg::ep

namespace pg {

bool
EPSolver::parseOptions(std::string& options)
{
    if (!configurable) return Solver::parseOptions(options);
    if (!options.empty() && options != "clipped" && options != "sign") return false;
    mode = options == "sign" ? Mode::Sign : Mode::Clipped;
    return true;
}

void
EPSolver::run()
{
    for (;;) {
        const ep::Arena graph(game, disabled);
        if (graph.size() == 0) return;
        const ep::EnergyMeasure measure(graph);
        std::vector<int> moves(graph.size());
        auto* log = trace ? &logger : nullptr;
        auto winners = ep::lift_measures(measure, moves, log);
        if (winners.empty()) winners = ep::improve_policy(measure, moves, log);
        if (winners.empty()) {
            winners = mode == Mode::Sign ? ep::sign_solution(measure, moves, log)
                                         : ep::solve_measures(measure, moves, log);
        }
        if (winners.size() != moves.size()) THROW_ERROR("EP failed to return a winning region");
        const auto undecided = std::count(winners.begin(), winners.end(), -1);
        if (undecided == graph.size()) THROW_ERROR("EP returned an empty dominion");
        for (int v = 0; v < graph.size(); ++v) {
            if (winners[v] < 0) continue;
            const int to = graph.owners[v] == winners[v] ? graph.vertices[moves[v]] : -1;
            solve(graph.vertices[v], winners[v], to);
        }
        if (undecided == 0) return;
        if (trace) logger << "EP certified " << graph.size() - undecided << " of " << graph.size()
                          << " vertices\n";
        const auto before = disabled.count();
        flush();
        if (trace) logger << "EP attractor added " << disabled.count() - before << '\n';
    }
}

} // namespace pg
