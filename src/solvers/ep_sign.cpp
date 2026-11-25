// SPDX-License-Identifier: Apache-2.0
#include "ep_arena.hpp"

#include <ostream>

namespace pg::ep {

namespace {

// On the Even region, lengths n*w[v]-1 have no negative cycle under a
// winning strategy. Distances to a zero-cost exit witness the positive
// price demands. Negating weights and prices gives Odd's negative demands.
std::vector<int> sign_regions(const Arena& graph, const Prices& weight, const Integer& radius,
                              std::vector<int>& moves, std::ostream* trace = nullptr)
{
    const int n = graph.size();
    std::vector<int> winners(n);
    if (trace) *trace << "EP searching for winner signs\n";
    ep::extremal_prices(n, radius,
        [&](const Prices& p, const Integer& unit, std::vector<int>& labels) {
            auto known = graph.supported_dominions(weight, p, p, unit, moves);
            if (!known.empty() && std::find(known.begin(), known.end(), -1) == known.end()) {
                winners = std::move(known);
                return true;
            }
            graph.greedy_moves(p, moves);
            labels.resize(n);
            for (int v = 0; v < n; ++v)
                labels[v] = weight[v] * unit + p[moves[v]] > p[v] ? 1 : -1;
            return false;
        },
        [&](const Prices& p, const Integer&) {
            graph.greedy_moves(p, moves);
            for (int v = 0; v < n; ++v) {
                if (p[v] == 0) THROW_ERROR("EP returned a zero winner price");
                winners[v] = p[v] > 0 ? 0 : 1;
            }
        });
    return winners;
}

// Keep a nonempty contiguous range of moves at just one vertex.
Arena restrict_moves(const Arena& graph, int v, std::size_t lo, std::size_t hi)
{
    Arena trial = graph;
    const auto first = trial.edges.begin() + trial.offsets[v];
    const auto count = trial.offsets[v + 1] - trial.offsets[v];
    trial.edges.erase(first + hi, first + count);
    trial.edges.erase(first, first + lo);
    for (auto i = std::size_t(v) + 1; i < trial.offsets.size(); ++i)
        trial.offsets[i] -= count - (hi - lo);
    return trial;
}

} // namespace

std::vector<int> sign_solution(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace)
{
    const Arena& graph = game.arena;
    Prices weight = game.weight;
    const Integer radius = graph.size() * (game.upper - game.lower) + 1;
    for (auto& w : weight) w *= graph.size();
    const auto winners = sign_regions(graph, weight, radius, moves, trace);
    if (graph.certifies(winners, moves)) return winners;

    // Positional determinacy supplies one policy winning from its whole region.
    // If restricting an owned winning vertex to the first half loses a start,
    // a policy that preserves every winner must use the other half.
    Arena current = graph;
    std::vector<int> trial_moves(graph.size());
    std::size_t decisions = 0;
    for (int v = 0; v < graph.size(); ++v) {
        if (graph.owners[v] != winners[v]) continue;
        const auto out = graph.successors(v);
        std::size_t lo = 0, hi = out.size();
        while (hi - lo > 1) {
            const auto mid = lo + (hi - lo) / 2;
            const auto trial = restrict_moves(current, v, lo, mid);
            ++decisions;
            if (sign_regions(trial, weight, radius, trial_moves) == winners) hi = mid;
            else lo = mid;
        }
        moves[v] = out[lo];
        if (out.size() > 1) current = restrict_moves(current, v, lo, hi);
    }
    if (!graph.certifies(winners, moves)) THROW_ERROR("EP winner prices failed strategy recovery");
    if (trace) *trace << "EP recovered the sign strategies with " << decisions
                      << (decisions == 1 ? " restricted game\n" : " restricted games\n");
    return winners;
}

} // namespace pg::ep
