// SPDX-License-Identifier: Apache-2.0
#include "ep_arena.hpp"

#include <ostream>
#include <utility>

namespace pg::ep {

namespace {
Integer assign_priority_weights(const Arena& graph, Prices& weights)
{
    Integer totals[2] = {0, 0};
    for (int v = graph.size() - 1; v >= 0; --v)
        weights[v] = next_priority_weight(graph.colors[v], totals);
    return totals[0] + totals[1];
}

// A cycle whose top rank is v lies in the SCC of v induced by ranks at
// or below it. Only earlier opposite-parity weights there can oppose v.
// A vertex on no cycle in its prefix can never be the top of a cycle.
void assign_cycle_weights(const Arena& arena, Prices& reward)
{
    const int n = arena.size();
    std::fill(reward.begin(), reward.end(), 0);
    std::vector<std::vector<int>> predecessors(n);
    for (int v = 0; v < n; ++v) {
        for (int to : arena.successors(v)) predecessors[to].push_back(v);
    }
    std::vector<int> forward(n, -1), backward(n, -1), stack;
    for (int top = n - 1; top >= 0; --top) {
        forward[top] = top;
        stack.push_back(top);
        bool on_cycle = false;
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            for (int to : arena.successors(v)) {
                if (to == top) on_cycle = true;
                if (to < top || forward[to] == top) continue;
                forward[to] = top;
                stack.push_back(to);
            }
        }
        if (!on_cycle) continue;
        Integer magnitude = 1;
        backward[top] = top;
        stack.push_back(top);
        while (!stack.empty()) {
            const int v = stack.back();
            stack.pop_back();
            for (int from : predecessors[v]) {
                if (forward[from] != top || backward[from] == top) continue;
                backward[from] = top;
                stack.push_back(from);
                if (arena.colors[from] != arena.colors[top]) magnitude += abs(reward[from]);
            }
        }
        reward[top] = arena.colors[top] ? -magnitude : magnitude;
    }
}

} // namespace

EnergyMeasure::EnergyMeasure(const Arena& arena) : arena(arena), weight(arena.size())
{
    const Integer n = arena.size();
    if (assign_priority_weights(arena, weight) > (n + 1) * (n + 1))
        assign_cycle_weights(arena, weight);
    // Any simple path gains at most the positive sum and loses at most the
    // negative magnitude. The two winners' credit ranges are one unit apart.
    for (const auto& w : weight) {
        if (w > 0) lower -= w;
        else upper -= w;
    }
}

Integer EnergyMeasure::lift(int v, const Integer& successor_score) const
{
    const Integer y = weight[v] + successor_score;
    return std::clamp(y, lower, upper);
}

Integer EnergyMeasure::lift(int v, const Prices& x) const
{
    return lift(v, x[arena.best_successor(v, x)]);
}

std::vector<int> EnergyMeasure::solution(const Prices& x, std::vector<int>& moves) const
{
    arena.greedy_moves(x, moves);
    std::vector<int> winners(arena.size());
    for (int v = 0; v < arena.size(); ++v) winners[v] = x[v] > 0 ? 0 : 1;
    if (!arena.certifies(winners, moves)) THROW_ERROR("EP energy measures failed parity verification");
    return winners;
}

namespace {
// With any boundary values pinned, the other equations have a unique integer
// solution. A positive maximum difference between two solutions must propagate
// unchanged along an unclipped extremal edge, and cannot end in a zero-sum
// cycle. The same argument bounds changes by the largest boundary change.
// Thus the scalar return map minus its trial coordinate is nonincreasing.
template<class Restore>
void solve_boundary(const EnergyMeasure& game, const std::vector<int>& boundary,
                    std::size_t index, Prices& x, const Restore& restore)
{
    if (index == boundary.size()) {
        restore();
        return;
    }
    const int v = boundary[index];
    Integer lo = game.lower, hi = game.upper;
    while (lo <= hi) {
        x[v] = lo + (hi - lo) / 2;
        solve_boundary(game, boundary, index + 1, x, restore);
        const Integer next = game.lift(v, x);
        if (next == x[v]) return;
        if (next > x[v]) lo = x[v] + 1;
        else hi = x[v] - 1;
    }
    THROW_ERROR("EP lost the integer boundary fixed point");
}

// The heads of DFS back edges meet every directed cycle. Deleting these
// vertices leaves an acyclic graph in reverse DFS finishing order.
void feedback_boundary(const Arena& arena, std::vector<int>& boundary, std::vector<int>& downstream_first)
{
    std::vector<bool> pinned(arena.size());
    std::vector<int> state(arena.size());
    std::vector<std::pair<int, std::size_t>> stack;
    for (int root = 0; root < arena.size(); ++root) {
        if (state[root]) continue;
        state[root] = 1;
        stack.emplace_back(root, 0);
        while (!stack.empty()) {
            const int v = stack.back().first;
            auto& edge = stack.back().second;
            const auto successors = arena.successors(v);
            if (edge == successors.size()) {
                state[v] = 2;
                downstream_first.push_back(v);
                stack.pop_back();
                continue;
            }
            const int to = successors[edge++];
            if (state[to] == 1) pinned[to] = true;
            if (state[to] != 0) continue;
            state[to] = 1;
            stack.emplace_back(to, 0);
        }
    }
    for (int v = 0; v < arena.size(); ++v) if (pinned[v]) boundary.push_back(v);
    downstream_first.erase(std::remove_if(downstream_first.begin(), downstream_first.end(),
        [&](int v) { return pinned[v]; }), downstream_first.end());
}

// Pin the smaller of a feedback set and the branching vertices of either
// owner. The remaining rows are acyclic or have only one player with a choice.
struct FreeRows
{
    std::vector<int> boundary, rows, saturated;
    int controller = -1;

    explicit FreeRows(const Arena& arena)
    {
        feedback_boundary(arena, boundary, rows);
        auto choices = arena.branching_vertices();
        for (int player : {0, 1}) {
            if (choices[player].size() < boundary.size()) {
                boundary = std::move(choices[player]);
                controller = 1 - player;
            }
        }
        if (controller >= 0)
            find_cycles(arena, [&](int v) { return arena.successors(v); });
    }

    FreeRows(const Arena& arena, int player, const std::vector<int>& policy) : controller(1 - player)
    {
        find_cycles(arena, [&](int v) { return arena.policy_edges(v, player, policy); });
    }

    template<class Edges>
    void find_cycles(const Arena& arena, const Edges& out)
    {
        // The controller can follow any free cycle. At the top of a cycle
        // of its own parity every nonempty prefix has that sign, and a full
        // lap has nonzero weight. Its exact score is the favorable wall.
        rows.clear();
        std::vector<int> seen(arena.size(), -1), stack;
        for (int v : boundary) seen[v] = -2;
        for (int v = 0; v < arena.size(); ++v) {
            if (seen[v] == -2) continue;
            if (arena.colors[v] == controller && arena.cycle_at_rank(v, out, seen, stack)) {
                saturated.push_back(v);
                seen[v] = -2;
            } else {
                rows.push_back(v);
            }
        }
    }

    template<class Update>
    void restore(Prices& x, const Integer& lower, const Integer& upper, const Update& update) const
    {
        if (controller < 0) {
            for (int v : rows) x[v] = update(v);
            return;
        }
        const Integer& favorable = controller == 0 ? upper : lower;
        const Integer& opposite = controller == 0 ? lower : upper;
        for (int v : saturated) x[v] = favorable;
        for (int v : rows) x[v] = opposite;
        // Each surviving cycle has the opposite sign. An optimal policy
        // on rows strictly off the starting wall therefore has no cycle.
        for (std::size_t pass = 0; pass < rows.size(); ++pass) {
            bool changed = false;
            for (int v : rows) {
                Integer next = update(v);
                if (next == x[v]) continue;
                x[v] = std::move(next);
                changed = true;
            }
            if (!changed) break;
        }
    }
};

} // namespace

// The signs point toward the unique solution on every coordinate attaining
// the largest absolute error. An empty result means the query is exact.
std::vector<int> EnergyMeasure::price_directions(const Prices& trial, const Integer& unit,
                                                  std::vector<int>& greedy) const
{
    const Integer floor = unit * lower, ceiling = unit * upper;
    Prices x = trial;
    for (auto& value : x) value = std::clamp(value, floor, ceiling);
    arena.greedy_moves(x, greedy);
    std::vector<int> labels(arena.size()), ready;
    std::vector<std::vector<int>> predecessors(arena.size());
    std::vector<std::size_t> remaining(arena.size());
    bool fixed = true;
    for (int v = 0; v < arena.size(); ++v) {
        const Integer raw = unit * weight[v] + x[greedy[v]];
        fixed &= trial[v] == std::clamp(raw, floor, ceiling);
        if (trial[v] != x[v]) labels[v] = trial[v] < x[v] ? 1 : -1;
        else if (x[v] != raw) labels[v] = raw > x[v] ? 1 : -1;
        if (labels[v] != 0) {
            ready.push_back(v);
        } else {
            for (int to : arena.successors(v)) {
                if (x[to] != x[greedy[v]]) continue;
                predecessors[to].push_back(v);
                ++remaining[v];
            }
        }
    }
    if (fixed) return {};
    // Tight equalities telescope around any cycle, so these edges form a DAG.
    // Even takes OR of the positive replies, Odd takes AND.
    for (std::size_t i = 0; i < ready.size(); ++i) {
        const int to = ready[i];
        for (int v : predecessors[to]) {
            if (labels[v] != 0) continue;
            if (--remaining[v] == 0 || labels[to] == (arena.owners[v] ? -1 : 1)) {
                labels[v] = labels[to];
                ready.push_back(v);
            }
        }
    }
    if (ready.size() != x.size()) THROW_ERROR("EP found a zero-weight cycle of tight edges");
    return labels;
}

// Ascend from the lower wall and descend from the upper wall. Each changing
// integer sweep consumes at least one unit of the sum of remaining distances.
// Use the full bound only when the width itself is polynomial in the arena.
std::vector<int> lift_measures(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace)
{
    const Arena& arena = game.arena;
    const Integer n = arena.size(), width = game.upper - game.lower;
    const Integer limit = width <= (n + 1) * (n + 1) + 1 ? Integer(n * width + 1) : n;
    Prices low(arena.size(), game.lower), high(arena.size(), game.upper);
    for (Integer round = 0; round <= limit; ++round) {
        auto winners = arena.supported_dominions(game.weight, low, high, 1, moves);
        if (!winners.empty()) {
            if (trace) *trace << "EP certified a dominion after " << round << " energy sweeps\n";
            return winners;
        }
        if (round == limit) break;
        bool lower_changed = false, upper_changed = false;
        for (int v = 0; v < arena.size(); ++v) {
            Integer l = game.lift(v, low), h = game.lift(v, high);
            lower_changed |= l != low[v];
            upper_changed |= h != high[v];
            low[v] = std::move(l);
            high[v] = std::move(h);
        }
        if (!lower_changed || !upper_changed) return game.solution(lower_changed ? high : low, moves);
    }
    return {};
}

std::vector<int> improve_policy(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace)
{
    const Arena& arena = game.arena;
    const auto branching = arena.branching_vertices();
    const int player = branching[1].size() < branching[0].size() ? 1 : 0;
    const auto& choices = branching[player];
    Prices x(arena.size());
    std::vector<int> policy(arena.size());
    arena.greedy_moves(x, policy);

    // Evaluate the initial policy and at most one improvement round per
    // branching vertex. Each evaluation keeps every opponent reply. Only
    // the certificate on the full arena can publish a winning region.
    const auto limit = choices.size() + 1;
    for (std::size_t evaluation = 1; evaluation <= limit; ++evaluation) {
        const FreeRows rows(arena, player, policy);
        rows.restore(x, game.lower, game.upper, [&](int v) {
            const int to = arena.owners[v] == player ? policy[v] : arena.best_successor(v, x);
            return game.lift(v, x[to]);
        });
        auto winners = arena.supported_dominions(game.weight, x, x, 1, moves);
        if (!winners.empty()) {
            if (trace) *trace << "EP certified a dominion after " << evaluation
                              << " policy evaluation" << (evaluation == 1 ? "\n" : "s\n");
            return winners;
        }
        if (evaluation == limit) break;
        for (int v : choices) {
            const int to = arena.best_successor(v, x);
            if (x[to] != x[policy[v]]) policy[v] = to;
        }
    }
    return {};
}

std::vector<int> solve_measures(const EnergyMeasure& game, std::vector<int>& moves, std::ostream* trace)
{
    const Arena& arena = game.arena;
    const FreeRows free(arena);
    const auto& boundary = free.boundary;
    if (trace) *trace << "EP solving " << boundary.size() << " boundary coordinates using "
                      << (free.controller < 0 ? "DAG" : "one-player") << " substitution\n";
    Prices x(arena.size());
    if (boundary.size() <= 2) {
        solve_boundary(game, boundary, 0, x, [&] {
            free.restore(x, game.lower, game.upper, [&](int v) { return game.lift(v, x); });
        });
        return game.solution(x, moves);
    }

    // A single hypothesis has demands +s_i and -s_i on both signs of every
    // boundary coordinate. An exact probe stops. At all other probes its
    // largest surplus is the positive maximum error. Free-row restoration
    // is nonexpansive, so the full directions label every boundary maximizer.
    auto directions = [&](const Prices& p, const Integer& unit) {
        const Integer floor = unit * game.lower, ceiling = unit * game.upper;
        for (std::size_t i = 0; i < boundary.size(); ++i)
            x[boundary[i]] = std::clamp(p[i], floor, ceiling);
        free.restore(x, floor, ceiling, [&](int v) -> Integer {
            Integer y = unit * game.weight[v] + x[arena.best_successor(v, x)];
            return std::clamp(y, floor, ceiling);
        });
        for (std::size_t i = 0; i < boundary.size(); ++i) x[boundary[i]] = p[i];
        return game.price_directions(x, unit, moves);
    };
    std::vector<int> winners;
    const Integer radius = std::max(Integer(-game.lower), game.upper);
    ep::extremal_prices(boundary.size(), radius,
        [&](const Prices& p, const Integer& unit, std::vector<int>& reply) {
            const auto labels = directions(p, unit);
            if (labels.empty()) {
                winners = game.solution(x, moves);
                return true;
            }
            winners = arena.supported_dominions(game.weight, x, x, unit, moves);
            if (!winners.empty()) return true;
            reply.clear();
            for (int v : boundary) reply.push_back(labels[v]);
            return false;
        },
        [&](const Prices& p, const Integer& unit) {
            if (!directions(p, unit).empty()) THROW_ERROR("EP prices failed the energy equations");
            winners = game.solution(x, moves);
        });
    return winners;
}

} // namespace pg::ep
