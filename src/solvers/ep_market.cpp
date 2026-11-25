// SPDX-License-Identifier: Apache-2.0

#include "ep_market.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>
#include <map>

namespace pg::ep {
namespace {

// Dinic's algorithm uses a number of augmentations polynomial in the graph
// size, independent of the magnitudes of the integer capacities.
class MinCut
{
    struct Edge {
        std::size_t to;
        std::size_t reverse;
        Integer capacity;
    };

    static constexpr auto unreachable = std::numeric_limits<std::size_t>::max();
    std::vector<std::vector<Edge>> edges;
    std::vector<std::size_t> level, next, queue, path;
    const std::size_t source, sink;

    bool build_levels()
    {
        std::fill(level.begin(), level.end(), unreachable);
        level[source] = 0;
        queue.assign(1, source);
        for (std::size_t i = 0; i < queue.size(); ++i) {
            const auto vertex = queue[i];
            for (const auto& edge : edges[vertex]) {
                if (edge.capacity == 0 || level[edge.to] != unreachable) continue;
                level[edge.to] = level[vertex] + 1;
                queue.push_back(edge.to);
            }
        }
        return level[sink] != unreachable;
    }

    Integer augment()
    {
        path.assign(1, source);
        while (!path.empty() && path.back() != sink) {
            const auto from = path.back();
            auto& i = next[from];
            while (i < edges[from].size()) {
                const auto& edge = edges[from][i];
                if (edge.capacity > 0 && level[edge.to] == level[from] + 1) break;
                ++i;
            }
            if (i < edges[from].size()) {
                path.push_back(edges[from][i].to);
            } else {
                path.pop_back();
                if (!path.empty()) ++next[path.back()];
            }
        }
        if (path.empty()) return 0;

        Integer sent = edges[source][next[source]].capacity;
        path.pop_back();
        for (auto from : path) sent = std::min(sent, edges[from][next[from]].capacity);
        for (auto from : path) {
            auto& edge = edges[from][next[from]];
            edge.capacity -= sent;
            edges[edge.to][edge.reverse].capacity += sent;
        }
        return sent;
    }

public:
    MinCut(std::size_t source, std::size_t sink)
        : edges(sink + 1), source(source), sink(sink)
    {}

    std::size_t add_vertex()
    {
        edges.emplace_back();
        return edges.size() - 1;
    }

    void add_edge(std::size_t from, std::size_t to, const Integer& capacity)
    {
        if (capacity == 0) return;
        edges[from].push_back({to, edges[to].size(), capacity});
        edges[to].push_back({from, edges[from].size() - 1, 0});
    }

    Integer solve()
    {
        level.resize(edges.size());
        next.resize(edges.size());
        Integer flow = 0;
        while (build_levels()) {
            std::fill(next.begin(), next.end(), 0);
            while (Integer sent = augment()) flow += sent;
        }
        return flow;
    }

    bool on_source_side(std::size_t vertex) const
    {
        return level[vertex] != unreachable;
    }
};

struct Charges {
    Integer hits_moved;
    Integer hits_unmoved;
};

using Coverage = std::map<GroupSlots, Charges>;

void add_group_coverage(Coverage& coverage, std::size_t group, const GroupSlots& slots,
                        const GroupCost& group_cost, const std::vector<Integer>& prices,
                        const std::vector<bool>& movable, const Integer& step,
                        bool increase, bool reflected)
{
    GroupSlots variables{};
    std::array<unsigned, 4> option_bits{};
    GroupPrices high, low;
    unsigned count = 0;
    // Reflection negates the price and swaps the two slots of each coordinate.
    const bool lowers_price = increase == reflected;
    for (std::size_t i = 0; i < slots.size(); ++i) {
        assert(slots[i] < prices.size());
        const auto slot = slots[i] ^ std::size_t(reflected);
        high[i] = prices[slot];
        if (reflected) high[i] = -high[i];
        low[i] = high[i];
        if (!movable[slot]) continue;

        // Options on the same slot move together.
        const auto end = variables.begin() + count;
        const auto variable = std::find(variables.begin(), end, slot);
        if (variable == end) variables[count++] = slot;
        option_bits[i] = 1u << (variable - variables.begin());
        if (lowers_price) low[i] -= step;
        else high[i] += step;
    }
    if (count == 0) return;

    // From the all-high prices, an expert's utility increases by the largest
    // improvement among the options given low prices. This is a sum of
    // nonnegative charges for hitting the sets at successive improvement levels.
    const unsigned full = (1u << count) - 1;
    std::array<Integer, 16> values, charges;
    for (unsigned mask = 0; mask <= full; ++mask) {
        GroupPrices candidate = high;
        for (std::size_t i = 0; i < slots.size(); ++i) {
            if (mask & option_bits[i]) candidate[i] = low[i];
        }
        values[mask] = group_cost(group, candidate);
    }

    // The utility lost by omitting a set of low-price slots is the sum of
    // charges whose entire support was omitted. Invert this subset sum.
    for (unsigned mask = 0; mask <= full; ++mask) charges[mask] = values[full] - values[full ^ mask];
    for (unsigned bit = 1; bit <= full; bit <<= 1) {
        for (unsigned mask = 1; mask <= full; ++mask) {
            if (mask & bit) charges[mask] -= charges[mask ^ bit];
        }
    }
    for (unsigned mask = 1; mask <= full; ++mask) {
        assert(charges[mask] >= 0);
        if (charges[mask] == 0) continue;
        GroupSlots key;
        key.fill(prices.size());
        unsigned size = 0;
        for (unsigned i = 0; i < count; ++i) {
            if (mask & (1u << i)) key[size++] = variables[i];
        }
        std::sort(key.begin(), key.end());
        auto& charge = coverage[key];
        if (lowers_price) charge.hits_moved += charges[mask];
        else charge.hits_unmoved += charges[mask];
    }
}

struct Move {
    Integer change;
    std::vector<bool> slots;
};

Move best_move(const std::vector<Integer>& prices, const std::vector<Integer>& limits,
               const Integer& step, bool increase, const std::vector<GroupSlots>& groups,
               const GroupCost& group_cost)
{
    const auto slot_count = prices.size();
    std::vector<bool> movable(slot_count);
    for (std::size_t i = 0; i < slot_count; ++i) {
        if (increase) movable[i] = prices[i] + step <= limits[i];
        else movable[i] = prices[i] - step >= limits[i];
    }

    Coverage coverage;
    for (std::size_t group = 0; group < groups.size(); ++group) {
        add_group_coverage(coverage, group, groups[group], group_cost,
                           prices, movable, step, increase, false);
        add_group_coverage(coverage, group, groups[group], group_cost,
                           prices, movable, step, increase, true);
    }

    // The empty move pays every charge for hitting unmoved slots. Its cost
    // bounds the minimum cut and lets us encode implications with finite arcs.
    Integer baseline = 0;
    for (const auto& entry : coverage) baseline += entry.second.hits_unmoved;
    if (baseline == 0) return {};
    const Integer implication = baseline + 1;
    const auto source = slot_count;
    const auto sink = slot_count + 1;
    MinCut cut(source, sink);
    for (const auto& [slots, charge] : coverage) {
        if (slots[1] == slot_count) {
            cut.add_edge(slots[0], sink, charge.hits_moved);
            cut.add_edge(source, slots[0], charge.hits_unmoved);
            continue;
        }
        // A moved slot lies on the source side of the cut.
        if (charge.hits_moved > 0) {
            const auto vertex = cut.add_vertex();
            cut.add_edge(vertex, sink, charge.hits_moved);
            for (auto slot : slots) {
                if (slot != slot_count) cut.add_edge(slot, vertex, implication);
            }
        }
        if (charge.hits_unmoved > 0) {
            const auto vertex = cut.add_vertex();
            cut.add_edge(source, vertex, charge.hits_unmoved);
            for (auto slot : slots) {
                if (slot != slot_count) cut.add_edge(vertex, slot, implication);
            }
        }
    }
    const Integer change = cut.solve() - baseline;
    if (change == 0) return {};
    Move move{change, std::vector<bool>(slot_count)};
    for (std::size_t i = 0; i < slot_count; ++i) move.slots[i] = cut.on_source_side(i);
    return move;
}

}

std::vector<Integer>
market_prices(std::size_t coordinates, const std::vector<GroupSlots>& groups,
              const Integer& bid_bound, const GroupCost& group_cost)
{
    assert(bid_bound >= 0);
    std::vector<Integer> result(coordinates);
    if (bid_bound == 0 || groups.empty()) return result;

    // Use independent signed-slot prices z. If G is the market objective on
    // these prices and J swaps and negates opposite slots, minimize G(z)+G(Jz).
    // Clamping a slot outside the bid bounds leaves one of the terms unchanged
    // because its price exceeds all bids, and raises the price in the other.
    const auto slot_count = 2 * coordinates;
    std::vector<Integer> prices(slot_count), lower(slot_count), upper(slot_count);
    const Integer negative_bound = -bid_bound;

    // The objective is L-natural convex. Rounding the midpoint prices up
    // rounds an expert's utility down, and vice versa, since bids are integral.
    // The rounded utilities sum to twice the exact midpoint utility. Convexity
    // bounds that sum by the sum at the two endpoints. Reflection swaps the
    // two roundings, so the same inequality holds for the reflected term.
    // Use Murota (2010), Section 4.2. The previous grid's optimum is within
    // slot_count steps per coordinate of some optimum on the finer grid.
    // Zero is optimal on the grid with spacing twice the first step.
    Integer step = Integer(1) << boost::multiprecision::msb(bid_bound);
    for (; step > 0; step >>= 1) {
        const Integer radius = step * slot_count;
        for (std::size_t i = 0; i < slot_count; ++i) {
            lower[i] = std::max(negative_bound, Integer(prices[i] - radius));
            upper[i] = std::min(bid_bound, Integer(prices[i] + radius));
        }
        // Steepest descent needs at most the sum of the positive and negative
        // sup distances to a minimizer (Murota and Shioura, Theorem 1.2).
        // From the center of this box each distance is at most slot_count.
        std::size_t steps = 0;
        for (;;) {
            const Move up = best_move(prices, upper, step, true, groups, group_cost);
            const Move down = best_move(prices, lower, step, false, groups, group_cost);
            const bool increase = up.change <= down.change;
            const Move& best = increase ? up : down;
            if (best.change >= 0) break;
            if (++steps > 2 * slot_count) throw std::logic_error("market descent exceeded its proximity bound");
            const Integer delta = increase ? step : Integer(-step);
            for (std::size_t i = 0; i < slot_count; ++i) {
                if (best.slots[i]) prices[i] += delta;
            }
        }
    }

    // The integer optimum minimizes the continuous objective too. The regions
    // where it is affine have integral vertices, since their constraints are
    // integer bounds on individual prices and pairwise differences. Averaging
    // z with its reflection yields a signed price vector of the same value.
    for (std::size_t i = 0; i < coordinates; ++i) result[i] = prices[2 * i] - prices[2 * i + 1];
    return result;
}

}
