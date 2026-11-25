// SPDX-License-Identifier: Apache-2.0
// Extremal prices, Section 6 of "Strongly polynomial linear programming"
// (OpenAI, November 24, 2026). Prices in the recursion share the denominator
// 2^depth_limit.
#include "ep_search.hpp"
#include "ep_experts.hpp"
#include "ep_market.hpp"

#include <algorithm>

namespace pg::ep {
namespace {

std::size_t depth_limit(std::size_t n)
{
    Integer numerator = 5, denominator = 4;
    std::size_t depth = 1;
    while (numerator <= 2 * Integer(n) * denominator) {
        numerator *= 25;
        denominator *= 16;
        depth += 2;
    }
    return depth;
}

Prices median(const History& samples)
{
    const std::size_t count = samples.size();
    Prices result(samples.front().size()), column(count);
    const auto middle = column.begin() + count / 2;
    for (std::size_t i = 0; i < result.size(); ++i) {
        for (std::size_t j = 0; j < count; ++j) column[j] = samples[j][i];
        std::sort(column.begin(), column.end());
        result[i] = *middle;
    }
    return result;
}

std::vector<std::size_t> required_agreements(std::size_t n, std::size_t depth)
{
    // A slot starts with mass at least 1/(2n). Every ancestor divides mass by
    // 8/5. A full boost also multiplies it by four (paper, Section 6.2).
    std::vector<std::size_t> counts(depth);
    Integer numerator = 1, denominator = 2 * Integer(n);
    for (std::size_t level = 0; level < depth; ++level) {
        Integer boosted = numerator;
        while (boosted <= denominator) {
            boosted <<= 2;
            ++counts[level];
        }
        numerator *= 5;
        denominator *= 8;
    }
    return counts;
}

bool certifies_median(const History& history, const Prices& candidate, std::size_t required)
{
    // If a signed slot were violated at the node's threshold, each ancestor
    // at or below its signed candidate price would have fully boosted it.
    // Enough such ancestors would give the slot mass greater than one.
    for (std::size_t v = 0; v < candidate.size(); ++v) {
        std::size_t below = 0, above = 0;
        for (const auto& ancestor : history) {
            below += ancestor[v] <= candidate[v];
            above += ancestor[v] >= candidate[v];
        }
        if (below < required || above < required) return false;
    }
    return true;
}

template<class Vector>
Vector signed_slots(const Vector& coordinates)
{
    Vector result;
    result.reserve(2 * coordinates.size());
    for (const auto& value : coordinates) {
        result.push_back(value);
        result.push_back(-value);
    }
    return result;
}

History signed_history(const History& history, std::size_t shift)
{
    History result;
    for (const auto& vector : history) {
        auto slots = signed_slots(vector);
        // All ancestors' prices are exactly representable on this grid.
        for (auto& value : slots) value >>= shift;
        result.push_back(std::move(slots));
    }
    return result;
}

class Search
{
public:
    Search(std::size_t n, const Integer& radius, detail::ProbeSchedule schedule,
           const Integer& initial_wealth, const Probe& query) :
        dimension(n), schedule(schedule), demand(radius),
        unit(Integer(1) << schedule.depth), agreements(required_agreements(n, schedule.depth)),
        probe(query), pool(n, radius, schedule.depth, initial_wealth)
    {
    }

    void run(const Integer& root_probes, const Decode& decode)
    {
        Prices result;
        if (visit(root_probes, true, result)) return;
        decode(result, unit);
    }

private:
    // True propagates an independently certified solution out of the tree.
    bool visit(const Integer& count, bool uniform, Prices& result)
    {
        const auto level = history.size();
        const auto depth = schedule.depth;
        if (level == depth) {
            result = median(history);
            return false;
        }
        if (schedule.use_early_leaves && agreements[level] <= level) {
            Prices candidate = median(history);
            if (certifies_median(history, candidate, agreements[level])) {
                result = std::move(candidate);
                return false;
            }
        }

        const Integer scale = Integer(1) << level;
        const auto ancestors = signed_history(history, depth - level);
        const auto finer_ancestors = signed_history(history, depth - level - 1);
        Integer bound = demand * 2 * scale;
        for (const auto& ancestor : ancestors) {
            for (const auto& value : ancestor) {
                const Integer bid_bound = 2 * demand * scale + abs(value);
                if (bound < bid_bound) bound = bid_bound;
            }
        }

        History children;
        for (Integer i = 0; i < count; ++i) {
            Prices y(dimension, 0);
            if (!uniform) {
                y = pool.minimizing_prices(bound, scale, ancestors);
            }
            // At the first probe along the initial descent the expert weights
            // are uniform and all ancestors are zero. Reversing every slot
            // leaves the market invariant, hence zero is an exact minimizer.
            pool.open_probe(signed_slots(y), 2 * scale, finer_ancestors);
            for (auto& value : y) value <<= depth - level - 1;
            history.push_back(y);
            Prices x;
            if (visit(schedule.probes, uniform, x)) return true;
            history.pop_back();
            std::vector<int> feedback;
            if (probe(x, unit, feedback)) return true;
            for (std::size_t v = 0; v < dimension; ++v) {
                if (x[v] != y[v]) feedback[v] = x[v] > y[v] ? 1 : -1;
            }
            pool.close_probe(signed_slots(feedback));
            children.push_back(std::move(x));
            uniform = false;
        }
        result = median(children);
        return false;
    }

    std::size_t dimension;
    detail::ProbeSchedule schedule;
    Integer demand, unit;
    std::vector<std::size_t> agreements;
    const Probe& probe;
    ExpertPool pool;
    History history;
};

} // namespace

void detail::price_tree(std::size_t dimension, const Integer& radius, ProbeSchedule schedule,
                        const Probe& label, const Decode& finish)
{
    const Integer experts = ExpertPool::count(dimension, radius, schedule.depth);
    // The bit length of the positive integer |E| equals ceil(log2(|E|+1)).
    const Integer root_probes = schedule.probes + schedule.root_factor * (Integer(boost::multiprecision::msb(experts)) + 1);
    // Q counts settlements across the whole tree. Rounding each expert up by
    // < 1/Q per settlement adds < |E| to total wealth (the actual increment
    // at depth l also has the factor c^l). Thus each root return is <= 2|E|.
    // The lower inequality (6.7) still holds. With c=2^-14 the bonus
    // ln(1+c*return/2) is below log2(1+|E|), as required by Lemma 6.4.
    Integer settlements = root_probes;
    for (std::size_t level = 1; level < schedule.depth; ++level) {
        settlements *= schedule.probes;
        settlements += root_probes;
    }
    Search(dimension, radius, schedule, settlements, label).run(root_probes, finish);
}

void extremal_prices(std::size_t dimension, const Integer& radius,
                     const Probe& probe, const Decode& decode)
{
    detail::price_tree(dimension, radius, {depth_limit(dimension), (1u << 19) + 1, 3072}, probe, decode);
}

} // namespace pg::ep
