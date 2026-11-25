// SPDX-License-Identifier: Apache-2.0

#include "solvers/ep_experts.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace {
using namespace pg::ep;

using Trial = std::array<int, 2>;
using TrialHistory = std::vector<Trial>;

History
grid_history(const TrialHistory& history)
{
    History result;
    for (const auto& prices : history) result.emplace_back(prices.begin(), prices.end());
    return result;
}

std::size_t
find_group(const ExpertPool& pool, const GroupSlots& slots)
{
    const auto& groups = pool.group_slots();
    const auto match = std::find(groups.begin(), groups.end(), slots);
    if (match == groups.end()) throw std::runtime_error("missing sorted expert group");
    return static_cast<std::size_t>(match - groups.begin());
}

class ExplicitPool
{
    struct Expert {
        GroupSlots slots;
        std::array<int, 4> bids;
        std::size_t group;
        int pivot_level = -1;
        std::size_t pivot_slot = 0;
        int pivot_bid = 0;

        std::pair<int, int> evaluate(const Trial& prices, int scale, const TrialHistory& history,
                                     const Trial& feedback) const
        {
            int threshold = 0;
            if (pivot_level >= 0) {
                if (std::size_t(pivot_level) >= history.size()) return {0, 0};
                threshold = scale * pivot_bid - history[std::size_t(pivot_level)][pivot_slot];
            }
            int utility = 0, gain = 0;
            for (unsigned i = 0; i < slots.size(); ++i) {
                const int value = scale * bids[i] - threshold - prices[slots[i]];
                if (value > utility) {
                    utility = value;
                    gain = feedback[slots[i]];
                } else if (value == utility) {
                    gain = std::min(gain, feedback[slots[i]]);
                }
            }
            return {utility, gain};
        }
    };

    struct Frame {
        Trial prices;
        int scale;
        TrialHistory history;
        std::vector<Integer> weights;
    };

    ExpertPool& pool;
    std::vector<Expert> experts;
    std::vector<Integer> weights;
    std::vector<Frame> frames;

public:
    ExplicitPool(ExpertPool& pool, const Integer& initial_wealth) : pool(pool)
    {
        // Two signed slots, three possible demands, and three pivot levels.
        for (int code = 0; code < 6*6*6*6; ++code) {
            Expert expert{};
            int choices = code;
            for (unsigned i = 0; i < expert.slots.size(); ++i) {
                const int choice = choices % 6;
                choices /= 6;
                expert.slots[i] = std::size_t(choice % 2);
                expert.bids[i] = choice / 2 - 1;
            }
            auto sorted = expert.slots;
            std::sort(sorted.begin(), sorted.end());
            expert.group = find_group(pool, sorted);
            experts.push_back(expert);
            for (expert.pivot_level = 0; expert.pivot_level < 3; ++expert.pivot_level)
                for (expert.pivot_slot = 0; expert.pivot_slot < 2; ++expert.pivot_slot)
                    for (expert.pivot_bid = -1; expert.pivot_bid <= 1; ++expert.pivot_bid)
                        experts.push_back(expert);
        }
        weights.assign(experts.size(), initial_wealth);
        if (ExpertPool::count(1, 1, 3) != experts.size()) throw std::runtime_error("incorrect ordered expert count");
    }

    void check(const Trial& prices, int scale, const TrialHistory& history) const
    {
        std::vector<Integer> costs(pool.group_slots().size());
        for (std::size_t i = 0; i < experts.size(); ++i) {
            const int utility = experts[i].evaluate(prices, scale, history, {}).first;
            costs[experts[i].group] += weights[i] * utility;
        }
        auto ancestors = grid_history(history);
        for (std::size_t group = 0; group < costs.size(); ++group) {
            GroupPrices selected;
            const auto& slots = pool.group_slots()[group];
            for (unsigned i = 0; i < slots.size(); ++i) selected[i] = prices[slots[i]];
            if (pool.group_objective(group, selected, scale, ancestors) != costs[group])
                throw std::runtime_error("compressed expert objective differs from exact enumeration");
        }
    }

    void open(const Trial& prices, int scale, const TrialHistory& history)
    {
        pool.open_probe(Prices(prices.begin(), prices.end()), scale, grid_history(history));
        frames.push_back({prices, scale, history, weights});
    }

    void close(const Trial& feedback)
    {
        pool.close_probe(std::vector<int>(feedback.begin(), feedback.end()));
        const auto& parent = frames.back();
        // Compute the ceiling separately for each ordered expert.
        for (std::size_t i = 0; i < experts.size(); ++i) {
            const int gain = experts[i].evaluate(parent.prices, parent.scale, parent.history, feedback).second;
            const Integer total = parent.weights[i] * (16383 * (16 + gain)) + 16 * weights[i];
            Integer quotient, remainder;
            boost::multiprecision::divide_qr(total, Integer(262144), quotient, remainder);
            weights[i] = quotient + (remainder != 0);
        }
        frames.pop_back();
    }
};

void
check_distinct_slots()
{
    ExpertPool pool(2, 1, 2);
    const History history{{-2, 1, 3, -3}};
    pool.open_probe(history.front(), 2, {});

    const GroupSlots distinct{0, 1, 2, 3};
    const GroupPrices prices{1, -3, 1, -2};
    auto group = find_group(pool, distinct);
    Integer expected = 0;
    for (int code = 0; code < 81; ++code) {
        int choices = code;
        Integer best = -10;
        for (unsigned i = 0; i < prices.size(); ++i) {
            int bid = choices % 3 - 1;
            choices /= 3;
            best = std::max(best, Integer(2 * bid - prices[i]));
        }
        expected += std::max(Integer(0), best);
        for (const auto& pivot_price : history.front())
            for (int bid = -1; bid <= 1; ++bid)
                expected += std::max(Integer(0), Integer(best - 2 * bid + pivot_price));
    }
    // The sorted tuple represents the 24 permutations of the distinct slots.
    if (pool.group_objective(group, prices, 2, history) != 24 * expected)
        throw std::runtime_error("five-guess expert group differs from exact enumeration");
}

void check_symmetric_market()
{
    ExpertPool pool(2, 1, 2, Integer(1) << 24);
    const Prices origin(4);
    const std::vector<int> same_feedback{1, -1, 1, -1};
    for (int i = 0; i < 32; ++i) {
        pool.open_probe(origin, 1, {});
        pool.close_probe(same_feedback);
    }
    const Prices reduced = pool.minimizing_prices(2, 1, {});
    if (reduced[0] != reduced[1] || reduced[0] <= 0)
        throw std::runtime_error("equivalent coordinates failed to share a nonzero price");
    const Prices signed_price{reduced[0], -reduced[0], reduced[1], -reduced[1]};
    const Integer minimum = pool.objective(signed_price, 2, {});
    for (int a = -4; a <= 4; ++a) {
        for (int b = -4; b <= 4; ++b) {
            if (pool.objective({a, -a, b, -b}, 2, {}) < minimum)
                throw std::runtime_error("symmetric price failed full half-grid comparison");
        }
    }

    // Merging positive slots 0 and 2 changes the tuple's permutation count
    // from 12 to 4. Its ordered-guess objective, including all four possible
    // pivot slots of the pending ancestor, must be identical.
    pool.open_probe(origin, 1, {});
    const History history{origin};
    const Integer original = pool.group_objective(find_group(pool, {0, 0, 2, 3}), {1, 1, 1, -2}, 2, history);
    const Integer repeated = pool.group_objective(find_group(pool, {0, 0, 0, 1}), {1, 1, 1, -2}, 2, history);
    if (original != 3 * repeated) throw std::runtime_error("class representative changed ordered expert wealth");
    pool.close_probe(same_feedback);

    ExpertPool half(3, 1, 1, 12345);
    for (int i = 0; i < 7; ++i) {
        half.open_probe(Prices(6), 1, {});
        half.close_probe({1, -1, 1, -1, 1, -1});
    }
    if (half.minimizing_prices(1, 1, {}) != Prices(3, 1))
        throw std::runtime_error("class-common market missed the half-integral optimum");
}

}

int test_ep_experts()
{
    const Integer initial_wealth = 12345;
    ExpertPool pool(1, 1, 3, initial_wealth);
    ExplicitPool brute(pool, initial_wealth);
    brute.check({0, 0}, 1, {});
    brute.open({0, 0}, 1, {});
    brute.check({-1, -7}, 4, {{0, 0}});
    brute.open({3, -2}, 2, {{0, 0}});
    brute.open({0, 3}, 4, {{0, 0}, {6, -4}});
    brute.close({1, -1});
    brute.check({-1, 9}, 4, {{0, 0}, {6, -4}});
    brute.close({-1, 1});
    brute.check({3, -21}, 16, {{0, 0}});
    brute.open({-4, 4}, 4, {{0, 0}});
    brute.check({-7, -12}, 8, {{0, 0}, {-8, 8}});
    brute.close({-1, 1});
    brute.close({1, -1});
    brute.check({-8, -8}, 2, {});

    check_distinct_slots();
    check_symmetric_market();

    Integer radius = Integer(1) << 200;
    ExpertPool large(1, radius, 5);
    Integer n = 2 * radius + 1, m = n - 1;
    // At these prices the utility is the maximum of four guesses in [1,n].
    // Summing the differences in successive fourth powers gives this total.
    Integer fourth_sum = m * (m+1) * (2*m+1) * (3*m*m+3*m-1) / 30;
    if (large.objective({-radius-1, -radius-1}, 1, {}) != 16 * (n*n*n*n*n - fourth_sum))
        throw std::runtime_error("expert objective differs from the closed-form large-radius sum");
    return 0;
}
