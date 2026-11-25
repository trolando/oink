// SPDX-License-Identifier: Apache-2.0

#include "solvers/ep_market.hpp"

#include <algorithm>
#include <stdexcept>

namespace {
using namespace pg::ep;

struct Expert {
    GroupPrices bids;
    Integer weight;
};

struct Group {
    GroupSlots slots;
    std::vector<Expert> experts;
};

Integer group_cost(const Group& group, const GroupPrices& prices, int scale)
{
    Integer total = 0;
    for (const auto& expert : group.experts) {
        Integer utility = 0;
        for (std::size_t i = 0; i < prices.size(); ++i) {
            utility = std::max(utility, Integer(scale * expert.bids[i] - prices[i]));
        }
        total += expert.weight * utility;
    }
    return total;
}

std::vector<Integer>
propose(std::size_t coordinates, const Integer& bound, const std::vector<Group>& groups)
{
    std::vector<GroupSlots> slots;
    for (const auto& group : groups) slots.push_back(group.slots);
    return market_prices(coordinates, slots, bound, [&](std::size_t group, const GroupPrices& prices) {
        for (std::size_t i = 1; i < prices.size(); ++i) {
            for (std::size_t j = 0; j < i; ++j) {
                if (slots[group][i] == slots[group][j] && prices[i] != prices[j])
                    throw std::runtime_error("market moved repeated occurrences of a slot separately");
            }
        }
        return group_cost(groups[group], prices, 1);
    });
}

Integer objective(const std::vector<Group>& groups, const Integer& a, const Integer& b)
{
    const GroupPrices signed_prices{a, -a, b, -b};
    Integer total = 0;
    for (const auto& group : groups) {
        GroupPrices prices;
        for (std::size_t i = 0; i < prices.size(); ++i) prices[i] = signed_prices[group.slots[i]];
        total += group_cost(group, prices, 2);
    }
    return total;
}

std::vector<Integer> check_brute(const std::vector<Group>& groups, int bound)
{
    const auto proposal = propose(2, bound, groups);
    const auto value = objective(groups, proposal[0], proposal[1]);
    for (int a = -2 * bound; a <= 2 * bound; ++a) {
        for (int b = -2 * bound; b <= 2 * bound; ++b) {
            if (objective(groups, a, b) < value)
                throw std::runtime_error("market proposal failed exhaustive half-grid comparison");
        }
    }
    return proposal;
}

}

int test_ep_market()
{
    // The two positive slots are tied at zero prices. A descent step must
    // raise both prices together to reduce the first group's utility.
    const std::vector<Group> coupled{
        {{0, 2, 0, 2}, {{{2, 2, -2, -2}, 7}}},
        {{1, 3, 1, 3}, {{{-2, -2, -2, -2}, 1}}},
    };
    if (check_brute(coupled, 2) != std::vector<Integer>{4, 4})
        throw std::runtime_error("market missed a simultaneous price move");

    check_brute({
        {{1, 3, 0, 2}, {{{0, 3, 1, -2}, 13}, {{1, -2, 2, 1}, 7}}},
        {{0, 2, 2, 3}, {{{2, -3, 3, -1}, 11}, {{-1, 3, -2, 2}, 5}}},
        {{3, 3, 1, 1}, {{{-2, 0, -1, 3}, 17}}},
    }, 3);

    if (propose(1, 3, {{{0, 1, 0, 1}, {{{3, 0, 1, -1}, 1}}}}) != std::vector<Integer>{3})
        throw std::runtime_error("market lost a half-unit price");

    const Integer bid = (Integer(1) << 192) + 1;
    const Integer weight = Integer(1) << 512;
    const std::vector<Group> precise{
        {{0, 0, 0, 0}, {{{bid, bid, bid, bid}, weight + 1}}},
        {{1, 1, 1, 1}, {{{bid, bid, bid, bid}, weight}}},
    };
    if (propose(1, bid, precise) != std::vector<Integer>{2 * bid})
        throw std::runtime_error("market lost precision in competing weights");

    return 0;
}
