// SPDX-License-Identifier: Apache-2.0

#ifndef OINK_EP_EXPERTS_HPP
#define OINK_EP_EXPERTS_HPP

#include <cstddef>
#include <memory>
#include <vector>

#include "ep_market.hpp"

namespace pg::ep {

using Prices = std::vector<Integer>;
using History = std::vector<Prices>;

class ExpertPool
{
public:
    ExpertPool(std::size_t coordinates, const Integer& radius, std::size_t max_depth,
               const Integer& initial_wealth = 1);
    ~ExpertPool();

    // Counts the ordered experts, including the currently inactive pivot choices.
    static Integer count(std::size_t coordinates, const Integer& radius, std::size_t max_depth);

    // Even slots are positive, and odd slots are negative. Prices and histories
    // use integer numerators over scale. All groups use a common wealth scale.
    const std::vector<GroupSlots>& group_slots() const;
    Integer group_objective(std::size_t group, const GroupPrices& prices,
                            const Integer& scale, const History& history) const;
    Integer objective(const Prices& prices, const Integer& scale, const History& history) const;
    Prices minimizing_prices(const Integer& bid_bound, const Integer& scale, const History& history) const;

    // Opening a probe keeps its parent share and leaves the child share in the current pool.
    // Closing the probe combines the child's wealth with the parent's updated share.
    void open_probe(const Prices& prices, const Integer& scale, const History& history);
    void close_probe(const std::vector<int>& feedback);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}

#endif
