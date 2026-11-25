// SPDX-License-Identifier: Apache-2.0

#ifndef OINK_EP_MARKET_HPP
#define OINK_EP_MARKET_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

#include <boost/multiprecision/cpp_int.hpp>

namespace pg::ep {

using Integer = boost::multiprecision::cpp_int;
using GroupSlots = std::array<std::size_t, 4>;
using GroupPrices = std::array<Integer, 4>;
using GroupCost = std::function<Integer(std::size_t, const GroupPrices&)>;

// Slot 2*i has price y[i], and slot 2*i+1 has price -y[i]. The callback sums
// the experts' utilities in one group with nonnegative weights. Each utility
// is the maximum of zero and the four differences between bids and prices.
// Bids are integral in the callback's units, and bid_bound bounds their
// absolute values. The returned coordinates are in half those units.
std::vector<Integer> market_prices(std::size_t coordinates, const std::vector<GroupSlots>& groups,
                                   const Integer& bid_bound, const GroupCost& group_cost);

}

#endif
