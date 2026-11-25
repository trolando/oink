// SPDX-License-Identifier: Apache-2.0
#ifndef OINK_EP_SEARCH_HPP
#define OINK_EP_SEARCH_HPP

#include <boost/multiprecision/cpp_int.hpp>
#include <cstddef>
#include <functional>
#include <vector>

namespace pg::ep {

using Integer = boost::multiprecision::cpp_int;
using Prices = std::vector<Integer>;

// The trial vector is point / unit. Return one sign per coordinate, or
// stop early if the application has independently certified its solution.
using Probe = std::function<bool(const Prices& point, const Integer& unit,
                                 std::vector<int>& signs)>;
using Decode = std::function<void(const Prices& point, const Integer& unit)>;

namespace detail {
// Each node makes an odd number of probes. The root gets an additional
// root_factor * ceil(log2(|experts|+1)) probes. Separating the schedule from
// the update allows the complete recursion to be checked on finite examples.
struct ProbeSchedule {
    std::size_t depth;
    unsigned probes;
    unsigned root_factor;
    bool use_early_leaves = true;
};
void price_tree(std::size_t dimension, const Integer& radius, ProbeSchedule schedule,
                const Probe& label, const Decode& finish);
}

// Theorem 2. For each consistent signed-slot hypothesis with integer demands
// in [-radius,radius], the returned prices meet every demand of that hypothesis.
void extremal_prices(std::size_t dimension, const Integer& radius,
                     const Probe& probe, const Decode& decode);

} // namespace pg::ep

#endif
