#pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace pg::detail {

// Both edge arrays use int offsets and contain one terminator per vertex.
constexpr std::size_t max_edge_slots = INT_MAX;

inline std::size_t
checked_edge_capacity(std::size_t capacity, std::size_t grown)
{
    if (capacity >= max_edge_slots) throw std::runtime_error("too many edges for int edge offsets");
    return std::min(grown, max_edge_slots);
}

inline int
pgsolver_vertex_capacity(uint64_t header)
{
    // The header can be the count or the highest ID. A valid game needs
    // at least one successor and a terminator per vertex.
    if (header > max_edge_slots / 2) throw std::runtime_error("number of nodes too high");
    return static_cast<int>(header) + 1;
}

}
