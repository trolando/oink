/*
 * Copyright 2026 Tom van Dijk
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef OINK_SUCCESSOR_CHECK_HPP
#define OINK_SUCCESSOR_CHECK_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace pg::detail {

constexpr int small_successor_list = 8;

inline int
repeated_short_successor(const int* successors, int count)
{
    for (int i = 1; i < count; i++) {
        for (int j = 0; j < i; j++) {
            if (successors[i] == successors[j]) return successors[i];
        }
    }
    return -1;
}

// Check each source once per stamp array. Targets must be valid vertex IDs.
// Short lists stay in cache and need no storage indexed by vertex ID.
inline int
repeated_successor(const int* successors, int count, int source, std::size_t vertices,
                   std::vector<int>& stamps)
{
    if (count <= small_successor_list) return repeated_short_successor(successors, count);
    if (stamps.empty()) stamps.assign(vertices, -1);
    for (int i = 0; i < count; i++) {
        const int target = successors[i];
        if (stamps[target] == source) return target;
        stamps[target] = source;
    }
    return -1;
}

// Pack 32 targets into each word to keep the fast parser's stamp table small.
// The high bits name the source and the low bits mark its targets.
inline int
repeated_successor(const int* successors, int count, int source, std::size_t vertices,
                   std::vector<uint64_t>& stamps)
{
    if (count <= small_successor_list) return repeated_short_successor(successors, count);
    if (stamps.empty()) stamps.resize((vertices + 31) / 32);
    const uint64_t stamp = uint64_t(source) << 32;
    for (int i = 0; i < count; i++) {
        const int target = successors[i];
        const uint64_t bit = uint64_t(1) << (target & 31);
        uint64_t& word = stamps[target >> 5];
        const uint64_t seen = (word & 0xffffffff00000000ULL) == stamp ? word : stamp;
        if (seen & bit) return target;
        word = seen | bit;
    }
    return -1;
}

inline std::string
repeated_edge_message(int source, int target)
{
    return "repeated edge from vertex " + std::to_string(source) + " to successor " + std::to_string(target);
}

}

#endif
