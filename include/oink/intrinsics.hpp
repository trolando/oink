/*
 * Copyright 2024 Tom van Dijk
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

#ifndef INTRINSICS_HPP
#define INTRINSICS_HPP

#include <cstdint>

/**
 * Portability wrappers for the bit-manipulation intrinsics used in hot loops.
 *
 * These are kept inline (header-only) on purpose: they each compile to a single
 * machine instruction and are called from performance-critical bitset code, so
 * out-of-lining them into a .cpp would prevent inlining. Non-GCC/Clang compilers
 * fall back to portable implementations.
 */

namespace pg {
namespace intrinsics {

/** Number of set bits in x. */
inline int popcount64(uint64_t x) noexcept
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int c = 0;
    while (x) { x &= x - 1; c++; }
    return c;
#endif
}

/** Number of trailing zero bits (index of the lowest set bit). Undefined for x == 0. */
inline int countr_zero64(uint64_t x) noexcept
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    while ((x & 1) == 0) { x >>= 1; n++; }
    return n;
#endif
}

/** Number of leading zero bits. Undefined for x == 0. */
inline int countl_zero64(uint64_t x) noexcept
{
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#else
    int n = 0;
    for (int i = 63; i >= 0; i--) {
        if (x & (uint64_t(1) << i)) break;
        n++;
    }
    return n;
#endif
}

}
}

#endif
