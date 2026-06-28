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

#ifndef MULTISTRATEGY_HPP
#define MULTISTRATEGY_HPP

#include <cstddef>

#include <oink/bitset.hpp>

namespace pg {

/**
 * An alternative to the single strategy in Solution: instead of one strategy
 * edge per vertex, MultiStrategy records a *set* of winning strategy edges per
 * vertex (the maximal permissive strategy).
 *
 * It is backed by a single bitset indexed by position in the game's outgoing
 * edge array. A game stores its edges consecutively in one array (see Game),
 * where the outgoing edges of a vertex <v> begin at firstout(v) and out-edge
 * <k> of <v> lives at absolute index firstout(v)+k. The bit at that index is
 * set iff out-edge <k> of <v> is a winning strategy move for <v>.
 *
 * This deliberately stores only the strategy edges; whether a vertex is solved
 * and who won it remain in the game's Solution. Indexing by edge-array position
 * keeps marking/testing an edge O(1) and makes the representation invariant
 * under vertex permutation (the edge blocks do not move; only firstout pointers
 * do).
 */
class MultiStrategy
{
public:
    explicit MultiStrategy(size_t edge_array_size = 0)
        : edges_(edge_array_size)
    { }

    /**
     * Resize to hold an edge array of <edge_array_size> entries.
     * Any new bits are zero.
     */
    void resize(size_t edge_array_size)
    {
        edges_.resize(edge_array_size);
    }

    /**
     * Clear all strategy edges.
     */
    void reset()
    {
        edges_.reset();
    }

    /**
     * Whether the multi-strategy has no storage (no edges representable).
     */
    [[nodiscard]] bool empty() const noexcept { return edges_.size() == 0; }

    /**
     * The number of edge slots (i.e. the size of the game's edge array).
     */
    [[nodiscard]] size_t size() const noexcept { return edges_.size(); }

    /**
     * Mark/unmark/test out-edge at absolute array index <idx> (= firstout(v)+k)
     * as a winning strategy move.
     */
    void add(size_t idx) { edges_.set(idx); }
    void remove(size_t idx) { edges_.reset(idx); }
    [[nodiscard]] bool has(size_t idx) const { return edges_.test(idx); }

    /**
     * Clear all strategy edges in the contiguous block [start, start+count),
     * i.e. one vertex's outgoing edge block, before recomputing it.
     */
    void clear_range(size_t start, size_t count)
    {
        for (size_t i = 0; i < count; i++) edges_.reset(start + i);
    }

    /**
     * Direct access to the underlying bitset, for cheap copy/swap.
     */
    [[nodiscard]] const bitset& bits() const noexcept { return edges_; }
    [[nodiscard]] bitset& bits() noexcept { return edges_; }

private:
    bitset edges_; // bit firstout(v)+k set => out-edge k of v is a strategy move
};

}

#endif
