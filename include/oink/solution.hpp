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

#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <algorithm>
#include <vector>

#include <oink/bitset.hpp>

namespace pg {

/**
 * Holds the mutable solver output for a parity game: for each vertex whether it
 * has been solved, its winner (0 for Even, 1 for Odd), and the strategy edge.
 *
 * This separates solver state from the (eventually immutable) Game structure.
 * The winner of an unsolved vertex is undefined; the strategy of a vertex is
 * the next vertex to play to, or -1 for none.
 */
class Solution
{
public:
    explicit Solution(int vertex_count = 0)
        : solved_(vertex_count), winner_(vertex_count), strategy_(vertex_count, -1)
    { }

    /**
     * Resize to <vertex_count> vertices. Existing entries are preserved; new
     * vertices are unsolved with strategy -1.
     */
    void resize(int vertex_count)
    {
        solved_.resize(vertex_count);
        winner_.resize(vertex_count);
        strategy_.resize(vertex_count, -1);
    }

    /**
     * Mark all vertices as unsolved and clear all strategies.
     */
    void reset()
    {
        solved_.reset();
        winner_.reset();
        std::fill(strategy_.begin(), strategy_.end(), -1);
    }

    [[nodiscard]] int vertex_count() const noexcept { return (int)strategy_.size(); }

    [[nodiscard]] bool is_solved(int v) const noexcept { return solved_[v]; }

    /** Raw winner bit (0 or 1); only meaningful when is_solved(v). */
    [[nodiscard]] int winner(int v) const noexcept { return winner_[v] ? 1 : 0; }

    [[nodiscard]] int strategy(int v) const noexcept { return strategy_[v]; }

    /** The set of solved vertices, for fast set operations. */
    [[nodiscard]] const bitset& solved() const noexcept { return solved_; }

    /** Direct access to the strategy array (for performance-critical solvers). */
    [[nodiscard]] const int* strategy_data() const noexcept { return strategy_.data(); }
    [[nodiscard]] int* strategy_data() noexcept { return strategy_.data(); }

    /**
     * Mark vertex <v> as solved, won by <winner> (0 or 1), with the given
     * <strategy> (the next vertex to play to, or -1 for none).
     */
    void solve(int v, int winner, int strategy) noexcept
    {
        solved_[v] = true;
        winner_[v] = winner;
        strategy_[v] = strategy;
    }

    void set_winner(int v, int winner) noexcept { winner_[v] = winner; }
    void set_strategy(int v, int strategy) noexcept { strategy_[v] = strategy; }

private:
    bitset solved_;             // set if vertex is solved
    bitset winner_;             // for solved vertices, 1 if won by Odd, else 0
    std::vector<int> strategy_; // strategy per vertex, or -1 for none
};

}

#endif
