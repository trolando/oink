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
#include <cstddef>
#include <vector>

#include <oink/bitset.hpp>

namespace pg {

class Game;

/**
 * Holds the mutable solver output for a parity game: for each vertex whether it
 * has been solved, its winner (0 for Even, 1 for Odd), and the strategy.
 *
 * A solution natively represents both single and multi strategies. The single
 * strategy (one move per vertex) is primary. A vertex won by its owner with
 * single strategy -1 (the "sentinel") instead takes its winning moves from the
 * multi-strategy: a set of outgoing edges (the maximal permissive strategy),
 * stored as a bitset indexed by position in the game's outgoing edge array
 * (out-edge k of v is at index firstout(v)+k). A single solver never leaves a
 * won vertex at -1, so -1 with winner==owner unambiguously means multi.
 * strategyTargets() applies this dispatch so both kinds read uniformly.
 *
 * Interpreting the multi-strategy needs the game's edge layout, so a Solution
 * that carries one references its Game (which must outlive it); the index-based
 * edge operations used in solver hot loops need no Game.
 */
class Solution
{
public:
    explicit Solution(int vertex_count = 0)
        : solved_(vertex_count), winner_(vertex_count), strategy_(vertex_count, -1)
    { }

    /**
     * Resize to <vertex_count> vertices. Existing entries are preserved; new
     * vertices are unsolved with strategy -1. (Does not touch the multi-strategy,
     * which is sized to the edge array; see init_multi.)
     */
    void resize(int vertex_count)
    {
        solved_.resize(vertex_count);
        winner_.resize(vertex_count);
        strategy_.resize(vertex_count, -1);
    }

    /**
     * Mark all vertices as unsolved and clear all strategies (single and multi).
     */
    void reset()
    {
        solved_.reset();
        winner_.reset();
        std::fill(strategy_.begin(), strategy_.end(), -1);
        edges_.reset();
        has_multi_ = false;
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
     * Mark vertex <v> as solved, won by <winner> (0 or 1), with the given single
     * <strategy> (the next vertex to play to, or -1 for none / multi sentinel).
     */
    void solve(int v, int winner, int strategy) noexcept
    {
        solved_[v] = true;
        winner_[v] = winner;
        strategy_[v] = strategy;
    }

    void set_winner(int v, int winner) noexcept { winner_[v] = winner; }
    void set_strategy(int v, int strategy) noexcept { strategy_[v] = strategy; }

    /** Swap the solution state of two vertices (used when permuting vertices). */
    void swap_vertices(int a, int b) noexcept
    {
        // The multi-strategy is indexed by edge-array position (not by vertex) and
        // is permutation-invariant, so only the per-vertex state is swapped here.
        bool sa = solved_[a]; solved_[a] = (bool)solved_[b]; solved_[b] = sa;
        bool wa = winner_[a]; winner_[a] = (bool)winner_[b]; winner_[b] = wa;
        std::swap(strategy_[a], strategy_[b]);
    }

    /**
     * Apply a vertex permutation (same convention as Game::permute) so the
     * solution stays consistent with a permuted game. Only the per-vertex state
     * and the single strategy targets are remapped.
     *
     * IMPORTANT, two preconditions:
     *  - <mapping> is NOT modified here, but Game::permute() destroys its mapping
     *    in place (it is a cycle-sort). So call solution.permute(mapping) BEFORE
     *    game.permute(mapping) (or pass each its own copy).
     *  - The multi-strategy (edges_) is indexed by edge-array position and is
     *    intentionally left untouched: it is correct only because the game is
     *    permuted with the SAME mapping right after (the bits are permutation-
     *    invariant — blocks don't move; firstout pointers and edge targets move
     *    with the vertices). A solution-only renumber (without permuting its
     *    game) leaves the multi / strategy_targets() path reading the OLD edge
     *    layout while the single strategy is in the new numbering — inconsistent.
     *    The single strategy is self-contained; the multi strategy is only
     *    meaningful relative to its game's current edge layout.
     */
    void permute(const int* mapping)
    {
        const int n = (int)strategy_.size();
        for (int i=0; i<n; i++) if (strategy_[i] != -1) strategy_[i] = mapping[strategy_[i]];
        std::vector<int> m(mapping, mapping + n);
        for (int i=0; i<n; i++) {
            while (m[i] != i) {
                int k = m[i]; m[i] = m[k]; m[k] = k;
                swap_vertices(i, k);
            }
        }
    }

    /** Cheap swap (pointer swaps only). */
    void swap(Solution& other) noexcept
    {
        solved_.swap(other.solved_);
        winner_.swap(other.winner_);
        strategy_.swap(other.strategy_);
        edges_.swap(other.edges_);
        std::swap(has_multi_, other.has_multi_);
    }

    /* --- multi-strategy --- */

    /**
     * Whether this solution carries a multi-strategy.
     */
    [[nodiscard]] bool has_multi() const noexcept { return has_multi_; }

    /**
     * Allocate (sized to the game's edge array, see Game::edgeArraySize) and clear
     * the multi-strategy.
     */
    void init_multi(std::size_t edge_array_size)
    {
        edges_.resize(edge_array_size);
        edges_.reset();
        has_multi_ = true;
    }

    /**
     * Edge operations by out-edge array index (idx = firstout(v)+k). These need no
     * game and are used in the solver hot loops. The edge-indexed multi-strategy
     * is meaningless without the game that defines the edge layout; the game-aware
     * reads are the free functions strategy_targets()/first_strategy_edge()/
     * has_strategy_edge_to() below (also reachable via Oink).
     */
    void add_edge_index(std::size_t idx) { edges_.set(idx); }
    void remove_edge_index(std::size_t idx) { edges_.reset(idx); }
    [[nodiscard]] bool has_edge_index(std::size_t idx) const { return edges_.test(idx); }
    void clear_edge_range(std::size_t start, std::size_t count)
    {
        for (std::size_t i = 0; i < count; i++) edges_.reset(start + i);
    }

private:
    bitset solved_;             // set if vertex is solved
    bitset winner_;             // for solved vertices, 1 if won by Odd, else 0
    std::vector<int> strategy_; // single strategy per vertex, or -1 (none / multi sentinel)

    bitset edges_;              // multi-strategy: set of strategy edges (empty unless has_multi_)
    bool has_multi_ = false;    // whether a multi-strategy is recorded
};

/* --- game-aware reads of a solution's strategy ---
 *
 * The multi-strategy is indexed by the game's outgoing edge array, so reading it
 * needs the game. These are free functions (over a game and a solution) rather
 * than members, so Solution stays a plain value type; they are also reachable
 * via Oink (oink.strategyTargets(v)) for callers that hold the solver.
 */

/**
 * Append the winning strategy moves of <v> to <out>: the single strategy if one
 * is set, otherwise (a vertex won by its owner with the -1 sentinel) the recorded
 * multi-strategy edges. A losing or unsolved vertex contributes nothing. The
 * uniform way to read a strategy, single or multi.
 */
void strategy_targets(const Game& game, const Solution& solution, int v, std::vector<int>& out);

/** First recorded strategy target of <v>, or -1 if none. */
[[nodiscard]] int first_strategy_edge(const Game& game, const Solution& solution, int v);

/** Whether the edge <v> -> <to> is a recorded strategy edge of <solution>. */
[[nodiscard]] bool has_strategy_edge_to(const Game& game, const Solution& solution, int v, int to);

/**
 * Permute a game and its solution together, consistently (the safe way to
 * renumber a solution). Equivalent to solution.permute(mapping) followed by
 * game.permute(mapping), but guarantees the required order: Game::permute()
 * destroys <mapping> in place, so the solution must be permuted first. Because
 * the multi-strategy is only meaningful relative to its game's edge layout
 * (see Solution::permute), prefer this over permuting either alone.
 */
void permute(Game& game, Solution& solution, int* mapping);

}

#endif
