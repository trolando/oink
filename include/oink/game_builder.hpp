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

#ifndef GAME_BUILDER_HPP
#define GAME_BUILDER_HPP

#include <string>
#include <vector>

#include <oink/bitset.hpp>
#include <oink/game.hpp>

namespace pg {

/**
 * Mutable construction API for parity games.
 *
 * A GameBuilder accumulates vertices (priority, owner, label) and edges, and
 * produces an immutable Game with build(). It is meant for parsing and
 * test/tool construction, not for solver hot loops, so it favours clarity over
 * raw performance.
 *
 * Vertices are numbered 0..vertex_count-1. Owner is 0 for Even, 1 for Odd.
 */
class GameBuilder
{
public:
    explicit GameBuilder(int vertex_count);

    void set_priority(int v, int priority);
    void set_owner(int v, int owner);
    void set_label(int v, std::string label);
    void add_edge(int from, int to);

    /**
     * Reduce the vertex count, dropping trailing vertices (and their data).
     * Used when the final vertex count is only known after parsing.
     */
    void truncate(int new_count);

    [[nodiscard]] int vertex_count() const noexcept { return vertex_count_; }

    /**
     * Construct a Game from the accumulated vertices and edges.
     * The builder remains usable afterwards.
     */
    Game build();

private:
    int vertex_count_;
    size_t edge_count_ = 0;
    std::vector<int> priority_;
    bitset owner_;
    std::vector<std::string> labels_;
    std::vector<std::vector<int>> outgoing_;
};

}

#endif
