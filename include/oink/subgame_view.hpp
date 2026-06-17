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

#ifndef SUBGAME_VIEW_HPP
#define SUBGAME_VIEW_HPP

#include <cstddef>

#include <oink/bitset.hpp>

namespace pg {

class Game;

/**
 * A non-owning view of a parity game restricted to a set of active vertices.
 *
 * Unlike Game::extract_subgame, it does not copy or reindex the game: it simply
 * pairs a game with the bitset of vertices currently in play, for read-only
 * subgame operations. It is a view and must not outlive the game or the bitset.
 */
class SubgameView
{
public:
    SubgameView(const Game& game, const bitset& active) noexcept
        : game_(&game), active_(&active) { }

    [[nodiscard]] const Game& game() const noexcept { return *game_; }
    [[nodiscard]] const bitset& active_vertices() const noexcept { return *active_; }

    [[nodiscard]] bool contains(int v) const noexcept { return (*active_)[v]; }
    [[nodiscard]] std::size_t size() const noexcept { return active_->count(); }

private:
    const Game* game_;
    const bitset* active_;
};

}

#endif
