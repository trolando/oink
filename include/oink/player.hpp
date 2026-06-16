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

#ifndef PLAYER_HPP
#define PLAYER_HPP

#include <cstdint>

namespace pg {

/**
 * The two players of a parity game. The underlying values match the historical
 * integer convention (0 for Even, 1 for Odd) so conversion at the boundary with
 * int-based solver code is cheap.
 */
enum class Player : uint8_t {
    Even = 0,
    Odd = 1
};

/** The opposing player. */
inline Player opponent(Player p) noexcept
{
    return p == Player::Even ? Player::Odd : Player::Even;
}

/** The integer index of a player (0 for Even, 1 for Odd). */
inline int player_index(Player p) noexcept
{
    return static_cast<int>(p);
}

/** Player from an integer (0 is Even, anything else is Odd). */
inline Player player_from_int(int v) noexcept
{
    return v ? Player::Odd : Player::Even;
}

}

#endif
