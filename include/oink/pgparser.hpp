/*
 * Copyright 2017-2024 Tom van Dijk
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

#ifndef PGPARSER_HPP
#define PGPARSER_HPP

#include <sstream>
#include <memory>
#include <oink/game.hpp>

namespace pg {

/**
 * Parse a .pg file
 */
class PGParser
{
public:
    PGParser() = delete;

    /**
     * Parse a pgsolver game.
     */
    static Game parse_pgsolver(std::istream &in, bool removeBadLoops);

    /**
     * Parse a pgsolver game. This version can handle priorities above INT_MAX, such as the mlsolver benchmarks.
     * It does this by renumbering (compressing) the priorities afterwards.
     */
    static Game parse_pgsolver_renumber(std::istream &in, bool removeBadLoops);
};

}

#endif 
