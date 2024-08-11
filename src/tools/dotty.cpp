/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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

#include <iostream>
#include <fstream>

#include "oink/game.hpp"
#include "oink/pgparser.hpp"

using namespace pg;

int
main(int argc, char **argv)
{
    Game game;
    if (argc >= 2) {
        std::ifstream file(argv[1]);
        game = PGParser::parse_pgsolver(file, false);
        file.close();
    } else {
        game = PGParser::parse_pgsolver(std::cin, false);
    }

    if (argc >= 3) {
        std::ofstream file(argv[2]);
        game.write_dot(file);
        file.close();
    } else {
        game.write_dot(std::cout);
    }
    return 0;
}
