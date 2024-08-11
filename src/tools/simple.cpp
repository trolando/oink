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

#include <fstream>

#include "oink/oink.hpp"
#include "oink/pgparser.hpp"

int
main(int argc, char** argv)
{
    /**
     * Read game from file or stdin
     */
    pg::Game pg;
    if (argc > 1) {
        std::ifstream file(argv[1], std::ios_base::binary);
        pg = pg::PGParser::parse_pgsolver(file, true);
        file.close();
    } else {
        pg = pg::PGParser::parse_pgsolver(std::cin, true);
    }

    /**
     * Solve with TL solver and no WCWC / SP solving
     */
    std::stringstream log;
    pg::Oink solver(pg, log);
    solver.setSolveSingle(false);
    solver.setRemoveWCWC(false);
    solver.setSolver("tl");
    solver.run();

    /**
     * Write solution to file or stdout
     */
    if (argc > 2) {
        std::ofstream file(argv[2]);
        pg.write_sol(file);
        file.close();
    } else {
        pg.write_sol(std::cout);
    }

    return 0;
}
