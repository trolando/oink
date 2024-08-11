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
#include <sys/time.h>

#include "oink/game.hpp"
#include "verifier.hpp"
#include "oink/pgparser.hpp"

using namespace std;
using namespace pg;

static double
wctime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + 1E-6 * time.tv_usec;
}

int
main(int argc, const char **argv)
{
    if (argc != 3) {
        std::cout << "Syntax: " << argv[0] << " pg_file sol_file" << std::endl;
        return -1;
    }

    try {
        std::ifstream inp(argv[1]);
        Game pg = PGParser::parse_pgsolver_renumber(inp, false);
        inp.close();
        std::cout << "game loaded." << std::endl;

        std::ifstream inpsol(argv[2]);
        pg.parse_solution(inpsol);
        inpsol.close();
        std::cout << "solution loaded." << std::endl;

        pg.sort();

        Verifier v(pg, std::cout);
        auto begin = wctime();
        v.verify(true);
        auto end = wctime();
        std::cout << "verified in " << (end - begin) << " sec." << std::endl;
        
        return 0;
    } catch (std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
        return -1;
    }
}
