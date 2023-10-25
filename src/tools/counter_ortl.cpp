/*
 * Copyright 2021 Tom van Dijk, Johannes Kepler University Linz
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

#include "oink/game.hpp"

using namespace pg;

int
main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    const int n = std::stoi(argv[1]);

    Game game(1 + 5*n);
    game.vec_init();

    game.init_vertex(0, 0, 1);              // root of the distracted tangle
    game.vec_add_edge(0, 0);                // self-loop

    for (int i=0; i<n; i++) {
        int c = 5*i+1;
        int d = 2*(n+i+1);
        game.init_vertex(c+0, 2*(n-i), 0);  // distracted vertex
        game.init_vertex(c+1, d, 0);        // distraction
        game.init_vertex(c+2, 1, 1);        // opponent tangle start
        game.init_vertex(c+3, 1, 0);        // opponent tangle end
        game.init_vertex(c+4, d+1, 0);      // attracting odd vertex
        game.vec_add_edge(c+0, 0);          // from distracted vertex to root
        game.vec_add_edge(0, c+0);          // from root to distracted vertex
        game.vec_add_edge(c+0, c+1);        // from distracted vertex to distraction
        game.vec_add_edge(c+1, c+2);        // from distraction to tangle
        game.vec_add_edge(c+2, c+3);        // tangle forward edge
        game.vec_add_edge(c+3, c+2);        // tangle backward edge
        game.vec_add_edge(c+3, c+4);        // tangle to attracting odd vertex
        game.vec_add_edge(c+4, c);          // attracting odd vertex to distracted tangle
    }

    game.vec_finish();
    game.sort();
    game.renumber();
    game.write_pgsolver(std::cout);
}
