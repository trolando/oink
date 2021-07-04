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

#include "game.hpp"

using namespace pg;

int
main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    const int n = std::stoi(argv[1]);  // index of game
    const int piv = 2+2*n;             // pivot priority

    Game game(8 + 4*n);

    game.vec_init();

    game.init_vertex(0, piv+6, 1);  // top vertex of dominion
    game.init_vertex(1, 0, 1);      // first [distracted] tangle
    game.init_vertex(2, piv+5, 1);  // barrier
    game.init_vertex(3, 2, 1);      // second [distracted] tangle
    game.init_vertex(4, 2, 1);      // second [distracted] tangle "self-loop" assist
    game.vec_add_edge(0, 3);        // top vertex -> second tangle
    game.vec_add_edge(1, 0);        // first tangle -> top vertex
    game.vec_add_edge(2, 1);        // barrier -> first tangle
    game.vec_add_edge(3, 2);        // second tangle -> barrier
    game.vec_add_edge(3, 4);        // second tangle -> self-loop assist
    game.vec_add_edge(4, 3);        // self-loop assist -> second tangle

    game.init_vertex(5, piv+2, 0);  // extra distracted vertex connected to first tangle
    game.init_vertex(6, piv+4, 0);  // distraction
    game.init_vertex(7, piv+5, 0);  // dominant
    game.vec_add_edge(1, 5);        // first tangle -> distracted vertex
    game.vec_add_edge(5, 1);        // distracted vertex -> first tangle
    game.vec_add_edge(5, 6);        // distracted vertex -> distraction
    game.vec_add_edge(6, 7);        // distraction -> dominant
    game.vec_add_edge(7, 5);        // dominant -> distracted vertex

    for (int i=0; i<n; i++) {
        int c = 8+4*i;
        int d = piv+8+2*i;
        int e = piv-2*i;
        game.init_vertex(c+0, e, 0);    // distracted vertex for first tangle
        game.init_vertex(c+1, e, 0);    // distracted vertex for second tangle
        game.init_vertex(c+2, d, 0);    // distraction
        game.init_vertex(c+3, d+1, 0);  // dominant

        game.vec_add_edge(1, c+0);      // first tangle -> distracted vertex 1
        game.vec_add_edge(c+0, 1);      // distracted vertex 1 -> first tangle
        game.vec_add_edge(3, c+1);      // second tangle -> distracted vertex 2
        game.vec_add_edge(c+1, 3);      // distracted vertex 2 -> second tangle
        game.vec_add_edge(c+0, c+2);    // distracted vertex 1 -> distraction
        game.vec_add_edge(c+1, c+2);    // distracted vertex 2 -> distraction
        game.vec_add_edge(c+2, c+3);    // distraction -> dominant
        game.vec_add_edge(c+3, c+1);    // dominant -> distracted vertex 2
    }

    game.vec_finish();

    //game.sort();  // optional sort vertices
    //game.renumber();  // optional renumber vertices for tight priorities

    game.write_pgsolver(std::cout);
}
