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

#include "game.hpp"

using namespace pg;

// The core family of difficult parity games by Benerecetti et al
int
main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    int n = 2*std::stoi(argv[1]);

    Game game(3+3*n);
    game.vec_init();

    /* create n+1 pieces */
    for (int i=0; i<=n; i++) {
        game.init_vertex(3*i+0, n+1+i, (i&1));
        game.init_vertex(3*i+1, i, (i&1));
        game.init_vertex(3*i+2, i, 1-(i&1));
        game.vec_add_edge(3*i+0, 3*i+1);
        game.vec_add_edge(3*i+1, 3*i+2);
        game.vec_add_edge(3*i+2, 3*i+1);
        game.vec_add_edge(3*i+2, 3*i+2);
    }
    
    /* connect the pieces */
    for (int i=0; i<n; i++) {
        game.vec_add_edge(3*i+2, 3*i+3);
        game.vec_add_edge(3*i+4, 3*i+0);
    }

    game.vec_finish();
    game.write_pgsolver(std::cout);
}
