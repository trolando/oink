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

// The "robust" SCC variation of the core family by Benerecetti et al
int
main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    int n = 2*std::stoi(argv[1]);

    // core: 3*(n+1)
    // extension: 3*(n*n+(n&1))/4 + n
    int n_positions = 3*(n+1) + n + (3*n*n+(n&1))/4;

    Game game(n_positions);
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

    /* create more connectors */
    int nxt = 3*n+3;
    for (int i=0; i<=n; i++) {
        for (int j=i+1; j<=n; j++) {
            int i_parity = 1-(i&1); // parity of c_i vertex
            int j_parity = 1-(j&1); // parity of c_j vertex
            int ci = 3*i+2;
            int cj = 3*j+2;
            if (i_parity == j_parity) {
                // same parity, only need 1 
                game.init_vertex(nxt, 0, 1-i_parity);
                game.vec_add_edge(ci, nxt);
                game.vec_add_edge(cj, nxt);
                game.vec_add_edge(nxt, ci);
                game.vec_add_edge(nxt, cj);
                nxt++;
            } else {
                // different parity, so we need two
                game.init_vertex(nxt,   0, 1-i_parity);
                game.init_vertex(nxt+1, 0, 1-j_parity);
                game.vec_add_edge(ci, nxt);
                game.vec_add_edge(nxt, ci);
                game.vec_add_edge(nxt, nxt+1);
                game.vec_add_edge(nxt+1, nxt);
                game.vec_add_edge(cj, nxt+1);
                game.vec_add_edge(nxt+1, cj);
                nxt+=2;
            }
        }        
    }

    game.vec_finish();
    game.write_pgsolver(std::cout);
}
