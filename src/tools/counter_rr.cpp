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

int
main(int argc, char** argv)
{
    if (argc != 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    int n = std::stoi(argv[1]);

    Game game(4*n);

    // generate the Nth game
    for (int i=0; i<n; i++) {
        game.initNode(4*i,   n-1-i, (n+i)&1);
        game.initNode(4*i+1, n-1-i, (n+i)&1);
        game.initNode(4*i+2, n+1+i, (n+i)&1);
        game.initNode(4*i+3, n+1+i, (n+i)&1);
        // add edges within the group
        game.addEdge(4*i,   4*i+1);
        game.addEdge(4*i+1, 4*i);
        game.addEdge(4*i+2, 4*i+3);
        game.addEdge(4*i+3, 4*i+2);
        game.addEdge(4*i,   4*i+2);
        game.addEdge(4*i+1, 4*i+3);
        // add edges to the next group
        if (i == (n-1)) {
            game.addEdge(4*i+2, 0);
            game.addEdge(4*i+3, 1);
        } else {
            game.addEdge(4*i+2, 4*i+4);
            game.addEdge(4*i+3, 4*i+5);
        }
    }

    game.write_pgsolver(std::cout);
}
