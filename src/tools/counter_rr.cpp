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
