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

    Game game(3*n+3);

    /* create 2n+1 pieces */
    for (int i=0; i<=n; i++) {
        game.initNode(3*i+0, i+2, (i&1));
        game.initNode(3*i+1, 1-(i&1), (i&1));
        game.initNode(3*i+2, 1-(i&1), 1-(i&1));
        game.addEdge(3*i+0, 3*i+1);
        game.addEdge(3*i+1, 3*i+2);
        game.addEdge(3*i+2, 3*i+1);
    }
    
    /* connect the pieces */
    for (int i=0; i<n; i++) {
        game.addEdge(3*i+0, 3*i+3);
        game.addEdge(3*i+1, 3*i+3);
        game.addEdge(3*i+5, 3*i+2);
    }

    game.write_pgsolver(std::cout);
}
