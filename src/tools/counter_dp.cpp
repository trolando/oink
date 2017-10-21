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

    Game game(4+n*4);

    /* create n+1 pieces */
    for (int i=0; i<=n; i++) {
        game.initNode(4*i+0, i, 1-(i&1));
        game.initNode(4*i+1, i, 1-(i&1));
        game.initNode(4*i+2, i, (i&1));
        game.initNode(4*i+3, i+3, (i&1));
        game.addEdge(4*i+0, 4*i+1);
        game.addEdge(4*i+1, 4*i+0);
        game.addEdge(4*i+1, 4*i+2);
        game.addEdge(4*i+2, 4*i+1);
        game.addEdge(4*i+3, 4*i+2);
    }
    
    /* connect the pieces */
    for (int i=0; i<n; i++) {
        game.addEdge(4*i+6, 4*i+3);
        game.addEdge(4*i+1, 4*i+7);
    }

    game.write_pgsolver(std::cout);
}
