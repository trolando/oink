#include <iostream>
#include <fstream>

#include "game.hpp"

using namespace pg;

int
main(int argc, char **argv)
{
    Game game;
    if (argc >= 2) {
        std::ifstream file(argv[1]);
        game.parse_pgsolver(file);
        file.close();
    } else {
        game.parse_pgsolver(std::cin);
    }
    if (argc >= 3) {
        std::ofstream file(argv[2]);
        game.write_dot(file);
        file.close();
    } else {
        game.write_dot(std::cout);
    }
    return 0;
}
