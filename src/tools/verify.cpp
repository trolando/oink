#include <iostream>
#include <fstream>
#include <sys/time.h>

#include "game.hpp"
#include "verifier.hpp"

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
        Game pg;

        std::ifstream inp(argv[1]);
        pg.parse_pgsolver(inp);
        inp.close();
        std::cout << "game loaded." << std::endl;

        std::ifstream inpsol(argv[2]);
        pg.parse_solution(inpsol);
        inpsol.close();
        std::cout << "solution loaded." << std::endl;

        pg.reindex();

        Verifier v(&pg, std::cout);
        auto begin = wctime();
        v.verify(true);
        auto end = wctime();
        std::cout << "verified in " << (end - begin) << " sec." << std::endl;
        
        return 0;
    } catch (const char *err) {
        std::cerr << err << std::endl;
        return -1;
    }
}
