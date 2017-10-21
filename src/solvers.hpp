#include <iostream>

#ifndef SOLVERS_HPP
#define SOLVERS_HPP

namespace pg {

enum Solvers
{
    NONE = 0,
    ZLK,
    PP,
    PPP,
    RR,
    DP,
    RRDP,
    PSI,
    SPM,
    TSPM,
    MSPM,
    QPT,
    TL
};

class Oink;
class Game;
class Solver;

Solvers solverToId(std::string s);

std::string solverToString(Solvers s);

bool solverIsParallel(Solvers s);

void listSolvers(std::ostream &out);

Solver *constructSolver(Solvers s, Oink *oink, Game *game, std::ostream &lgr);

}

#endif 
