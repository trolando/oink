#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "game.hpp"
#include "oink.hpp"
#include <signal.h>

#define LOGIC_ERROR { printf("\033[1;7mlogic error %s:%d!\033[m\n", __FILE__, __LINE__); raise(SIGABRT); }

namespace pg {

class Solver
{
public:
    Solver(Oink *oink, Game *game, std::ostream &lgr = std::cout) :
        oink(oink), game(game), logger(lgr),
        n_nodes(game->n_nodes), priority(game->priority), owner(game->owner),
        out(game->out), in(game->in), disabled(game->disabled),
        outa(oink->outa), ina(oink->ina), outs(oink->outs), ins(oink->ins) { }
    virtual ~Solver() { }
    virtual void run() = 0;

    void setTrace(int value) { trace = value; }

protected:
    Oink *oink;
    Game *game;
    std::ostream &logger;
    int trace = 0;

    const int n_nodes;
    const int * const priority;
    const int * const owner;
    const std::vector<int> * const out;
    const std::vector<int> * const in;
    const int * const disabled;

    const int* const outa;
    const int* const ina;
    const int* const outs;
    const int* const ins;
};

}

#endif 
