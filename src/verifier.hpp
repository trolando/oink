#ifndef VERIFIER_HPP
#define VERIFIER_HPP

#include <istream>
#include "game.hpp"

namespace pg {

class Verifier
{
public:
    Verifier(Game* pg, std::ostream &lgr) : game(pg), logger(lgr) { }

    void verify(bool fullgame=true, bool even=true, bool odd=true);

    int n_strategies = 0;

protected:
    Game *game;
    std::ostream &logger;
};

}

#endif 
