#ifndef APT_HPP
#define APT_HPP

#include <queue>

#include "solver.hpp"

namespace pg {

class APTSolver : public Solver
{
public:
    APTSolver(Oink *oink, Game *game, std::ostream &lgr);
    virtual ~APTSolver();

    virtual void run();
};

}

#endif 
