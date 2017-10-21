#ifndef PPP_HPPP
#define PPP_HPPP

#include <queue>

#include "pp.hpp"

namespace pg {

class PPPSolver : public PPSolver
{
public:
    PPPSolver(Oink *oink, Game *game, std::ostream &logger);

    virtual void run();

protected:
    int reset0, reset1;
};

}

#endif 
