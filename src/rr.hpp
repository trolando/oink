#ifndef RR_HPP
#define RR_HPP

#include <queue>

#include "pp.hpp"

namespace pg {

class RRSolver : public PPSolver
{
public:
    RRSolver(Oink *oink, Game *game, std::ostream &lgr);

    virtual void run();

protected:
    virtual bool checkRegion(int priority);
};

}

#endif 
