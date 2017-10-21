#ifndef RRDP_HPP
#define RRDP_HPP

#include <queue>
#include <set>

#include "rr.hpp"

namespace pg {

class RRDPSolver : public RRSolver
{
public:
    RRDPSolver(Oink *oink, Game *game, std::ostream &lgr);

    virtual void run();

    int delayed = 0;

protected:
    int *region_;

    int getRegionStatus(int index, int priority);
};

}

#endif 
