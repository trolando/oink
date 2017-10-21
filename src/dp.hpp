#ifndef DP_HPP
#define DP_HPP

#include <queue>
#include <set>

#include "ppp.hpp"

namespace pg {

class DPSolver : public PPPSolver
{
public:
    DPSolver(Oink *oink, Game *game, std::ostream &logger);

    virtual void run();

    int delayed = 0;

protected:
    int *region_;

    virtual int getRegionStatus(int index, int priority);
};

}

#endif 
