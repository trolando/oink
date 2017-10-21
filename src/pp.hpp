#ifndef PP_HPP
#define PP_HPP

#include <queue>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class PPSolver : public Solver
{
public:
    PPSolver(Oink *oink, Game *game, std::ostream &lgr = std::cout);
    virtual ~PPSolver();
    virtual void run();

    int promotions;

protected:
    int *inverse;
    int max_prio;

    std::vector<int> *regions;
    int *region;
    int *strategy;
    //std::set<std::pair<int,int>> seen;

    virtual void attract(int prio, std::queue<int> queue=std::queue<int>());
    virtual void promote(int from, int to);
    virtual void resetRegion(int priority);
    virtual bool setupRegion(int index, int priority, bool mustReset);
    virtual void setDominion(int priority);
    virtual int getRegionStatus(int index, int priority);
    virtual void reportRegion(int priority);
    virtual void printState();
};

}

#endif 
