#ifndef QPT_HPP
#define QPT_HPP

#include <queue>

#include "solver.hpp"

namespace pg {

class QPTSolver : public Solver
{
public:
    QPTSolver(Oink *oink, Game *game, std::ostream &lgr);
    virtual ~QPTSolver();

    virtual void run();

    int iterations = 0;
    int lift_attempt = 0;
    int lift_count = 0;

protected:
    int *pm_nodes;
    int *strategy;
    int k;

    void print_state(std::vector<int> *choices);
    bool try_lift(int node, std::vector<int> &vec);
    bool lift(int node);
    bool liftR(int node, int target);
};

}

#endif 
