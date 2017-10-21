#ifndef ZLK_HPP
#define ZLK_HPP

#include <queue>

#include "solver.hpp"
#include "lace.h"

namespace pg {

class ZLKSolver : public Solver
{
public:
    ZLKSolver(Oink *oink, Game *game, std::ostream &lgr);
    virtual ~ZLKSolver();

    virtual void run();

    int iterations;

protected:
    int *inverse;
    int max_prio;

    int *region;
    int *winning;
    int *strategy;
    int *outcount;

    int attractExt(int i, int r, std::vector<int> *R);
    int attractLosing(int i, int r, std::vector<int> *S, std::vector<int> *R);

    friend void attractParT_WORK(WorkerP*, Task*, int, int, int, ZLKSolver*);
    friend int attractPar_WORK(WorkerP*, Task*, int, int, std::vector<int>*, ZLKSolver*);
    friend void updateOutcount_WORK(WorkerP*, Task*, int, int, ZLKSolver*);
};

}

#endif 
