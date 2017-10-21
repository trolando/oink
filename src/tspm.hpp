#ifndef TSPM_HPP
#define TSPM_HPP

#include <queue>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class TSPMSolver : public Solver
{
public:
    TSPMSolver(Oink *oink, Game *game, std::ostream &lgr = std::cout);
    virtual ~TSPMSolver();

    virtual void run();

    int64_t lift_attempt = 0;
    int64_t lift_count = 0;

protected:
    int *pms;
    int *tmp, *best;
    int *strategy;
    int *counts;
    int64_t k;

    std::deque<int> todo;
    int *dirty;
    int *unstable;

    bool canlift(int node, int pl);
    bool lift(int node, int target);
    bool pm_less(int *a, int *b, int d, int pl);
    void pm_copy(int *dst, int *src, int pl);
    void pm_stream(std::ostream &out, int *pm);
    void Prog(int *dst, int *src, int d, int pl);
    void update(int pl);

    void todo_push(int node) {
        if (dirty[node]) return;
        todo.push_back(node);
        dirty[node] = 1;
#ifndef NDEBUG
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
#endif
    }

    int todo_pop() {
        int node = todo.front();
        todo.pop_front();
        dirty[node] = 0;
        if (trace >= 2) logger << "pop() => " << node << std::endl;
        return node;
    }
};

}

#endif 
