#ifndef MSPM_HPP
#define MSPM_HPP

#include <queue>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class MSPMSolver : public Solver
{
public:
    MSPMSolver(Oink *oink, Game *game, std::ostream &lgr = std::cout);
    virtual ~MSPMSolver();

    virtual void run();

    int iterations = 0;
    int lift_attempt = 0;
    int lift_count = 0;

protected:
    int *pms;
    int *tmp, *best;
    int *strategy;
    int *counts;
    int *cover;
    int coverdepth;
    int k;

    std::deque<int> todo;
    std::vector<int> looping;
    int *dirty;

    bool lift(int node, int target);
    bool pm_less(int *a, int *b, int d, int pl);
    void pm_copy(int *dst, int *src, int pl);
    void pm_stream(std::ostream &out, int *pm);
    void Prog(int *dst, int *src, int d, int pl);

    void solve(int node, int str);
    void coverlower(int node, int k);
    void uncover(int k);

    void todo_push(int node) {
        if (dirty[node]) return;
        todo.push_back(node);
        dirty[node] = 1;
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
    }

    int todo_pop() {
        int node = todo.front();
        dirty[node] = 0;
        todo.pop_front();
        if (trace >= 2) logger << "pop() => " << node << std::endl;
        return node;
    }
};

}

#endif 
