/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef SPM_HPP
#define SPM_HPP

#include <queue>

#include "oink/solver.hpp"

namespace pg {

class SPMSolver : public Solver
{
public:
    SPMSolver(Oink& oink, Game& game);
    virtual ~SPMSolver();

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
    int pm_cycles(int *a, int pl);
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
