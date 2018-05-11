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

#ifndef QPT_HPP
#define QPT_HPP

#include <boost/dynamic_bitset.hpp>

#include "solver.hpp"
#include "uintqueue.hpp"

namespace pg {

class QPTSolver : public Solver
{
public:
    QPTSolver(Oink *oink, Game *game);
    virtual ~QPTSolver();

    virtual void run();

    int iterations = 0;
    int lift_attempt = 0;
    int lift_count = 0;

protected:
    int *pm_nodes;
    int *strategy;
    int k, k0, k1;
    int max, max_even, max_odd;
    unsigned long goal, goal0, goal1;

    uintqueue todo;
    boost::dynamic_bitset<unsigned long long> dirty;

    bool lift(int v, int target, int pl);
    void liftloop(int pl);

    void todo_push(int node) {
        if (dirty[node]) return;
        todo.push(node);
        dirty[node] = 1;
#ifndef NDEBUG
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
#endif
    }

    int todo_pop() {
        int node = todo.pop();
        dirty[node] = 0;
#ifndef NDEBUG
        if (trace >= 2) logger << "pop() => " << node << std::endl;
#endif
        return node;
    }
};

}

#endif 
