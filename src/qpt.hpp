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

#include "oink/solver.hpp"

namespace pg {

class QPTSolver : public Solver
{
public:
    QPTSolver(Oink& oink, Game& game);
    virtual ~QPTSolver();

    virtual void run();

protected:
    int *pm_nodes;
    int *strategy;

    int pl;              // current measures of <pl> parity
    int k;               // current measures: <k> components
    int max, maxo;       // highest priority for current player / opponent
    unsigned long goal;  // number of vertices of <pl> parity

    long lift_count;
    long lift_attempt;

    uintqueue todo;
    bitset dirty;

    bool bounded = false;

    bool lift(int v, int target);
    void liftloop();
    void updateState(unsigned long &_n0, unsigned long &_n1, int &_max0, int &_max1, int &_k0, int &_k1);

    void todo_push(int node) {
        if (dirty[node]) return;
        todo.push(node);
        dirty[node] = 1;
#ifndef NDEBUG
        if (trace >= 3) logger << "push(" << node << ")" << std::endl;
#endif
    }

    int todo_pop() {
        int node = todo.pop();
        dirty[node] = 0;
#ifndef NDEBUG
        if (trace >= 3) logger << "pop() => " << node << std::endl;
#endif
        return node;
    }
};

class BoundedQPTSolver : public QPTSolver
{
public:
    BoundedQPTSolver(Oink& oink, Game& game) : QPTSolver(oink, game) { bounded = true; }
    virtual ~BoundedQPTSolver() { }
};

}

#endif 
