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
