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

#ifndef PP_HPP
#define PP_HPP

#include <queue>

#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class PPSolver : public Solver
{
public:
    PPSolver(Oink *oink, Game *game);
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
