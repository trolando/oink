/*
 * Copyright 2017-2020 Tom van Dijk, Johannes Kepler University Linz
 * Copyright 2019-2020 Ruben Lapauw, KU Leuven
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

#ifndef PPJ_HPP
#define PPJ_HPP

#include <queue>
#include "pp.hpp"
#include "oink.hpp"
#include "solver.hpp"

namespace pg {

class PPJSolver : public PPSolver
{
public:
    PPJSolver(Oink *oink, Game *game);
    virtual void run();

protected:
// Todo: possibly a large memory consumption.
    std::vector<int> *waiting;
    int waitingPriority;
    std::vector<int> lost;
    bitset is_lost;

    virtual void unattracted(int node);
    virtual void endAttracting(int prio);
    virtual bool setupRegion(int index, int priority, bool mustReset);

    void setWaiting(int node, int priority);
    escape bestEscape(int node);
};

}

#endif
