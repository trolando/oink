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

#ifndef RR_HPP
#define RR_HPP

#include <queue>

#include "solvers/pp.hpp"

namespace pg {

class RRSolver : public PPSolver
{
public:
    RRSolver(Oink& oink, Game& game);

    virtual void run();

protected:
    virtual bool checkRegion(int priority);
};

}

#endif 
