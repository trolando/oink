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

#ifndef PPP_HPPP
#define PPP_HPPP

#include <queue>

#include "pp.hpp"

namespace pg {

class PPPSolver : public PPSolver
{
public:
    PPPSolver(Oink& oink, Game& game);

    virtual void run();

protected:
    int reset0, reset1;
};

}

#endif 
