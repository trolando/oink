/*
 * Copyright 2020 Tom van Dijk
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

#ifndef FPJ_HPP
#define FPJ_HPP

#include "oink/solver.hpp"
#include "lace.h"

namespace pg {

class FPJSolver : public Solver
{
public:
    FPJSolver(Oink& oink, Game& game);
    virtual ~FPJSolver();

    unsigned long long iterations = 0;
    bool greedy = false;

    void runSeqGreedy(void);
    void runSeq(void);

    virtual void run();
};

class FPJGSolver : public FPJSolver
{
public:
    FPJGSolver(Oink& oink, Game& game) : FPJSolver(oink, game) { greedy = true; }
    virtual ~FPJGSolver() { }
};

}

#endif 
