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

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "game.hpp"
#include "oink.hpp"
#include <signal.h>

#define LOGIC_ERROR { printf("\033[1;7mlogic error %s:%d!\033[m\n", __FILE__, __LINE__); raise(SIGABRT); }

namespace pg {

class Solver
{
public:
    Solver(Oink *oink, Game *game, std::ostream &lgr = std::cout) :
        oink(oink), game(game), logger(lgr),
        n_nodes(game->n_nodes), priority(game->priority), owner(game->owner),
        out(game->out), in(game->in), disabled(game->disabled),
        outa(oink->outa), ina(oink->ina), outs(oink->outs), ins(oink->ins) { }
    virtual ~Solver() { }
    virtual void run() = 0;

    void setTrace(int value) { trace = value; }

protected:
    Oink *oink;
    Game *game;
    std::ostream &logger;
    int trace = 0;

    const int n_nodes;
    const int * const priority;
    const int * const owner;
    const std::vector<int> * const out;
    const std::vector<int> * const in;
    const int * const disabled;

    const int* const outa;
    const int* const ina;
    const int* const outs;
    const int* const ins;
};

}

#endif 
