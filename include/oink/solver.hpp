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

#include "oink/game.hpp"
#include "oink/oink.hpp"
#include "oink/error.hpp"

namespace pg {

class Solver
{
public:
    Solver(Oink* oink, Game* game);
    virtual ~Solver() { }

    /**
     * Run the solver.
     */
    virtual void run() = 0;

    /**
     * Returns true if the solver always solves all enabled vertices
     * before leaving run().
     */
    virtual bool isFullSolver() { return true; }

    /**
     * Set solver options (via -c "...")
     */
    virtual bool parseOptions(std::string&) { return true; }

protected:
    Game* game;
    std::ostream &logger;
    int trace = 0;

    const bitset &disabled; // TODO make function...

    inline long nodecount() { return game->nodecount(); }
    inline long edgecount() { return game->edgecount(); }

    inline int priority(const int vertex) { return game->priority(vertex); }
    inline int owner(const int vertex) { return game->owner(vertex); }
    inline const int* outs(const int vertex) { return game->outedges() + game->firstout(vertex); }
    inline const int* ins(const int vertex) { return game->inedges() + game->firstin(vertex); }
    inline Game::_label_vertex label_vertex(int v) { return game->label_vertex(v); }

    void solve(int node, int winner, int strategy) { oink->solve(node, winner, strategy); }
    void flush(void) { oink->flush(); }

private:
    Oink* oink;
};

}

#endif 
