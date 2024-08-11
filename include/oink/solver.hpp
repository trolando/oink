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

/**
 * Base class for the parity game solvers.
 * A derived class should implement the run method, and solve all vertices of the game.
 * TODO: instead of having 'disabled' vertices, the constructor should simply get
 *       a 'subgame' bitset that identifies the part of the game that we need to solve.
 */
class Solver
{
public:
    Solver(Oink& oink, Game& game);
    virtual ~Solver() = default;

    /**
     * Run the solver.
     */
    virtual void run() = 0;

    /**
     * Set solver options (via -c "...")
     */
    virtual bool parseOptions(std::string&) { return true; }

protected:
    Game& game;
    std::ostream &logger;
    int trace = 0;

    const bitset &disabled; // TODO change into subgame
    //inline bool disabled(int vertex) { return oink.disabled[vertex]; }

    [[nodiscard]] long nodecount() const { return game.nodecount(); }
    [[nodiscard]] long edgecount() const { return game.edgecount(); }
    [[nodiscard]] int priority(int vertex) const { return game.priority(vertex); }
    [[nodiscard]] int owner(int vertex) const { return game.owner(vertex); }
    [[nodiscard]] const int* outs(int vertex) const { return game.outedges() + game.firstout(vertex); }
    [[nodiscard]] const int* ins(int vertex) const { return game.inedges() + game.firstin(vertex); }
    [[nodiscard]] Game::_label_vertex label_vertex(int v) const { return game.label_vertex(v); }

    void solve(int node, int winner, int strategy) { oink.solve(node, winner, strategy); }
    void flush() { oink.flush(); }

private:
    Oink& oink;
};

}

#endif 
