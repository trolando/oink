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

#include <cstdint>
#include <optional>

#include "oink/game.hpp"
#include "oink/oink.hpp"
#include "oink/error.hpp"

namespace pg {

/**
 * Explicit configuration for a solver. Bundles the solver-facing settings that
 * were previously read individually from Oink, as a step towards making solver
 * dependencies explicit instead of reaching into Oink for everything.
 */
struct SolverConfig
{
    int workers = -1;             // requested workers for parallel solvers (-1 sequential, 0 autodetect)
    int trace = 0;                // verbosity (0 normal, 1 trace, 2 debug)
    std::optional<uint64_t> seed; // optional seed for randomized solvers
};

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
    SolverConfig config;
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

    /** The solution being built (owned by Oink). */
    [[nodiscard]] Solution& solution() { return oink.solution(); }
    [[nodiscard]] const Solution& solution() const { return oink.solution(); }

    /**
     * Access to the solution being built (owned by Oink). Solvers read these as
     * working memory during solving and finalize via solve().
     */
    [[nodiscard]] bool isSolved(int v) const { return oink.solution().is_solved(v); }
    [[nodiscard]] long count_unsolved() const { return game.nodecount() - (long)oink.solution().solved().count(); }
    [[nodiscard]] int getWinner(int v) const { const Solution& s = oink.solution(); return s.is_solved(v) ? s.winner(v) : -1; }
    [[nodiscard]] int getStrategy(int v) const { return oink.solution().strategy(v); }
    [[nodiscard]] int* getStrategy() { return oink.solution().strategy_data(); }

    /** Multi-strategy access (used by fpim/fpjm). */
    [[nodiscard]] bool hasMultiStrategy() const { return oink.solution().has_multi(); }
    void initMultiStrategy() { oink.solution().init_multi(game.edgeArraySize()); }
    void addStrategyEdge(int v, int k) { oink.solution().add_edge_index((std::size_t)game.firstout(v) + k); }
    void clearStrategyEdges(int v) { oink.solution().clear_edge_range(game.firstout(v), game.outcount(v)); }
    [[nodiscard]] bool isStrategyEdgeIndex(int idx) const { return oink.solution().has_edge_index(idx); }
    void removeStrategyEdgeIndex(int idx) { oink.solution().remove_edge_index(idx); }

private:
    Oink& oink;
};

}

#endif 
