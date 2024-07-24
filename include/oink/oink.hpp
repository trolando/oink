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

#ifndef OINK_HPP
#define OINK_HPP

#include <iostream>
#include <vector>
#include <optional>

#include "oink/error.hpp"
#include "oink/game.hpp"
#include "oink/uintqueue.hpp"
#include "oink/solvers.hpp"

namespace pg {

class Solver;

class Oink
{
public:
    /**
     * Construct a solver for the given game, sending log output to <out>.
     */
    Oink(Game &game, std::ostream &out=std::cout);

    virtual ~Oink();

    /**
     * After configuring Oink, use run() to run the solver.
     */
    void run();

    /**
     * Instruct Oink to use the given solver.
     */
    void setSolver(std::string solver);

    /**
     * Instruct Oink to inflate as a preprocessing step.
     */
    void setInflate() { inflate = true; renumber = false; compress = false; }

    /**
     * Instruct Oink to compress as a preprocessing step.
     */
    void setCompress() { inflate = false; compress = true; renumber = false; }

    /**
     * Instruct Oink to renumber as a preprocessing step.
     */
    void setRenumber() { inflate = false; compress = false; renumber = true; }

    /**
     * Instruct Oink whether to remove self-loops as a preprocessing step. (Default true)
     */
    void setRemoveLoops(bool val) { removeLoops = val; }

    /**
     * Instruct Oink whether to solve winner-controlled winning cycles as a preprocessing step. (Default true)
     */
    void setRemoveWCWC(bool val) { removeWCWC = val; }

    /**
     * Instruct Oink whether to solve single parity games. (Default true)
     */
    void setSolveSingle(bool val) { solveSingle = val; }

    /**
     * Instruct Oink whether solve per bottom SCC. (Default false)
     */
    void setBottomSCC(bool val) { bottomSCC = val; }

    /**
     * Set the number of workers for parallel solvers (psi and zielonka).
     * -1 for sequential code, 0 for autodetect.
     */
    void setWorkers(int count) { workers = count; }

    /**
     * Set verbosity level (0 = normal, 1 = trace, 2 = debug)
     */
    void setTrace(int level) { trace = level; }

    /**
     * Set options for the solver
     */
    void setSolverOptions(std::string options) { this->options = options; }

    /**
     * Mark node <node> as won by <winner> with strategy <strategy>.
     * (Set <strategy> to -1 for no strategy.)
     * After marking nodes, call flush().
     */
    void solve(int node, int winner, int strategy);

    /**
     * After marking nodes as solved using solve(), flush attracts to the solved dominion.
     */
    void flush(void);

protected:
    /**
     * Solve winner-controlled winning cycles.
     * Returns number of cycles solved.
     */
    int solveTrivialCycles(void);

    /**
     * Resolve self-loops.
     * Returns number of resolved self-loops.
     */
    int solveSelfloops(void);

    /**
     * Solve single parity games.
     * Returns true if the game was solved as a single parity game.
     */
    bool solveSingleParity(void);

    /**
     * Find a bottom SCC starting from the first non-disabled node.
     * Avoids "disabled" nodes.
     * (if nonempty is set, only obtain a non-empty bottom SCC.)
     */
    void getBottomSCC(std::vector<int> &scc, bool nonempty=false);

    /**
     * Find a bottom SCC starting from the given node.
     * Avoids "disabled" nodes.
     * (if nonempty is set, only obtain a non-empty bottom SCC.)
     */
    void getBottomSCC(int start, std::vector<int> &scc, bool nonempty=false);

    /**
     * Tarjan's SCC algorithm, modified to only compute the bottom SCC and avoid disabled nodes.
     */
    void tarjan(int start_node, std::vector<int> &res, bool nonempty);

    /**
     * Run the solver in a loop until the game is solved.
     */
    void solveLoop(void);
    friend void _solve_loop(Oink*); // access point from a Lace worker

    Game *game;              // game being solved
    std::ostream &logger;    // logger for trace/debug messages
    std::optional<std::string> solver; // which solver to use
    int workers = -1;        // number of workers, 0 = autodetect, -1 = use non parallel
    int trace = 0;           // verbosity (0 for normal, 1 for trace, 2 for debug)
    bool inflate = false;    // inflate the game before solving
    bool compress = false;   // compress the game before solving
    bool renumber = false;   // renumber the game before solving (removes gaps)
    bool removeLoops = true; // resolve self-loops before solving
    bool removeWCWC = true;  // solve winner-controlled winning cycles before solving
    bool solveSingle = true; // solve games with only 1 parity
    bool bottomSCC = false;  // solve per bottom SCC
    std::string options = "";// options for the solver

    uintqueue todo;          // internal queue for solved nodes for flushing
    int *outcount;           // number of unsolved outgoing edges per node (for fast attraction)
    bitset disabled;         // which vertices are disabled

    friend class pg::Solver; // to allow access to edges
};

}

#endif 
