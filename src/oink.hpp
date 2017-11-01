#ifndef OINK_HPP
#define OINK_HPP

#include <iostream>
#include <vector>
#include "game.hpp"

namespace pg {

class Solver;

class Oink
{
public:
    Oink(Game &game, std::ostream &out=std::cout) : game(&game),  logger(out) { }
    ~Oink() {}

    /**
     * After configuring Oink, use run() to run the solver.
     */
    void run();

    /**
     * Instruct Oink to use the given solver.
     */
    void setSolver(int solverid);
    void setSolver(std::string label);

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
     * Solve node <node> as won by <winner> with strategy <strategy>.
     * (Set <strategy> to -1 for no strategy.
     */
    void solve(int node, int winner, int strategy);

    /**
     * After marking nodes as solved using solve(), use flush to attract to the solved dominion.
     */
    void flush();

protected:
    /**
     * Solve winner-controlled winning cycles.
     */
    int solveTrivialCycles(void);

    /**
     * Resolve self-loops.
     */
    int solveSelfloops(void);

    /**
     * During flush(), attract to solved node <i>.
     */
    void attractDominion(int i);

    Game *game;              // game being solved
    std::ostream &logger;    // logger for trace/debug messages
    int solver = -1;         // which solver to use
    int workers = -1;        // number of workers, 0 = autodetect, -1 = use non parallel
    int trace = 0;           // verbosity (0 for normal, 1 for trace, 2 for debug)
    bool inflate = false;    // inflate the game before solving
    bool compress = false;   // compress the game before solving
    bool renumber = false;   // renumber the game before solving (removes gaps)
    bool removeLoops = true; // resolve self-loops before solving
    bool removeWCWC = true;  // solve winner-controlled winning cycles before solving
    bool solveSingle = true; // solve games with only 1 parity
    bool bottomSCC = false;  // solve per bottom SCC

    std::vector<int> todo;   // internal queue for solved nodes for flushing
    int *outcount;           // number of unsolved outgoing edges per node (for fast attraction)
    int *outa;               // index array for outgoing edges
    int *ina;                // index array for incoming edges
    int *outs;               // all outgoing edges
    int *ins;                // all incoming edges

    friend class pg::Solver; // to allow access to edges
};

}

#endif 
