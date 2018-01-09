#ifndef GAME_HPP
#define GAME_HPP

#include <sstream>
#include <vector>
#include <map>

namespace pg {

class Game
{
public:
    /**
     * Construct a new (uninitialized) parity game.
     */
    Game();

    /**
     * Construct a new parity game for <count> nodes.
     */
    Game(int count);

    /**
     * Construct a copy of an existing parity game.
     */
    Game(const Game& other);

    /**
     * Deconstructor.
     */
    ~Game();

    /**
     * Initialize the game with <count> uninitialized nodes.
     */
    void initGame(int count);

    /**
     * Initialize a node <node> with given <priority>, <owner> and <label>.
     * Also resets disabled, dominion and strategy.
     */
    void initNode(int node, int priority, int owner, std::string label="");

    /**
     * Add an edge from <from> to <to>.
     * Returns true if the edge was added or false if it already existed.
     */
    bool addEdge(int from, int to);

    /**
     * Parse a pgsolver game and return the number of edges.
     */
    size_t parse_pgsolver(std::istream &inp);

    /**
     * Parse a [full or partial] pgsolver solution.
     */
    void parse_solution(std::istream &in);

    /**
     * Write the game in pgsolver format.
     */
    void write_pgsolver(std::ostream &out);

    /**
     * Write the game as a DOT graph.
     */
    void write_dot(std::ostream &out);

    /**
     * Write the solution in pgsolver format.
     */
    void write_sol(std::ostream &cout);

    /**
     * Sort the nodes in order of priority (low to high).
     * Afterwards, <mapping> is such that node <i> was originally at <mapping[i]>.
     */
    void reindex(int *mapping = NULL);

    /**
     * Apply a permutation, moving node <i> to position <mapping[i]>.
     * This reverses a reindex operation.
     */
    void permute(int *mapping); // undo reindex

    /**
     * Reassign priorities such that every node has a unique priority.
     * (Assumes reindex() has been called earlier.)
     * Returns number of distinct priorities.
     */
    int inflate(void);

    /**
     * Reassign priorities such that no priority is skipped. ("compression")
     * (Assumes reindex() has been called earlier.)
     * Returns number of distinct priorities.
     */
    int compress(void);

    /**
     * Reassign priorities in order, but do not inflate or compress.
     * (Assumes reindex() has been called earlier.)
     * Returns number of distinct priorities.
     */
    int renumber(void);

    /**
     * Swap players (priorities and ownership)
     * (Assumes reindex() has been called earlier.)
     */
    void evenodd(void);

    /**
     * Change a max game into a min game and vice versa.
     * (Assumes reindex() has been called earlier.)
     */
    void minmax(void);

    /**
     * Count and return the number of edges.
     */
    size_t edgecount();

    /**
     * Returns whether every node has dominion 0 or 1.
     */
    bool solved();

    /**
     * Count and return how many nodes have dominion -1.
     */
    int countUnsolved();

    /**
     * Create a new Game of the subgame of the nodes given in <selection>.
     * (We don't check if the relation is total.)
     * Afterwards, <mapping> contains for index i the corresponding index in the full game.
     * (mapping must be int[selection.size()])
     */
    Game *extract_subgame(std::vector<int> &selection, int *mapping=NULL);

    /**
     * Set disabled to 0 for all nodes in <selection>, to 1 otherwise.
     */
    void restrict(std::vector<int> &selection);

    /**
     * Reset <dominion>, <strategy> and <disabled>.
     */
    void reset();

    /**
     * Copy the game.
     */
    Game &operator=(const Game &other);

    int n_nodes;           // number of nodes
    int *priority;         // priority of each node
    int *owner;            // owner of each node
    std::vector<int> *out; // outgoing edges
    std::vector<int> *in;  // incoming edges
    std::string *label;    // (optional) node labels
    int *dominion;         // winner of node (or -1 if unknown)
    int *strategy;         // strategy for the winner

    int *disabled;         // if the node is disabled (for restricting to subgames)
};

}

#endif 
