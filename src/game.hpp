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

#ifndef GAME_HPP
#define GAME_HPP

#include <sstream>
#include <vector>
#include <map>

#include <boost/dynamic_bitset.hpp>

namespace pg {

typedef boost::dynamic_bitset<unsigned long long> bitset;

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
     * Parse a pgsolver game.
     */
    Game(std::istream &inp);

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
     */
    void initNode(int node, int priority, int owner, std::string label="");

    /**
     * Add an edge from <from> to <to>.
     * Returns true if the edge was added or false if it already existed.
     */
    bool addEdge(int from, int to);

    /**
     * Remove an edge from <from> to <to>.
     * Returns true if the edge was removed or false if it did not exist.
     */
    bool removeEdge(int from, int to);

    /**
     * Parse a pgsolver game.
     */
    void parse_pgsolver(std::istream &in);

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
     * Reindex the game if it was not yet reindexed.
     */
    void reindex_once(void);

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
     * Return the number of vertices.
     */
    inline size_t nodecount() { return n_nodes; }

    /**
     * Count and return the number of edges.
     */
    inline size_t edgecount() { edge_recount(); return n_edges; }

    /**
     * Returns whether every node has dominion 0 or 1.
     */
    inline bool gameSolved() { return solved.count() == nodecount(); }

    /**
     * Count and return how many nodes have dominion -1.
     */
    inline int countUnsolved() { return n_nodes - solved.count(); }

    inline void edge_recount()
    {
        n_edges = 0;
        for (int n=0; n<n_nodes; n++) n_edges += out[n].size();
    }

    /**
     * Create a new Game of the subgame of the nodes given in <selection>.
     * (We don't check if the relation is total.)
     * Afterwards, <mapping> contains for index i the corresponding index in the full game.
     * (mapping must be int[selection.size()])
     */
    Game *extract_subgame(std::vector<int> &selection, int *mapping=NULL);

    /**
     * Reset <solved>, <winner> and <strategy>.
     */
    void reset();

    /**
     * Copy the game.
     */
    Game &operator=(const Game &other);

    /**
     * Swap with other game.
     */
    void swap(Game& other);

    /**
     * Game fields
     */

    int n_nodes;           // number of nodes
    int n_edges;           // number of edges
    int *priority;         // priority of each node
    bitset owner;          // owner of each node (1 for odd, 0 for even)
    std::string *label;    // (optional) node labels
    std::vector<int> *out; // outgoing edges
    std::vector<int> *in;  // incoming edges

    bitset solved;         // set true if node solved
    bitset winner;         // for solved vertices, set 1 if won by 1, else 0
    int *strategy;         // strategy for winning vertices

    bool reindexed;        // records if the game was reindexed (before solving)

};

}

#endif 
