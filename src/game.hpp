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

#include <cassert>
#include <sstream>
#include <vector>
#include <map>
#include <random>

#include <bitset.hpp>

namespace pg {

/**
 * The main class holding a Parity game.
 *
 * Internally, two ways exist to represent the edges.
 * - as vector: each vertex has a std::vector of integers.
 * - as array: we store all edges in the game in a single array, seperated by a -1 value.
 *
 * All methods use the array form, except add_edge and remove_edge.
 *
 * Usage scenario for making a game:
 * - initialize with constructor Game(count) or using init_game(count)
 * - now only the vector representation is available
 * - use init_vertex, add_edge to make the game
 * - then use build_arrays() to convert vector to array representation
 *
 * Usage scenario for changing a game:
 * - use build_vectors() to convert array to vector representation
 * - use add_edge and remove_edge to change the edges
 * - use rebuild_arrays() to convert vector to array representation
 *
 * Notice that build_arrays() only converts vector to array representation if none exists yet.
 * So use rebuild_arrays() to update the array representation after a change.
 */

class Game
{
public:
    /**
     * Construct a new (uninitialized) parity game.
     */
    Game();

    /**
     * Construct a new parity game for <count> vertices.
     */
    Game(int count);

    /**
     * Construct a copy of an existing parity game.
     */
    Game(const Game& other);

    /**
     * Construct a subgame of a game.
     */
    Game(const Game& other, bitset mask);

    /**
     * Parse a pgsolver game.
     */
    Game(std::istream &inp, bool removeBadLoops=true);

    /**
     * Deconstructor.
     */
    ~Game();

    /**
     * Initialize the game with <count> uninitialized vertices.
     * After init_game(x), the state of the object is as after constructor Game(x).
     */
    void init_game(int count);

     /**
      * Create random game with <n> vertices.
      * - maximum priority <maxP>
      * - allow self-loops
      * - each vertex minimum 1 random edge
      * - then generate at most <maxE> more edges
      */
    void init_random_game(int n, int maxP, int maxE);
    inline void set_random_seed(unsigned int seed) { generator.seed(seed); }

    /**
     * Initialize a vertex <v> with given <priority>, <owner> and <label>.
     */
    void init_vertex(int v, int priority, int owner, std::string label="");

    /**
     * Change the priority of a vertex.
     * At the end of this method, is_ordered is updated.
     */
    void set_priority(int vertex, int priority);

    /**
     * Change the owner of a vertex.
     * (zero is owner Even; non-zero is owner Odd)
     */
    void set_owner(int vertex, int owner);

    /**
     * Change the label of a vertex.
     */
    void set_label(int vertex, std::string label);

    /**
     * Add an edge from <from> to <to>.
     * Returns true if the edge was added or false if it already existed.
     *
     * Important: this manipulates the vector representation only.
     * After adding/removing edges, call rebuild_arrays().
     */
    bool add_edge(int from, int to);

    /**
     * Remove an edge from <from> to <to>.
     * Returns true if the edge was removed or false if it did not exist.
     *
     * Important: this manipulates the vector representation only.
     * After adding/removing edges, call rebuild_arrays().
     */
    bool remove_edge(int from, int to);

    /**
     * Check if the edge from <from> to <to> exists, using the vector representation.
     */
    bool has_edge_vec(int from, int to);

    /**
     * Check if a certain edge exists.
     */
    bool has_edge(int from, int to);

    /**
     * Get the index of an edge in the edge array (or -1 if not found)
     * NOTE: uses the edge arrays.
     */
    int find_edge(int from, int to);

    /**
     * Parse a pgsolver game.
     * After parse_pgsolver(x), the state of the object is as after constructor Game(x).
     */
    void parse_pgsolver(std::istream &in, bool removeBadLoops=true);

    /**
     * Parse a [full or partial] pgsolver solution.
     */
    void parse_solution(std::istream &in);

    /**
     * Write the game in pgsolver format to the stream <out>.
     */
    void write_pgsolver(std::ostream &out);

    /**
     * Write the game as a DOT graph to the stream <out>.
     */
    void write_dot(std::ostream &out);

    /**
     * Write the (partial) solution in pgsolver format to the stream <out>.
     */
    void write_sol(std::ostream &out);

    /**
     * Sort the vertices in order of priority (low to high).
     * If <mapping> is given as an int array of size vertexcount(),
     * then it can be used with permute to reverse the procedure.
     */
    void sort(int *mapping = NULL);

    /**
     * Ensure that vertices are ordered by priority.
     */
    inline void ensure_sorted(void) { sort(NULL); }

    /**
     * Apply a permutation, moving each vertex <i> to position <mapping[i]>.
     * Afterwards, <is_ordered> is updated.
     */
    void permute(int *mapping); // undo reindex

    /**
     * Reassign priorities such that every vertex has a unique priority.
     * Returns number of distinct priorities.
     * (Only valid if vertices are ordered.)
     */
    int inflate(void);

    /**
     * Reassign priorities such that no priority is skipped. ("compression")
     * Returns number of distinct priorities.
     * (Only valid if vertices are ordered.)
     */
    int compress(void);

    /**
     * Reassign priorities in order, but do not inflate or compress.
     * Returns number of distinct priorities.
     * (Only valid if vertices are ordered.)
     */
    int renumber(void);

    /**
     * Swap players (priorities and ownership) and renumbers on the fly.
     * (Only valid if vertices are ordered.)
     */
    void evenodd(void);

    /**
     * Change a max game into a min game and vice versa. Renumbers on the fly.
     * (Only valid if vertices are ordered.)
     */
    void minmax(void);

    /**
     * Fix the game, no more changes.
     */
    void build_arrays(void); // build the out/in arrays that are not build
    void build_vectors(void); // if we never built the vector, do it now
    void rebuild_arrays(void); // force rebuilding from vectors
    void rebuild_vectors(void); // force rebuilding from arrays
    bool test_consistency(void); // test if the two representations agree

    /**
     * Return the number of vertices.
     */
    inline long vertexcount() const { return n_vertices; }
    inline long nodecount() const { return n_vertices; }

    /**
     * Count and return the number of edges.
     */
    inline long edgecount() const { return n_edges; }

    /**
     * Returns whether every vertex has dominion 0 or 1.
     */
    inline bool game_solved() const { return (unsigned)vertexcount() == solved.count(); }

    /**
     * Count and return how many vertices have dominion -1.
     */
    inline long count_unsolved() const { return vertexcount() - solved.count(); }

    /**
     * Create a new Game of the subgame of the vertices given in <selection>.
     */
    Game *extract_subgame(std::vector<int> &selection);

    /**
     * Reset <solved>, <winner> and <strategy>.
     */
    void reset_solution();

    /**
     * Copy solution (<other> must be a subgame and have the same number of vertices)
     */
    void copy_solution(Game &other);

    /**
     * Copy the game.
     */
    Game &operator=(const Game &other);

    /**
     * Swap with other game.
     */
    void swap(Game& other);

    /**
     * Get the priority of a vertex
     */
    inline int priority(const int vertex) const
    {
        return _priority[vertex];
    }

    /**
     * Get the owner of a vertex
     */
    inline int owner(const int vertex) const
    {
        return _owner[vertex];
    }

    /**
     * Get the "real" label of a vertex
     */
    inline std::string& rawlabel(const int vertex) const
    {
        return _label[vertex];
    }

    /**
     * Bunch of methods to *read* the ingoing/outgoing edges.
     */

    inline const int* outedges() const
    {
        return _outedges;
    }

    inline const int* inedges() const
    {
        return _inedges;
    }

    inline int firstout(const int vertex) const
    {
        return _firstouts[vertex];
    }

    inline int firstin(const int vertex) const
    {
        return _firstins[vertex];
    }

    inline long outcount(const int vertex) const
    {
        return _outcount[vertex];
    }

    inline long incount(const int vertex) const
    {
        return _incount[vertex];
    }

    inline const int *outs(const int vertex) const
    {
        return outedges() + firstout(vertex);
    }

    inline const int *ins(const int vertex) const
    {
        return inedges() + firstin(vertex);
    }

    /*inline auto out_iter(const int vertex) const
    {
        return out[vertex].cbegin();
    }*/

    inline const std::vector<int> outvec(const int vertex) const
    {
        return out[vertex];
    }

    class _label_vertex;

    _label_vertex label_vertex(int v)
    {
        return _label_vertex(*this, v);
    }

    _label_vertex label(int v)
    {
        return _label_vertex(*this, v);
    }

    /**
     * Game fields
     */

private:
    int n_vertices;        // number of vertices
    int n_edges;           // number of edges
    int *_priority;        // priority of each vertex
    bitset _owner;         // owner of each vertex (1 for odd, 0 for even)
    std::string *_label;   // (optional) vertex labels

    std::vector<int> *out; // outgoing edges as vector

    int *_outedges;        // outgoing edges as array
    int *_firstouts;       // first outgoing edge of each vertex
    int *_outcount;        // outgoing edge count of each vertex

    int *_inedges;
    int *_firstins;
    int *_incount;

    bool is_ordered;       // records if the game is in-order
    bool is_mmap;          // outgoing edges allocated as virtual memory
    size_t oe_allocated;   // amount of memory allocated as virtual memory

public:
    bitset solved;         // set true if vertex solved
    bitset winner;         // for solved vertices, set 1 if won by 1, else 0
    int *strategy;         // strategy for winning vertices

    class _label_vertex
    {
        protected:
            _label_vertex(Game &g, int v) : g(g), v(v) { }
            friend class Game;
        public:
            friend std::ostream& operator<<(std::ostream& out, const _label_vertex &lv) {
                if (lv.v == -1) {
                    out << "-1";
                } else {
                    std::string& l = lv.g.rawlabel(lv.v);
                    if (l.empty()) out << lv.v << "/" << lv.g.priority(lv.v);
                    else out << l;
                }
                return out;
            }
        protected:
            Game &g;
            int v;
    };

private:
    void unsafe_permute(int *mapping); // apply a reordering
    void test_vectors_built(void);
    void test_arrays_built(void);
    
    static std::random_device rd;
    std::mt19937 generator;
    inline int rng(int low, int high) { return std::uniform_int_distribution<int>(low, high)(generator); }
};

}

#endif 
