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
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <bitset.hpp>

namespace pg {

/**
 * The main class holding a Parity game.
 *
 * Internally, the implementation leverages virtual memory vs real memory.
 * We overallocate virtual memory, but use less actual memory.
 *
 * Edges are stored consecutively in a single array, ending with th evalue -1.
 * For game solving, the build_in_array method creates the reverse array.
 *
 * There is a special vector representation of edges to allow for "random-order" modification
 * or game-building:
 * - vec_init() initializes the vectors with the current edges
 * - vec_add_edge, vec_remove_edge, vec_has_edge to manipulate the edges
 * - vec_finish() rebuilds the array representation
 *
 * Usage scenario for making a game (random order):
 * - initialize with constructor Game(count) or using init_game(count)
 * - use vec_init
 * - use init_vertex, vec_add_edge to make the game
 * - use vec_finish
 * - use sort, renumber, write_pgsolver, et
 *
 * Usage scenario for making a game (streaming):
 * - initialize with some number of vertices (e.g. 1000)
 * - use init_vertex to initialize the vertex
 * - use e_start to start adding successors of the vertex
 * - use e_add to add each successor
 * - use e_finish to finish adding successors of the vertex
 * - use v_sizeup (or v_resize) to increase the number of vertices if needed
 * - use v_resize to set the final number of vertices
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
    Game(int count, int ecount = -1);

    /**
     * Construct a deep clone of an existing parity game.
     * Does not clone the vector representation.
     * Does not clone the <in> array.
     */
    Game(const Game& other);

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
    void init_random_game(int n, long maxP, long maxE);
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
     * For vector-based (random order) manipulation of edges
     */
    void vec_init(void); // initialize for manipulating vertices using the vectors
    void vec_finish(void); // finalize
    bool vec_add_edge(int from, int to); // return true if changed
    bool vec_remove_edge(int from, int to); // return true if changed
    bool vec_has_edge(int from, int to); // return true if exists

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
     * (re)build the <in> array for PG solving.
     */
    void build_in_array(bool rebuild=false);

    /**
     * Dynamic size methods.
     */
    void e_sizeup(void);
    void v_sizeup(void);
    void v_resize(size_t newsize); // does not do cleaning up..?

    /**
     * Edge making helpers.
     */
    void e_start(int source);
    void e_add(int source, int target);
    void e_finish(void);

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
    Game *extract_subgame(bitset mask);

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
    inline std::string* rawlabel(const int vertex) const
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

    inline const std::vector<int> outvec(const int vertex) const
    {
        return _outvec[vertex];
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
    long n_vertices;        // number of vertices
    long n_edges;           // number of edges
    int *_priority;        // priority of each vertex
    bitset _owner;         // owner of each vertex (1 for odd, 0 for even)
    std::string **_label;  // (optional) vertex labels

    int *_outedges;        // outgoing edges as array
    int *_firstouts;       // first outgoing edge of each vertex
    int *_outcount;        // outgoing edge count of each vertex

    int *_inedges;         // incoming edges as array
    int *_firstins;        // first incoming edge of each vertex
    int *_incount;         // incoming edge count of each vertex

    std::vector<int> *_outvec; // outgoing edges as vector

    bool is_ordered;       // records if the game is in-order
    size_t v_allocated;    // number of vertices allocated as virtual memory
    size_t e_allocated;    // number of edges allocated as virtual memory
    size_t e_size;         // number of entries used in edge array

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
                    std::string* l = lv.g.rawlabel(lv.v);
                    if (l == NULL or l->empty()) out << lv.v << "/" << lv.g.priority(lv.v);
                    else out << *l;
                }
                return out;
            }
        protected:
            Game &g;
            int v;
    };

private:
    void unsafe_permute(int *mapping); // apply a reordering
    
    boost::random::mt19937 generator;
    inline long rng(long low, long high) { return boost::random::uniform_int_distribution<> (low, high)(generator); }
};

}

#endif 
