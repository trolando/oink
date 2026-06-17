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

#include <algorithm>
#include <cassert>
#include <iostream>
#include <ctime>

#include "oink/game.hpp"
#include "oink/subgame_view.hpp"

#define USE_MMAP 1

using namespace std;

namespace pg {

Game::Game() : _owner(0), solution_(0)
{
    n_vertices = 0;
    n_edges = 0;
    is_ordered = true;
    v_allocated = 0;
    e_allocated = 0;
    e_size = 0;
    set_random_seed(static_cast<unsigned int>(std::time(0)));
}

Game::Game(size_t nv, size_t ne, std::vector<int>& priorities, bitset& owners, std::vector<std::vector<int>>& edges, std::vector<std::string*>& labels) : Game(nv, ne)
{
    n_vertices = nv;
    n_edges = ne;

    // copy priorities
    std::move(priorities.begin(), priorities.end(), _priority.begin());

    // copy owners
    _owner = owners;

    // copy edges
    int e = 0;
    for (auto v=0; v<n_vertices; v++) {
        const auto f = e;
        _firstouts[v] = e;
        for (auto & to : edges[v]) {
            _outedges[e++] = to;
        }
        _outcount[v] = e-f;
        _outedges[e++] = -1;
    }
    e_size = e;

    // copy labels
    for (auto v=0; v<n_vertices; v++) {
        if (labels[v] != nullptr) _label[v] = *labels[v];
    }

    // check if ordered
    is_ordered = true;
    for (int i=1; i<n_vertices; i++) {
        if (_priority[i-1] > _priority[i]) {
            is_ordered = false;
            break;
        }
    }
}

Game::~Game()
{
    // all storage is RAII-managed (std::vector / bitset / Solution)
}

Game::Game(int vcount, int ecount) : _owner(vcount), solution_(vcount)
{
    assert(vcount > 0);
    if (ecount == -1) ecount = size_t(4) * vcount; // reasonable default outdegree

    n_vertices = vcount;
    n_edges = 0;

    v_allocated = vcount;
    e_allocated = vcount+ecount+1;  // extra space for -1
    e_size = 0;

    // vertex arrays sized to capacity, value-initialized to 0
    _priority.resize(v_allocated);
    _label.resize(v_allocated);
    _firstouts.resize(v_allocated);
    _outcount.resize(v_allocated);
    _outedges.resize(e_allocated);

    is_ordered = true;

    _outedges[0] = -1;
    e_size++;

    set_random_seed(static_cast<unsigned int>(std::time(0)));
}

/**
 * Make a deep clone of the given game <other>.
 * Does not clone the vector representation.
 * Does not clone the <in> array.
 */
Game::Game(const Game& other) : Game(other.n_vertices, other.e_size)
{
    n_edges = other.n_edges;

    std::copy_n(other._priority.data(), n_vertices, _priority.data());
    _owner = other._owner;
    for (int i=0; i<n_vertices; i++) {
        _label[i] = other._label[i];
    }

    // clone the edge out ARRAY
    e_size = other.e_size;
    std::copy_n(other._outedges.data(), e_size, _outedges.data());
    std::copy_n(other._firstouts.data(), n_vertices, _firstouts.data());
    std::copy_n(other._outcount.data(), n_vertices, _outcount.data());

    // copy inedges
    if (!other._inedges.empty()) {
         size_t len = n_vertices + n_edges;
         _inedges.resize(len);
         _firstins.resize(n_vertices);
         _incount.resize(n_vertices);
         std::copy_n(other._inedges.data(), len, _inedges.data());
         std::copy_n(other._firstins.data(), n_vertices, _firstins.data());
         std::copy_n(other._incount.data(), n_vertices, _incount.data());
    }

    is_ordered = other.is_ordered;

    solution_ = other.solution_;

    set_random_seed(static_cast<unsigned int>(std::time(0)));
}

Game::Game(Game&& other) noexcept : Game()
{
    swap(other);
}

void
Game::init_game(int count)
{
    Game g(count);
    swap(g);
}

/**
 * Create random game with <n> vertices.
 * - maximum priority <maxP>
 * - allow self-loops
 * - each vertex minimum 1 random edge
 * - then generate at most <maxE> more edges
 */
void
Game::init_random_game(int n, long maxP, long maxE)
{
    // reset/initialize game with vector representation
    init_game(n);

    // First initialize all vertices, and give each vertex one random successor
    vec_init();

    for (int i=0; i<n; i++) {
        // initialize vertex i with random priority and random owner
        init_vertex(i, rng(0, maxP), rng(0, 1));
        // add 1 random edge (including self-loops)
        vec_add_edge(i, rng(0, n-1));
    }

    // Then add more edges randomly, at most maxE extra edges (random)
    int sources[n], source_count = n;
    for (int i=0; i<n; i++) sources[i] = i;

    // This is optimized for SPARSE random graphs!!

    for (maxE = rng(0, maxE); maxE != 0; maxE--) {
        if (source_count == 0) break;
        // select a random source vertex
        int src_idx = rng(0, source_count-1);
        int from = sources[src_idx];
        // select a random target vertex
        auto to = rng(0, n-1);
        if (vec_add_edge(from, to)) {
            if (outvec(from).size() == (unsigned)n) {
                // last target, so remove from sources
                sources[src_idx] = sources[--source_count];
            }
        } else {
            maxE++;
        }
    }

    vec_finish();
}

void
Game::init_vertex(int v, int priority, int owner, std::string label)
{
    assert(v >= 0);
    while (v >= n_vertices) v_sizeup();
    set_priority(v, priority);
    set_owner(v, owner);
    set_label(v, label);
    solution_.set_strategy(v, -1); // initialize strategy
}

void
Game::set_priority(int node, int priority)
{
    _priority[node] = priority;
    if (is_ordered) {
        if (node > 0 and _priority[node-1] > _priority[node]) is_ordered = false;
        else if (node < (n_vertices-1) and _priority[node] > _priority[node+1]) is_ordered = false;
    }
}

void
Game::set_owner(int node, int owner)
{
    this->_owner[node] = owner ? 1 : 0;
}

void
Game::set_label(int node, std::string label)
{
    _label[node] = std::move(label);
}


/**
 * Vector stuff
 */

void
Game::vec_init(void)
{
    _outvec.assign(n_vertices, {});

    // copy current edges to vectors
    for (int v=0; v<n_vertices; v++) {
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            _outvec[v].push_back(*curedge);
        }
    }
}

void
Game::vec_finish(void)
{
    e_size = 0;
    n_edges = 0;
    for (int v=0; v<n_vertices; v++) {
        e_start(v);
        for (int to : _outvec[v]) e_add(v, to);
        e_finish();
    }
    _outvec.clear();
    _outvec.shrink_to_fit();
}

bool
Game::vec_add_edge(int from, int to)
{
    assert(from >= 0 and from < n_vertices);
    assert(to >= 0 and to < n_vertices);

    if (!vec_has_edge(from, to)) {
        _outvec[from].push_back(to);
        return true;
    } else {
        return false;
    }
}

bool
Game::vec_remove_edge(int from, int to)
{
    assert(from >= 0 and from < n_vertices);
    assert(to >= 0 and to < n_vertices);

    if (vec_has_edge(from, to)) {
        auto &o = _outvec[from];
        o.erase(std::remove(o.begin(), o.end(), to), o.end());
        return true;
    } else {
        return false;
    }
}

bool
Game::vec_has_edge(int from, int to)
{
    return std::find(_outvec[from].begin(), _outvec[from].end(), to) != _outvec[from].end();
}

bool
Game::has_edge(int from, int to) const
{
    return find_edge(from, to) != -1;
}

int
Game::find_edge(int from, int to) const
{
    for (int idx = _firstouts[from]; _outedges[idx] != -1; idx++) {
        if (_outedges[idx] == to) return idx;
    }
    return -1;
}

void
Game::parse_solution(std::istream &in)
{
    string line;
    while (getline(in, line)) {
        stringstream ss(line);
        string token;

        // ignore empty line
        if (!(ss >> token)) continue;

        // ignore line with "paritysol"
        if (token == "paritysol") continue;

        // get node
        int ident = stoi(token);
        if (ident < 0 || ident >= n_vertices) {
            throw std::runtime_error("node index out of bounds");
        }

        if (solution_.is_solved(ident)) throw std::runtime_error("node already solved");

        // parse winner
        int w;
        if (!(ss >> w)) throw std::runtime_error("missing winner");
        if (w!= 0 && w!= 1) throw std::runtime_error("invalid winner");

        // parse strategy
        int str = -1;
        if (w == _owner[ident]) {
            if (!(ss >> str)) throw std::runtime_error("missing strategy for winning node");
            // if (!has_edge(ident, str)) throw std::runtime_error("strategy not successor of node");
            // actually this is already checked by the verifier
        }

        solution_.solve(ident, w, str);
    }
}

void
Game::write_pgsolver(std::ostream &os)
{
    // print banner
    os << "parity " << n_vertices << ";" << std::endl;

    // print vertices
    for (int i=0; i<n_vertices; i++) {
        os << i << " " << priority(i) << " " << owner(i) << " ";
        bool first = true;
        for (auto curedge = outs(i); *curedge != -1; curedge++) {
            if (first) first = false;
            else os << ",";
            os << *curedge;
        }
        if (!_label[i].empty()) os << " \"" << _label[i] << "\"";
        os << ";" << std::endl;
    }
}

void
Game::write_dot(std::ostream &out)
{
    out << "digraph G {" << std::endl;
    for (int i=0; i<n_vertices; i++) {
        out << i << " [ shape=\"" << (owner(i) ? "box" : "diamond")
            << "\", label=\"" << priority(i) << "\"];" << std::endl;
        for (auto curedge = outs(i); *curedge != -1; curedge++) {
            out << i << " -> " << (*curedge) << ";" << std::endl;
        }
    }
    out << "}" << std::endl;
}

/**
 * Write a (partial) solution to the stream <out>.
 */
void
Game::write_sol(std::ostream &out)
{
    // print banner
    out << "paritysol " << solution_.solved().count() << ";" << std::endl;

    // print solution
    for (int i=0; i<n_vertices; i++) {
        if (solution_.is_solved(i)) {
            out << i << " " << (solution_.winner(i) ? "1" : "0");
            if (solution_.winner(i) == _owner[i] and solution_.strategy(i) != -1) out << " " << solution_.strategy(i);
            out << ";" << std::endl;
        }
    }
}

/**
 * Sort all vertices by priority.
 */
void
Game::sort(int *mapping)
{
    if (is_ordered) {
        // already ordered, only update mapping if given
        if (mapping != NULL)  {
            for (int i=0; i<n_vertices; i++) mapping[i] = i;
        }
    } else if (mapping == NULL) {
        // no mapping given, so we allocate one and then free it afterwards
        mapping = new int[n_vertices];
        sort(mapping);
        delete[] mapping;
    } else {
        // initialize mapping
        for (int i=0; i<n_vertices; i++) mapping[i] = i;

        // sort the mapping
        std::sort(mapping, mapping+n_vertices, [&](const int &a, const int &b) { return (unsigned int)priority(a)<(unsigned int)priority(b); });

        // now mapping stores the reorder, all we need to do now is reorder in-place
        int *inverse = new int[n_vertices];
        for (int i=0; i<n_vertices; i++) inverse[mapping[i]] = i;

        // apply the permutation
        unsafe_permute(inverse);

        // free used memory
        delete[] inverse;

        // record that the vertices are now ordered
        is_ordered = true;
    }
}

/**
 * The "safe" permute: apply the mapping, then update is_ordered.
 */
void
Game::permute(int *mapping)
{
    unsafe_permute(mapping);

    // check if ordered...
    is_ordered = true;
    for (int i=1; i<n_vertices; i++) {
        if (_priority[i-1] > _priority[i]) {
            is_ordered = false;
            break;
        }
    }
}

/**
 * Apply permutation (only to arrays)
 */
void
Game::unsafe_permute(int *mapping)
{
    // first update vectors and arrays and the strategies
    int* strat = solution_.strategy_data();
    for (int i=0; i<n_vertices; i++) {
        if (strat[i] != -1) strat[i] = mapping[strat[i]];
    }
    unsigned long len = n_vertices + n_edges;
    for (unsigned long i=0; i<len; i++) {
        if (_outedges[i] != -1) _outedges[i] = mapping[_outedges[i]];
    }
    if (!_inedges.empty()) {
        for (unsigned long i=0; i<len; i++) {
            if (_inedges[i] != -1) _inedges[i] = mapping[_inedges[i]];
        }
    }
    // swap nodes until done
    for (int i=0; i<n_vertices; i++) {
        // this is basically a loop, that swaps mapping[i] and i, until mapping[i] equals i.
        while (mapping[i] != i) {
            int k = mapping[i];
            mapping[i] = mapping[k];
            mapping[k] = k;
            // swap i and k
            std::swap(_priority[i], _priority[k]);
            { bool b = _owner[k]; _owner[k] = _owner[i]; _owner[i] = b; }
            std::swap(_label[i], _label[k]);
            // swap out array
            std::swap(_firstouts[i], _firstouts[k]);
            std::swap(_outcount[i], _outcount[k]);
            // swap in array
            if (!_inedges.empty()) {
                std::swap(_firstins[i], _firstins[k]);
                std::swap(_incount[i], _incount[k]);
            }
            // swap solution
            solution_.swap_vertices(i, k);
        }
    }
}

int
Game::inflate()
{
    assert(is_ordered);

    if (n_vertices == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1;
    for (int i=0; i<n_vertices; i++) {
        const int p_mod_i = _priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        else { prio += 2; d++; }
        _priority[i] = prio;
    }

    return d;
}

int
Game::compress()
{
    assert(is_ordered);

    if (n_vertices == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1;
    for (int i=0; i<n_vertices; i++) {
        const int p_mod_i = _priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        _priority[i] = prio;
    }

    return d;
}

int
Game::renumber()
{
    assert(is_ordered);

    if (n_vertices == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=0; i<n_vertices; i++) {
        const int p_mod_i = _priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        else if (last != _priority[i]) { prio += 2; d++; }
        last = _priority[i];
        _priority[i] = prio;
    }

    return d;
}

void
Game::evenodd()
{
    assert(is_ordered);

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=0; i<n_vertices; i++) {
        const int d = _priority[i]+1;
        _owner[i] = 1-_owner[i];

        const int p_mod_i = d&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) prio += 1;
        else if (last != d) prio += 2;
        last = d;
        _priority[i] = prio;
    }
}

void
Game::minmax()
{
    assert(is_ordered);

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=n_vertices-1; i>=0; i--) {
        const int p_mod_i = _priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) prio += 1;
        else if (last != _priority[i]) prio += 2;
        last = _priority[i];
        _priority[i] = prio;
    }
}

std::unique_ptr<Game>
Game::extract_subgame(const std::vector<int>& selection)
{
    // translate to a bitmask
    bitset sel(n_vertices);
    for (int i : selection) sel[i] = true;
    return extract_subgame(sel);
}

std::unique_ptr<Game>
Game::extract_subgame(const bitset& mask)
{
    std::vector<int> mapping;
    return extract_subgame(mask, mapping);
}

std::unique_ptr<Game>
Game::extract_subgame(const bitset& mask, std::vector<int>& subgame_to_game)
{
    // check if there are any dead ends (not allowed)
    // also count the number of edges in the subgame

    SubgameView view(*this, mask);
    int nv = view.size();
    int ne = 0;
    for (int v=0; v<n_vertices; v++) {
        if (view.contains(v)) {
            bool bad = true;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                if (view.contains(*curedge)) {
                    bad = false;
                    ne++;
                }
            }
            if (bad) {
                std::cerr << "no successor for vertex " << label_vertex(v) << " in extract_subgame!" << std::endl;
                std::cerr << "successors not in subgame:";
                for (auto curedge = outs(v); *curedge != -1; curedge++) {
                    std::cerr << " " << label_vertex(*curedge);
                }
                std::cerr << std::endl;
                abort(); // bad!
            }
        }
    }

    auto res = std::make_unique<Game>(nv, ne);

    // create mapping from game to subgame
    // also count the number of vertices in the subgame
    std::map<int, int> game_to_subgame;
    subgame_to_game.assign(n_vertices, 0);

    int vertices = 0;
    for (int v=0; v<n_vertices; v++) {
        if (!view.contains(v)) continue;

        // update game_to_subgame and subgame_to_game
        int w = vertices++;
        subgame_to_game[w] = v;
        game_to_subgame[v] = w;

        // initialize most stuff (except edges)
        res->init_vertex(w, _priority[v], _owner[v], _label[v]);
    }

    // now add all edges

    for (int w=0; w<nv; w++) {
        int v = subgame_to_game[w];
        res->e_start(w);
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            if (view.contains(*curedge)) res->e_add(w, game_to_subgame[*curedge]);
        }
        res->e_finish();
    }

    // TODO: fix is_ordered?

    return res;
}

Game&
Game::operator=(const Game &other)
{
    Game g(other);
    swap(g);
    return *this;
}

void
Game::swap(Game &other)
{
    std::swap(n_vertices, other.n_vertices);
    std::swap(n_edges, other.n_edges);
    std::swap(_priority, other._priority);
    std::swap(_owner, other._owner);
    std::swap(_label, other._label);
    std::swap(_outvec, other._outvec);
    std::swap(_outedges, other._outedges);
    std::swap(_firstouts, other._firstouts);
    std::swap(_outcount, other._outcount);
    std::swap(_inedges, other._inedges);
    std::swap(_firstins, other._firstins);
    std::swap(_incount, other._incount);
    solution_.swap(other.solution_);
    std::swap(is_ordered, other.is_ordered);
    std::swap(v_allocated, other.v_allocated);
    std::swap(e_allocated, other.e_allocated);
    std::swap(e_size, other.e_size);
}

void
Game::reset_solution()
{
    solution_.reset();
}

void
Game::copy_solution(Game &other)
{
    solution_ = other.solution_;
}

void 
Game::e_sizeup(void)
{
    e_allocated += e_allocated/2 + 1; // +1 so it always grows (e.g. from size 1)
    _outedges.resize(e_allocated);
}

void
Game::v_sizeup(void)
{
    v_allocated += v_allocated/2 + 1; // +1 so it always grows (e.g. from size 1)
    n_vertices = v_allocated;
    // resize grows the arrays, value-initializing the new tail to 0
    _priority.resize(v_allocated);
    _firstouts.resize(v_allocated);
    _outcount.resize(v_allocated);
    _label.resize(v_allocated);
    _owner.resize(v_allocated);
    solution_.resize(v_allocated);
}

void
Game::v_resize(size_t newsize)
{
    while (newsize > v_allocated) v_sizeup();
    n_vertices = newsize;
    _owner.resize(n_vertices);
    solution_.resize(n_vertices);
}

void
Game::e_start(int source)
{
    _firstouts[source] = e_size;
    _outcount[source] = 0;
}

void
Game::e_add(int source, int target)
{
    if (e_size == e_allocated) e_sizeup();
    _outedges[e_size++] = target;
    _outcount[source] += 1;
    n_edges++;
}

void
Game::e_finish(void)
{
    if (e_size == e_allocated) e_sizeup();
    _outedges[e_size++] = -1;
}

void
Game::build_in_array(bool rebuild)
{
    if (!_inedges.empty() and !rebuild) return;

    _inedges.assign(e_size, 0);
    _firstins.assign(n_vertices, 0);

    // set incount of each vertex

    _incount.assign(n_vertices, 0);
    for (int v=0; v<n_vertices; v++) {
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            _incount[to]++;
        }
    }

    unsigned long pos = 0;
    for (int v=0; v<n_vertices; v++) {
        _firstins[v] = pos+_incount[v]; // start at end!!
        _inedges[_firstins[v]] = -1;
        pos += (_incount[v] + 1);
    }

    for (int v=0; v<n_vertices; v++) {
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            int to = *curedge;
            _inedges[--_firstins[to]] = v;
        }
    }
}

}
