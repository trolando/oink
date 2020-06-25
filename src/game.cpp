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
#include <cstring> // memset
#include <iostream>
#include <sys/mman.h> // for mmap, munmap

#include "game.hpp"

#define USE_MMAP 1

using namespace std;

namespace pg {

Game::Game() : _owner(0), solved(0), winner(0)
{
    n_vertices = 0;
    n_edges = 0;
    _priority = NULL;
    _label = NULL;
    _outvec = NULL;
    _outedges = NULL;
    _firstouts = NULL;
    _outcount = NULL;
    _inedges = NULL;
    _firstins = NULL;
    _incount = NULL;
    is_ordered = true;
    v_allocated = 0;
    e_allocated = 0;
    e_size = 0;
    strategy = NULL;
    set_random_seed(rd());
}

Game::~Game()
{
    for (int i=0; i<n_vertices; i++) {
        if (_label[i]) delete _label[i];
    }

    munmap(_priority, sizeof(int[v_allocated]));
    munmap(_label, sizeof(int[v_allocated]));
    munmap(strategy, sizeof(int[v_allocated]));
    munmap(_firstouts, sizeof(int[v_allocated]));
    munmap(_outcount, sizeof(int[v_allocated]));
    munmap(_outedges, sizeof(int[e_allocated]));

    if (_outvec != NULL) {
        delete[] _outvec;
    }

    if (_inedges != NULL) {
        delete[] _inedges;
        delete[] _firstins;
        delete[] _incount;
    }
}

Game::Game(int vcount, int ecount) : _owner(vcount, true), solved(vcount, true), winner(vcount, true)
{
    assert(vcount > 0);
    if (ecount == -1) ecount = size_t(4) * vcount; // reasonable default outdegree

    n_vertices = vcount;
    n_edges = 0;

    v_allocated = vcount;
    e_allocated = vcount+ecount+1;  // extra space for -1
    e_size = 0;

    _priority = (int*)mmap(0, sizeof(int[v_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    _label = (string**)mmap(0, sizeof(string*[v_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    strategy = (int*)mmap(0, sizeof(int[v_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    _firstouts = (int*)mmap(0, sizeof(int[v_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    _outcount = (int*)mmap(0, sizeof(int[v_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    _outedges = (int*)mmap(0, sizeof(int[e_allocated]), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (_priority == (int*)MAP_FAILED) abort();
    if (_label == (string**)MAP_FAILED) abort();
    if (strategy == (int*)MAP_FAILED) abort();
    if (_firstouts == (int*)MAP_FAILED) abort();
    if (_outcount == (int*)MAP_FAILED) abort();
    if (_outedges == (int*)MAP_FAILED) abort();

    _outvec = NULL;
    _inedges = NULL;
    _firstins = NULL;
    _incount = NULL;
    is_ordered = true;

    _outedges[0] = -1;
    e_size++;

    set_random_seed(rd());
}

/**
 * Make a deep clone of the given game <other>.
 * Does not clone the vector representation.
 * Does not clone the <in> array.
 */
Game::Game(const Game& other) : Game(other.n_vertices, other.e_size)
{
    n_edges = other.n_edges;

    memcpy(_priority, other._priority, sizeof(int[n_vertices]));
    _owner = other._owner;
    for (int i=0; i<n_vertices; i++) {
        if (other._label[i]) _label[i] = new std::string(*other._label[i]);
    }

    // clone the edge out ARRAY
    e_size = other.e_size;
    memcpy(_outedges, other._outedges, sizeof(int[e_size]));
    memcpy(_firstouts, other._firstouts, sizeof(int[n_vertices]));
    memcpy(_outcount, other._outcount, sizeof(int[n_vertices]));

    // copy inedges
    if (other._inedges != NULL) {
         size_t len = n_vertices + n_edges;
         _inedges = new int[len];
         _firstins = new int[n_vertices];
         _incount = new int[n_vertices];
         memcpy(_inedges, other._inedges, sizeof(int[len]));
         memcpy(_firstins, other._firstins, sizeof(int[n_vertices]));
         memcpy(_incount, other._incount, sizeof(int[n_vertices]));
    }

    is_ordered = other.is_ordered;

    solved = other.solved;
    winner = other.winner;
    memcpy(strategy, other.strategy, sizeof(int[n_vertices]));

    set_random_seed(rd());
}

/**
 * Helper functions for parsing a PGSolver format parity game.
 */

static void
skip_whitespace(std::streambuf *rd)
{
    // read whitespace
    while (true) {
        int ch;
        if ((ch=rd->sbumpc()) == EOF) return;
        if (ch != ' ' and ch != '\n' and ch != '\t' and ch != '\r') break;
    }
    rd->sungetc();
}

static bool
read_uint64(std::streambuf *rd, uint64_t *res)
{
    uint64_t r = 0;
    int ch;
    if ((ch=rd->sbumpc()) == EOF) return false;
    if (ch < '0' or ch > '9') { rd->sungetc(); return false; }
    while (true) {
        r = (10*r)+(ch-'0');
        if ((ch=rd->sbumpc()) == EOF) break;
        if (ch < '0' or ch > '9') { rd->sungetc(); break; }
    }
    *res = r;
    return true;
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
Game::init_random_game(int n, int maxP, int maxE)
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
    if (this->_label[node]) delete this->_label[node];
    if (label != "") this->_label[node] = new std::string(label);
    else this->_label[node] = 0;
}


/**
 * Vector stuff
 */

void
Game::vec_init(void)
{
    if (_outvec != NULL) delete[] _outvec;
    _outvec = new std::vector<int>[n_vertices];

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
    delete[] _outvec;
    _outvec = NULL;
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
Game::has_edge(int from, int to)
{
    return find_edge(from, to) != -1;
}

int
Game::find_edge(int from, int to)
{
    for (int idx = _firstouts[from]; _outedges[idx] != -1; idx++) {
        if (_outedges[idx] == to) return idx;
    }
    return -1;
}

void
Game::parse_pgsolver(std::istream &inp, bool removeBadLoops)
{
    std::streambuf *rd = inp.rdbuf();

    char buf[64];
    uint64_t n;
    char ch;

    /**
     * Read header line...
     * "parity" <number of nodes> ;
     */

    inp.read(buf, 6);
    if (!inp) throw "expecting parity game specification";
    if (strncmp(buf, "parity", 6) != 0) throw "expecting parity game specification";

    skip_whitespace(rd);
    if (!read_uint64(rd, &n)) throw "missing number of nodes";

    skip_whitespace(rd);
    while ((inp >> ch) and ch != ';') continue;
    if (ch != ';') throw "missing ';'";

    /**
     * Construct game...
     */

    Game res(n+1, true);
    swap(res);

    /**
     * Parse the nodes...
     */

    int node_count = 0; // number of read nodes
    solved.set(); // we use solved to store whether a node has not yet been read

    while (node_count < n_vertices) {
        uint64_t id;
        skip_whitespace(rd);
        if (!read_uint64(rd, &id)) {
            // we expect maybe one more node...
            if (node_count == n_vertices-1) {
                v_resize(node_count);
                // ignore rest, they can be bigger
                break;
            }
            throw "unable to read id";
        }
        if (id >= (unsigned)n_vertices) throw "invalid id";

        if (!solved[id]) throw "duplicate id";
        solved[id] = false;
        node_count++;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw "missing priority";
        _priority[id] = n;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw "missing owner";

        if (n == 0) { /* nothing */ }
        else if (n == 1) { _owner[id] = true; }
        else { throw "invalid owner"; }

        e_start(id);

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            skip_whitespace(rd);
            if (!read_uint64(rd, &n)) throw "missing successor";
            if (n >= (uint64_t)n_vertices) {
                std::cout << "id " << id << " with successor " << n << std::endl;
                throw "invalid successor";
            }

            if (id == n and removeBadLoops and (_owner[id] != (_priority[id]&1))) {
                has_self = true;
            } else {
                // add edge to the vector
                e_add(id, n);
                count++;
            }

            char ch;
            skip_whitespace(rd);
            if (!(inp >> ch)) throw "missing ; to end line";
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                _label[id] = new std::string();
                while (true) {
                    inp >> ch;
                    if (ch == '\"') break;
                    *_label[id] += ch;
                }
                // now read ;
                skip_whitespace(rd);
                if (!(inp >> ch) or ch != ';') throw "missing ; to end line";
            }
            break;
        }

        if (has_self and count == 0) {
            e_add(id, id);
            count++;
        }

        e_finish();
    }

    if (solved.any()) {
        std::cout << "count : " << solved.count() << std::endl;
        throw "missing nodes";
    }

    // check if ordered...
    is_ordered = true;
    for (int i=1; i<n_vertices; i++) {
        if (_priority[i-1] > _priority[i]) {
            is_ordered = false;
            break;
        }
    }

    // ensure strategy empty
    std::fill(strategy, strategy+n_vertices, static_cast<int>(~0));
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
            throw "node index out of bounds";
        }

        if (solved[ident]) throw "node already solved";

        // parse winner
        int w;
        if (!(ss >> w)) throw "missing winner";
        if (w!= 0 && w!= 1) throw "invalid winner";

        // set winner
        solved[ident] = true;
        winner[ident] = w;

        // parse strategy
        if (w == _owner[ident]) {
            int str;
            if (!(ss >> str)) throw "missing strategy for winning node";
            // if (!has_edge(ident, str)) throw "strategy not successor of node";
            // actually this is already checked by the verifier
            strategy[ident] = str;
        } else {
            strategy[ident] = -1;
        }
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
        if (_label[i] != 0 and !_label[i]->empty()) os << " \"" << *_label[i] << "\"";
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
    out << "paritysol " << solved.count() << ";" << std::endl;

    // print solution
    for (int i=0; i<n_vertices; i++) {
        if (solved[i]) {
            out << i << " " << (winner[i] ? "1" : "0");
            if (winner[i] == _owner[i] and strategy[i] != -1) out << " " << strategy[i];
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
    for (int i=0; i<n_vertices; i++) {
        if (strategy[i] != -1) strategy[i] = mapping[strategy[i]];
    }
    unsigned long len = n_vertices + n_edges;
    for (unsigned long i=0; i<len; i++) {
        if (_outedges[i] != -1) _outedges[i] = mapping[_outedges[i]];
    }
    if (_inedges != NULL) {
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
            if (_inedges != NULL) {
                std::swap(_firstins[i], _firstins[k]);
                std::swap(_incount[i], _incount[k]);
            }
            // swap solution
            { bool b = solved[k]; solved[k] = solved[i]; solved[i] = b; }
            { bool b = winner[k]; winner[k] = winner[i]; winner[i] = b; }
            std::swap(strategy[i], strategy[k]);
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

Game *
Game::extract_subgame(std::vector<int> &selection)
{
    // translate to a bitmask
    bitset sel(n_vertices);
    for (int i : selection) sel[i] = true;
    return extract_subgame(sel);
}

Game *
Game::extract_subgame(bitset mask)
{
    // check if there are any dead ends (not allowed)
    // also count the number of edges in the subgame

    int nv = mask.count();
    int ne = 0;
    for (int v=0; v<n_vertices; v++) {
        if (mask[v]) {
            bool bad = true;
            for (auto curedge = outs(v); *curedge != -1; curedge++) {
                if (mask[*curedge]) {
                    bad = false;
                    ne++;
                }
            }
            if (bad) abort(); // bad!
        }
    }

    Game *res = new Game(nv, ne);

    // create mapping from game to subgame
    // our v is their invmap[v]
    // their w is our mapping[w]
    // also count the number of vertices in the subgame
    int *mapping = new int[n_vertices];
    int *invmap = new int[nv];

    int vertices = 0;
    for (int v=0; v<n_vertices; v++) {
        if (!mask[v]) continue;

        // update mapping and invmapping
        int w = vertices++;
        invmap[w] = v;
        mapping[v] = w;

        // initialize most stuff (except edges)
        if (_label[v] != 0) res->init_vertex(w, _priority[v], _owner[v], *_label[v]);
        else res->init_vertex(w, _priority[v], _owner[v], "");
    }

    // now add all edges

    for (int w=0; w<nv; w++) {
        int v = invmap[w];
        res->e_start(v);
        for (auto curedge = outs(v); *curedge != -1; curedge++) {
            if (mask[*curedge]) res->e_add(v, mapping[*curedge]);
        }
        res->e_finish();
    }

    delete[] mapping;
    delete[] invmap;

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
    std::swap(solved, other.solved);
    std::swap(winner, other.winner);
    std::swap(strategy, other.strategy);
    std::swap(is_ordered, other.is_ordered);
    std::swap(v_allocated, other.v_allocated);
    std::swap(e_allocated, other.e_allocated);
    std::swap(e_size, other.e_size);
}

void
Game::reset_solution()
{
    solved.reset();
    winner.reset();
    memset(strategy, -1, sizeof(int[n_vertices]));
}

void
Game::copy_solution(Game &other)
{
    solved = other.solved;
    winner = other.winner;
    memcpy(strategy, other.strategy, sizeof(int[n_vertices]));
}

void 
Game::e_sizeup(void)
{
    size_t old = e_allocated;
    e_allocated += e_allocated/2;
    _outedges = (int*)mremap(_outedges, sizeof(int[old]), sizeof(int[e_allocated]), MREMAP_MAYMOVE);
    if (_outedges == (int*)MAP_FAILED) abort();
}

void
Game::v_sizeup(void)
{
    size_t old = v_allocated;
    v_allocated += v_allocated/2;
    n_vertices = v_allocated;
    _priority = (int*)mremap(_priority, sizeof(int[old]), sizeof(int[v_allocated]), MREMAP_MAYMOVE);
    strategy = (int*)mremap(strategy, sizeof(int[old]), sizeof(int[v_allocated]), MREMAP_MAYMOVE);
    _firstouts = (int*)mremap(_firstouts, sizeof(int[old]), sizeof(int[v_allocated]), MREMAP_MAYMOVE);
    _outcount = (int*)mremap(_outcount, sizeof(int[old]), sizeof(int[v_allocated]), MREMAP_MAYMOVE);
    _label = (string**)mremap(_label, sizeof(string*[old]), sizeof(string*[v_allocated]), MREMAP_MAYMOVE);
    if (_priority == (int*)MAP_FAILED) abort();
    if (_label == (string**)MAP_FAILED) abort();
    if (strategy == (int*)MAP_FAILED) abort();
    if (_firstouts == (int*)MAP_FAILED) abort();
    if (_outcount == (int*)MAP_FAILED) abort();
    _owner.resize(v_allocated);
    solved.resize(v_allocated);
    winner.resize(v_allocated);
}

void
Game::v_resize(size_t newsize)
{
    while (newsize > v_allocated) v_sizeup();
    n_vertices = newsize;
    _owner.resize(n_vertices);
    solved.resize(n_vertices);
    winner.resize(n_vertices);
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
    if (_inedges != NULL) {
        if (rebuild) {
            delete[] _inedges;
            delete[] _firstins;
            delete[] _incount;
        } else {
            return;
        }
    }

    _inedges = new int[e_size];
    _firstins = new int[n_vertices];
    _incount = new int[n_vertices];

    // set incount of each vertex

    memset(_incount, 0, sizeof(int[n_vertices]));
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

std::random_device Game::rd;

}
