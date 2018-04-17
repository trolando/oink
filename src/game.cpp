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

#include <cassert>
#include <cstring> // memset
#include <iostream>
#include <algorithm>

#include "game.hpp"

using namespace std;

namespace pg {

Game::Game() :
    n_nodes(0), n_edges(0),
    priority(new int[n_nodes]), owner(n_nodes), label(new std::string[n_nodes]),
    out(new std::vector<int>[n_nodes]), in(new std::vector<int>[n_nodes]),
    solved(n_nodes), winner(n_nodes), strategy(new int[n_nodes]), reindexed(false)
{
}

Game::Game(int count) :
    n_nodes(count), n_edges(0),
    priority(new int[n_nodes]), owner(n_nodes), label(new std::string[n_nodes]),
    out(new std::vector<int>[n_nodes]), in(new std::vector<int>[n_nodes]),
    solved(n_nodes), winner(n_nodes), strategy(new int[n_nodes]), reindexed(false)
{
    assert(count > 0);
    memset(strategy, -1, sizeof(int[n_nodes]));
}

Game::Game(const Game& other) : Game(other.n_nodes)
{
    n_edges = other.n_edges;
    memcpy(priority, other.priority, sizeof(int[n_nodes]));
    owner = other.owner;
    for (int i=0; i<n_nodes; i++) label[i] = other.label[i];
    for (int i=0; i<n_nodes; i++) out[i] = other.out[i];
    for (int i=0; i<n_nodes; i++) in[i] = other.in[i];

    reindexed = other.reindexed;

    solved = other.solved;
    winner = other.winner;
    memcpy(strategy, other.strategy, sizeof(int[n_nodes]));
}

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

Game::Game(istream &inp)
{
    std::streambuf *rd = inp.rdbuf();

    char buf[64];
    uint64_t n;
    char ch;

    /**
     * Read header line...
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

    n_nodes = n+1; // plus 1, in case this parity game encodes "max id" instead of "n_nodes"
    n_edges = 0;

    priority = new int[n_nodes];
    owner.resize(n_nodes);
    label = new std::string[n_nodes];
    out = new std::vector<int>[n_nodes];
    in = new std::vector<int>[n_nodes];

    reindexed = false;

    solved.resize(n_nodes);
    winner.resize(n_nodes);
    strategy = new int[n_nodes];

    memset(strategy, -1, sizeof(int[n_nodes]));
    solved.set(); // we use solved for temporary storage

    /**
     * Read nodes...
     */

    int node_count = 0; // number of read nodes

    while (node_count < n_nodes) {
        uint64_t id;
        skip_whitespace(rd);
        if (!read_uint64(rd, &id)) {
            if (node_count == n_nodes-1) {
                n_nodes--;
                owner.resize(n_nodes);
                solved.resize(n_nodes);
                winner.resize(n_nodes);
                // ignore rest, they can be bigger
                break;
            }
            throw "unable to read id";
        }
        if (id >= (unsigned)n_nodes) throw "invalid id";

        if (!solved[id]) throw "duplicate id";
        solved[id] = false;
        node_count++;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw "missing priority";
        priority[id] = n;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw "missing owner";

        if (n == 0) { /* nothing */ }
        else if (n == 1) { owner[id] = true; }
        else { throw "invalid owner"; }

        // parse successors and optional label
        for (;;) {
            skip_whitespace(rd);
            if (!read_uint64(rd, &n)) throw "missing successor";
            if (n >= (uint64_t)n_nodes) {
                std::cout << "id " << id << " with successor " << n << std::endl;
                throw "invalid successor";}

            out[id].push_back(n);
            in[n].push_back(id);
            n_edges++;

            char ch;
            skip_whitespace(rd);
            if (!(inp >> ch)) throw "missing ; to end line";
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                while (true) {
                    inp >> ch;
                    if (ch == '\"') break;
                    label[id] += ch;
                }
                // now read ;
                skip_whitespace(rd);
                if (!(inp >> ch) or ch != ';') throw "missing ; to end line";
            }
            else label[id] = "";
            break;
        }
    }

    if (solved.any()) throw "missing nodes";
}

Game::~Game()
{
    delete[] priority;
    delete[] out;
    delete[] in;
    delete[] label;
    delete[] strategy;
}

void
Game::initGame(int count)
{
    Game g(count);
    swap(g);
}

void
Game::initNode(int node, int priority, int owner, std::string label)
{
    assert(node >= 0 and node < n_nodes);
    assert(owner == 0 or owner == 1);

    this->priority[node] = priority;
    this->owner[node] = owner;
    this->label[node] = label;
}

bool
Game::addEdge(int from, int to)
{
    assert(from >= 0 and from < n_nodes);
    assert(to >= 0 and to < n_nodes);

    if (std::find(out[from].begin(), out[from].end(), to) == out[from].end()) {
        out[from].push_back(to);
        in[to].push_back(from);
        n_edges++;
        return true;
    } else {
        return false;
    }
}

bool
Game::removeEdge(int from, int to)
{
    assert(from >= 0 and from < n_nodes);
    assert(to >= 0 and to < n_nodes);

    auto &o = out[from];
    auto &i = in[to];

    auto pre = o.size();
    o.erase(std::remove(o.begin(), o.end(), to), o.end());
    if (o.size() == pre) return false;
    i.erase(std::remove(i.begin(), i.end(), from), i.end());
    return true;
}

void
Game::parse_pgsolver(std::istream &in)
{
    Game g(in);
    swap(g);
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
        if (ident < 0 || ident >= n_nodes) {
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
        if (w == owner[ident]) {
            int str;
            if (!(ss >> str)) throw "missing strategy for winning node";

            bool done = false;
            for (auto o : out[ident]) {
                if (o == str) {
                    strategy[ident] = o;
                    done = true;
                    break;
                }
            }

            if (!done) throw "strategy not successor of node";
        }
    }
}

void
Game::write_pgsolver(std::ostream &os)
{
    // print banner
    os << "parity " << n_nodes << ";" << endl;

    // print nodes
    for (int i=0; i<n_nodes; i++) {
        os << i << " " << priority[i] << " " << owner[i];
        for (unsigned j=0; j<out[i].size(); j++) {
            if (j == 0) os << " ";
            else os << ",";
            os << out[i][j];
        }
        if (label[i] != "") os << " \"" << label[i] << "\"";
        os << ";" << endl;
    }
}

void
Game::write_dot(std::ostream &out)
{
    out << "digraph G {" << endl;
    for (int i=0; i<n_nodes; i++) {
        out << i << " [ shape=\"" << (owner[i] ? "box" : "diamond")
            << "\", label=\"" << priority[i] << "\"];" << endl;
        for (auto j : this->out[i]) {
            out << i << " -> " << j << ";" << endl;
        }
    }
    out << "}" << endl;
}

void
Game::write_sol(std::ostream &out)
{
    // print banner
    out << "paritysol " << solved.count() << ";" << endl;

    // print solution
    for (int i=0; i<n_nodes; i++) {
        if (solved[i]) {
            out << i << " " << (winner[i] ? "1" : "0");
            if (strategy[i] != -1) out << " " << strategy[i];
            out << ";" << endl;
        }
    }
}

void
Game::reindex(int *mapping)
{
    int *index = mapping == NULL ? new int[n_nodes] : mapping;
    for (int i=0; i<n_nodes; i++) index[i] = i;
    std::sort(index, index+n_nodes, [&](const int &a, const int &b) { return (unsigned int)priority[a]<(unsigned int)priority[b]; });
    // now index stores the reorder, all we need to do now is reorder in-place
    int *inv = new int[n_nodes];
    for (int i=0; i<n_nodes; i++) inv[index[i]] = i;
    if (index != mapping) delete[] index;
    // apply the permutation
    permute(inv);
    delete[] inv;
    reindexed = true;
}

void
Game::reindex_once(void)
{
    if (!reindexed) reindex(NULL);
}

void
Game::permute(int *mapping)
{
    // first apply reorder to "in" and "out"
    for (int i=0; i<n_nodes; i++) {
        for (auto it = in[i].begin(); it != in[i].end(); it++) *it = mapping[*it];
        for (auto it = out[i].begin(); it != out[i].end(); it++) *it = mapping[*it];
        if (strategy[i] != -1) strategy[i] = mapping[strategy[i]];
    }
    // now swap nodes until done
    for (int i=0; i<n_nodes; i++) {
        for (;;) {
            int k = mapping[i];
            if (k == i) break;
            // swap i and mapping[i]
            std::swap(priority[i], priority[k]);
            { bool b = owner[k]; owner[k] = owner[i]; owner[i] = b; }
            std::swap(in[i], in[k]);
            std::swap(out[i], out[k]);
            std::swap(label[i], label[k]);
            { bool b = solved[k]; solved[k] = solved[i]; solved[i] = b; }
            { bool b = winner[k]; winner[k] = winner[i]; winner[i] = b; }
            std::swap(strategy[i], strategy[k]);
            mapping[i] = mapping[k];
            mapping[k] = k;
        }
    }
}

int
Game::inflate()
{
    // assumption: reindex has been called first!!!
    if (n_nodes == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1;
    for (int i=0; i<n_nodes; i++) {
        const int p_mod_i = priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        else { prio += 2; d++; }
        priority[i] = prio;
    }

    return d;
}

int
Game::compress()
{
    // assumption: reindex has been called first!!!
    if (n_nodes == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1;
    for (int i=0; i<n_nodes; i++) {
        const int p_mod_i = priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        priority[i] = prio;
    }

    return d;
}

int
Game::renumber()
{
    // assumption: reindex has been called first!!!
    if (n_nodes == 0) return 0;
    int d = 1;

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=0; i<n_nodes; i++) {
        const int p_mod_i = priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) { prio += 1; d++; }
        else if (last != priority[i]) { prio += 2; d++; }
        last = priority[i];
        priority[i] = prio;
    }

    return d;
}

void
Game::evenodd()
{
    // assumption: reindex has been called first!!!

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=0; i<n_nodes; i++) {
        const int d = priority[i]+1;
        owner[i] = 1-owner[i];

        const int p_mod_i = d&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) prio += 1;
        else if (last != d) prio += 2;
        last = d;
        priority[i] = prio;
    }
}

void
Game::minmax()
{
    // assumption: reindex has been called first!!!

    // reassign priorities and reindex nodes
    int prio = -1, last = -1;
    for (int i=n_nodes-1; i>=0; i--) {
        const int p_mod_i = priority[i]&1;
        if (prio == -1) prio = p_mod_i;
        else if (p_mod_i != prio%2) prio += 1;
        else if (last != priority[i]) prio += 2;
        last = priority[i];
        priority[i] = prio;
    }
}

Game *
Game::extract_subgame(std::vector<int> &selection, int *mapping)
{
    // order selection
    sort(selection.begin(), selection.end());

    // analyse subgame to remove vertices without outgoing edges
    int *inv = new int[n_nodes];
    for (int i=0; i<n_nodes; i++) inv[i] = -1;

    // set inv[i] := 0 for every selected vertex
    std::vector<int> q;
    for (int i : selection) {
        assert(i >= 0 and i < n_nodes);
        assert(inv[i] == -1);
        inv[i] = 0;
    }

    // count number of successors in selection for every vertex
    for (int i : selection) {
        int count = 0;
        for (int j : out[i]) if (inv[j] != -1) count++;
        inv[i] = count;
        if (count == 0) q.push_back(i);
    }

    // handle dead ends (turn into -1)
    while (!q.empty()) {
        int i = q.back();
        q.pop_back();
        inv[i] = -1;
        for (int j : in[i]) {
            if (inv[j] != -1 and --inv[j] == 0) q.push_back(j);
        }
    }

    // temporary helper mapping
    int *map = mapping != NULL ? mapping : new int[selection.size()];
    int counter = 0;
    for (int i : selection) {
        if (inv[i] == -1) continue;
        inv[i] = counter;
        map[counter++] = i;
    }

    // construct subgame
    Game* res = new Game(counter);

    // create nodes of subgame
    for (int i=0; i<counter; i++) {
        int k = map[i];
        res->priority[i] = priority[k];
        res->owner[i] = owner[k];
        res->label[i] = label[k];
        res->solved[i] = solved[k];
        res->winner[i] = winner[k];
        res->strategy[i] = strategy[k];
        for (auto j : in[k]) {
            if (inv[j] != -1) res->in[i].push_back(inv[j]);
        }
        for (auto j : out[k]) {
            if (inv[j] != -1) res->out[i].push_back(inv[j]);
        }
    }

    if (mapping == NULL) delete[] map;
    delete[] inv;
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
    std::swap(n_nodes, other.n_nodes);
    std::swap(n_edges, other.n_edges);
    std::swap(priority, other.priority);
    std::swap(owner, other.owner);
    std::swap(label, other.label);
    std::swap(out, other.out);
    std::swap(in, other.in);
    std::swap(solved, other.solved);
    std::swap(winner, other.winner);
    std::swap(strategy, other.strategy);
}

void
Game::reset()
{
    solved.reset();
    winner.reset();
    memset(strategy, -1, sizeof(int[n_nodes]));
}

}
