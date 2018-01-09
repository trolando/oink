#include <cassert>
#include <iostream>
#include <algorithm>

#include "game.hpp"

using namespace std;

namespace pg {

Game::Game()
{
    n_nodes = 0;
    priority = NULL;
    owner = NULL;
    out = NULL;
    in = NULL;
    label = NULL;
    dominion = NULL;
    strategy = NULL;
    disabled = NULL;
}

Game::Game(int count)
{
    n_nodes = count < 0 ? 0 : count;
    if (n_nodes > 0) {
        priority = new int[n_nodes+1];
        owner = new int[n_nodes+1];
        label = new std::string[n_nodes+1];
        out = new std::vector<int>[n_nodes+1];
        in = new std::vector<int>[n_nodes+1];
        dominion = new int[n_nodes+1];
        strategy = new int[n_nodes+1];
        disabled = new int[n_nodes+1];
    } else {
        priority = NULL;
        owner = NULL;
        out = NULL;
        in = NULL;
        label = NULL;
        dominion = NULL;
        strategy = NULL;
        disabled = NULL;
    }
}

Game::Game(const Game& other)
{
    n_nodes = other.n_nodes;
    if (n_nodes > 0) {
        priority = new int[n_nodes];
        owner = new int[n_nodes];
        label = new std::string[n_nodes];
        out = new std::vector<int>[n_nodes];
        in = new std::vector<int>[n_nodes];
        dominion = new int[n_nodes];
        strategy = new int[n_nodes];
        disabled = new int[n_nodes];
    }
    for (int i=0; i<n_nodes; i++) {
        priority[i] = other.priority[i];
        owner[i] = other.owner[i];
        label[i] = other.label[i];
        dominion[i] = other.dominion[i];
        strategy[i] = other.strategy[i];
        disabled[i] = other.disabled[i];
        out[i] = other.out[i];
        in[i] = other.in[i];
    }
}

Game::~Game()
{
    if (n_nodes > 0) {
        delete[] priority;
        delete[] owner;
        delete[] out;
        delete[] in;
        delete[] label;
        delete[] dominion;
        delete[] strategy;
        delete[] disabled;
        n_nodes = 0;
    }
}

void
Game::initGame(int count)
{
    if (n_nodes > 0) {
        delete[] priority;
        delete[] owner;
        delete[] out;
        delete[] in;
        delete[] label;
        delete[] dominion;
        delete[] strategy;
        delete[] disabled;
    }

    n_nodes = count < 0 ? 0 : count;
    if (n_nodes > 0) {
        priority = new int[n_nodes+1];
        owner = new int[n_nodes+1];
        label = new std::string[n_nodes+1];
        out = new std::vector<int>[n_nodes+1];
        in = new std::vector<int>[n_nodes+1];
        dominion = new int[n_nodes+1];
        strategy = new int[n_nodes+1];
        disabled = new int[n_nodes+1];
    } else {
        priority = NULL;
        owner = NULL;
        out = NULL;
        in = NULL;
        label = NULL;
        dominion = NULL;
        strategy = NULL;
        disabled = NULL;
    }
}

void
Game::initNode(int node, int priority, int owner, std::string label)
{
    assert(node < n_nodes and node >= 0);

    this->priority[node] = priority;
    this->owner[node] = owner;
    this->label[node] = label;
    this->dominion[node] = -1;
    this->strategy[node] = -1;
    this->disabled[node] = 0;
}

bool
Game::addEdge(int from, int to)
{
    assert(from >= 0 and from < n_nodes);
    assert(to >= 0 and to < n_nodes);

    if (std::find(out[from].begin(), out[from].end(), to) == out[from].end()) {
        out[from].push_back(to);
        in[to].push_back(from);
        return true;
    } else {
        return false;
    }
}

size_t
Game::parse_pgsolver(istream &inp)
{
    string line;
    while (getline(inp, line)) {
        stringstream ss(line);
        string token;

        // ignore empty line
        if (!(ss >> token)) continue;

        // process line with "parity"
        if (token != "parity") throw "expecting parity game specification";
        if (!(ss >> n_nodes)) throw "missing number of nodes";
        break;
    }

    initGame(n_nodes);

    for (int i=0; i<=n_nodes; i++) dominion[i] = -1;
    for (int i=0; i<=n_nodes; i++) strategy[i] = -1;
    for (int i=0; i<=n_nodes; i++) disabled[i] = 1; // !

    size_t node_count = 0;
    size_t edge_count = 0;

    while (getline(inp, line)) {
        stringstream ss(line);
        string token;

        // ignore empty line
        if (!(ss >> token)) continue;

        // parse id, priority, owner
        int id;
        try {
            id = stoi(token);
        } catch (const std::invalid_argument) {
            // ignore lines starting with a non number
            continue;
        }
        if (id == n_nodes) n_nodes++; // work-around
        if (id >= n_nodes || id < 0) {
            throw "invalid id";
        }

        if (disabled[id] != 1) {
            throw "duplicate id";
        }

        disabled[id] = 0;
        node_count++;

        if (!(ss >> priority[id])) {
            throw "missing priority";
        }

        if (!(ss >> owner[id])) {
            throw "missing owner";
        }

        if (owner[id] != 0 && owner[id] != 1) {
            throw "invalid owner";
        }

        // parse successors and optional label
        for (;;) {
            int to;
            if (!(ss >> to)) throw "missing successor";

            if (to == n_nodes) n_nodes++; // work-around
            if (to >= n_nodes) throw "invalid successor";

            out[id].push_back(to);
            in[to].push_back(id);
            edge_count++;

            char ch;
            if (!(ss >> ch)) throw "missing ; to end line";

            if (ch == ',') continue;
            if (ch == '\"') getline(ss, label[id], '\"');
            else label[id] = "";
            break;
        }
    }

    if ((int)node_count != n_nodes) throw "missing nodes";
    for (int i=0; i<n_nodes; i++) if (disabled[i]) throw "missing nodes";

    return edge_count;
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

        if (dominion[ident] != -1) throw "node appears twice in solution";

        // parse winner
        int winner;
        if (!(ss >> winner)) throw "missing winner";
        if (winner != 0 && winner != 1) throw "invalid winner";

        // set winner
        dominion[ident] = winner;

        // parse strategy
        if (winner == owner[ident]) {
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
    // count size of solution
    int size = 0;
    for (int i=0; i<n_nodes; i++) if (dominion[i] != -1) size++;

    // print banner
    out << "paritysol " << size << ";" << endl;

    // print solution
    for (int i=0; i<n_nodes; i++) {
        if (dominion[i] != -1) {
            out << i << " " << dominion[i];
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
            std::swap(owner[i], owner[k]);
            std::swap(in[i], in[k]);
            std::swap(out[i], out[k]);
            std::swap(label[i], label[k]);
            std::swap(dominion[i], dominion[k]);
            std::swap(strategy[i], strategy[k]);
            std::swap(disabled[i], disabled[k]);
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

size_t
Game::edgecount()
{
    size_t res = 0;
    for (int i=0; i<n_nodes; i++) res += out[i].size();
    return res;
}

bool
Game::solved()
{
    for (int i=0; i<n_nodes; i++) {
        if (dominion[i] == -1) return false;
    }
    return true;
}

int
Game::countUnsolved()
{
    int count = 0;
    for (int i=0; i<n_nodes; i++) {
        if (dominion[i] == -1) count++;
    }
    return count;
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
        res->dominion[i] = dominion[k];
        res->strategy[i] = strategy[k];
        res->disabled[i] = disabled[k];
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

void 
Game::restrict(std::vector<int> &selection)
{
    for (int i=0; i<n_nodes; i++) disabled[i] = 1;
    for (int i : selection) disabled[i] = 0;
}

Game&
Game::operator=(const Game &other)
{
    initGame(other.n_nodes);
    for (int i=0; i<n_nodes; i++) {
        priority[i] = other.priority[i];
        owner[i] = other.owner[i];
        label[i] = other.label[i];
        dominion[i] = other.dominion[i];
        strategy[i] = other.strategy[i];
        disabled[i] = other.disabled[i];
        out[i] = other.out[i];
        in[i] = other.in[i];
    }
    return *this;
}

void
Game::reset()
{
    for (int i=0; i<n_nodes; i++) {
        dominion[i] = -1;
        strategy[i] = -1;
        disabled[i] = 0;
    }
}

}
