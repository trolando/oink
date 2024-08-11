/*
 * Copyright 2024 Tom van Dijk
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

#include <stdexcept>
#include <map>
#include <set>
#include <vector>
#include <boost/container/flat_map.hpp>
#include "oink/pgparser.hpp"
#include "printf.hpp"

namespace pg {

/**
 * Helper functions for parsing a PGSolver format parity game.
 */

static void
skip_whitespace(std::streambuf *rd)
{
    // read whitespace
    while (true) {
        int ch;
        if ((ch=rd->sgetc()) == EOF) return;
        if (ch != ' ' and ch != '\n' and ch != '\t' and ch != '\r') break;
        rd->sbumpc();
    }
}

static void
skip_line(std::streambuf *rd)
{
    while (true) {
        int ch;
        if ((ch=rd->sbumpc()) == EOF) return;
        if (ch == '\n' or ch == '\r') return;
    }
}

static bool
read_uint64(std::streambuf *rd, uint64_t *res)
{
    uint64_t r = 0;
    int ch;
    if ((ch=rd->sgetc()) == EOF) return false;
    if (ch < '0' or ch > '9') { return false; }
    while (true) {
        rd->sbumpc();
        r = (10*r)+(ch-'0');
        if ((ch=rd->sgetc()) == EOF) break;
        if (ch < '0' or ch > '9') { break; }
    }
    *res = r;
    return true;
}

Game
PGParser::parse_pgsolver(std::istream &inp, bool removeBadLoops)
{
    std::streambuf *rd = inp.rdbuf();

    char buf[64];
    uint64_t n;

    /**
     * Read header line...
     * "parity" <number of nodes> ;
     */

    inp.read(buf, 6);
    if (!inp) throw std::runtime_error("expecting parity game specification");
    if (strncmp(buf, "parity", 6) != 0) throw std::runtime_error("expecting parity game specification");

    skip_whitespace(rd);
    if (!read_uint64(rd, &n)) throw std::runtime_error("missing number of nodes");

    skip_whitespace(rd);
    if (rd->sbumpc() != ';') throw std::runtime_error("missing ';'");

    // check if next token is 'start'
    {
        skip_whitespace(rd);
        int ch = rd->sgetc();
        if (ch == 's') skip_line(rd);
    }

    /**
     * Construct game...
     */

    Game res(n+1, true);

    /**
     * Parse the nodes...
     */

    int node_count = 0; // number of read nodes
    res.solved.set(); // we use solved to store whether a node has not yet been read

    while (node_count < res.n_vertices) {
        uint64_t id;
        skip_whitespace(rd);
        if (!read_uint64(rd, &id)) {
            // we expect maybe one more node...
            if (node_count == res.n_vertices-1) {
                res.v_resize(node_count);
                // ignore rest, they can be bigger
                break;
            }
            throw std::runtime_error("unable to read id");
        }
        if (id >= (unsigned)res.n_vertices) throw std::runtime_error("invalid id");

        if (!res.solved[id]) throw std::runtime_error("duplicate id");
        res.solved[id] = false;
        node_count++;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing priority");
        if (n > INT_MAX) throw std::runtime_error("priority too high"); // don't be ridiculous
        res._priority[id] = n;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing owner");

        if (n == 0) { /* nothing */ }
        else if (n == 1) { res._owner[id] = true; }
        else { throw std::runtime_error("invalid owner"); }

        res._label[id] = 0;
        res.e_start(id);

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            skip_whitespace(rd);
            if (!read_uint64(rd, &n)) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)res.n_vertices) {
                std::cout << "id " << id << " with successor " << n << std::endl;
                throw std::runtime_error("invalid successor");
            }

            if (id == n and removeBadLoops and (res._owner[id] != (res._priority[id]&1))) {
                has_self = true;
            } else {
                // add edge to the vector
                res.e_add(id, n);
                count++;
            }

            char ch;
            skip_whitespace(rd);
            if (!(inp >> ch)) throw std::runtime_error("missing ; to end line");
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                res._label[id] = new std::string();
                while (true) {
                    inp >> ch;
                    if (ch == '\"') break;
                    *res._label[id] += ch;
                }
                // now read ;
                skip_whitespace(rd);
                if (!(inp >> ch) or ch != ';') throw std::runtime_error("missing ; to end line");
            }
            break;
        }

        if (has_self and count == 0) {
            res.e_add(id, id);
            count++;
        }

        res.e_finish();
    }

    if (res.solved.any()) {
        std::cout << "count : " << res.solved.count() << std::endl;
        throw std::runtime_error("missing nodes");
    }

    // check if ordered...
    res.is_ordered = true;
    for (int i=1; i<res.n_vertices; i++) {
        if (res._priority[i-1] > res._priority[i]) {
            res.is_ordered = false;
            break;
        }
    }

    // ensure strategy empty
    std::fill(res.strategy, res.strategy+res.n_vertices, static_cast<int>(~0));

    return res;
}

/*
 * this version does not need friend class access, but it can be 10% slower (like the _compress version)
 *
Game
PGParser::parse_pgsolver(std::istream &inp, bool removeBadLoops)
{
    std::streambuf *rd = inp.rdbuf();

    char buf[64];
    uint64_t n;

    inp.read(buf, 6);
    if (!inp) throw std::runtime_error("expecting parity game specification");
    if (strncmp(buf, "parity", 6) != 0) throw std::runtime_error("expecting parity game specification");

    skip_whitespace(rd);
    if (!read_uint64(rd, &n)) throw std::runtime_error("missing number of nodes");

    skip_whitespace(rd);
    if (rd->sbumpc() != ';') throw std::runtime_error("missing ';'");

    // check if next token is 'start'
    // that means "start" <initial vertex id> ";"
    // which we ignore...
    {
        skip_whitespace(rd);
        int ch = rd->sgetc();
        if (ch == 's') skip_line(rd);
    }

    // the given number is either the number of vertices, or the highest vertex
    // so, we expect n or n+1 nodes
    unsigned long n_vertices = n+1;

    // initialize variables
    bitset seen(n_vertices);
    std::vector<int> priority(n_vertices);
    bitset owner(n_vertices);
    std::vector<std::vector<int>> edges(n_vertices);
    std::vector<std::string*> labels(n_vertices);

    size_t node_count = 0; // number of read nodes
    size_t edge_count = 0; // number of read edges

    while (node_count < n_vertices) {
        uint64_t id;
        skip_whitespace(rd);
        // TODO: assume size is n, not n+1
        if (!read_uint64(rd, &id)) {
            // we expect maybe one more node...
            if (node_count == n_vertices - 1 && !seen[n_vertices - 1]) {
                n_vertices -= 1;
                owner.resize(n_vertices);
                edges.resize(n_vertices);
                labels.resize(n_vertices);
                // ignore rest, they can be bigger
                break;
            } else {
                throw std::runtime_error("unable to read id");
            }
        }
        if (id >= n_vertices) throw std::runtime_error("invalid id (too high)");
        if (seen[id]) throw std::runtime_error("duplicate id");
        seen[id] = true;
        node_count++;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing priority");
        if (n > INT_MAX) throw std::runtime_error("priority too high"); // don't be ridiculous
        priority[id] = n;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing owner");

        if (n == 1) owner[id] = true;
        else if (n != 0) throw std::runtime_error("invalid owner (must be 0 or 1)");

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            skip_whitespace(rd);
            if (!read_uint64(rd, &n)) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)n_vertices) {
                std::stringstream err;
                err << "invalid successor for id " << id;
                throw std::runtime_error(err.str());
            }

            if (id == n and removeBadLoops and (owner[id] != (priority[id]&1))) {
                has_self = true;
            } else {
                // add edge to the vector
                edges[id].push_back((int) n);
                count++;
            }

            char ch;
            skip_whitespace(rd);
            if (!(inp >> ch)) throw std::runtime_error("missing ; to end line");
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                std::string label;
                while (true) {
                    inp >> ch;
                    if (ch == '\"') break;
                    label += ch;
                }
                // now read ;
                skip_whitespace(rd);
                if (!(inp >> ch) or ch != ';') throw std::runtime_error("missing ; to end line");
                labels[id] = new std::string(label);
            }
            break;
        }

        if (has_self and count == 0) {
            // we must keep it
            edges[id].push_back((int)id);
            count++;
        }
        edge_count += count;
    }

    if (!seen.all()) {
        std::stringstream err;
        err << "missing nodes, seen only " << seen.count();
        throw std::runtime_error(err.str());
    }

    return { node_count, edge_count, priority, owner, edges, labels };
}
*/

Game
PGParser::parse_pgsolver_renumber(std::istream &in, bool removeBadLoops)
{
    std::streambuf *rd = in.rdbuf();

    char buf[64];
    uint64_t n;

    /**
     * Read header line... one of the following:
     * "parity" <number of vertices> ;
     * "parity" <highest vertex id> ;
     */

    in.read(buf, 6);
    if (!in) throw std::runtime_error("expecting parity game specification");
    if (strncmp(buf, "parity", 6) != 0) throw std::runtime_error("expecting parity game specification");

    skip_whitespace(rd);
    if (!read_uint64(rd, &n)) throw std::runtime_error("missing number of nodes");

    skip_whitespace(rd);
    if (rd->sbumpc() != ';') throw std::runtime_error("missing ';'");

    // check if next token is 'start'
    // that means "start" <initial vertex id> ";"
    // which we ignore...
    {
        skip_whitespace(rd);
        int ch = rd->sgetc();
        if (ch == 's') skip_line(rd);
    }

    /**
     * Construct game...
     */

    // the given number is either the number of vertices, or the highest vertex
    // so, we expect n or n+1 nodes
    unsigned long n_vertices = n+1;

    // initialize variables
    bitset seen(n_vertices);
    std::vector<uint64_t> priority(n_vertices);
    bitset owner(n_vertices);
    std::vector<std::vector<int>> edges(n_vertices);
    std::vector<std::string*> labels(n_vertices);

    /**
     * Parse the nodes...
     */

    size_t node_count = 0; // number of read nodes
    size_t edge_count = 0; // number of read edges

    while (node_count < n_vertices) {
        uint64_t id;
        skip_whitespace(rd);
        // TODO: assume size is n, not n+1
        if (!read_uint64(rd, &id)) {
            // we expect maybe one more node...
            if (node_count == n_vertices - 1 && !seen[n_vertices - 1]) {
                n_vertices -= 1;
                seen.resize(n_vertices);
                owner.resize(n_vertices);
                edges.resize(n_vertices);
                labels.resize(n_vertices);
                // ignore rest, they can be bigger
                break;
            } else {
                throw std::runtime_error("unable to read id");
            }
        }
        if (id >= n_vertices) throw std::runtime_error("invalid id (too high)");
        if (seen[id]) throw std::runtime_error("duplicate id");
        seen[id] = true;
        node_count++;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing priority");
        priority[id] = n;

        skip_whitespace(rd);
        if (!read_uint64(rd, &n)) throw std::runtime_error("missing owner");

        if (n == 1) owner[id] = true;
        else if (n != 0) throw std::runtime_error("invalid owner (must be 0 or 1)");

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            skip_whitespace(rd);
            if (!read_uint64(rd, &n)) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)n_vertices) {
                std::stringstream err;
                err << "invalid successor for id " << id;
                throw std::runtime_error(err.str());
            }

            if (id == n and removeBadLoops and (owner[id] != (priority[id]&1))) {
                has_self = true;
            } else {
                // add edge to the vector
                edges[id].push_back((int) n);
                count++;
            }

            char ch;
            skip_whitespace(rd);
            if (!(in >> ch)) throw std::runtime_error("missing ; to end line");
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                std::string label;
                while (true) {
                    in >> ch;
                    if (ch == '\"') break;
                    label += ch;
                }
                // now read ;
                skip_whitespace(rd);
                if (!(in >> ch) or ch != ';') throw std::runtime_error("missing ; to end line");
                labels[id] = new std::string(label);
            }
            break;
        }

        if (has_self and count == 0) {
            // we must keep it
            edges[id].push_back((int)id);
            count++;
        }
        edge_count += count;
    }

    if (!seen.all()) {
        std::stringstream err;
        err << "missing nodes, seen only " << seen.count();
        throw std::runtime_error(err.str());
    }

    // we now need to fix the priorities first...
    boost::container::flat_map<uint64_t, int> map;
    for (const auto& entry : priority) {
        map[entry] = -1;
    }

    int counter = 0;
    uint64_t previous = 0;
    for (auto& pair : map) {
        if (previous != pair.first) {
            counter++;
            if ((counter&1) != (pair.first&1)) counter++;
            previous = pair.first;
        }
        pair.second = counter;
    }

    std::vector<int> int_priorities(node_count);
    for (unsigned v=0; v<node_count; v++) {
        int_priorities[v] = map[priority[v]];
    }

    return { node_count, edge_count, int_priorities, owner, edges, labels };
}

}
