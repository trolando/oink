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

#include <algorithm>
#include <climits>
#include <cstdio>
#include <istream>
#include <stdexcept>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "oink/pgparser.hpp"
#include "oink/game_builder.hpp"
#include "oink/player.hpp"
#include "printf.hpp"

namespace pg {

/**
 * Helper for parsing a PGSolver format parity game: reads the stream in large blocks and
 * scans them with a pointer, which is much faster than reading characters through the
 * stream (formatted extraction such as "in >> ch" costs a sentry and a locale per call).
 * It reads ahead, so the stream is left after the game, at an unspecified position.
 */
namespace {

class Reader
{
public:
    explicit Reader(std::istream &in) : rd(in.rdbuf()), buf(1 << 20) { }

    /** The next character, or EOF. */
    inline int peek()
    {
        if (p == end and !fill()) return EOF;
        return (unsigned char)*p;
    }

    /** Read the next character, or EOF. */
    inline int get()
    {
        if (p == end and !fill()) return EOF;
        return (unsigned char)*p++;
    }

    inline void skip_whitespace()
    {
        for (;;) {
            if (p == end and !fill()) return;
            const char ch = *p;
            if (ch != ' ' and ch != '\n' and ch != '\t' and ch != '\r') return;
            p++;
        }
    }

    inline void skip_line()
    {
        for (;;) {
            const int ch = get();
            if (ch == EOF or ch == '\n' or ch == '\r') return;
        }
    }

    inline bool read_uint64(uint64_t *res)
    {
        int ch = peek();
        if (ch < '0' or ch > '9') return false;
        uint64_t r = 0;
        do {
            r = (10*r)+(ch-'0');
            p++;
            ch = peek();
        } while (ch >= '0' and ch <= '9');
        *res = r;
        return true;
    }

    /** Skip whitespace, then read the next character (EOF at the end). */
    inline int next_token()
    {
        skip_whitespace();
        return get();
    }

    /** Read a quoted label (after the opening quote); whitespace is left out, as before. */
    std::string read_label()
    {
        std::string label;
        for (;;) {
            const int ch = next_token();
            if (ch == EOF) throw std::runtime_error("unterminated label");
            if (ch == '\"') return label;
            label += (char)ch;
        }
    }

private:
    bool fill()
    {
        const std::streamsize k = rd->sgetn(buf.data(), (std::streamsize)buf.size());
        p = buf.data();
        end = p + (k > 0 ? k : 0);
        return k > 0;
    }

    std::streambuf *rd;
    std::vector<char> buf;
    const char *p = nullptr;
    const char *end = nullptr;
};

/** Read the "parity" keyword that starts a game. */
void
read_parity_keyword(Reader &rd)
{
    for (const char *kw = "parity"; *kw; kw++) {
        if (rd.get() != *kw) throw std::runtime_error("expecting parity game specification");
    }
}

}

Game
PGParser::parse_pgsolver(std::istream &inp, bool removeBadLoops)
{
    Reader rd(inp);

    uint64_t n;

    /**
     * Read header line...
     * "parity" <number of nodes> ;
     */

    read_parity_keyword(rd);

    rd.skip_whitespace();
    if (!rd.read_uint64(&n)) throw std::runtime_error("missing number of nodes");

    rd.skip_whitespace();
    if (rd.get() != ';') throw std::runtime_error("missing ';'");

    // check if next token is 'start'
    {
        rd.skip_whitespace();
        int ch = rd.peek();
        if (ch == 's') rd.skip_line();
    }

    /**
     * Construct game...
     */

    Game res(n+1, true);

    /**
     * Parse the nodes...
     */

    int node_count = 0; // number of read nodes
    bitset seen(res.n_vertices); // tracks which nodes have been read

    while (node_count < res.n_vertices) {
        uint64_t id;
        rd.skip_whitespace();
        if (!rd.read_uint64(&id)) {
            // we expect maybe one more node...
            if (node_count == res.n_vertices-1) {
                res.v_resize(node_count);
                seen.resize(node_count);
                // ignore rest, they can be bigger
                break;
            }
            throw std::runtime_error("unable to read id");
        }
        if (id >= (unsigned)res.n_vertices) throw std::runtime_error("invalid id");

        if (seen[id]) throw std::runtime_error("duplicate id");
        seen[id] = true;
        node_count++;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n)) throw std::runtime_error("missing priority");
        if (n > INT_MAX) throw std::runtime_error("priority too high"); // don't be ridiculous
        res._priority[id] = n;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n)) throw std::runtime_error("missing owner");

        if (n == 0) { /* nothing */ }
        else if (n == 1) { res._owner[id] = true; }
        else { throw std::runtime_error("invalid owner"); }

        res.e_start(id);

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            rd.skip_whitespace();
            if (!rd.read_uint64(&n)) throw std::runtime_error("missing successor");
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

            const int ch = rd.next_token();
            if (ch == EOF) throw std::runtime_error("missing ; to end line");
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                res._label[id] += rd.read_label();
                // now read ;
                if (rd.next_token() != ';') throw std::runtime_error("missing ; to end line");
            }
            break;
        }

        if (has_self and count == 0) {
            res.e_add(id, id);
            count++;
        }

        res.e_finish();
    }

    if (!seen.all()) {
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
    Reader rd(in);

    uint64_t n;

    /**
     * Read header line... one of the following:
     * "parity" <number of vertices> ;
     * "parity" <highest vertex id> ;
     */

    read_parity_keyword(rd);

    rd.skip_whitespace();
    if (!rd.read_uint64(&n)) throw std::runtime_error("missing number of nodes");

    rd.skip_whitespace();
    if (rd.get() != ';') throw std::runtime_error("missing ';'");

    // check if next token is 'start'
    // that means "start" <initial vertex id> ";"
    // which we ignore...
    {
        rd.skip_whitespace();
        int ch = rd.peek();
        if (ch == 's') rd.skip_line();
    }

    /**
     * Construct game...
     */

    // the given number is either the number of vertices, or the highest vertex
    // so, we expect n or n+1 nodes
    unsigned long n_vertices = n+1;

    // accumulate vertices and edges directly into a builder; priorities are kept
    // as uint64 here so they can be renumbered (compressed) before building
    GameBuilder builder(n_vertices);
    bitset seen(n_vertices);
    std::vector<uint64_t> priority(n_vertices);

    /**
     * Parse the nodes...
     */

    size_t node_count = 0; // number of read nodes

    while (node_count < n_vertices) {
        uint64_t id;
        rd.skip_whitespace();
        // TODO: assume size is n, not n+1
        if (!rd.read_uint64(&id)) {
            // we expect maybe one more node...
            if (node_count == n_vertices - 1 && !seen[n_vertices - 1]) {
                n_vertices -= 1;
                seen.resize(n_vertices);
                builder.truncate(n_vertices);
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

        rd.skip_whitespace();
        if (!rd.read_uint64(&n)) throw std::runtime_error("missing priority");
        priority[id] = n;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n)) throw std::runtime_error("missing owner");
        if (n != 0 and n != 1) throw std::runtime_error("invalid owner (must be 0 or 1)");
        const int cur_owner = (int) n;
        builder.set_owner((int)id, player_from_int(cur_owner));

        bool has_self = false;
        int count = 0;

        // parse successors and optional label
        for (;;) {
            rd.skip_whitespace();
            if (!rd.read_uint64(&n)) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)n_vertices) {
                std::stringstream err;
                err << "invalid successor for id " << id;
                throw std::runtime_error(err.str());
            }

            if (id == n and removeBadLoops and (cur_owner != (int)(priority[id]&1))) {
                has_self = true;
            } else {
                builder.add_edge((int)id, (int)n);
                count++;
            }

            const int ch = rd.next_token();
            if (ch == EOF) throw std::runtime_error("missing ; to end line");
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of line
            if (ch == '\"') {
                builder.set_label((int)id, rd.read_label());
                // now read ;
                if (rd.next_token() != ';') throw std::runtime_error("missing ; to end line");
            }
            break;
        }

        if (has_self and count == 0) {
            // we must keep it
            builder.add_edge((int)id, (int)id);
            count++;
        }
    }

    if (!seen.all()) {
        std::stringstream err;
        err << "missing nodes, seen only " << seen.count();
        throw std::runtime_error(err.str());
    }

    // we now need to fix the priorities first...
    // Compress priorities: assign each distinct priority a dense value that
    // preserves order and parity. This needs the distinct priorities in sorted
    // order, so we work on a sorted, de-duplicated vector (the parser is on the
    // hot path for huge inputs, so we avoid node-based maps).
    std::vector<uint64_t> distinct(priority);
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());

    std::vector<int> compressed(distinct.size());
    int counter = 0;
    uint64_t previous = 0;
    for (size_t i=0; i<distinct.size(); i++) {
        const uint64_t p = distinct[i];
        if (previous != p) {
            counter++;
            if ((counter&1) != (int)(p&1)) counter++;
            previous = p;
        }
        compressed[i] = counter;
    }

    for (unsigned v=0; v<node_count; v++) {
        auto it = std::lower_bound(distinct.begin(), distinct.end(), priority[v]);
        builder.set_priority(v, compressed[it - distinct.begin()]);
    }

    return builder.build();
}

}
