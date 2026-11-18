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
#include "game_limits.hpp"
#include "printf.hpp"
#include "successor_check.hpp"

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
    explicit Reader(std::istream &in) : rd(in.rdbuf()), buf((1 << 20) + 1) { }

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

    inline bool read_uint64(uint64_t *res, const char *field)
    {
        if (p == end and !fill()) return false;
        const char *start = p;
        unsigned digit = *p - '0';
        if (digit > 9) return false;
        uint64_t r = 0;
        // A non-digit at end lets the digit loop scan each block without
        // checking the buffer boundary after every digit.
        while (digit <= 9) {
            r = (10*r) + digit;
            digit = *++p - '0';
        }
        // Up to 19 decimal digits always fit. Check longer tokens and
        // tokens that might continue in the next block out of line.
        if (__builtin_expect(p == end or p - start >= 20, false)) r = read_long_uint64(start, r, field);
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
    [[gnu::noinline]] uint64_t read_long_uint64(const char *start, uint64_t r, const char *field)
    {
        if (p - start >= 20) {
            // Re-read a possibly wrapped prefix. Allow arbitrary leading zeros.
            p = start;
            r = 0;
        }
        for (;;) {
            const unsigned digit = peek() - '0';
            if (digit > 9) return r;
            if (r > UINT64_MAX / 10 or (r == UINT64_MAX / 10 and digit > UINT64_MAX % 10)) {
                throw std::runtime_error(std::string(field) + " exceeds UINT64_MAX");
            }
            r = 10*r + digit;
            ++p;
        }
    }

    bool fill()
    {
        const std::streamsize k = rd->sgetn(buf.data(), (std::streamsize)buf.size() - 1);
        p = buf.data();
        end = p + (k > 0 ? k : 0);
        buf[end - p] = '\0';
        return k > 0;
    }

    std::streambuf *rd;
    std::vector<char> buf;
    const char *p = nullptr;
    const char *end = nullptr;
};

int
read_header(Reader &rd)
{
    for (const char *kw = "parity"; *kw; kw++) {
        if (rd.get() != *kw) throw std::runtime_error("expecting parity game specification");
    }
    uint64_t n;
    rd.skip_whitespace();
    if (!rd.read_uint64(&n, "number of nodes")) throw std::runtime_error("missing number of nodes");
    if (rd.next_token() != ';') throw std::runtime_error("missing ';' after parity header");
    return detail::pgsolver_vertex_capacity(n);
}

int
read_start(Reader &rd, int vertices)
{
    rd.skip_whitespace();
    if (rd.peek() != 's') return -1;
    for (const char *kw = "start"; *kw; kw++) {
        if (rd.get() != *kw) throw std::runtime_error("expecting start statement");
    }
    uint64_t start;
    rd.skip_whitespace();
    if (!rd.read_uint64(&start, "start vertex ID")) throw std::runtime_error("missing start vertex ID");
    if (start >= static_cast<uint64_t>(vertices)) throw std::runtime_error("invalid start vertex ID");
    if (rd.next_token() != ';') throw std::runtime_error("missing ';' after start statement");
    return static_cast<int>(start);
}

[[noreturn]] void
invalid_successor(uint64_t source, uint64_t target)
{
    throw std::runtime_error("invalid successor " + std::to_string(target) +
                             " for vertex " + std::to_string(source));
}

void
check_omitted_vertex(const Game& game, int omitted, int start)
{
    // Scan the packed array first. Look up the source only when reporting an error.
    const auto* edges = game.outedges();
    const auto* end = edges + game.edgeArraySize();
    if (std::find(edges, end, omitted) != end) {
        for (int v = 0; v < omitted; v++) {
            if (game.has_edge(v, omitted)) invalid_successor(v, omitted);
        }
    }
    if (start == omitted) throw std::runtime_error("invalid start vertex ID");
}

// Keep validation out of line so the renumbering parser's hot path stays small.
[[gnu::noinline]] void
check_successors(const int* successors, int count, int source, int vertices,
                 std::vector<int>& stamps)
{
    const int repeated = detail::repeated_successor(successors, count, source, vertices, stamps);
    if (repeated != -1) throw std::runtime_error(detail::repeated_edge_message(source, repeated));
}

void
check_successors(const int* successors, int count, int source, int vertices,
                 std::vector<uint64_t>& stamps)
{
    const int repeated = detail::repeated_successor(successors, count, source, vertices, stamps);
    if (repeated != -1) throw std::runtime_error(detail::repeated_edge_message(source, repeated));
}

int
initial_edge_count(uint64_t vertices)
{
    // Use Game's usual four edges per vertex unless the initial array would
    // overflow its int size calculation. Large games grow the array as needed.
    return vertices <= (INT_MAX - 1) / 5 ? -1 : 1;
}

}

Game
PGParser::parse_pgsolver(std::istream &inp, bool removeBadLoops)
{
    Reader rd(inp);
    const int vertices = read_header(rd);
    const int start = read_start(rd, vertices);
    uint64_t n;

    // Skip the first few allocations while keeping the same growth sequence
    // for larger edge arrays. The initial array includes the terminators.
    size_t capacity = static_cast<size_t>(vertices) + 2;
    for (int growth = 0; growth < 3; growth++) capacity = Game::grown_edge_capacity(capacity);
    Game res(vertices, capacity <= detail::max_edge_slots ? capacity - vertices - 1 : 1);

    /**
     * Parse the nodes...
     */

    int node_count = 0; // number of read nodes
    bitset seen(res.n_vertices); // tracks which nodes have been read
    std::vector<uint64_t> stamps;

    while (node_count < res.n_vertices) {
        uint64_t id;
        rd.skip_whitespace();
        if (!rd.read_uint64(&id, "vertex ID")) {
            // we expect maybe one more node...
            if (rd.peek() == EOF and node_count == vertices - 1 and !seen[vertices - 1]) {
                check_omitted_vertex(res, node_count, start);
                // The omitted vertex has no row, and no edge refers to it.
                res.n_vertices = node_count;
                res._owner.resize(node_count);
                seen.resize(node_count);
                break;
            }
            throw std::runtime_error("unable to read id");
        }
        if (id >= (unsigned)res.n_vertices) throw std::runtime_error("invalid id");

        if (seen[id]) throw std::runtime_error("duplicate id");
        seen[id] = true;
        node_count++;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n, "priority")) throw std::runtime_error("missing priority");
        if (n > INT_MAX) throw std::runtime_error("priority too high"); // don't be ridiculous
        res._priority[id] = n;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n, "owner")) throw std::runtime_error("missing owner");

        if (n == 0) { /* nothing */ }
        else if (n == 1) { res._owner[id] = true; }
        else { throw std::runtime_error("invalid owner"); }

        const size_t first = res.e_size;
        res._firstouts[id] = first;

        const bool remove_self = removeBadLoops and res._owner[id] != (res._priority[id]&1);
        bool has_self = false;

        // parse successors and optional label
        for (;;) {
            rd.skip_whitespace();
            if (!rd.read_uint64(&n, "successor")) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)res.n_vertices) invalid_successor(id, n);

            if (remove_self and id == n) {
                if (has_self) throw std::runtime_error(detail::repeated_edge_message(id, n));
                has_self = true;
            } else {
                res.append_out(n);
            }

            const int ch = rd.next_token();
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of record
            if (ch == '\"') {
                res._label[id] += rd.read_label();
                if (rd.next_token() == ';') break;
            }
            throw std::runtime_error("missing ; to end line");
        }

        if (has_self and res.e_size == first) res.append_out(id);
        const int count = res.e_size - first;
        if (count > 1) check_successors(res.outs(id), count, id, res.n_vertices, stamps);
        res._outcount[id] = count;
        res.n_edges += count;
        res.append_out(-1);
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
    const int n_vertices = read_header(rd);
    const int start = read_start(rd, n_vertices);
    uint64_t n;

    // Store each outgoing list directly. The priorities are kept separately
    // until they have been compressed to fit in a Game.
    Game res(n_vertices, initial_edge_count(n_vertices));
    res.e_size = 0;
    bitset seen(n_vertices);
    std::vector<uint64_t> priority(n_vertices);
    std::vector<int> stamps;

    /**
     * Parse the nodes...
     */

    int node_count = 0; // number of read nodes
    bool in_vertex_order = true;

    while (node_count < n_vertices) {
        uint64_t id;
        rd.skip_whitespace();
        if (!rd.read_uint64(&id, "vertex ID")) {
            if (rd.peek() == EOF and node_count == n_vertices - 1 and !seen[n_vertices - 1]) {
                check_omitted_vertex(res, node_count, start);
                // The omitted vertex has no row, and no edge refers to it.
                res.n_vertices = node_count;
                res._owner.resize(node_count);
                seen.resize(node_count);
                break;
            }
            throw std::runtime_error("unable to read id");
        }
        if (id >= (uint64_t)n_vertices) throw std::runtime_error("invalid id (too high)");
        if (seen[id]) throw std::runtime_error("duplicate id");
        seen[id] = true;
        in_vertex_order &= id == node_count;
        node_count++;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n, "priority")) throw std::runtime_error("missing priority");
        priority[id] = n;

        rd.skip_whitespace();
        if (!rd.read_uint64(&n, "owner")) throw std::runtime_error("missing owner");
        if (n != 0 and n != 1) throw std::runtime_error("invalid owner (must be 0 or 1)");
        res._owner[id] = n != 0;

        const size_t first = res.e_size;
        res._firstouts[id] = first;
        const bool remove_self = removeBadLoops and res._owner[id] != (priority[id]&1);
        bool has_self = false;

        // parse successors and optional label
        for (;;) {
            rd.skip_whitespace();
            if (!rd.read_uint64(&n, "successor")) throw std::runtime_error("missing successor");
            if (n >= (uint64_t)n_vertices) invalid_successor(id, n);

            if (remove_self and id == n) {
                if (has_self) throw std::runtime_error(detail::repeated_edge_message(id, n));
                has_self = true;
            } else {
                res.append_out(n);
            }

            const int ch = rd.next_token();
            if (ch == ',') continue; // next successor
            if (ch == ';') break; // end of record
            if (ch == '\"') {
                res._label[id] = rd.read_label();
                if (rd.next_token() == ';') break;
            }
            throw std::runtime_error("missing ; to end line");
        }

        if (has_self and res.e_size == first) res.append_out(id);
        const int count = res.e_size - first;
        if (count > 1) check_successors(res.outs(id), count, id, res.n_vertices, stamps);
        res._outcount[id] = count;
        res.n_edges += count;
        res.append_out(-1);
    }

    if (!seen.all()) {
        std::stringstream err;
        err << "missing nodes, seen only " << seen.count();
        throw std::runtime_error(err.str());
    }

    if (node_count == 0) return Game();

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
        res._priority[v] = compressed[it - distinct.begin()];
    }

    // Preserve the renumbering parser's layout, with outgoing lists in vertex
    // order. The common case already has that layout.
    if (!in_vertex_order) {
        std::vector<int> edges(res.e_size);
        size_t pos = 0;
        for (int v = 0; v < res.n_vertices; v++) {
            const int count = res._outcount[v] + 1;
            std::copy_n(res.outs(v), count, edges.data() + pos);
            res._firstouts[v] = pos;
            pos += count;
        }
        res._outedges.swap(edges);
        res.e_allocated = res.e_size;
    }

    res.is_ordered = std::is_sorted(res._priority.begin(), res._priority.begin() + node_count);
    return res;
}

}
