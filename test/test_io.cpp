/*
 * Copyright 2026 Tom van Dijk
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

/**
 * Tests for open_input: decompressing .bz2 input (in parallel for large files, which have
 * many blocks), with several streams, trailing data, and corrupt or truncated files. And
 * for writing and parsing solutions, with single and multi-strategies.
 */

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include "oink/game.hpp"
#include "oink/game_builder.hpp"
#include "oink/io.hpp"
#include "oink/pgparser.hpp"
#include "oink/solution.hpp"
#include "verifier.hpp"

#ifdef OINK_HAVE_BZIP2
#include <bzlib.h>
#include <vector>
#include "../src/bzip2_parallel.hpp" // internal, to check that the parallel path is taken
#endif

using namespace pg;

static int failures = 0;

static void
check(bool ok, const std::string& name)
{
    std::cout << (ok ? "ok:   " : "FAIL: ") << name << std::endl;
    if (!ok) failures++;
}

#ifdef OINK_HAVE_BZIP2

/** Pseudo-random bytes (hardly compressible, so the compressed file has many blocks). */
static std::string
random_bytes(size_t n, uint64_t seed)
{
    std::string s(n, '\0');
    for (size_t i = 0; i < n; i++) {
        seed = seed * 6364136223846793005ULL + 1442695040888963407ULL;
        s[i] = (char)(seed >> 56);
    }
    return s;
}

static std::string
compress(const std::string& data, int level)
{
    unsigned int size = (unsigned int)(data.size() + data.size() / 100 + 600);
    std::string out(size, '\0');
    const int ret = BZ2_bzBuffToBuffCompress(&out[0], &size, const_cast<char*>(data.data()), (unsigned int)data.size(), level, 0, 0);
    if (ret != BZ_OK) throw std::runtime_error("cannot compress test data");
    out.resize(size);
    return out;
}

static std::string
write_temp(const std::string& name, const std::string& contents)
{
    const std::string path = std::string("oink_test_io_") + name + ".bz2";
    std::ofstream f(path, std::ios::binary);
    f << contents;
    return path;
}

static std::string
read_back(const std::string& path)
{
    auto in = open_input(path);
    return std::string((std::istreambuf_iterator<char>(*in)), std::istreambuf_iterator<char>());
}

/** Decompress with the parallel decompressor directly; false if it declines. */
static bool
parallel(const std::string& compressed, std::string& out)
{
    const std::vector<unsigned char> data(compressed.begin(), compressed.end());
    std::vector<std::string> chunks;
    if (!bzip2_decompress_parallel(data, chunks, 4)) return false;
    out.clear();
    for (const auto& c : chunks) out += c;
    return true;
}

static bool
throws(const std::string& path)
{
    try {
        read_back(path);
    } catch (std::runtime_error&) {
        return true;
    }
    return false;
}

static void
test_bzip2()
{
    // 2.5 MB at level 1 (100k blocks): about 25 blocks, over 1 MiB compressed
    const std::string a = random_bytes(2500000, 1);
    const std::string b = random_bytes(1300000, 2) + std::string(200000, 'x');
    const std::string za = compress(a, 1);
    const std::string zb = compress(b, 9);
    check(za.size() >= (1 << 20), "test data is large enough for parallel decompression");

    std::string out;
    check(parallel(za, out) and out == a, "parallel: one stream of many blocks");
    check(parallel(za + zb, out) and out == a + b, "parallel: two streams");
    check(!parallel(za + "not bzip2 data", out), "parallel: declines trailing data (sequential then)");

    std::string p = write_temp("single", za);
    check(read_back(p) == a, "one stream of many blocks");
    std::remove(p.c_str());

    p = write_temp("multi", za + zb);
    check(read_back(p) == a + b, "two concatenated streams (levels 1 and 9)");
    std::remove(p.c_str());

    p = write_temp("trailing", za + "not bzip2 data");
    check(read_back(p) == a, "data after the stream is ignored");
    std::remove(p.c_str());

    std::string corrupt = za;
    corrupt[corrupt.size() / 2] ^= 0x55;
    p = write_temp("corrupt", corrupt);
    check(throws(p), "a corrupt block is an error");
    std::remove(p.c_str());

    p = write_temp("truncated", za.substr(0, za.size() / 2));
    check(throws(p), "a truncated stream is an error");
    std::remove(p.c_str());

    const std::string small = "parity 1;\n0 0 0 1;\n1 1 1 0;\n";
    p = write_temp("small", compress(small, 9));
    check(read_back(p) == small, "a small file (decompressed sequentially)");
    Game game = PGParser::parse_pgsolver_renumber(*open_input(p), false);
    check(game.edgecount() == 2, "parse a valid bzip2 game");
    std::remove(p.c_str());

    const auto duplicate = compress("parity 1;\n0 1 0 1,0,0;\n1 0 0 0;\n", 9);
    for (const auto& contents : {duplicate, duplicate + za}) {
        // Appending the many-block stream also exercises parallel open_input.
        p = write_temp("repeated_edge", contents);
        bool rejected = false;
        try {
            PGParser::parse_pgsolver_renumber(*open_input(p), true);
        } catch (const std::runtime_error& e) {
            rejected = std::string(e.what()) == "repeated edge from vertex 0 to successor 0";
        }
        check(rejected, "reject a repeated edge in a bzip2 game");
        std::remove(p.c_str());
    }
}

#endif

/** Parse <text> as a solution of <game>; false if parsing fails. */
static bool
parse(const Game& game, const std::string& text, Solution& solution)
{
    solution = Solution(game.nodecount());
    std::istringstream in(text);
    try {
        parse_solution(game, solution, in);
    } catch (std::runtime_error&) {
        return false;
    }
    return true;
}

/** Whether the verifier accepts <solution> of the (sorted) <game>. */
static bool
verifies(const Game& game, const Solution& solution)
{
    std::ostringstream log;
    Verifier verifier(game, solution, log);
    try {
        verifier.verify(true, true, true);
    } catch (std::runtime_error&) {
        return false;
    }
    return true;
}

static void
test_solution_format()
{
    // Even wins everywhere: 0 may move to 1 or 2, which loop; 3 (Odd) must move to 0
    // (priorities in order, so the vertices keep their numbers)
    GameBuilder builder(4);
    const int prio[] = {2, 2, 2, 3};
    for (int v = 0; v < 4; v++) builder.set_priority(v, prio[v]);
    builder.set_owner(3, Player::Odd);
    builder.add_edge(0, 1);
    builder.add_edge(0, 2);
    builder.add_edge(1, 1);
    builder.add_edge(2, 2);
    builder.add_edge(3, 0);
    Game game = builder.build();
    game.ensure_sorted(); // (already sorted)

    // a multi-strategy is written with all moves
    Solution multi(4);
    multi.init_multi(game.edgeArraySize());
    multi.solve(0, 0, -1);
    multi.add_edge_index(game.firstout(0));
    multi.add_edge_index(game.firstout(0) + 1);
    multi.solve(1, 0, 1); // (single and multi mix)
    multi.solve(2, 0, -1);
    multi.add_edge_index(game.firstout(2));
    multi.solve(3, 0, -1);
    std::ostringstream out;
    write_solution(game, multi, out);
    check(out.str() == "paritysol 4;\n0 0 1,2;\n1 0 1;\n2 0 2;\n3 0;\n", "write a multi-strategy");

    // and read back
    Solution read;
    check(parse(game, out.str(), read), "parse a multi-strategy");
    check(read.has_multi() and read.strategy(0) == -1 and has_strategy_edge_to(game, read, 0, 1) and has_strategy_edge_to(game, read, 0, 2),
          "both moves of the multi-strategy");
    check(read.strategy(1) == 1 and read.strategy(2) == 2 and read.is_solved(3) and read.winner(3) == 0, "the single moves");
    check(verifies(game, read), "the parsed multi-strategy verifies");
    std::ostringstream again;
    write_solution(game, read, again);
    check(again.str() == out.str(), "write it again the same");

    // single moves only: no multi-strategy
    check(parse(game, "paritysol 4;\n0 0 2;\n1 0 1;\n2 0 2;\n3 0;\n", read) and !read.has_multi() and read.strategy(0) == 2, "single moves");
    check(parse(game, "0 0 2, 1 ;\n", read) and has_strategy_edge_to(game, read, 0, 1) and has_strategy_edge_to(game, read, 0, 2), "spaces around the comma");
    check(!parse(game, "0 0 1,3;\n", read), "a move that is not an edge");
    check(!parse(game, "0 0 1,;\n", read), "a comma without a move");
    check(!parse(game, "0 0;\n", read), "no move for the owner");
    check(parse(game, "0 0 1,1,2;\n", read) and has_strategy_edge_to(game, read, 0, 1) and has_strategy_edge_to(game, read, 0, 2), "a duplicate move counts once");
    std::ostringstream dup;
    write_solution(game, read, dup);
    check(dup.str() == "paritysol 1;\n0 0 1,2;\n", "and is written once");
    check(parse(game, "0 0 1\n1 0 1 ;\n3 0\r\n", read) and read.strategy(0) == 1 and read.is_solved(3), "the semicolon is optional");

    // malformed lines are errors, not cut short
    check(!parse(game, "0 0 1,2x,0;\n", read), "a move with trailing letters");
    check(!parse(game, "0 0 1,2 0;\n", read), "a move without a comma");
    check(!parse(game, "0 0 1;junk\n", read), "text after the semicolon");
    check(!parse(game, "0x 0 1;\n", read), "a vertex with trailing letters");
    check(!parse(game, "x 0;\n", read), "a vertex that is not a number");
    check(!parse(game, "3 0 0;\n", read), "a move of the losing player");
    check(!parse(game, "3 0 0,1;\n", read), "moves of the losing player");
}

int
test_io()
{
#ifdef OINK_HAVE_BZIP2
    test_bzip2();
#else
    std::cout << "bzip2 support not built, nothing to test" << std::endl;
#endif
    test_solution_format();
    if (failures) {
        std::cout << failures << " test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all tests passed" << std::endl;
    return 0;
}
