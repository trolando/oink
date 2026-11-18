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

/**
 * Check PGSolver syntax and the games produced at input boundaries.
 */

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "game_limits.hpp"
#include "oink/io.hpp"
#include "oink/pgparser.hpp"

using namespace pg;

static int failures = 0;
using Parser = Game (*)(std::istream&, bool);

/** A malformed input must throw std::runtime_error. */
static void
expect_reject(const std::string& name, const std::string& input,
              Parser parse = PGParser::parse_pgsolver, bool remove_loops = false,
              const std::string& message = "")
{
    std::istringstream in(input);
    try {
        parse(in, remove_loops);
        std::cerr << "FAIL [" << name << "]: expected rejection, but parsing succeeded" << std::endl;
        failures++;
    } catch (const std::runtime_error& e) {
        if (!message.empty() and e.what() != message) {
            std::cerr << "FAIL [" << name << "]: expected " << message << ", got: " << e.what() << std::endl;
            failures++;
        }
    } catch (const std::exception& e) {
        std::cerr << "FAIL [" << name << "]: expected std::runtime_error, got: " << e.what() << std::endl;
        failures++;
    }
}

/** A valid input must parse without throwing. */
static void
expect_accept(const std::string& name, const std::string& input, Parser parse = PGParser::parse_pgsolver)
{
    std::istringstream in(input);
    try {
        parse(in, false);
    } catch (const std::exception& e) {
        std::cerr << "FAIL [" << name << "]: expected acceptance, got: " << e.what() << std::endl;
        failures++;
    }
}

static void
test_repeated_edges()
{
    const std::string header = "parity 6;\r\nstart 2;\r\n\r\n";
    const std::string rest = "0 0 0 0;\n1 0 0 1;\n3 0 0 3;\n4 0 0 4;\n5 0 0 5;\n";
    const std::string long_header = "parity 12;\n";
    const std::string long_rest = "10 0 0 0,1,2,3,4,5,6,7,8;\n"
        "9 0 0 9;\n8 0 0 8;\n7 0 0 7;\n6 0 0 6;\n5 0 0 5;\n4 0 0 4;\n"
        "3 0 0 3;\n2 0 0 2;\n1 0 0 1;\n0 0 0 0;\n";
    for (Parser parse : {PGParser::parse_pgsolver, PGParser::parse_pgsolver_renumber}) {
        for (bool remove_loops : {false, true}) {
            expect_reject("adjacent repeated successor", header + "2 3 1 5,0,0;\n" + rest,
                          parse, remove_loops, "repeated edge from vertex 2 to successor 0");
            expect_reject("nonadjacent repeated successor", header + "2 3 1 5,0,3,0;\n" + rest,
                          parse, remove_loops, "repeated edge from vertex 2 to successor 0");
            expect_reject("repeated self-loop loss", header + "2 3 0 2,1,2;\n" + rest,
                          parse, remove_loops, "repeated edge from vertex 2 to successor 2");
            expect_reject("repeated self-loop win", header + "2 3 1 2,1,2;\n" + rest,
                          parse, remove_loops, "repeated edge from vertex 2 to successor 2");
        }
        expect_reject("CR-only line endings", "parity 1;\rstart 0;\r0 0 0 0,0;",
                      parse, false, "repeated edge from vertex 0 to successor 0");
        expect_reject("repeated successor in a long list", long_header + "11 0 0 0,1,2,3,4,5,6,7,8,0;\n" + long_rest,
                      parse, false, "repeated edge from vertex 11 to successor 0");
        expect_accept("different vertices can share successors", long_header + "11 0 0 0,1,2,3,4,5,6,7,8;\n" + long_rest, parse);
        std::string bucket_rest;
        for (int v = 1; v <= 33; v++) bucket_rest += std::to_string(v) + " 0 0 0;\n";
        expect_accept("wide list with distinct targets", "parity 33;\n0 0 0 0,32,1,33,2,3,4,5,6;\n" + bucket_rest, parse);
        expect_reject("wide list with repeated target 32", "parity 33;\n0 0 0 0,32,1,33,2,3,4,5,6,32;\n" + bucket_rest,
                      parse, false, "repeated edge from vertex 0 to successor 32");
        // Split the second occurrence of 10 between the reader's blocks.
        const std::string prefix = long_header + "11 0 0 10, 0,";
        const std::string split = prefix + std::string((1 << 20) - prefix.size() - 1, ' ');
        expect_reject("duplicate split across input blocks", split + "10;\n" + long_rest,
                      parse, false, "repeated edge from vertex 11 to successor 10");

        std::istringstream in(header + "2 3 0 5,2,0;\n" + rest);
        Game g = parse(in, true);
        if (g.outcount(2) != 2 or g.outs(2)[0] != 5 or g.outs(2)[1] != 0) {
            std::cerr << "FAIL [losing self-loop]: loop removal changed the other successors" << std::endl;
            failures++;
        }
    }
    expect_reject("renumbering preserves input vertex IDs",
                  "parity 2;\n2 4294967297 0 1,1;\n0 0 0 0;\n1 0 0 1;\n",
                  PGParser::parse_pgsolver_renumber, false,
                  "repeated edge from vertex 2 to successor 1");
}

static void
test_bundled_games()
{
    namespace fs = std::filesystem;
    const fs::path tests = OINK_TEST_DATA_DIR;
    std::size_t parsed = 0;
    for (const auto& root : {tests, tests.parent_path() / "examples"}) {
        for (const auto& entry : fs::recursive_directory_iterator(root)) {
            if (!entry.is_regular_file()) continue;
            const auto& path = entry.path();
            const std::string name = path.filename().string();
            const bool bzip2 = path.extension() == ".bz2" and path.stem().extension() == ".pg";
            const bool fixture = root == tests and path.extension().empty() and name.compare(0, 2, "vb") == 0;
#ifndef OINK_HAVE_BZIP2
            if (bzip2) continue;
#endif
            if (!bzip2 and !fixture and path.extension() != ".pg") continue;
            try {
                for (Parser parse : {PGParser::parse_pgsolver, PGParser::parse_pgsolver_renumber}) {
                    parse(*open_input(path.string()), false);
                }
                parsed++;
            } catch (const std::exception& e) {
                std::cerr << "FAIL [" << path << "]: " << e.what() << std::endl;
                failures++;
            }
        }
    }
    std::cout << "Loaded " << parsed << " bundled games with both parsers" << std::endl;
}

struct ExpectedVertex {
    int priority;
    int owner;
    const char* label;
    std::vector<int> successors;
};

static void
expect_game(const char* name, std::istream& in, Parser parse,
            const std::vector<ExpectedVertex>& expected, bool remove_loops = false)
{
    try {
        const Game game = parse(in, remove_loops);
        if (game.nodecount() != static_cast<long>(expected.size())) {
            throw std::runtime_error("wrong vertex count");
        }
        long edges = 0;
        for (int v = 0; v < game.nodecount(); v++) {
            const auto& vertex = expected[v];
            const std::string* label = game.rawlabel(v);
            const auto successors = game.out_edges(v);
            if (game.priority(v) != vertex.priority or game.owner(v) != vertex.owner
                or (vertex.label ? !label or *label != vertex.label : label != nullptr)
                or !std::equal(successors.begin(), successors.end(),
                               vertex.successors.begin(), vertex.successors.end())) {
                throw std::runtime_error("wrong data at vertex " + std::to_string(v));
            }
            edges += vertex.successors.size();
        }
        if (game.edgecount() != edges) throw std::runtime_error("wrong edge count");
    } catch (const std::exception& e) {
        std::cerr << "FAIL [" << name << "]: " << e.what() << std::endl;
        failures++;
    }
}

static std::vector<ExpectedVertex>
expected_records(Parser parse)
{
    const bool renumber = parse == PGParser::parse_pgsolver_renumber;
    return {{0, 0, nullptr, {2, 1}},
            {renumber ? 1 : 3, 0, "node1", {0}},
            {renumber ? 3 : 5, 1, "2,;name", {0, 2}}};
}

int
test_parser_input(bool renumber, const std::string& path)
{
    const Parser parse = renumber ? PGParser::parse_pgsolver_renumber : PGParser::parse_pgsolver;
    if (path == "-") {
        expect_game("standard input", std::cin, parse, expected_records(parse));
    } else {
        auto in = open_input(path);
        expect_game(path.c_str(), *in, parse, expected_records(parse));
    }
    return failures ? 1 : 0;
}

static void
test_reader_boundaries()
{
    constexpr size_t block_size = 1 << 20; // Reader reads 1 MiB at a time.
    std::string tail = ";\n";
    std::vector<ExpectedVertex> boundary = {{0, 0, nullptr, {1, 10}}};
    for (int v = 1; v <= 10; v++) {
        tail += std::to_string(v) + " 0 0 0;\n";
        boundary.push_back({0, 0, nullptr, {0}});
    }

    for (Parser parse : {PGParser::parse_pgsolver, PGParser::parse_pgsolver_renumber}) {
        for (const char* header : {"parity 2;\r\n", "parity 3;\r\n"}) {
            std::istringstream in(std::string(header) + "\r\nstart 2;\r\n\r\n"
                "2 5 1 0,2 \"2, ; name\";\r\n0 0 0 2,1;\r\n1 3 0 0 \"node 1\";");
            expect_game("header, CRLF, start, labels and EOF after semicolon",
                        in, parse, expected_records(parse));
        }

        const std::string prefix = "parity 10;\n0 0 0 1,";
        // These positions put the block boundary after, inside, and before 10.
        for (size_t digits_before_boundary : {2, 1, 0}) {
            const std::string padding(block_size - prefix.size() - digits_before_boundary, ' ');
            std::istringstream in(prefix + padding + "10" + tail);
            expect_game("successor at a block boundary", in, parse, boundary);
        }

        expect_reject("missing final semicolon", "parity 0;\n0 0 0 0",
                      parse, false, "missing ; to end line");
        const std::string record = "0 0 0 0";
        const std::string header = "parity 0;\n";
        const std::string padding(block_size - header.size() - record.size(), ' ');
        expect_reject("missing semicolon at block boundary", header + padding + record,
                      parse, false, "missing ; to end line");

        std::istringstream empty("parity 0;");
        expect_game("zero vertices", empty, parse, {});

        std::istringstream loops("parity 1;\n0 0 1 0,1;\n1 0 1 1;");
        expect_game("losing loops with and without an alternative", loops, parse,
                    {{0, 1, nullptr, {1}}, {0, 1, nullptr, {1}}}, true);
    }

    std::istringstream large_priority("parity 1;\n1 4294967297 1 0,1;\n0 0 0 1;");
    expect_game("renumber a large priority without changing vertex IDs", large_priority,
                PGParser::parse_pgsolver_renumber, {{0, 0, nullptr, {1}}, {1, 1, nullptr, {0, 1}}});
}

static void
test_numeric_limits()
{
    constexpr auto max_slots = detail::max_edge_slots;
    const auto max_header = max_slots / 2;
    // Test the same capacity checks used by the parsers without allocating
    // billions of vertices or edges.
    if (detail::pgsolver_vertex_capacity(max_header) != static_cast<int>(max_header) + 1
        or detail::checked_edge_capacity(max_slots - 1, max_slots + 1) != max_slots) {
        std::cerr << "FAIL [parser capacity]: largest representable sizes should fit" << std::endl;
        failures++;
    }
    try {
        detail::checked_edge_capacity(max_slots, max_slots + 1);
        std::cerr << "FAIL [edge slots]: expected rejection of one more slot" << std::endl;
        failures++;
    } catch (const std::runtime_error& e) {
        if (std::string(e.what()) != "too many edges for int edge offsets") {
            std::cerr << "FAIL [edge slots]: " << e.what() << std::endl;
            failures++;
        }
    }

    constexpr char max_u64[] = "18446744073709551615";
    constexpr char overflow[] = "18446744073709551616";
    constexpr size_t block_size = 1 << 20;
    const std::string prefix = "parity 1;\n0 ";
    const std::string tail = " 0 1;\n1 1 1 0;\n";
    const std::string int_max = std::to_string(INT_MAX);
    const std::string above_int_max = std::to_string(uint64_t(INT_MAX) + 1);
    for (Parser parse : {PGParser::parse_pgsolver, PGParser::parse_pgsolver_renumber}) {
        const bool renumber = parse == PGParser::parse_pgsolver_renumber;
        expect_reject("priority overflow changes the winner", prefix + overflow + tail,
                      parse, false, "priority exceeds UINT64_MAX");
        expect_reject("successor overflow wraps to vertex zero",
                      std::string("parity 0;\n0 0 0 ") + overflow + ";\n",
                      parse, false, "successor exceeds UINT64_MAX");
        expect_reject("vertex ID overflow wraps to vertex zero",
                      std::string("parity 0;\n") + overflow + " 0 0 0;",
                      parse, false, "vertex ID exceeds UINT64_MAX");
        expect_reject("owner overflow wraps to player zero",
                      std::string("parity 0;\n0 0 ") + overflow + " 0;",
                      parse, false, "owner exceeds UINT64_MAX");
        expect_reject("header overflow wraps to zero", std::string("parity ") + overflow + ";\n",
                      parse, false, "number of nodes exceeds UINT64_MAX");
        expect_reject("header narrows to zero", "parity 4294967296;\n0 0 0 0;\n",
                      parse, false, "number of nodes too high");
        expect_reject("smallest unrepresentable header",
                      "parity " + std::to_string(max_header + 1) + ";",
                      parse, false, "number of nodes too high");
        expect_reject("maximum UINT64 header", std::string("parity ") + max_u64 + ";",
                      parse, false, "number of nodes too high");
        expect_reject("smallest invalid ID with parity 1", "parity 1;\n2 0 0 0;",
                      parse, false, renumber ? "invalid id (too high)" : "invalid id");
        expect_reject("smallest invalid successor with parity 1", "parity 1;\n0 0 0 2;",
                      parse, false, "invalid successor 2 for vertex 0");
        expect_reject("smallest invalid owner", "parity 0;\n0 0 2 0;",
                      parse, false, renumber ? "invalid owner (must be 0 or 1)" : "invalid owner");

        std::istringstream game("parity 1;\n0 0 0 1;\n1 " + int_max + " 1 1;");
        expect_game("maximum ID and successor for this header, maximum owner and int priority",
                    game, parse, {{0, 0, nullptr, {1}}, {renumber ? 1 : INT_MAX, 1, nullptr, {1}}});
        if (!renumber) {
            expect_reject("smallest priority above INT_MAX", "parity 0;\n0 " + above_int_max + " 0 0;",
                          parse, false, "priority too high");
            expect_reject("maximum UINT64 priority in the fast parser", prefix + max_u64 + tail,
                          parse, false, "priority too high");
        }

        for (size_t digits_before_boundary : {0, 1, 19, 20}) {
            const auto padded = prefix + std::string(block_size - prefix.size() - digits_before_boundary, ' ');
            expect_reject("overflow at a reader block boundary", padded + overflow + tail,
                          parse, false, "priority exceeds UINT64_MAX");
            if (renumber) {
                std::istringstream maximum(padded + max_u64 + tail);
                expect_game("maximum UINT64 priority at a reader block boundary", maximum, parse,
                            {{3, 0, nullptr, {1}}, {1, 1, nullptr, {0}}});
            }
        }

        const auto zeros = std::string(block_size + 20, '0');
        std::istringstream leading_zeros("parity 0;\n0 " + zeros + "1 0 0;");
        expect_game("arbitrary leading zeros across reader blocks", leading_zeros, parse,
                    {{1, 0, nullptr, {0}}});
        expect_reject("overflow after arbitrary leading zeros", prefix + zeros + overflow + tail,
                      parse, false, "priority exceeds UINT64_MAX");
    }
}

static void
test_record_endings()
{
    for (Parser parse : {PGParser::parse_pgsolver, PGParser::parse_pgsolver_renumber}) {
        expect_reject("wrong vertex terminator", "parity 0;\n0 0 0 0!\n",
                      parse, false, "missing ; to end line");
        expect_reject("wrong labeled vertex terminator", "parity 0;\n0 0 0 0 \"node\"!\n",
                      parse, false, "missing ; to end line");
        expect_reject("negative optional vertex ID", "parity 1;\n0 0 0 0;\n-1 1 1 0;\n",
                      parse, false, "unable to read id");
        expect_reject("misspelled start", "parity 0;\nsomething 0;\n0 1 0 0;",
                      parse, false, "expecting start statement");
        expect_reject("missing start terminator", "parity 0;\nstart 0 0 1 0 0;",
                      parse, false, "missing ';' after start statement");
        expect_reject("out-of-range start", "parity 0;\nstart 1; 0 1 0 0;",
                      parse, false, "invalid start vertex ID");
        expect_reject("start at omitted vertex", "parity 1;\nstart 1; 0 1 0 0;",
                      parse, false, "invalid start vertex ID");
        expect_reject("overflowing start ID", "parity 0;\nstart 18446744073709551616; 0 1 0 0;",
                      parse, false, "start vertex ID exceeds UINT64_MAX");
        expect_reject("successor of missing last vertex", "parity 1;\n0 0 0 1;\n",
                      parse, false, "invalid successor 1 for vertex 0");
        expect_reject("missing last vertex with records out of order",
                      "parity 2;\n1 5 0 0,2;\n0 0 0 0;\n",
                      parse, false, "invalid successor 2 for vertex 1");

        std::istringstream same_line("parity 0;\nstart 0; 0 1 0 0;\n");
        expect_game("start followed by a vertex on the same line", same_line, parse,
                    {{1, 0, nullptr, {0}}});
    }
}

int
test_parser()
{
    // Positive controls.
    expect_accept("valid game", "parity 2;\n0 0 0 0,1;\n1 1 1 1;\n");
    expect_accept("valid with label", "parity 2;\n0 0 0 0,1 \"a\";\n1 1 1 1;\n");

    // Labels round-trip through the parser into Game storage.
    {
        std::istringstream in("parity 2;\n0 0 0 0,1 \"hello\";\n1 1 1 1;\n");
        Game g = PGParser::parse_pgsolver(in, false);
        const std::string* l0 = g.rawlabel(0);
        const std::string* l1 = g.rawlabel(1);
        if (l0 == nullptr or *l0 != "hello") {
            std::cerr << "FAIL [label storage]: vertex 0 label not stored correctly" << std::endl;
            failures++;
        }
        if (l1 != nullptr) {
            std::cerr << "FAIL [label storage]: vertex 1 should have no label" << std::endl;
            failures++;
        }
    }

    // Malformed input.
    expect_reject("empty file", "");
    expect_reject("missing header", "notparity 2;\n0 0 0 1;\n");
    expect_reject("missing semicolon after count", "parity 2\n0 0 0 1;\n");
    expect_reject("invalid vertex id (too high)", "parity 1;\n5 0 0 0;\n");
    expect_reject("invalid owner", "parity 1;\n0 0 2 0;\n");
    expect_reject("invalid priority (above INT_MAX)", "parity 1;\n0 9999999999 0 0;\n");
    expect_reject("missing successor", "parity 1;\n0 0 0;\n");
    expect_reject("malformed edge list (no terminator)", "parity 1;\n0 0 0 0");
    expect_reject("out-of-range successor", "parity 1;\n0 0 0 5;\n");
    expect_reject("duplicate vertex id", "parity 1;\n0 0 0 0;\n0 0 0 0;\n");
    expect_reject("missing nodes", "parity 3;\n0 0 0 0;\n");

    test_reader_boundaries();
    test_numeric_limits();
    test_record_endings();
    test_repeated_edges();
    test_bundled_games();

    if (failures) {
        std::cerr << failures << " parser test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all parser tests passed" << std::endl;
    return 0;
}
