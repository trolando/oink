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
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
expect_accept(const std::string& name, const std::string& input)
{
    std::istringstream in(input);
    try {
        PGParser::parse_pgsolver(in, false);
    } catch (const std::exception& e) {
        std::cerr << "FAIL [" << name << "]: expected acceptance, got: " << e.what() << std::endl;
        failures++;
    }
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

    if (failures) {
        std::cerr << failures << " parser test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all parser tests passed" << std::endl;
    return 0;
}
