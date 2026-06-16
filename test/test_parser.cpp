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
 * Unit tests for the pgsolver parser, focusing on malformed input.
 * Each malformed case must be rejected with a std::runtime_error; valid input
 * must parse successfully.
 */

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "oink/pgparser.hpp"

using namespace pg;

static int failures = 0;

/** A malformed input must throw std::runtime_error. */
static void
expect_reject(const std::string& name, const std::string& input)
{
    std::istringstream in(input);
    try {
        PGParser::parse_pgsolver(in, false);
        std::cerr << "FAIL [" << name << "]: expected rejection, but parsing succeeded" << std::endl;
        failures++;
    } catch (const std::runtime_error&) {
        // expected
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

int
main()
{
    // Positive controls.
    expect_accept("valid game", "parity 2;\n0 0 0 0,1;\n1 1 1 1;\n");
    expect_accept("valid with label", "parity 2;\n0 0 0 0,1 \"a\";\n1 1 1 1;\n");

    // Malformed input (see roadmap Phase 1, task 5).
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

    if (failures) {
        std::cerr << failures << " parser test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all parser tests passed" << std::endl;
    return 0;
}
