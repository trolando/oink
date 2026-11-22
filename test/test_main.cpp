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
 * Runs the unit tests: all of them, or the one given on the command line.
 */

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int test_measure_ops();
int test_measure_codes();
int test_parser();
int test_parser_input(bool renumber, const std::string& path);
int test_io();
int test_solution();
int test_pma(const char* fixture);
int test_ppj();
int test_zlkj();
int test_region_renumber();
int test_region_renumber_files(const std::vector<std::string>& paths);
int test_psi();
int test_lace_restart();
int test_fpj_batch();
int test_fpi_frontier();
int test_fpi_parallel();
int test_tlq_cache();
int test_tlj_cache();
int test_gpm_cache();
int test_pmtl_cache();
int test_in_edge_map();
int test_builder();
int test_game();
int test_verifier();
int test_terminal();
int test_chooser();
int test_bench();
int test_table();

namespace {

struct Test {
    const char* name;
    int (*run)();
};

const Test tests[] = {
    {"measure-ops", test_measure_ops},
    {"measure-codes", test_measure_codes},
    {"parser", test_parser},
    {"io", test_io},
    {"solution", test_solution},
    {"pma", [] { return test_pma(OINK_TEST_DATA_DIR "/vb232"); }},
    {"ppj", test_ppj},
    {"zlkj", test_zlkj},
    {"region-renumber", test_region_renumber},
    {"psi", test_psi},
    {"lace-restart", test_lace_restart},
    {"fpj-batch", test_fpj_batch},
    {"fpi-frontier", test_fpi_frontier},
    {"fpi-parallel", test_fpi_parallel},
    {"tlq-cache", test_tlq_cache},
    {"tlj-cache", test_tlj_cache},
    {"gpm-cache", test_gpm_cache},
    {"pmtl-cache", test_pmtl_cache},
    {"in-edge-map", test_in_edge_map},
    {"builder", test_builder},
    {"game", test_game},
    {"verifier", test_verifier},
    {"terminal", test_terminal},
    {"chooser", test_chooser},
    {"bench", test_bench},
    {"table", test_table},
};

int
run(const Test& test)
{
    std::cout << "[" << test.name << "]" << std::endl;
    try {
        return test.run();
    } catch (const std::exception& error) {
        std::cerr << test.name << ": " << error.what() << std::endl;
        return 1;
    }
}

}

int
main(int argc, char** argv)
{
    if (argc == 1) {
        int failures = 0;
        for (const auto& test : tests) failures += run(test) != 0;
        return failures ? 1 : 0;
    }
    if (argc == 2) {
        const std::string name = argv[1];
        for (const auto& test : tests) {
            if (name == test.name) return run(test);
        }
        if (name == "--list") {
            for (const auto& test : tests) std::cout << test.name << std::endl;
            return 0;
        }
    }
    if (argc == 3) {
        const std::string mode = argv[1];
        if (mode == "parser-input-fast" or mode == "parser-input-renumber") {
            return test_parser_input(mode == "parser-input-renumber", argv[2]);
        }
    }
    if (argc >= 3 and std::string(argv[1]) == "region-renumber-files") {
        // the region-renumber test on the given games (see test/test_region_renumber.cmake)
        std::cout << "[region-renumber-files]" << std::endl;
        try {
            return test_region_renumber_files(std::vector<std::string>(argv+2, argv+argc));
        } catch (const std::exception& error) {
            std::cerr << "region-renumber-files: " << error.what() << std::endl;
            return 1;
        }
    }
    std::cerr << "Usage: " << argv[0] << " [--list|TEST]" << std::endl;
    return 2;
}
