/*
 * Copyright 2020-2026 Tom van Dijk
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

#include "external_solver.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "oink/io.hpp"
#include "oink/solvers.hpp"

namespace fs = std::filesystem;

// Generate a unique, non-existing path in the system temp directory.
// Replaces boost::filesystem::unique_path() (no std::filesystem equivalent).
static fs::path
unique_temp_path()
{
    static std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist;
    const auto dir = fs::temp_directory_path();
    for (;;) {
        std::ostringstream name;
        name << "oink-" << std::hex << dist(gen);
        auto p = dir / name.str();
        if (!fs::exists(p)) return p;
    }
}

// Run a shell command and stream its combined stdout+stderr to `logger`.
// Portable replacement for boost::process (popen / _popen on Windows).
// Returns the exit status of the command, or -1 if it could not be run.
static int
run_command_logged(const std::string& cmd, std::ostream& logger)
{
    const std::string full = cmd + " 2>&1";
#ifdef _WIN32
    FILE* pipe = _popen(full.c_str(), "r");
#else
    FILE* pipe = popen(full.c_str(), "r");
#endif
    if (pipe == nullptr) {
        logger << "failed to run command: " << cmd << std::endl;
        return -1;
    }
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), pipe) != nullptr) logger << buf;
#ifdef _WIN32
    return _pclose(pipe);
#else
    const int status = pclose(pipe);
    if (status < 0) return -1;
    if (WIFEXITED(status)) return WEXITSTATUS(status);
    return WIFSIGNALED(status) ? 128+WTERMSIG(status) : status;
#endif
}

static void
replace_all(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Move past the replaced substring
    }
}

void
ExternalSolver::run()
{
    // TODO: handle partial files (where some nodes are disabled / only solving a subgame)

    // Obtain temporary filenames
    auto path1 = unique_temp_path();
    auto path2 = unique_temp_path();

    // Write parity game to temporary file
    std::ofstream file(path1);
    game.write_pgsolver(file);
    file.close();

    // Replace %I with input filename and %O with output filename
    std::string command = cmd;
    replace_all(command, "%I", path1.string());
    replace_all(command, "%O", path2.string());

    // Run the external program (no timeout), logging its output
    const int status = run_command_logged(command, logger);
    if (status != 0) {
        fs::remove(path1);
        fs::remove(path2);
        throw std::runtime_error("external solver failed with exit status " + std::to_string(status));
    }

    // Read solution from temporary file into the (Oink-owned) solution
    std::ifstream sol(path2);
    pg::parse_solution(game, solution(), sol);
    sol.close();

    // Delete temporary files
    fs::remove(path1);
    fs::remove(path2);
}

std::string
register_external_solver(const std::string& spec)
{
    // Parse "name:command"; fall back to using the whole string as both if no colon
    const size_t colon = spec.find(':');
    const std::string name = (colon != std::string::npos) ? spec.substr(0, colon) : spec;
    const std::string cmd = (colon != std::string::npos) ? spec.substr(colon + 1) : spec;
    pg::Solvers::add(name, "external tool", 0, [cmd](pg::Oink& oink, pg::Game& game) {
        return std::make_unique<ExternalSolver>(oink, game, cmd);
    });
    return name;
}
