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
 * The command line of bench_solvers.
 */

#ifndef OINK_BENCH_OPTIONS_HPP
#define OINK_BENCH_OPTIONS_HPP

#include <ostream>
#include <string>
#include <vector>

#include "bench.hpp"

namespace bench {

/**
 * Parse the command line into <options> and the roster <solvers>; external
 * solvers (-e) are registered in the solver registry. Returns false after
 * printing the help to <out>. Reports implicit solver exclusions to <notices>
 * if supplied. Throws std::runtime_error for a bad command line.
 */
bool parse_command_line(int argc, const char* const* argv, Options& options, std::vector<SolverRun>& solvers,
                        std::ostream& out, std::ostream* notices = nullptr);

}

#endif
