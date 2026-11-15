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

/**
 * External solvers for test_solvers and bench_solvers: a command that reads a
 * parity game from the file %I and writes its solution to the file %O.
 */

#ifndef OINK_TEST_EXTERNAL_SOLVER_HPP
#define OINK_TEST_EXTERNAL_SOLVER_HPP

#include <string>

#include "oink/solver.hpp"

class ExternalSolver : public pg::Solver
{
public:
    ExternalSolver(pg::Oink& oink, pg::Game& game, std::string cmd) : Solver(oink, game), cmd(cmd) { }

    void run() override;

private:
    std::string cmd;
};

/**
 * Register the external solver given as "name:command" (or a command alone,
 * which is then also its name) in the solver registry, and return its name.
 */
std::string register_external_solver(const std::string& spec);

#endif
