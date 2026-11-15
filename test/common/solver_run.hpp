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

#ifndef OINK_TEST_SOLVER_RUN_HPP
#define OINK_TEST_SOLVER_RUN_HPP

#include <string>

/**
 * A solver to test or benchmark: its id, its configuration (see
 * Solver::parseOptions; empty for the defaults), and the name it is reported
 * under, which is the id, or id[configuration] for a configured variant.
 */
struct SolverRun
{
    std::string id;
    std::string config;
    std::string name;

    SolverRun() = default;
    SolverRun(const std::string& id, const std::string& config) : id(id), config(config), name(name_of(id, config)) { }

    static std::string name_of(const std::string& id, const std::string& config)
    {
        return config.empty() ? id : id + "[" + config + "]";
    }
};

#endif
