/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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

#include <iostream>
#include <functional>
#include <vector>
#include <map>
#include <set>
#include <memory>

#ifndef SOLVERS_HPP
#define SOLVERS_HPP

namespace pg {

class Oink;
class Game;
class Solver;

class Solvers final
{
public:
    using SolverConstructor = std::function<std::unique_ptr<Solver>(Oink&, Game&)>;

    ~Solvers() = default;
    Solvers(const Solvers&) = delete;
    Solvers(const Solvers&&) = delete;
    Solvers& operator=(const Solvers&) = delete;
    Solvers& operator=(const Solvers&&) = delete;

    /**
     * Get number of solvers
     */
    static unsigned count()
    {
        return instance().solvers.size();
    }

    /**
     * Get description of solver <id>
     */
    static std::string desc(const std::string& id)
    {
        return instance().solvers[id].description;
    }

    /**
     * Get whether solver <id> can be run in parallel
     */
    static bool isParallel(const std::string& id)
    {
        return instance().solvers[id].isParallel;
    }

    /**
     * Construct solver with the given parameters
     */
    static std::unique_ptr<Solver> construct(const std::string& id, Oink& oink, Game& game);

    /**
     * Write a formatted list of all solvers to the given ostream
     */
    static void list(std::ostream &out);

    /**
     * Add a solver to the set of solvers
     */
    static void add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor);

    static std::set<std::string> getSolverIDs()
    {
        std::set<std::string> ids;
        for (const auto& entry : instance().solvers) {
            ids.insert(entry.first);
        }
        return ids;
    }

private:
    struct SolverInfo {
        std::string description;
        bool isParallel;
        SolverConstructor constructor;
    };

    std::map<std::string, SolverInfo> solvers;

    Solvers();

    static Solvers& instance() {
        static Solvers instance;
        return instance;
    }

    /**
     * Add a solver to the set of solvers
     */
    void _add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor)
    {
        solvers[id] = {description, isParallel, constructor};
    }
};

}

#endif
