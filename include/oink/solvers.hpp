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

    enum class Family {
        TangleLearning, ProgressMeasures, FixpointIteration, PriorityPromotion,
        StrategyImprovement, Zielonka, Other
    };
    static Family family(const std::string& id)
    {
        const auto it = instance().solvers.find(id);
        return it == instance().solvers.end() ? Family::Other : it->second.family;
    }
    static std::string familyName(Family family);
    static std::vector<std::string> getSolverIDsByFamily();

    /** Configuration syntax for this solver; empty when it accepts no options. */
    static std::string optionsHelp(const std::string& id)
    {
        const auto it = instance().solvers.find(id);
        return it == instance().solvers.end() ? std::string() : it->second.optionsHelp;
    }

    /** Validate with the solver parser, without solving or modifying an input game. */
    static bool validateOptions(const std::string& id, const std::string& options, std::ostream& errors);

    static bool hasOptions(const std::string& id) { return !optionsHelp(id).empty(); }

    /** Experimental entries may be hidden by solver-selection interfaces. */
    static bool isExperimental(const std::string& id)
    {
        const auto it = instance().solvers.find(id);
        return it != instance().solvers.end() && it->second.experimental;
    }

    /**
     * Get whether solver <id> can be run in parallel
     */
    static bool isParallel(const std::string& id)
    {
        return instance().solvers[id].isParallel;
    }

    /**
     * Get the number of (unsolved) vertices from which parallel solver <id> runs in parallel;
     * on smaller games starting Lace costs more than it saves, so they are solved sequentially
     * (see Oink::setForceParallel)
     */
    static long parallelMinVertices(const std::string& id)
    {
        auto it = instance().solvers.find(id);
        return it == instance().solvers.end() ? 0 : it->second.parallelMinVertices;
    }

    /**
     * Get the solver that solves a game instead of parallel solver <id> when it runs
     * sequentially, or "" if <id> itself does (it then runs its sequential code)
     */
    static std::string sequentialFallback(const std::string& id)
    {
        auto it = instance().solvers.find(id);
        return it == instance().solvers.end() ? std::string() : it->second.sequentialFallback;
    }

    /**
     * Get whether solver <id> uses the game's in-edges (Game::ins()); true if unknown
     */
    static bool needsInEdges(const std::string& id)
    {
        auto it = instance().solvers.find(id);
        return it == instance().solvers.end() or it->second.needsInEdges;
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
    static void add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor, const std::string& optionsHelp, bool experimental = false);

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
        bool experimental = false;
        bool needsInEdges = true;
        long parallelMinVertices = 100000; // (only for parallel solvers)
        std::string sequentialFallback = "";
        std::string optionsHelp = "";
        Family family = Family::Other;
        int familyOrder = 0;
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
    void _add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor, bool experimental = false)
    {
        solvers[id] = {description, isParallel, constructor, experimental};
    }
};

}

#endif
