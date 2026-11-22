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
 * Unit tests for the renumbering of the region ids in ZLK, ZLKJ, ZLKW and ZLKWJ (see renumber()
 * in their sources). With a tiny limit, the solvers renumber often, also with a deep stack; they
 * must give the same solution and the same counts (the "solved with" line of the log) as with
 * the normal limit, where nothing is renumbered.
 */

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "oink/io.hpp"
#include "oink/oink.hpp"
#include "oink/pgparser.hpp"
#include "oink/solvers.hpp"
#include "solvers/zlk.hpp"
#include "solvers/zlkj.hpp"
#include "solvers/zlkw.hpp"
#include "solvers/zlkwj.hpp"
#include "verifier.hpp"

using namespace pg;
namespace fs = std::filesystem;

namespace {

int limit = 16;          // the region limit of the solvers "tiny-..."
uint64_t renumbered = 0; // how often they renumbered the region ids

/**
 * The solver <S> with region limit <limit>; adds its number of renumberings to <renumbered>.
 */
template <class S>
class Tiny : public S
{
public:
    Tiny(Oink& oink, Game& game) : S(oink, game) { this->region_limit = limit; }

    void run() override
    {
        S::run();
        renumbered += this->renumbered;
    }
};

const char* const solvers[] = {"zlk", "zlkj", "zlkw", "zlkwj"};

void
add_tiny_solvers()
{
    static bool added = false;
    if (added) return;
    added = true;
    // experimental, so that solver choosers hide them
    Solvers::add("tiny-zlk", "ZLK with a tiny region limit", true,
        [] (Oink& oink, Game& game) { return std::make_unique<Tiny<ZLKSolver>>(oink, game); }, "", true);
    Solvers::add("tiny-zlkj", "ZLKJ with a tiny region limit", false,
        [] (Oink& oink, Game& game) { return std::make_unique<Tiny<ZLKJSolver>>(oink, game); }, "", true);
    Solvers::add("tiny-zlkw", "ZLKW with a tiny region limit", false,
        [] (Oink& oink, Game& game) { return std::make_unique<Tiny<ZLKWSolver>>(oink, game); }, "", true);
    Solvers::add("tiny-zlkwj", "ZLKWJ with a tiny region limit", false,
        [] (Oink& oink, Game& game) { return std::make_unique<Tiny<ZLKWJSolver>>(oink, game); }, "", true);
}

struct Result {
    std::vector<int> winner, strategy;
    std::string summary; // the "solved with" line of the log
};

/**
 * Solve a copy of <game> with solver <id> and <workers> (-1: sequential), without
 * preprocessing, and verify the solution.
 */
Result
solve(const Game& game, const std::string& id, int workers)
{
    Game copy(game);
    copy.ensure_sorted();
    std::ostringstream log;
    Oink controller(copy, log, true);
    controller.setRemoveLoops(false);
    controller.setRemoveWCWC(false);
    controller.setSolveSingle(false);
    controller.setWorkers(workers);
    controller.setSolver(id);
    controller.run();
    Verifier(copy, controller.solution(), log).verify(true, true, true);

    Result res;
    for (int v = 0; v < copy.nodecount(); v++) {
        res.winner.push_back(controller.solution().winner(v));
        res.strategy.push_back(controller.solution().strategy(v));
    }
    std::istringstream lines(log.str());
    for (std::string line; std::getline(lines, line); ) {
        const size_t at = line.find("solved with");
        if (at != std::string::npos) res.summary += line.substr(at) + "\n";
    }
    if (res.summary.empty()) throw std::runtime_error(id + " logged no summary");
    return res;
}

/**
 * Solve <game> with each solver, with the normal limit and with each of the tiny <limits>, and
 * compare. The sequential solvers must give the same strategies and summaries; ZLK in parallel
 * the same winners (its strategies and iterations depend on the schedule: the parallel
 * attractor may end a run of heads at a vertex that it attracts at the same time).
 */
void
check(const Game& game, const std::string& name, const std::vector<int>& limits)
{
    for (const char* id : solvers) {
        const Result normal = solve(game, id, -1);
        for (int l : limits) {
            limit = l;
            const Result tiny = solve(game, std::string("tiny-") + id, -1);
            if (tiny.winner != normal.winner or tiny.strategy != normal.strategy or tiny.summary != normal.summary) {
                throw std::runtime_error(std::string(id) + " with region limit " + std::to_string(l) +
                    " differs on " + name + ":\n" + normal.summary + tiny.summary);
            }
        }
    }
    const Result normal = solve(game, "zlk", 2);
    for (int l : limits) {
        limit = l;
        const Result tiny = solve(game, "tiny-zlk", 2);
        if (tiny.winner != normal.winner) {
            throw std::runtime_error("parallel zlk with region limit " + std::to_string(l) + " differs on " + name);
        }
    }
}

/**
 * A random game as made by rngame: <n> vertices with priorities 0..<maxp>, and <mind>..<maxd>
 * distinct successors each (self-loops allowed).
 */
Game
make_rngame(std::mt19937& rng, int n, int maxp, int mind, int maxd)
{
    auto pick = [&](int low, int high) { return std::uniform_int_distribution<int>(low, high)(rng); };
    Game game(n);
    game.vec_init();
    for (int v = 0; v < n; v++) {
        const int priority = pick(0, maxp);
        const int owner = pick(0, 1);
        game.init_vertex(v, priority, owner);
        const int degree = pick(mind, maxd);
        for (int e = 0; e < degree; ) {
            if (game.vec_add_edge(v, pick(0, n-1))) e++;
        }
    }
    game.vec_finish();
    return game;
}

/**
 * The core family of difficult parity games by Benerecetti et al, as made by counter_core <n>:
 * Zielonka's algorithm takes exponentially many steps, with a deep stack.
 */
Game
make_counter_core(int n)
{
    n = 2*n;
    Game game(3+3*n);
    game.vec_init();
    for (int i = 0; i <= n; i++) {
        game.init_vertex(3*i+0, n+1+i, i&1);
        game.init_vertex(3*i+1, i, i&1);
        game.init_vertex(3*i+2, i, 1-(i&1));
        game.vec_add_edge(3*i+0, 3*i+1);
        game.vec_add_edge(3*i+1, 3*i+2);
        game.vec_add_edge(3*i+2, 3*i+1);
        game.vec_add_edge(3*i+2, 3*i+2);
    }
    for (int i = 0; i < n; i++) {
        game.vec_add_edge(3*i+2, 3*i+3);
        game.vec_add_edge(3*i+4, 3*i+0);
    }
    game.vec_finish();
    return game;
}

Game
read_game(const fs::path& path)
{
    auto in = open_input(path.string());
    return PGParser::parse_pgsolver_renumber(*in, false);
}

/**
 * Check all games in the files or directories <paths>.
 */
void
check_files(const std::vector<fs::path>& paths, const std::vector<int>& limits)
{
    std::vector<fs::path> files;
    for (const auto& p : paths) {
        if (fs::is_directory(p)) {
            for (const auto& entry : fs::directory_iterator(p)) {
                if (entry.is_regular_file() and entry.path().extension() != ".sh") files.push_back(entry.path());
            }
        } else {
            files.push_back(p);
        }
    }
    std::sort(files.begin(), files.end());
    if (files.empty()) throw std::runtime_error("no games found");
    for (const auto& f : files) check(read_game(f), f.filename().string(), limits);
}

void
report(const std::string& what, uint64_t before)
{
    if (renumbered == before) throw std::runtime_error("the region ids were never renumbered on " + what);
    std::cout << "Same solutions with tiny region limits on " << what << " (" << (renumbered - before) << " renumberings)." << std::endl;
}

}

int
test_region_renumber()
{
    add_tiny_solvers();

    // random games as made by test_solvers: dense (up to 4 edges per vertex) and sparse
    uint64_t before = renumbered;
    const std::pair<int, int> kinds[] = {{20, 80}, {30, 120}, {50, 200}, {30, 45}};
    for (const auto& kind : kinds) {
        const int size = kind.first, maxe = kind.second;
        for (unsigned seed = 1; seed <= 40; seed++) {
            Game game;
            game.set_random_seed(seed);
            game.init_random_game(size, size, maxe - size);
            check(game, "random game size " + std::to_string(size) + " maxe " + std::to_string(maxe) + " seed " + std::to_string(seed), {1, 16, 64});
        }
    }
    report("random games", before);

    // sparse random games with many priorities, as rngame 300 300 1 2
    before = renumbered;
    std::mt19937 rng(42);
    for (int k = 0; k < 20; k++) check(make_rngame(rng, 300, 300, 1, 2), "rngame 300 300 1 2 #" + std::to_string(k), {1, 16, 64});
    report("sparse random games with many priorities", before);

    // a counter family, with thousands of steps and a deep stack (see also test_region_renumber.cmake)
    before = renumbered;
    for (int n : {4, 6, 8}) check(make_counter_core(n), "counter_core " + std::to_string(n), {1, 16, 64});
    report("counter_core", before);

    // the small games in tests/ and some examples with many steps
    before = renumbered;
    check_files({OINK_TEST_DATA_DIR}, {16});
    report("the games in tests", before);

    before = renumbered;
    std::vector<fs::path> examples;
    for (const char* name : {"ltl2dpa03", "ltl2dpa12", "ltl2dpa13", "ltl2dpa21", "ltl2dba08", "EscalatorSmart", "amba_decomposed_arbiter"}) {
        examples.push_back(fs::path(OINK_TEST_DATA_DIR) / ".." / "examples" / (std::string(name) + ".tlsf.ehoa.pg"));
    }
    check_files(examples, {1, 16});
    report("examples", before);
    return 0;
}

/**
 * The same for the games in the files or directories <paths> (see test_region_renumber.cmake).
 */
int
test_region_renumber_files(const std::vector<std::string>& paths)
{
    add_tiny_solvers();
    const uint64_t before = renumbered;
    check_files(std::vector<fs::path>(paths.begin(), paths.end()), {1, 16, 64});
    report(std::to_string(paths.size()) + " paths", before);
    return 0;
}
