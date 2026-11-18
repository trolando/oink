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
 * Unit tests for the Verifier.
 *
 * The point of these tests is to be confident the verifier never accepts a bad
 * strategy: for each failure mode we hand-build a small game with a deliberately
 * invalid solution and assert that verify() throws. We also assert that valid
 * solutions (single- and multi-strategy) pass. Multi-strategy gets special
 * attention: the verifier must check *every* recorded strategy edge, so a
 * permissive strategy that includes one losing move (e.g. an odd self-loop) must
 * be rejected even though a valid move also exists.
 *
 * Finally, on random games with random strategies in the winning regions, the
 * verifier must agree with the definition: the loser wins a play that follows
 * the strategy if and only if a cycle has a highest priority of its parity.
 */

#include <algorithm>
#include <array>
#include <initializer_list>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

#include "oink/game.hpp"
#include "oink/game_builder.hpp"
#include "oink/oink.hpp"
#include "verifier.hpp"

using namespace pg;

static int failures = 0;

/**
 * A vertex description for building a small test game.
 */
struct V {
    int prio;
    Player owner;
    std::vector<int> edges;
};

/**
 * Build a game from vertex descriptions. Priorities must be given in
 * non-decreasing order so that sorting the game preserves the vertex ids.
 */
static Game
make_game(std::initializer_list<V> verts)
{
    GameBuilder b((int)verts.size());
    int i = 0;
    for (const auto& v : verts) {
        b.set_priority(i, v.prio);
        b.set_owner(i, v.owner);
        for (int to : v.edges) b.add_edge(i, to);
        i++;
    }
    Game g = b.build();
    g.ensure_sorted(); // no-op: priorities are non-decreasing
    return g;
}

/** Add the strategy edge v->to to the (already initialized) multi-strategy. */
static void
add_strat(Game& g, Solution& s, int v, int to)
{
    const int* o = g.outs(v);
    for (int k = 0; o[k] != -1; k++) {
        if (o[k] == to) { s.add_edge_index(g.firstout(v) + k); return; }
    }
    std::cerr << "test bug: no edge " << v << "->" << to << std::endl;
    failures++;
}

/** Assert that verifying <g>/<s> throws (a runtime_error whose message contains <want>). */
static void
expect_reject(const char* name, Game& g, const Solution& s, const std::string& want,
              bool full = true, bool even = true, bool odd = true)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier ver(g, s, log);
    try {
        ver.verify(full, even, odd);
    } catch (std::runtime_error& e) {
        const std::string msg = e.what();
        if (msg.find(want) == std::string::npos) {
            std::cerr << "FAIL [" << name << "] threw \"" << msg
                      << "\" but expected to contain \"" << want << "\"" << std::endl;
            failures++;
        }
        return;
    }
    std::cerr << "FAIL [" << name << "] accepted a bad solution (expected reject \""
              << want << "\")" << std::endl;
    failures++;
}

/** Assert that verifying <g>/<s> succeeds (no throw). */
static void
expect_accept(const char* name, Game& g, const Solution& s,
              bool full = true, bool even = true, bool odd = true)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier ver(g, s, log);
    try {
        ver.verify(full, even, odd);
    } catch (std::runtime_error& e) {
        std::cerr << "FAIL [" << name << "] rejected a good solution: " << e.what() << std::endl;
        failures++;
    }
}

/* -------------------------------------------------------------------------- *
 * Random games, compared with the definition
 * -------------------------------------------------------------------------- */

/**
 * The moves of <v> in the solution graph: the strategy of a vertex won by its owner (the
 * single move if there is one, otherwise the moves of the multi-strategy), and all moves of
 * a vertex won by the other player.
 */
static std::vector<int>
solution_moves(const Game& g, const Solution& s, int v)
{
    std::vector<int> res;
    const bool own = s.winner(v) == g.owner(v);
    if (own and s.strategy(v) != -1) {
        res.push_back(s.strategy(v));
        return res;
    }
    const int* edges = g.outedges();
    for (int e = g.firstout(v); edges[e] != -1; e++) {
        if (!own or s.has_edge_index(e)) res.push_back(edges[e]);
    }
    return res;
}

/**
 * Whether the loser wins a play that follows the strategy in a checked region, by the
 * definition: in the solution graph, a vertex <w> of the region with a priority <p> of
 * the parity of the loser is on a cycle through vertices with a priority of at most <p>.
 * This takes quadratic time.
 */
static bool
loser_wins(const Game& g, const Solution& s, bool even, bool odd)
{
    const int n = g.vertexcount();
    auto checked = [&](int v) {
        return s.is_solved(v) and (s.winner(v) == 0 ? even : odd);
    };
    for (int w = 0; w < n; w++) {
        if (!checked(w)) continue;
        const int p = g.priority(w);
        if ((p&1) == s.winner(w)) continue;
        std::vector<char> seen(n, 0);
        std::vector<int> todo = solution_moves(g, s, w);
        while (!todo.empty()) {
            const int x = todo.back();
            todo.pop_back();
            if (!checked(x) or s.winner(x) != s.winner(w) or g.priority(x) > p) continue;
            if (x == w) return true;
            if (seen[x]) continue;
            seen[x] = 1;
            for (int y : solution_moves(g, s, x)) todo.push_back(y);
        }
    }
    return false;
}

/**
 * A random game with <n> vertices, priorities below <prios> and 1 to <degree> distinct
 * successors per vertex. The vertices are not sorted by priority.
 */
static Game
random_game(std::mt19937& rng, int n, int prios, int degree)
{
    auto below = [&](int k) { return std::uniform_int_distribution<int>(0, k-1)(rng); };
    GameBuilder b(n);
    std::vector<int> targets(n);
    for (int v = 0; v < n; v++) {
        b.set_priority(v, below(prios));
        b.set_owner(v, below(2) ? Player::Odd : Player::Even);
        for (int i = 0; i < n; i++) targets[i] = i;
        std::shuffle(targets.begin(), targets.end(), rng);
        const int count = 1 + below(std::min(n, degree));
        for (int i = 0; i < count; i++) b.add_edge(v, targets[i]);
    }
    return b.build();
}

/**
 * Compare the verifier with the definition on one solution; returns whether the loser wins.
 */
static bool
compare(const char* name, Game& g, const Solution& s, bool full, bool even, bool odd)
{
    const bool expect = loser_wins(g, s, even, odd);
    std::stringstream log;
    bool rejected = false;
    try {
        Verifier(g, s, log).verify(full, even, odd);
    } catch (std::runtime_error& e) {
        if (std::string(e.what()) != "loser can win") {
            std::cerr << "FAIL [" << name << "] unexpected error: " << e.what() << std::endl;
            failures++;
            return expect;
        }
        rejected = true;
    }
    if (rejected != expect) {
        std::cerr << "FAIL [" << name << "] the verifier " << (rejected ? "rejected" : "accepted")
                  << " a solution that the definition " << (expect ? "rejects" : "accepts") << std::endl;
        g.write_pgsolver(std::cerr);
        failures++;
    }
    return expect;
}

/**
 * Random games, solved by ZLK for the winning regions, with random strategies in these
 * regions: a mix of the strategy of ZLK, random single moves and random multi-strategies.
 * Also random games with random claimed regions, changed until they are closed (so that
 * only the cycles decide). For each, the full solution, the even and the odd region alone,
 * and partial solutions without the vertices that reach a random vertex in the solution
 * graph.
 */
static void
test_random_solutions()
{
    std::mt19937 rng(20261116);
    auto below = [&](int k) { return std::uniform_int_distribution<int>(0, k-1)(rng); };
    // per kind of check, how often the definition accepts [0] and rejects [1]
    std::map<std::string, std::array<long, 2>> outcomes;

    for (int game_index = 0; game_index < 4000; game_index++) {
        const bool claimed = game_index >= 3000;
        const int n = 1 + below(game_index % 3000 < 2000 ? 12 : 40);
        const int choices[] = {2, 3, 5, n+1, 2*n};
        const int prios = choices[below(5)];
        Game g = random_game(rng, n, prios, 1 + below(4));
        std::vector<int> winner(n), strategy(n, -1);

        if (!claimed) {
            // solve a sorted copy; mapping[i] is the vertex of <g> at position i of the copy
            Game sorted(g);
            std::vector<int> mapping(n);
            sorted.sort(mapping.data());
            std::ostringstream log;
            Oink controller(sorted, log, true);
            controller.setSolver("zlk");
            controller.run();
            const Solution& zs = controller.solution();
            for (int i = 0; i < n; i++) {
                if (sorted.priority(i) != g.priority(mapping[i]) or sorted.owner(i) != g.owner(mapping[i])) {
                    std::cerr << "FAIL [random] test bug: unexpected mapping of the sorted game" << std::endl;
                    failures++;
                    return;
                }
                winner[mapping[i]] = zs.winner(i);
                if (zs.strategy(i) != -1) strategy[mapping[i]] = mapping[zs.strategy(i)];
            }
        } else {
            // random claims, changed until each region is closed: a vertex won by its owner
            // has a move in the region, a vertex won by the opponent only has such moves
            for (int v = 0; v < n; v++) winner[v] = below(2);
            bool closed = false;
            for (int iteration = 0; iteration < 100 and !closed; iteration++) {
                closed = true;
                for (int v = 0; v < n; v++) {
                    bool some = false, all = true;
                    for (int to : g.out_edges(v)) {
                        if (winner[to] == winner[v]) some = true;
                        else all = false;
                    }
                    if (g.owner(v) == winner[v] ? !some : !all) {
                        winner[v] = 1 - winner[v];
                        closed = false;
                    }
                }
            }
            if (!closed) continue;
        }
        const std::string kind = claimed ? "claimed/" : "random/";
        auto check = [&](const std::string& what, const Solution& s, bool full, bool even, bool odd) {
            const std::string name = kind + what;
            outcomes[name][compare(name.c_str(), g, s, full, even, odd) ? 1 : 0]++;
        };

        const int* edges = g.outedges();
        for (int round = 0; round < 6; round++) {
            // keep the strategy of ZLK with probability <keep> in 10, else choose at random
            const int keep = round == 0 ? 10 : below(10);
            const bool multi = round > 0 and below(2);
            std::vector<std::vector<int>> moves(n); // the chosen edge indices
            for (int v = 0; v < n; v++) {
                if (g.owner(v) != winner[v]) continue;
                std::vector<int> good;
                for (int e = g.firstout(v); edges[e] != -1; e++) {
                    if (winner[edges[e]] == winner[v]) good.push_back(e);
                }
                if (good.empty()) {
                    std::cerr << "FAIL [random] test bug: no move in the winning region" << std::endl;
                    failures++;
                    return;
                }
                if (below(10) < keep and strategy[v] != -1) {
                    for (int e : good) if (edges[e] == strategy[v]) moves[v].push_back(e);
                } else if (multi and below(2)) {
                    for (int e : good) if (below(2)) moves[v].push_back(e);
                    if (moves[v].empty()) moves[v].push_back(good[below((int)good.size())]);
                } else {
                    moves[v].push_back(good[below((int)good.size())]);
                }
            }

            // the solution without the vertices in <removed>
            auto make = [&](const std::vector<char>& removed) {
                Solution s(n);
                if (multi) s.init_multi(g.edgeArraySize());
                for (int v = 0; v < n; v++) {
                    if (removed[v]) continue;
                    if (moves[v].size() == 1 and (!multi or below(2))) {
                        s.solve(v, winner[v], edges[moves[v][0]]);
                    } else {
                        s.solve(v, winner[v], -1);
                        for (int e : moves[v]) s.add_edge_index(e);
                    }
                }
                return s;
            };

            const std::vector<char> none(n, 0);
            const Solution s = make(none);
            check("full", s, true, true, true);
            check("even", s, true, true, false);
            check("odd", s, true, false, true);

            // remove the vertices that reach a random vertex in the solution graph
            std::vector<char> removed(n, 0);
            removed[below(n)] = 1;
            for (bool changed = true; changed; ) {
                changed = false;
                for (int v = 0; v < n; v++) {
                    if (removed[v]) continue;
                    for (int to : solution_moves(g, s, v)) {
                        if (removed[to]) {
                            removed[v] = 1;
                            changed = true;
                            break;
                        }
                    }
                }
            }
            if (std::count(removed.begin(), removed.end(), 1) == n) continue;
            const Solution part = make(removed);
            check("partial", part, false, true, true);
            check("partial even", part, false, true, false);
            check("partial odd", part, false, false, true);
        }
    }

    // the test only means something if each kind of check has both outcomes often
    for (const char* kind : {"random/", "claimed/"}) {
        for (const char* what : {"full", "even", "odd", "partial", "partial even", "partial odd"}) {
            const std::string name = std::string(kind) + what;
            const auto& count = outcomes[name];
            if (count[0] < 500 or count[1] < 500) {
                std::cerr << "FAIL [" << name << "] too few accepted (" << count[0]
                          << ") or rejected (" << count[1] << ") solutions" << std::endl;
                failures++;
            }
        }
    }
}

int
test_verifier()
{
    /* ---------------------------------------------------------------------- *
     * Single-strategy failure modes
     * ---------------------------------------------------------------------- */

    // GOOD: two even vertices in an even cycle, each with a valid strategy.
    {
        Game g = make_game({
            {0, Player::Even, {1}},
            {0, Player::Even, {0}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 1);
        s.solve(1, 0, 0);
        expect_accept("single/good even cycle", g, s);
    }

    // Each invocation counts the strategies it checks, even on a reused verifier.
    {
        Game g = make_game({
            {0, Player::Even, {0, 1}},
            {0, Player::Even, {0}},
            {1, Player::Odd, {2}},
        });
        Solution s(g.nodecount());
        s.init_multi(g.edgeArraySize());
        s.solve(0, 0, -1);
        add_strat(g, s, 0, 0);
        add_strat(g, s, 0, 1);
        s.solve(1, 0, 0);
        s.solve(2, 1, 2);
        std::stringstream log;
        Verifier ver(g, s, log);
        ver.checkClosure();
        const int closure = ver.numberOfStrategies();
        ver.verify();
        const int full = ver.numberOfStrategies();
        ver.verify(true, false, true);
        if (closure != 3 or full != 3 or ver.numberOfStrategies() != 1) {
            std::cerr << "FAIL [verifier reuse]: wrong strategy count" << std::endl;
            failures++;
        }
    }

    // BAD: a winning (owner == winner) vertex with no strategy.
    {
        Game g = make_game({
            {0, Player::Even, {1}},
            {0, Player::Even, {0}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1); // winner Even == owner, but strategy -1
        s.solve(1, 0, 0);
        expect_reject("single/winner no strategy", g, s, "no strategy");
    }

    // BAD: strategy that is not an outgoing edge.
    {
        Game g = make_game({
            {0, Player::Even, {0}}, // only a self-loop
            {0, Player::Even, {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 1); // 1 is not a successor of 0
        s.solve(1, 0, 1);
        expect_reject("single/strategy not a move", g, s, "not a valid move");
    }

    // BAD: strategy that leaves the dominion (points into the other player's region).
    {
        Game g = make_game({
            {0, Player::Even, {0, 1}}, // even self-loop, plus edge to the odd region
            {1, Player::Odd,  {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 1); // Even claims to win 0 but moves to Odd-won 1
        s.solve(1, 1, 1);
        expect_reject("single/strategy leaves dominion", g, s, "leaves dominion");
    }

    // BAD: a losing vertex that nonetheless has a strategy recorded.
    {
        Game g = make_game({
            {1, Player::Even, {1}}, // Even-owned but forced into the odd region
            {1, Player::Odd,  {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 1, -1); // 0 is won by Odd (owner Even is the loser)
        s.solve(1, 1, 1);
        s.strategy_data()[0] = 1; // force a strategy on the losing vertex
        expect_reject("single/loser has strategy", g, s, "losing vertex has strategy");
    }

    // BAD: a losing vertex that can escape to a vertex won by the other player.
    {
        Game g = make_game({
            {1, Player::Even, {1, 2}}, // claimed Odd-won, but has an edge to Even-won 2
            {1, Player::Odd,  {1}},
            {2, Player::Even, {2}},
        });
        Solution s(g.nodecount());
        s.solve(0, 1, -1); // claim 0 won by Odd
        s.solve(1, 1, 1);
        s.solve(2, 0, 2);
        expect_reject("single/loser can escape", g, s, "loser can escape");
    }

    // BAD: an odd self-loop claimed as an Even win. Passes the local checks (the
    // self-loop stays in the "dominion") but the SCC has odd top priority, so the
    // loser actually wins the cycle.
    {
        Game g = make_game({
            {1, Player::Even, {0}}, // odd priority, self-loop
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 0); // Even claims to win via the odd self-loop
        expect_reject("single/loser can win (odd self-loop)", g, s, "loser can win");
    }

    /* ---------------------------------------------------------------------- *
     * Multi-strategy failure modes
     * ---------------------------------------------------------------------- */

    // GOOD: a vertex with two valid winning moves recorded.
    {
        Game g = make_game({
            {2, Player::Even, {1, 2}},
            {2, Player::Even, {1}},
            {2, Player::Even, {2}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1); // -1 sentinel: moves are in the multi-strategy
        s.solve(1, 0, -1);
        s.solve(2, 0, -1);
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 1); add_strat(g, s, 0, 2); // both moves are winning
        add_strat(g, s, 1, 1);
        add_strat(g, s, 2, 2);
        expect_accept("multi/good two moves", g, s);
    }

    // BAD: one of the recorded strategy edges leaves the dominion.
    {
        Game g = make_game({
            {1, Player::Odd,  {0}},    // Odd-won self-loop
            {2, Player::Even, {1}},    // Even-won self-loop
            {2, Player::Even, {1, 0}}, // claimed Even, but one move goes to Odd-won 0
        });
        Solution s(g.nodecount());
        s.solve(0, 1, -1);
        s.solve(1, 0, -1);
        s.solve(2, 0, -1);
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 0);
        add_strat(g, s, 1, 1);
        add_strat(g, s, 2, 1); add_strat(g, s, 2, 0); // 2->0 leaves the dominion
        expect_reject("multi/edge leaves dominion", g, s, "leaves dominion");
    }

    // The self-loop trap, as a multi-strategy. Same game solved two ways:
    //   v0 (odd, Even-owned) has a self-loop AND an escape to the even region.
    // The escape (0->1) is winning; the self-loop (0->0) is a losing move.
    auto trap_game = [] {
        return make_game({
            {1, Player::Even, {0, 1}}, // odd priority: self-loop 0->0 and escape 0->1
            {2, Player::Even, {1}},    // even self-loop, Even-won
        });
    };

    // BAD: the multi-strategy includes the losing self-loop. Both recorded edges
    // stay in the dominion (0 is Even-won), so the local check passes, but the
    // self-loop forms an odd cycle the loser wins -> must be rejected.
    {
        Game g = trap_game();
        Solution s(g.nodecount());
        s.solve(0, 0, -1);
        s.solve(1, 0, -1);
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 0); add_strat(g, s, 0, 1); // self-loop included: BAD
        add_strat(g, s, 1, 1);
        expect_reject("multi/self-loop trap rejected", g, s, "loser can win");
    }

    // GOOD: same game, but the multi-strategy records only the winning escape.
    {
        Game g = trap_game();
        Solution s(g.nodecount());
        s.solve(0, 0, -1);
        s.solve(1, 0, -1);
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 1); // only the escape: GOOD
        add_strat(g, s, 1, 1);
        expect_accept("multi/escape only accepted", g, s);
    }

    // Mixed: single and multi vertices coexist in one solution. The dispatch is
    // per vertex via the sentinel: vertex 0 (-1) reads the multi-strategy, vertex
    // 1 (strategy set) reads the single strategy. Both valid -> accept.
    {
        Game g = make_game({
            {2, Player::Even, {0}},
            {2, Player::Even, {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1);  // multi
        s.solve(1, 0, 1);   // single (e.g. produced by a preprocessor/attractor)
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 0);
        expect_accept("multi/mixed single and multi ok", g, s);
    }

    // Single and multi on one cycle: 0 moves to 1 (single), 1 moves to 0 and 2 (multi). The
    // cycle 0 -> 1 -> 0 has highest priority 1, so Odd wins it; without the move 1 -> 0, the
    // only cycle is the even self-loop of 2.
    for (bool back : {true, false}) {
        Game g = make_game({
            {0, Player::Even, {1}},
            {1, Player::Even, {0, 2}},
            {2, Player::Even, {2}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 1);  // single
        s.solve(1, 0, -1); // multi
        s.solve(2, 0, 2);  // single
        s.init_multi(g.edgeArraySize());
        if (back) add_strat(g, s, 1, 0);
        add_strat(g, s, 1, 2);
        if (back) expect_reject("multi/mixed single and multi on an odd cycle", g, s, "loser can win");
        else expect_accept("multi/mixed single and multi on an even cycle", g, s);
    }

    // A vertex won by its owner with the -1 sentinel but no recorded multi edge
    // must be rejected (it claims a win with no move at all).
    {
        Game g = make_game({
            {2, Player::Even, {0}},
            {2, Player::Even, {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1);
        s.solve(1, 0, -1); // winning, sentinel, but no multi edge recorded
        s.init_multi(g.edgeArraySize());
        add_strat(g, s, 0, 0);
        expect_reject("multi/sentinel without edge rejected", g, s, "no strategy");
    }

    // BAD: the odd self-loop of 0 is in an SCC with the higher even priority of 1, so the
    // loser wins only a cycle below the highest priority of the SCC.
    {
        Game g = make_game({
            {1, Player::Odd,  {0, 1}}, // odd self-loop, and a move to 1
            {2, Player::Even, {0}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1);
        s.solve(1, 0, 0);
        expect_reject("single/loser can win below a higher even priority", g, s, "loser can win");
    }

    // BAD: the odd cycle 0 <-> 1 is in an SCC with the higher even priority of 2.
    {
        Game g = make_game({
            {0, Player::Odd,  {1, 2}},
            {1, Player::Odd,  {0, 2}},
            {2, Player::Even, {0}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, -1);
        s.solve(1, 0, -1);
        s.solve(2, 0, 0);
        expect_reject("single/odd cycle below a higher even priority", g, s, "loser can win");
    }

    // GOOD: a path of even vertices E0..E20 (vertex 2i, priority 2i, E0 with a self-loop) in
    // both directions, and odd vertices O1..O20 (vertex 2i-1, priority 2i-1) with only the
    // edges Ei <-> Oi. Every cycle through Oi passes Ei, so Even wins everywhere.
    {
        const int k = 20;
        GameBuilder b(2*k+1);
        for (int v = 0; v <= 2*k; v++) {
            b.set_priority(v, v);
            b.set_owner(v, Player::Odd);
        }
        b.add_edge(0, 0);
        for (int i = 1; i <= k; i++) {
            b.add_edge(2*i-2, 2*i);
            b.add_edge(2*i, 2*i-2);
            b.add_edge(2*i, 2*i-1);
            b.add_edge(2*i-1, 2*i);
        }
        Game g = b.build();
        Solution s(g.nodecount());
        for (int v = 0; v <= 2*k; v++) s.solve(v, 0, -1);
        expect_accept("single/odd cycles only through higher even priorities", g, s);
    }

    // GOOD: the game of counter_ortl 3000, with 3000 odd priorities that Even wins and a vertex
    // with 3001 moves, solved by TL. This took seconds while the search scanned the moves of a
    // vertex from the start each time it returned to it.
    {
        const int n = 3000;
        GameBuilder b(1 + 5*n);
        b.set_priority(0, 0);
        b.set_owner(0, Player::Odd);
        b.add_edge(0, 0);
        for (int i = 0; i < n; i++) {
            const int c = 5*i+1;
            const int d = 2*(n+i+1);
            const int prios[] = {2*(n-i), d, 1, 1, d+1};
            const Player owners[] = {Player::Even, Player::Even, Player::Odd, Player::Even, Player::Even};
            for (int k = 0; k < 5; k++) {
                b.set_priority(c+k, prios[k]);
                b.set_owner(c+k, owners[k]);
            }
            b.add_edge(c, 0);
            b.add_edge(0, c);
            b.add_edge(c, c+1);
            b.add_edge(c+1, c+2);
            b.add_edge(c+2, c+3);
            b.add_edge(c+3, c+2);
            b.add_edge(c+3, c+4);
            b.add_edge(c+4, c);
        }
        Game g = b.build();
        g.ensure_sorted();
        std::ostringstream log;
        Oink controller(g, log, true);
        controller.setSolver("tl");
        controller.run();
        expect_accept("single/many odd priorities and a vertex with many moves", g, controller.solution());
    }

    /* ---------------------------------------------------------------------- *
     * Partial solutions and the even and odd flags
     * ---------------------------------------------------------------------- */

    // An Even-won self-loop and an Odd-won vertex without a strategy. Checking only the even
    // region accepts; checking the odd region rejects.
    {
        Game g = make_game({
            {0, Player::Even, {0}},
            {1, Player::Odd,  {1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 0);
        s.solve(1, 1, -1);
        expect_accept("flags/even only", g, s, true, true, false);
        expect_reject("flags/odd only", g, s, "no strategy", true, false, true);
    }

    // An unsolved vertex: a full solution is rejected, a partial one is accepted.
    {
        Game g = make_game({
            {0, Player::Even, {0}},
            {1, Player::Even, {0, 1}},
        });
        Solution s(g.nodecount());
        s.solve(0, 0, 0);
        expect_reject("flags/full solution", g, s, "not every vertex is won");
        expect_accept("flags/partial solution", g, s, false, true, true);
    }

    test_random_solutions();

    if (failures) {
        std::cerr << failures << " verifier test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all verifier tests passed" << std::endl;
    return 0;
}
