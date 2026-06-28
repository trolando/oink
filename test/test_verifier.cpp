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
 */

#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

#include "oink/game.hpp"
#include "oink/game_builder.hpp"
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
 * non-decreasing order so the game is already sorted (the verifier requires a
 * sorted game) and vertex ids are preserved.
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
add_strat(Game& g, int v, int to)
{
    const int* o = g.outs(v);
    for (int k = 0; o[k] != -1; k++) {
        if (o[k] == to) { g.addStrategyEdge(v, k); return; }
    }
    std::cerr << "test bug: no edge " << v << "->" << to << std::endl;
    failures++;
}

/** Assert that verifying <g> throws (a runtime_error whose message contains <want>). */
static void
expect_reject(const char* name, Game& g, const std::string& want)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier ver(g, log);
    try {
        ver.verify(true, true, true);
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

/** Assert that verifying <g> succeeds (no throw). */
static void
expect_accept(const char* name, Game& g)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier ver(g, log);
    try {
        ver.verify(true, true, true);
    } catch (std::runtime_error& e) {
        std::cerr << "FAIL [" << name << "] rejected a good solution: " << e.what() << std::endl;
        failures++;
    }
}

int
main()
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
        g.solve(0, 0, 1);
        g.solve(1, 0, 0);
        expect_accept("single/good even cycle", g);
    }

    // BAD: a winning (owner == winner) vertex with no strategy.
    {
        Game g = make_game({
            {0, Player::Even, {1}},
            {0, Player::Even, {0}},
        });
        g.solve(0, 0, -1); // winner Even == owner, but strategy -1
        g.solve(1, 0, 0);
        expect_reject("single/winner no strategy", g, "no strategy");
    }

    // BAD: strategy that is not an outgoing edge.
    {
        Game g = make_game({
            {0, Player::Even, {0}}, // only a self-loop
            {0, Player::Even, {1}},
        });
        g.solve(0, 0, 1); // 1 is not a successor of 0
        g.solve(1, 0, 1);
        expect_reject("single/strategy not a move", g, "not a valid move");
    }

    // BAD: strategy that leaves the dominion (points into the other player's region).
    {
        Game g = make_game({
            {0, Player::Even, {0, 1}}, // even self-loop, plus edge to the odd region
            {1, Player::Odd,  {1}},
        });
        g.solve(0, 0, 1); // Even claims to win 0 but moves to Odd-won 1
        g.solve(1, 1, 1);
        expect_reject("single/strategy leaves dominion", g, "leaves dominion");
    }

    // BAD: a losing vertex that nonetheless has a strategy recorded.
    {
        Game g = make_game({
            {1, Player::Even, {1}}, // Even-owned but forced into the odd region
            {1, Player::Odd,  {1}},
        });
        g.solve(0, 1, -1); // 0 is won by Odd (owner Even is the loser)
        g.solve(1, 1, 1);
        g.getStrategy()[0] = 1; // force a strategy on the losing vertex
        expect_reject("single/loser has strategy", g, "losing vertex has strategy");
    }

    // BAD: a losing vertex that can escape to a vertex won by the other player.
    {
        Game g = make_game({
            {1, Player::Even, {1, 2}}, // claimed Odd-won, but has an edge to Even-won 2
            {1, Player::Odd,  {1}},
            {2, Player::Even, {2}},
        });
        g.solve(0, 1, -1); // claim 0 won by Odd
        g.solve(1, 1, 1);
        g.solve(2, 0, 2);
        expect_reject("single/loser can escape", g, "loser can escape");
    }

    // BAD: an odd self-loop claimed as an Even win. Passes the local checks (the
    // self-loop stays in the "dominion") but the SCC has odd top priority, so the
    // loser actually wins the cycle.
    {
        Game g = make_game({
            {1, Player::Even, {0}}, // odd priority, self-loop
        });
        g.solve(0, 0, 0); // Even claims to win via the odd self-loop
        expect_reject("single/loser can win (odd self-loop)", g, "loser can win");
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
        g.solve(0, 0, 1);
        g.solve(1, 0, 1);
        g.solve(2, 0, 2);
        g.initMultiStrategy();
        add_strat(g, 0, 1); add_strat(g, 0, 2); // both moves are winning
        add_strat(g, 1, 1);
        add_strat(g, 2, 2);
        expect_accept("multi/good two moves", g);
    }

    // BAD: one of the recorded strategy edges leaves the dominion.
    {
        Game g = make_game({
            {1, Player::Odd,  {0}},    // Odd-won self-loop
            {2, Player::Even, {1}},    // Even-won self-loop
            {2, Player::Even, {1, 0}}, // claimed Even, but one move goes to Odd-won 0
        });
        g.solve(0, 1, 0);
        g.solve(1, 0, 1);
        g.solve(2, 0, 1);
        g.initMultiStrategy();
        add_strat(g, 0, 0);
        add_strat(g, 1, 1);
        add_strat(g, 2, 1); add_strat(g, 2, 0); // 2->0 leaves the dominion
        expect_reject("multi/edge leaves dominion", g, "leaves dominion");
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
        g.solve(0, 0, 1);
        g.solve(1, 0, 1);
        g.initMultiStrategy();
        add_strat(g, 0, 0); add_strat(g, 0, 1); // self-loop included: BAD
        add_strat(g, 1, 1);
        expect_reject("multi/self-loop trap rejected", g, "loser can win");
    }

    // GOOD: same game, but the multi-strategy records only the winning escape.
    {
        Game g = trap_game();
        g.solve(0, 0, 1);
        g.solve(1, 0, 1);
        g.initMultiStrategy();
        add_strat(g, 0, 1); // only the escape: GOOD
        add_strat(g, 1, 1);
        expect_accept("multi/escape only accepted", g);
    }

    // Fallback: with a multi-strategy present, a winning vertex that has no
    // recorded multi edge (e.g. solved by a preprocessor) is checked against its
    // single strategy. A valid single strategy should still pass.
    {
        Game g = make_game({
            {2, Player::Even, {0}},
            {2, Player::Even, {1}},
        });
        g.solve(0, 0, 0);
        g.solve(1, 0, 1);
        g.initMultiStrategy();
        add_strat(g, 0, 0); // only vertex 0 has a multi edge; vertex 1 falls back
        expect_accept("multi/fallback to single ok", g);
    }

    // Fallback: a vertex with no multi edge and an invalid single strategy (-1)
    // must still be rejected.
    {
        Game g = make_game({
            {2, Player::Even, {0}},
            {2, Player::Even, {1}},
        });
        g.solve(0, 0, 0);
        g.solve(1, 0, -1); // winning but no strategy
        g.initMultiStrategy();
        add_strat(g, 0, 0);
        expect_reject("multi/fallback no strategy", g, "no strategy");
    }

    if (failures) {
        std::cerr << failures << " verifier test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all verifier tests passed" << std::endl;
    return 0;
}
