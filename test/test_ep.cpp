// SPDX-License-Identifier: Apache-2.0
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "solvers/ep_search.hpp"
#include "solvers/ep_arena.hpp"
#include "solvers/energy_credit.hpp"
#include "verifier.hpp"

int test_ep_experts();
int test_ep_market();

namespace {
using namespace pg;
using namespace pg::ep;

std::vector<int> solve_game(Game game, const std::string& solver,
                            const std::string& options = {}, bool scc = false, int solved = -1,
                            std::string* output = nullptr)
{
    game.ensure_sorted();
    std::ostringstream log;
    Oink oink(game, log, false);
    oink.setRemoveLoops(false);
    oink.setRemoveWCWC(false);
    oink.setSolveSingle(false);
    oink.setBottomSCC(scc);
    oink.setWorkers(-1);
    if (output) oink.setTrace(1);
    oink.setSolver(solver);
    oink.setSolverOptions(options);
    if (solved >= 0) oink.solution().solve(solved, 0, solved);
    oink.run();
    Verifier(game, oink.solution(), log).verify(true, true, true);
    if (output) *output = log.str();
    std::vector<int> winners;
    for (int i = 0; i < game.nodecount(); ++i) winners.push_back(oink.solution().winner(i));
    return winners;
}

void exhaustive_games()
{
    for (int first = 0; first < 4; ++first) {
        for (int second = 0; second < 4; ++second) {
            for (int owners = 0; owners < 4; ++owners) {
                for (int mask0 = 1; mask0 <= 3; ++mask0) {
                    for (int mask1 = 1; mask1 <= 3; ++mask1) {
                        Game game(2);
                        for (int v = 0; v < 2; ++v) {
                            game.init_vertex(v, v == 0 ? first : second, (owners >> v) & 1);
                            game.e_start(v);
                            const int mask = v == 0 ? mask0 : mask1;
                            for (int w = 0; w < 2; ++w) if (mask & (1 << w)) game.e_add(v, w);
                            game.e_finish();
                        }
                        const auto expected = solve_game(game, "zlk");
                        for (const char* solver : {"ep-clipped", "ep-sign"}) {
                            if (solve_game(game, solver) != expected)
                                throw std::runtime_error("EP disagrees with Zielonka");
                        }
                        for (const char* solver : {"gpm", "pma", "pmtl"}) {
                            if (solve_game(game, solver, "energy") != expected)
                                throw std::runtime_error("energy progress measure solver disagrees with Zielonka");
                        }
                    }
                }
            }
        }
    }
}

void partial_game()
{
    Game game(5);
    const std::vector<std::vector<int>> edges{{1, 2}, {0}, {2, 3}, {3}, {4}};
    for (int v = 0; v < 5; ++v) {
        game.init_vertex(v, v, v & 1);
        game.e_start(v);
        for (int w : edges[v]) game.e_add(v, w);
        game.e_finish();
    }
    const auto expected = solve_game(game, "zlk");
    for (bool scc : {false, true}) {
        if (solve_game(game, "ep", "", scc, 4) != expected)
            throw std::runtime_error("EP failed on a partial solution");
        for (const char* solver : {"gpm", "pma", "pmtl"}) {
            if (solve_game(game, solver, "energy", scc, 4) != expected)
                throw std::runtime_error("energy progress measures failed on a partial solution");
        }
    }
}

void signed_credits()
{
    Game game(3);
    for (int v = 0; v < 3; ++v) game.init_vertex(v, v, v & 1);
    EnergyCredits pm(game, 0);
    // The signed weights are +1, -2, +3 and the largest finite credit is 4.
    const Integer at_head = pm.progress(0, 2);
    if (at_head != 3 || pm.progress(at_head, 1) != 1)
        throw std::runtime_error("a lower odd priority should spend credit after the head");
    if (pm.progress(pm.progress(1, 1), 2) != 3 || pm.progress(1, 2) != 4)
        throw std::runtime_error("a lower odd priority should spend credit before the head");
    pm.raise(0, pm.progress(2, 2));
    pm.see(0, 1);
    if (!pm.is_top(0)) throw std::runtime_error("energy measure Top must absorb a debit");

    Game many(181);
    for (int v = 0; v < 181; ++v) many.init_vertex(v, v, v & 1);
    EnergyCredits wide(many, 0);
    wide.see(0, 180);
    const auto digits = wide.store(0);
    wide.load(1, digits);
    if (digits.size() <= 8 || wide.is_top(0) || !wide.eq(0, 1))
        throw std::runtime_error("energy credit lost precision on serialization");
}

void credit_handoff()
{
    Game game(5);
    const int highest = std::numeric_limits<int>::max();
    const int priorities[] = {0, 1, highest-3, highest-2, highest};
    const int owners[] = {0, 0, 1, 0, 1};
    const std::vector<std::vector<int>> edges{{1, 4, 0}, {2}, {3, 0}, {4, 3}, {0, 3}};
    for (int v = 0; v < 5; ++v) {
        game.init_vertex(v, priorities[v], owners[v]);
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    const std::vector<int> expected{0, 1, 1, 1, 1};
    std::string trace;
    if (solve_game(game, "gpa-energy", "no-shortcuts", false, -1, &trace) != expected ||
        trace.find("handing off 4 vertices after 5 paired credit rounds") == std::string::npos)
        throw std::runtime_error("PMA did not preserve the solved self-loop through its EP handoff");
    if (solve_game(game, "pmtl-energy", "no-shortcuts", false, -1, &trace) != expected ||
        trace.find("energy tangle") == std::string::npos)
        throw std::runtime_error("PMTL did not use a learned tangle to raise a credit");
    if (solve_game(game, "pma", "energy,up") != expected ||
        solve_game(game, "pmtl-energy", "up") != expected)
        throw std::runtime_error("energy measure up pass failed on sparse priorities");
    if (solve_game(game, "gpm-energy", "no-shortcuts", false, -1, &trace) != expected ||
        trace.find("credit wave limit reached") == std::string::npos)
        throw std::runtime_error("GPM did not switch to EP at its wave limit");
    for (const char* solver : {"gpm-energy", "gpa-energy", "pmtl-energy"}) {
        std::ostringstream errors;
        if (!Solvers::validateOptions(solver, "epm,no-shortcuts", errors) ||
            Solvers::validateOptions(solver, "small", errors))
            throw std::runtime_error("a named energy variant must retain its energy measure");
    }
}

void open_tangle()
{
    // At zero price, the negative weight at 2 prevents a certificate for
    // 2->0. Attraction to 0 resolves 2 and the opponent exit from 3.
    Game game(4);
    const std::vector<int> priorities{2, 1, 1, 0}, owners{0, 1, 0, 1};
    const std::vector<std::vector<int>> edges{{0}, {1}, {1, 0}, {3, 2}};
    for (int v = 0; v < 4; ++v) {
        game.init_vertex(v, priorities[v], owners[v]);
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    auto zero_price_dominions = [&](const bitset& disabled) {
        const ep::Arena arena(game, disabled);
        std::vector<int> moves(arena.size());
        const ep::EnergyMeasure measure(arena);
        const Prices zero(arena.size());
        auto seeds = arena.supported_dominions(measure.weight, zero, zero, 1, moves);
        std::vector<int> result(4, -1);
        for (int v = 0; v < arena.size(); ++v) result[arena.vertices[v]] = seeds.at(v);
        return result;
    };
    bitset removed(4);
    if (zero_price_dominions(removed) != std::vector<int>{0, 1, -1, -1})
        throw std::runtime_error("EP accepted a tangle with an unresolved escape");
    for (int v : {0, 1, 2}) removed[v] = true;
    if (zero_price_dominions(removed) != std::vector<int>{-1, -1, -1, 0})
        throw std::runtime_error("EP missed a tangle whose escape is solved");
}

void supported_move_choice()
{
    // At constant prices 0->1 is a greedy tie but loses for Even. The
    // supported safety subgame contains {0,2} and supplies the move 0->2.
    const std::vector<std::vector<int>> edges{{1, 2}, {1}, {2, 0}};
    for (int dual : {0, 1}) {
        Game game(3);
        for (int v = 0; v < 3; ++v) {
            game.init_vertex(v, (v == 1) + dual, (v != 0) ^ dual);
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        const ep::Arena arena(game, bitset(3));
        const ep::EnergyMeasure measure(arena);
        const Prices zero(3);
        std::vector<int> moves(3);
        const auto winners = arena.supported_dominions(measure.weight, zero, zero, 1, moves);
        if (winners.size() != 3 || !arena.certifies(winners, moves))
            throw std::runtime_error("invalid supported safety strategy");
        for (int v = 0; v < 3; ++v) {
            if (winners[v] != ((arena.vertices[v] == 1) ^ dual))
                throw std::runtime_error("EP missed a supported move");
            if (arena.vertices[v] == 0 && arena.vertices[moves[v]] != 2)
                throw std::runtime_error("EP chose the unsafe tied successor");
        }
    }
}

void partial_attractor_game()
{
    // Two excursions return to 2. The chooser at 27 must escape its decoy
    // via the Even sink at 0. Its attraction resolves the exit of the first
    // even tangle at 29, which in turn resolves the exit of the one at 30.
    std::vector<std::vector<int>> edges{{0}, {1}, {4, 3}, {5}, {6}, {7},
        {8}, {9}, {10}, {11}, {12}, {14}, {13}, {2}, {13}};
    std::vector<int> owners{0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0};
    for (int v = 15; v <= 26; ++v) {
        edges[2].push_back(v);
        edges.push_back({1, 2});
        owners.push_back(1);
    }
    edges.insert(edges.end(), {{28, 0}, {1, 2}, {29, 27}, {30, 31}, {32, 29}, {31}});
    owners.insert(owners.end(), {0, 1, 1, 1, 0, 0});
    std::vector<int> priorities;
    for (int v = 0; v <= 28; ++v) priorities.push_back(v);
    priorities.insert(priorities.end(), {0, 0, 31, 30});
    for (int dual : {0, 1}) {
        Game game(33);
        for (int v = 0; v < 33; ++v) {
            game.init_vertex(v, priorities[v] + dual, owners[v] ^ dual);
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        if (solve_game(game, "ep") != solve_game(game, "zlk"))
            throw std::runtime_error("EP lost a strategy at a dominion boundary");
    }
}

void endpoint_games()
{
    // Uniform priorities exercise each saturated endpoint.
    for (int player : {0, 1}) {
        Game game(3);
        for (int v = 0; v < 3; ++v) {
            game.init_vertex(v, 2*v + player, v % 2);
            game.e_start(v);
            game.e_add(v, (v + 1) % 3);
            game.e_add(v, (v + 2) % 3);
            game.e_finish();
        }
        if (solve_game(game, "ep", "clipped") != std::vector<int>(3, player))
            throw std::runtime_error("EP failed to recognize an exact endpoint");
    }
}

void clipped_mixed_game()
{
    // Both constant trial points select 0->1, which allows an odd cycle.
    // The direct equation must recognize the winning alternative 0->2.
    Game game(3);
    for (int v = 0; v < 3; ++v) {
        game.init_vertex(v, v + 2, v != 0);
        game.e_start(v);
        if (v == 0) {
            game.e_add(v, 1);
            game.e_add(v, 2);
        } else game.e_add(v, 0);
        game.e_finish();
    }
    if (solve_game(game, "ep", "clipped") != std::vector<int>(3, 0))
        throw std::runtime_error("EP failed the mixed clipped game");
}

void transient_priorities()
{
    // High, alternating transient priorities do not change the loop winners.
    const int n = 29;
    std::vector<int> priorities{0, 1, 1, 26, 1}, owners{0, 1, 1, 0, 0};
    std::vector<std::vector<int>> edges{{0}, {1}, {2}, {1}, {3, 0}};
    for (int v = 5; v < n; ++v) {
        priorities.push_back(v - 3);
        owners.push_back((v - 3) % 2);
        edges.push_back({2});
    }
    for (int copies : {1, 2}) {
        for (int dual : {0, 1}) {
            Game game(copies * n);
            for (int base = 0; base < copies * n; base += n) {
                for (int v = 0; v < n; ++v) {
                    game.init_vertex(base + v, priorities[v] + dual, owners[v] ^ dual);
                    game.e_start(base + v);
                    for (int to : edges[v]) game.e_add(base + v, base + to);
                    game.e_finish();
                }
            }
            if (solve_game(game, "ep", "clipped") != solve_game(game, "zlk"))
                throw std::runtime_error("EP failed to decode the slow return game");
            if (solve_game(game, "ep", "sign") != solve_game(game, "zlk"))
                throw std::runtime_error("EP failed on transient priorities with cycle weights");
        }
    }
}

// Exercise the exact backend even on games solved by the initial lifts.
void check_measure_search(Game game, std::string_view branch)
{
    const auto expected = solve_game(game, "zlk");
    game.ensure_sorted();
    const ep::Arena arena(game, bitset(game.nodecount()));
    const ep::EnergyMeasure measure(arena);
    std::vector<int> moves(arena.size());
    std::ostringstream trace;
    const auto result = ep::solve_measures(measure, moves, &trace);
    if (result.empty() || trace.str().find(branch) == std::string::npos)
        throw std::runtime_error("EP missed its reduced search");
    if (!arena.certifies(result, moves)) throw std::runtime_error("EP returned an invalid dominion");
    for (int v = 0; v < arena.size(); ++v) {
        if (result[v] >= 0 && result[v] != expected[arena.vertices[v]])
            throw std::runtime_error("EP reduced search gave a wrong strategy");
    }
}

void reduced_clipped_game()
{
    // Only two vertices of each owner branch. The other rows may contain cycles.
    const std::vector<int> priorities{2, 3, 4, 5, 7, 8, 10};
    const std::vector<int> owners{1, 0, 0, 1, 0, 1, 0};
    const std::vector<std::vector<int>> edges{{3}, {2}, {0, 2, 6}, {1, 3, 6}, {6}, {4, 3}, {6, 3, 5}};
    for (int dual : {0, 1}) {
        Game game(7);
        for (int v = 0; v < 7; ++v) {
            game.init_vertex(v, priorities[v] + dual, owners[v] ^ dual);
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        check_measure_search(game, "2 boundary coordinates using one-player substitution");
    }
}

void price_profile_game()
{
    // Four vertices of either owner branch, and all eight have self-loops.
    // The first production price probe can publish a certified region.
    for (int dual : {0, 1}) {
        Game game(8);
        for (int v = 0; v < 8; ++v) {
            game.init_vertex(v, v + dual, (v + dual) % 2);
            game.e_start(v);
            game.e_add(v, v);
            game.e_add(v, (v + 1) % 8);
            game.e_finish();
        }
        check_measure_search(game, "4 boundary coordinates");
    }
}

void energy_strategies()
{
    const std::vector<int> priorities{3, 1, 3, 4, 1, 1}, owners{0, 1, 0, 1, 0, 0};
    const std::vector<std::vector<int>> edges{{0, 3}, {5, 0}, {5, 1}, {0, 2}, {4, 2}, {5, 3}};
    for (int dual : {0, 1}) {
        Game game(6);
        for (int v = 0; v < 6; ++v) {
            game.init_vertex(v, priorities[v] + dual, owners[v] ^ dual);
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        if (solve_game(game, "ep", "sign") != solve_game(game, "zlk"))
            throw std::runtime_error("EP lost an energy strategy on a mixed game");
    }
}

void bounded_energy_sweeps()
{
    // The length-10 counter first admits its three-vertex Even dominion
    // at sweep 90. The large-width front end must hand off after n=16.
    const std::vector<int> priorities{15,14,12,11,10,9,8,7,6,5,4,3,2,1,1,0};
    const std::vector<int> owners{0,1,1,0,1,0,1,0,1,0,1,0,0,0,0,1};
    const std::vector<std::vector<int>> edges{{12}, {0,1,2,13}, {2,12,1,3}, {3,12,2,4},
        {4,12,3,5}, {5,12,4,6}, {6,12,5,7}, {7,12,6,8}, {8,12,7,9}, {9,12,8,10},
        {10,12,9,11}, {11,12,10}, {1,15}, {13,0}, {14,12}, {15,14}};
    Game game(16);
    for (int v = 0; v < 16; ++v) {
        game.init_vertex(v, priorities[v], owners[v]);
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    game.ensure_sorted();
    const ep::Arena arena(game, bitset(16));
    const ep::EnergyMeasure measure(arena);
    std::vector<int> moves(16);
    if (measure.upper - measure.lower != 564 || !ep::lift_measures(measure, moves, nullptr).empty())
        throw std::runtime_error("EP failed to bound local energy iteration");
}

void exact_price_labels()
{
    // Check all maximizing slots against an independently iterated fixed
    // point, including exact queries, rational prices, ties, and both walls.
    const std::vector<std::vector<int>> edges{{0, 1, 2}, {0, 1, 2}, {0, 1}};
    for (int colors = 0; colors < 8; ++colors) {
        for (int owners = 0; owners < 8; ++owners) {
            Game game(3);
            for (int v = 0; v < 3; ++v) {
                game.init_vertex(v, 2*v + ((colors >> v) & 1), (owners >> v) & 1);
                game.e_start(v);
                for (int to : edges[v]) game.e_add(v, to);
                game.e_finish();
            }
            const ep::Arena arena(game, bitset(3));
            const ep::EnergyMeasure measure(arena);
            Prices scores(3, measure.lower);
            for (Integer pass = 0; pass <= 3*(measure.upper - measure.lower); ++pass) {
                Prices before = scores;
                for (int v = 0; v < 3; ++v) scores[v] = measure.lift(v, scores);
                if (scores == before) break;
            }
            for (int v = 0; v < 3; ++v)
                if (measure.lift(v, scores) != scores[v]) throw std::runtime_error("unsettled credit");
            std::vector<int> moves(3);
            const auto winners = measure.solution(scores, moves);
            const auto expected = solve_game(game, "zlk");
            for (int v = 0; v < 3; ++v)
                if (expected[arena.vertices[v]] != winners[v]) throw std::runtime_error("wrong energy cut");
            for (const Integer& center : {measure.lower, Integer(0), measure.upper}) {
                for (int digits = 0; digits < 125; ++digits) {
                    Prices trial(3);
                    Integer distance = 0;
                    int code = digits;
                    for (int v = 0; v < 3; ++v) {
                        trial[v] = 2 * (center == 0 ? scores[v] : center) + code % 5 - 2;
                        code /= 5;
                        distance = std::max(distance, Integer(abs(2*scores[v] - trial[v])));
                    }
                    const auto labels = measure.price_directions(trial, 2, moves);
                    if (labels.empty() != (distance == 0)) throw std::runtime_error("wrong exact-price exit");
                    if (labels.empty()) continue;
                    for (int v = 0; v < 3; ++v) {
                        const Integer error = 2*scores[v] - trial[v];
                        if (abs(error) == distance && labels[v] != (error > 0 ? 1 : -1))
                            throw std::runtime_error("wrong direction on a maximal price error");
                    }
                }
            }
        }
    }
}

void price_recursion()
{
    // Exercise the same market, wealth and median code as the large schedule.
    // The theorem uses larger probe counts. Here we can check the whole tree.
    std::size_t visits = 0;
    bool decoded = false;
    detail::price_tree(1, 1, {3, 3, 0, false},
        [&](const Prices& x, const Integer& unit, std::vector<int>& labels) {
            if (visits++ == 0 && x[0] != 0) throw std::runtime_error("initial market is asymmetric");
            labels = {x[0] < unit ? 1 : -1};
            return false;
        },
        [&](const Prices& x, const Integer& unit) {
            decoded = true;
            if (x.size() != 1 || unit != 8) throw std::runtime_error("incorrect price units");
        });
    if (!decoded || visits != 39) throw std::runtime_error("incomplete price recursion");

    // With two signed slots one ancestor suffices for the mass proof:
    // (4*5/8)/2 > 1. The three first-level medians replace their subtrees.
    visits = 0;
    detail::price_tree(1, 1, {3, 3, 0},
        [&](const Prices& x, const Integer& unit, std::vector<int>& labels) {
            ++visits;
            labels = {x[0] < unit ? 1 : -1};
            return false;
        },
        [](const Prices&, const Integer&) {});
    if (visits != 3) throw std::runtime_error("certified price subtrees were traversed");

    // At depth one each query is a proposal and the output is their median.
    // Persistent feedback in either direction must actually move the market.
    for (int target : {-1, 1}) {
        Prices proposals;
        detail::price_tree(1, 1, {1, 19, 0},
            [&](const Prices& x, const Integer& unit, std::vector<int>& labels) {
                proposals.push_back(x[0]);
                labels = {x[0] < target * unit ? 1 : -1};
                return false;
            },
            [&](const Prices& x, const Integer& unit) {
                std::sort(proposals.begin(), proposals.end());
                if (unit != 2 || proposals.size() != 19 || x[0] != proposals[9] || target*x[0] <= 0)
                    throw std::runtime_error("price feedback or median is incorrect");
            });
    }
}

} // namespace

int test_ep()
{
    if (test_ep_experts() || test_ep_market()) return 1;
    price_recursion();
    exact_price_labels();
    bounded_energy_sweeps();
    signed_credits();
    credit_handoff();
    exhaustive_games();
    endpoint_games();
    clipped_mixed_game();
    transient_priorities();
    reduced_clipped_game();
    price_profile_game();
    energy_strategies();
    supported_move_choice();
    open_tangle();
    partial_attractor_game();
    partial_game();
    return 0;
}
