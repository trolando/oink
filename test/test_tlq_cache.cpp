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
 * Unit tests for TLQ: the completed-tangle markers, and chains of dominions.
 */

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "oink/oink.hpp"
#include "solvers/tlq.hpp"
#include "verifier.hpp"

using namespace pg;

/**
 * Attracts the tangle {0,1} with escapes 2 and 3 to a target that holds both
 * escapes, with an attraction epoch that wraps on the first attraction.
 */
class Probe : public TLQSolver
{
public:
    Probe(Oink& oink, Game& game) : TLQSolver(oink, game)
    {
        const int n = game.nodecount();
        G.resize(n);
        G.flip();
        Processed.resize(n);
        Q.resize(n);
        str = new int[n];
        std::fill(str, str+n, -1);
        escs = new unsigned[n];
        std::fill(escs, escs+n, 99u);
        esc_epoch = new unsigned[n];
        std::fill(esc_epoch, esc_epoch+n, 1u);
        tangleto = {2, 3};
        addBacklinks(0);
        tangleto.clear();
        tvpos.push_back(tvdata.size());
        tvdata.insert(tvdata.end(), {0, 1, 1, -1, -1});
        toutpos.push_back(toutdata.size());
        toutdata.insert(toutdata.end(), {2, 3, -1});
        tpr.push_back(0);
        tdone.push_back(1);
        attract_epoch = std::numeric_limits<unsigned>::max();
    }

    ~Probe()
    {
        delete[] str;
        delete[] escs;
        delete[] esc_epoch;
    }

    void check()
    {
        for (unsigned expected : {1u, 2u}) {
            bitset domain(G);
            bitset target(4);
            target[2] = true;
            target[3] = true;
            Q.push(2);
            Q.push(3);
            attract(0, target, domain);
            if (target.count() != 4 || str[0] != 1 || attract_epoch != expected || tdone[0] != expected) {
                throw std::runtime_error("completed-tangle cache failed across epoch transition");
            }
        }
    }
};

/**
 * Odd heads lead down a chain of even self-loops. Removing one even dominion
 * invalidates the previous odd attractor, exposing the next loop. Three links
 * already exceed the revisits in the initial precision tree.
 */
static void
check_dominion_restarts()
{
    for (int links : {1, 2, 3, 4, 7, 8, 15, 16, 31, 32}) {
        for (int winner : {0, 1}) {
            for (int head_owner : {0, 1}) {
                const int n = links*2+1;
                Game game(n);
                for (int v = 0; v < n; ++v) {
                    game.init_vertex(v, (v <= links ? 0 : 1)+winner, v <= links ? 1-winner : head_owner);
                    game.e_start(v);
                    if (v <= links) {
                        game.e_add(v, v);
                        if (v > 0) game.e_add(v, links+v);
                    } else {
                        game.e_add(v, v-links-1);
                    }
                    game.e_finish();
                }
                game.ensure_sorted();
                std::ostringstream log;
                Oink controller(game, log, true);
                controller.setRemoveLoops(false);
                controller.setRemoveWCWC(false);
                controller.setSolveSingle(false);
                controller.setWorkers(-1);
                controller.setSolver("tlq");
                controller.run();
                Verifier(game, controller.solution(), log).verify(true, true, true);
                for (int v = 0; v < n; ++v) {
                    if (controller.solution().winner(v) != winner) {
                        throw std::runtime_error("TLQ retained a stale region after removing dominions");
                    }
                }
            }
        }
    }
}

int
test_tlq_cache()
{
    Game game(4);
    for (int v = 0; v < 4; v++) game.init_vertex(v, v < 2 ? 0 : 2, v == 1 ? 1 : 0);
    const std::vector<std::vector<int>> edges = {{1}, {0, 2, 3}, {2}, {3}};
    for (int v = 0; v < 4; v++) {
        game.e_start(v);
        for (int to : edges[v]) game.e_add(v, to);
        game.e_finish();
    }
    game.ensure_sorted();
    std::ostringstream log;
    Oink controller(game, log, true);
    Probe(controller, game).check();
    std::cout << "Completed-tangle markers survive epoch wrap and repeated callbacks." << std::endl;
    check_dominion_restarts();
    std::cout << "Dominion removal revisits stale regions for both players." << std::endl;
    return 0;
}
