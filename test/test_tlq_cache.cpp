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
 * Unit tests for the completed-tangle markers of TLQ.
 */

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "oink/oink.hpp"
#include "solvers/tlq.hpp"

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
        tin = new std::vector<int>[n];
        has_tin.resize(n);
        for (int v : {2, 3}) {
            tin[v].push_back(0);
            has_tin[v] = true;
        }
        tv.push_back(new int[5]{0, 1, 1, -1, -1});
        tout.push_back(new int[3]{2, 3, -1});
        tpr.push_back(0);
        tdone.push_back(1);
        attract_epoch = std::numeric_limits<unsigned>::max();
    }

    ~Probe()
    {
        delete[] str;
        delete[] escs;
        delete[] esc_epoch;
        delete[] tin;
        for (auto p : tv) delete[] p;
        for (auto p : tout) delete[] p;
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
    return 0;
}
