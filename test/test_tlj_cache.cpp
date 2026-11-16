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
 * Unit tests for the blocking escapes that TLJ remembers when attracting tangles.
 */

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "oink/oink.hpp"
#include "solvers/tlj.hpp"

using namespace pg;

static void
require(bool condition)
{
    if (!condition) throw std::runtime_error("TLJ blocker or queued-escape regression");
}

/**
 * Exercise the attraction primitive independently of the full search. The tangle
 * is the 0<->1 cycle, with escapes 2 and 3. Region 4 grows towards its head 4.
 */
class CacheProbe : public TLJSolver
{
    std::vector<int> assignments, moves;
    int p, player;

    void rebuild()
    {
        std::fill(rcount.begin(), rcount.end(), 0);
        nonempty[0].reset();
        nonempty[1].reset();
        anyregion.reset();
        Z.reset();
        for (int v = 0; v < 7; ++v) {
            if (region[v] < 0) continue;
            ++rcount[region[v]];
            nonempty[region[v]&1].set(region[v]);
            anyregion.set(region[v]);
            if (region[v] == p) Z.set(v);
        }
    }

public:
    CacheProbe(Oink& controller, Game& game, int invert)
        : TLJSolver(controller, game), assignments(7), moves(7, -1), p(4+invert), player(invert)
    {
        region = assignments.data();
        str = moves.data();
        const int initial[] = {0, 0, 3, 3, 4, 5, 6};
        for (int v = 0; v < 7; ++v) region[v] = initial[v]+invert;
        rprio = assignments;
        G.resize(7);
        G.flip();
        Z.resize(7);
        C.resize(7);
        tangled.resize(7);
        rcount.resize(8);
        nonempty[0].resize(8);
        nonempty[1].resize(8);
        anyregion.resize(8);
        zlist.resize(7);
        flipped.resize(7);
        pull_lo = p;
        tpr = {invert};
        tcnt = {-1};
        tcpos = {0};
        tblock = {-1};
        tvdata = {0, 1, 1, 0, -1};
        tvpos = {0};
        toutdata = {2, 3, -1};
        toutpos = {0};
        zlist.push(4);
        rebuild();
    }

    bool attempt(unsigned position, int threshold = -1)
    {
        return attractTangleInto(0, position, threshold < 0 ? p : threshold, player);
    }

    void next_visit()
    {
        tcnt[0] = -1;
        ttouched.clear();
        zlist.clear();
        zlist.push(4);
    }

    void place(int vertex, int r, bool queued = false)
    {
        region[vertex] = r;
        rebuild();
        if (queued) zlist.push(vertex);
    }

    void join(int vertex) { place(vertex, p, true); }
    void remove(int vertex)
    {
        G.reset(vertex);
        place(vertex, -1);
    }
    int count() const { return tcnt[0]; }
    unsigned cutoff() const { return tcpos[0]; }
    bool touched_counter() const { return !ttouched.empty(); }
    bool deleted() const { return tpr[0] == -1; }

    void check_attracted() const
    {
        require(region[0] == p && region[1] == p && str[0] == 1 && str[1] == 0 && tangled[1]);
    }
};

int
test_tlj_cache()
{
    for (int invert : {0, 1}) {
        Game game(7);
        const int priorities[] = {0, 0, 1, 3, 4, 5, 6};
        const int owners[] = {0, 1, 0, 0, 0, 1, 0};
        for (int v = 6; v >= 0; --v) game.init_vertex(v, priorities[v]+invert, owners[v]^invert);
        const std::vector<std::vector<int>> edges = {{1}, {0, 2, 3}, {4}, {4}, {4}, {5}, {6}};
        for (int v = 0; v < 7; ++v) {
            game.e_start(v);
            for (int to : edges[v]) game.e_add(v, to);
            game.e_finish();
        }
        game.ensure_sorted();
        std::ostringstream log;
        Oink controller(game, log, true);
        {
            CacheProbe s(controller, game, invert);
            require(!s.attempt(1) && s.count() == 2);
            s.next_visit();
            require(!s.attempt(1) && s.count() == -1 && !s.touched_counter());
            s.join(2);
            require(!s.attempt(2) && s.count() == 1 && s.cutoff() == 2);
            // Escape 2 was assigned before initialization; its later callback
            // must not consume the count belonging to escape 3.
            require(!s.attempt(1) && s.count() == 1);
            s.join(3);
            require(s.attempt(2) && s.count() == -2);
            s.check_attracted();
        }
        {
            CacheProbe s(controller, game, invert);
            require(!s.attempt(1));
            s.next_visit();
            s.join(2);
            s.join(3);
            // Both escapes are assigned, but escape 3 is still queued. Early
            // initialization may attract now, and later callbacks are harmless.
            require(s.attempt(1) && s.count() == -2);
            s.check_attracted();
            require(!s.attempt(2) && s.count() == -2);
        }
        {
            CacheProbe s(controller, game, invert);
            require(!s.attempt(1));
            s.next_visit();
            s.place(2, 5+invert);
            require(!s.attempt(1) && s.count() == -1);
            s.place(2, -1);
            require(!s.attempt(1) && s.count() == -1);
            s.place(2, 4+invert);
            require(!s.attempt(1, 6+invert) && s.count() == -1);
        }
        {
            CacheProbe s(controller, game, invert);
            require(!s.attempt(1));
            s.next_visit();
            s.remove(2);
            s.join(3);
            require(s.attempt(1) && s.count() == -2);
            s.check_attracted();
        }
        {
            CacheProbe s(controller, game, invert);
            require(!s.attempt(1));
            s.next_visit();
            s.join(2);
            s.join(3);
            s.remove(0);
            require(!s.attempt(1) && s.deleted());
        }
    }
    std::cout << "TLJ blockers respect live escapes, region changes, removed members and queued callbacks." << std::endl;
    return 0;
}
