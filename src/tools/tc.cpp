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
#include <memory>

#include "oink/game.hpp"

#define DOUBLEDISTRACTION 0

using namespace std;
using namespace pg;

template<typename ... Args>
std::string string_format(const std::string& format, Args ... args)
{
    size_t size = snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
    std::unique_ptr<char[]> buf(new char[ size ]);
    snprintf(buf.get(), size, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

static int *in0, *in1; // indices of each bit's input/distraction
static int *hi0, *hi1; // indices of each bit's high vertex
static int n, m;

/**
 * Make bit <i> at base index <c> with given input/high priorities
 */
void
makeBit(Game &game, const int c, int hipr, int inpr, int lopr, int i)
{
    const int pl = hipr&1;
    int* _in = pl == 0 ? in1 : in0; // select _in of opponent

    /**
     * Create "core" of the bit.
     */

    const char *plch = pl ? "Odd-" : "Even-";

    game.init_vertex(c,   hipr,   pl,   string_format("%s%d-H", plch, i));  // high
    game.init_vertex(c+1, lopr,   1-pl, string_format("%s%d-T", plch, i));  // low (tangle)
    game.init_vertex(c+2, inpr,   1-pl, string_format("%s%d-L", plch, i));  // gate (distraction/input)
#if DOUBLEDISTRACTION
    game.init_vertex(c+3, inpr,   1-pl, string_format("%s%d-LL", plch, i));  // gate (distraction/input)
#endif

    game.vec_add_edge(c+1, c);   // from tangle to top
#if DOUBLEDISTRACTION
    game.vec_add_edge(c+2, c+3); // from input/distraction to tangle
    game.vec_add_edge(c+3, c+1); // from input/distraction to tangle
    game.vec_add_edge(c+1, c+4); // from tangle to first connector
#else
    game.vec_add_edge(c+2, c+1); // from input/distraction to tangle
    game.vec_add_edge(c+1, c+3); // from tangle to first connector
#endif

    /**
     * Create connectors (to higher bits)
     */

    for (int j=0; j<i; j++) {
        const int d = c + 3 + DOUBLEDISTRACTION + 3*j;

        game.init_vertex(d  , lopr-1, pl,   string_format("%s%d-S-%d", plch, i, j)); // selector
        game.init_vertex(d+1, lopr-1, 1-pl, string_format("%s%d-A-%d", plch, i, j)); // exit one (even)
        game.init_vertex(d+2, lopr-1, 1-pl, string_format("%s%d-B-%d", plch, i, j)); // exit two (odd)

        game.vec_add_edge(d, d+1);  // s to one
        game.vec_add_edge(d, d+2);  // s to two

        game.vec_add_edge(d+1, d+3); // one to next selector
        game.vec_add_edge(d+2, d+3); // two to next selector

        int *ina = pl ? in1 : in0;
        int *inb = pl ? in0 : in1;
        game.vec_add_edge(d+1, ina[j]); // one to even input of bit <j>
        game.vec_add_edge(d+2, inb[j]); // two to odd input of bit <j>

        // for (int j = i+pl; j<n; j++) game.vec_add_edge(d, _in[j]);
    }

    /**
     * Connect from <high> either to the next <input>
     */

#if 0
    game.vec_add_edge(c, (pl == 0 ? in0 : in1)[(i+n-1)%n]-1);
#else
    game.vec_add_edge(c, (pl == 0 ? in0 : in1)[(i+n-1)%n]);
#endif

    // even (lower or same bit)
    // odd (only to lower bits)
    int x = c+3+DOUBLEDISTRACTION+3*i;
    game.init_vertex(x, lopr-1, pl, string_format("%s%d-Z", plch, i));  // Z (distracted)
    game.vec_add_edge(x, c+1); // from Z to tangle
    for (int j = i+1-pl; j<n; j++) game.vec_add_edge(x, _in[j]);
}

int
main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Syntax: " << argv[0] << " N" << std::endl;
        return -1;
    }

    /**
     * Get n and allocate data structures
     */

    n = std::stoi(argv[1]);
    m = 3*n; // deceleration lane length
    in0 = new int[n];
    in1 = new int[n];
    hi0 = new int[n];
    hi1 = new int[n];

    /**
     * Compute number of vertices
     * - every bit has 3 vertices (low, high, input) plus 2 (to lower other, to lower own)
     * - every Nth bit has 3*N additional vertices connected to earlier counters
     */

    int size = 0;

    for (int i=0; i<n; i++) {
        hi0[i] = size;
        in0[i] = size + 2; // set even gate / distraction
        size += 4+3*i+DOUBLEDISTRACTION;
        hi1[i] = size;
        in1[i] = size + 2; // set odd gate / distraction
        size += 4+3*i+DOUBLEDISTRACTION;
    }

    /**
     * Create game...
     */

    Game game(size);
    game.vec_init();

    /**
     * Create the two counters
     * - the even bit is higher than the odd bit (hipr)
     * - lopr is only necessary to beat inflated oink DP (pgsolver DP not vulnerable) 
     * - strict separation between two bits priorities similarly
     */

    int lo = 0;
#if 0
    int toppo = (n+3)*n+6*n; // N*N + 9N actually
    for (int i=0; i<n; i++) {
        makeBit(game, hi0[i], toppo,   toppo-1, lo += 2*(n+1-i), i);
        makeBit(game, hi1[i], toppo-3, toppo-4, lo + 1, i);
        toppo -= 6;
    }
#else
    int toppo = (n+3)*n+6*n; // N*N + 9N actually
    for (int i=0; i<n; i++) {
        // makeBit(game, hi0[i], toppo,   toppo-2*n-1, lo += 2*(n+1-i), i);
        makeBit(game, hi0[i], toppo-2, toppo-2*n-3, lo + 2, i);
        makeBit(game, hi1[i], toppo-1, toppo-2*n-2, lo + 1, i);
        toppo -= 2;
    }
#endif

    game.vec_finish();
    game.sort();
    game.renumber();
    game.write_pgsolver(std::cout);

    delete[] in0;
    delete[] in1;
    delete[] hi0;
    delete[] hi1;
}
