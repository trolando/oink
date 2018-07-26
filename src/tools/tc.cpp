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

#include "game.hpp"

#define COT 0

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
    const int sela = c+3;

    /**
     * Create "core" of the bit.
     */

    game.initNode(c,   hipr,   pl,   string_format("h%d(%d)", pl, i));  // high
    game.initNode(c+1, lopr,   1-pl, string_format("l%d(%d)", pl, i));  // low (tangle)
    game.initNode(c+2, inpr,   1-pl, string_format("i%d(%d)", pl, i));  // gate (distraction/input)
#if COT
    game.initNode(c+3, 0,      1-pl, string_format("z%d(%d)", pl, i));  // sela (distracted)
#else
    game.initNode(c+3, 0,      pl,   string_format("z%d(%d)", pl, i));  // sela (distracted)
#endif

    game.addEdge(c+1, c);   // from low to high
    game.addEdge(c+2, c+1); // from gate to low
    game.addEdge(c+3, c+1); // from sela to low

    /**
     * Connect to first connector
     */

    if (i == 0) game.addEdge(c+1, c+3); // from <low> to <sela> (no connectors)
    else game.addEdge(c+1, c+4);        // from <low> to <sel0> (has connectors)

    /**
     * Create connectors (to higher bits)
     */

    for (int j=0; j<i; j++) {
        const int d = c + 4 + 3*j;

        game.initNode(d  , 0, pl,   string_format("s%d(%d,%d)", pl, i, j)); // selector
        game.initNode(d+1, 0, 1-pl, string_format("a%d(%d,%d)", pl, i, j)); // exit one (even)
        game.initNode(d+2, 0, 1-pl, string_format("b%d(%d,%d)", pl, i, j)); // exit two (odd)

        game.addEdge(d, d+1);  // s to one
        game.addEdge(d, d+2);  // s to two

        if (j+1 != i) {
            // not the last, connect to next
            game.addEdge(d+1, d+3); // one to next selector
            game.addEdge(d+2, d+3); // two to next selector
        } else {
            // last, connect to start
            game.addEdge(d+1, sela); // one to <sela>
            game.addEdge(d+2, sela); // two to <sela>
        }

        int *ina = pl ? in1 : in0;
        int *inb = pl ? in0 : in1;
        game.addEdge(d+1, ina[j]); // one to even input of bit <j>
        game.addEdge(d+2, inb[j]); // two to odd input of bit <j>
    }

    /**
     * Connect from <high> either to the next <input>
     */

    game.addEdge(c, (pl == 0 ? in0 : in1)[(i+n-1)%n]);

    /**
     * Connect from <sela> to lower bits.
     */
    int* _in = pl == 0 ? in1 : in0;
    for (int j = i+pl; j<n; j++) {
        // even (lower or same bit)
        // odd (only to lower bits)
#if COT
        const int e = c + 4 + 3*i + (j-i-pl);
        // const int pr = 2+pl + 2*(n-j+i+pl);
        const int pr = lopr - 2*(n-j+i+pl);
        game.initNode(e, pr, pl, string_format("d(%d,%d,%d)",pl,i,j)); // todo fix priority
        game.addEdge(e, sela);
        game.addEdge(sela, e);
        game.addEdge(e, _in[j]);
#else
        game.addEdge(sela, _in[j]);
#endif
    }
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
        in0[i] = size + 2; // set even gate
        size += (4+3*i) + (COT ? (n-i) : 0);
        hi1[i] = size;
        in1[i] = size + 2; // set odd gate
        size += (4+3*i) + (COT ? (n-i-1) : 0);
    }

    /**
     * Create game...
     */

    Game game(size);

    /**
     * Create the two counters
     * - the even bit is higher than the odd bit (hipr)
     * - lopr is only necessary to beat inflated oink DP (pgsolver DP not vulnerable) 
     * - strict separation between two bits priorities similarly
     */

    int lo = 0;
    int toppo = (n+3)*n+6*n; // N*N + 9N actually
#if 0
    for (int i=0; i<n; i++) {
        makeBit(game, hi0[i], toppo,   toppo-1, lo += 2*(n+1-i), i);
        makeBit(game, hi1[i], toppo-3, toppo-4, lo + 1, i);
        toppo -= 6;
    }
#else
    for (int i=0; i<n; i++) {
        makeBit(game, hi0[i], toppo,   toppo-2*n-1, lo += 2*(n+1-i), i);
        makeBit(game, hi1[i], toppo-1, toppo-2*n-2, lo + 1, i);
        toppo -= 2;
    }
#endif

    game.reindex();
    game.renumber();
    game.write_pgsolver(std::cout);

    delete[] in0;
    delete[] in1;
}
