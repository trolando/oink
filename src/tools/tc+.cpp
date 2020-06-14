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

#include "game.hpp"

#define DOUBLEDISTRACTION 1

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
    int *inmy = pl ? in1 : in0;
    int *inop = pl ? in0 : in1;

    /**
     * Create "core" of the bit.
     */

    string bitid_str = string_format(pl?"Odd-%d":"Even-%d", i);
    const char *bitid = bitid_str.c_str();

    game.init_vertex(c,   hipr, pl,   string_format("%s-H", bitid));  // high
    game.init_vertex(c+1, inpr, 1-pl, string_format("%s-I", bitid));  // gate (distraction/input)
    game.init_vertex(c+2, lopr, 1-pl, string_format("%s-T", bitid));  // low (tangle)
#if DOUBLEDISTRACTION
    game.init_vertex(c+3, inpr, 1-pl, string_format("%s-II", bitid));  // gate (distraction/input)
    game.vec_add_edge(c+1, c+3); // connect I -> T
    game.vec_add_edge(c+3, c+2); // connect I -> T
    game.vec_add_edge(c+2, c+4); // connect T -> first S
#else
    game.vec_add_edge(c+1, c+2); // connect I -> T
    game.vec_add_edge(c+2, c+2); // connect T -> first S
#endif
    game.vec_add_edge(c+2, c);   // connect T -> H
    game.vec_add_edge(c, inmy[i == 0 ? n-1 : i-1]); // connect H -> I{(j-1) mod n}

    /**
     * Create connectors (to higher bits)
     */

    for (int j=0; j<i; j++) {
        const int d = c + 3 + DOUBLEDISTRACTION + 3*j;

        game.init_vertex(d  , lopr-1, pl,   string_format("%s-S-%d", bitid, j)); // selector
        game.init_vertex(d+1, lopr-1, 1-pl, string_format("%s-A-%d", bitid, j)); // exit one (even)
        game.init_vertex(d+2, lopr-1, 1-pl, string_format("%s-B-%d", bitid, j)); // exit two (odd)

        const int next_s = (j+1 == i) ? c+2 : d+3;

        game.vec_add_edge(d, d+1);       // connect Sj -> Aj
        game.vec_add_edge(d, d+2);       // connect Sj -> Bj
        game.vec_add_edge(d+1, next_s);  // connect Aj -> S{j+1}
        game.vec_add_edge(d+2, next_s);  // connect Bj -> S{j+1}
        game.vec_add_edge(d+1, inmy[j]); // connect Aj -> <my> Ij
        game.vec_add_edge(d+2, inop[j]); // connect Bj -> <their> Ij
    }

#if 0
#if 1
    game.init_vertex(c+3+3*i, lopr+2, pl, string_format("%s-D-%d", bitid, i));
    game.vec_add_edge(c+3+3*i, c+2);
#else
    if (pl == 0) {
        game.init_vertex(c+3+3*i, lopr+2, pl, string_format("d(%d,%d,%d)",pl,i,i));
        game.vec_add_edge(c+3+3*i, c+2);
        // we want Odd to play first for DP,
        // but this breaks symmetry for RTL
        game.vec_add_edge(c+3+3*i, inop[i]);
    } else {
        game.init_vertex(c+3+3*i, lopr+2, pl, string_format("d(%d,%d,%d)",pl,i,i));
        game.vec_add_edge(c+3+3*i, c+2);
    }
#endif
#endif

    /**
     * Connect to distractions, i.e., lower bits.
     */

    int x = i+1;
    for (int j = x; j<n; j++) {
        const int dv = c + 3 + DOUBLEDISTRACTION + 3*i + (j-x);
        const int pr = lopr + 2*(j-i+1);
        game.init_vertex(dv, pr, pl, string_format("%s-D-%d", bitid, j)); // todo fix priority
        game.vec_add_edge(dv, c+2); // dv -> T
        game.vec_add_edge(c+2, dv); // T -> dv
        game.vec_add_edge(dv, inop[j]); // dv -> ...-I
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
        in0[i] = size + 1; // set even gate // distraction
        size += DOUBLEDISTRACTION + 3 + 3*i + n-i-1;
        hi1[i] = size;
        in1[i] = size + 1; // set odd gate // distraction
        size += DOUBLEDISTRACTION + 3 + 3*i + n-i-1;
    }

    /**
     * Create game...
     */

    Game game(size);
    game.vec_init();

    /**
     * Create the two counters
     * - the even bit is higher than the odd bit (hipr)
     * - strict separation between two bits priorities similarly
     */

    int toppo = (n+3)*n+6*n; // N*N + 9N actually
    for (int i=0; i<n; i++) {
        int lo = 2;// + (2*n-i+1)*i + 2*i;
        makeBit(game, hi0[i], toppo,   toppo-2*n-1, lo, i);
        makeBit(game, hi1[i], toppo-1, toppo-2*n-2, lo - 1, i);
        toppo -= 2;
    }

    game.vec_finish();
    int mapping[game.nodecount()];
    game.sort(mapping);
    game.renumber();
    game.permute(mapping);
    game.write_pgsolver(std::cout);

    delete[] in0;
    delete[] in1;
}
