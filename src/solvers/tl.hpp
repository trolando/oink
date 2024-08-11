/*
 * Copyright 2020 Tom van Dijk
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

#ifndef TL_HPP
#define TL_HPP

#include "oink/solver.hpp"

namespace pg {

class TLSolver : public Solver
{
public:
    TLSolver(Oink& oink, Game& game);
    virtual ~TLSolver();

    virtual void run();

protected:
    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    int steps = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle
    std::vector<unsigned int> tescs; // number of escapes of a tangle

    uintqueue Q; // main queue when attracting vertices

    int *str; // stores currently assigned strategy of each vertex
    unsigned int *escs; // stores remaining escapes of each vertex

    uintqueue pea_state; // v,i,...
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset R; // remaining region (in tl)
    bitset Z; // current region (in tl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset V; // heads 1
    bitset W; // heads 2

    bitset Even; // accumulate even regions
    bitset Odd; // accumulate odd regions

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, int maxpr);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, int maxpr);

    bool tl(void);
    bool extractTangles(int i, bitset &R, int *str);
};

}

#endif
