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

#ifndef DFTL_HPP
#define DFTL_HPP

#include "solver.hpp"

namespace pg {

class DFTLSolver : public Solver
{
public:
    DFTLSolver(Oink *oink, Game *game);
    virtual ~DFTLSolver();

    virtual void run();

protected:
    bool prepartition = false; // occasionally may result in MORE iterations
    bool interleaved = true;

    int iterations = 0;
    int dominions = 0;
    int tangles = 0;
    int steps = 0;

    std::vector<int*> tout; // for each tangle
    std::vector<int> *tin; // for each normal vertex
    std::vector<int*> tv; // the tangle (vertex-strategy pairs)
    std::vector<int> tpr; // priority of a tangle

    uintqueue Q; // main queue when attracting vertices
    std::vector<int> dom_vector; // stores the solved vertices (in dominions)

    int *str; // stores currently assigned strategy of each vertex

    uintqueue pea_vS; // vS
    uintqueue pea_iS; // iS
    uintqueue pea_S;  // S
    unsigned int* pea_vidx;    // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    uintqueue tangleto; // stores the vertices the tangle can escape to
    bitset escapes; // which escapes we considered

    bitset Z; // current region (in sptl)
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    bitset CurG; // in search

    bitset V;
    bitset W;

    bitset Parity; // vertices of parity 0/1
    bitset Even; // for partition
    bitset Odd; // for partition

    inline void attractVertices(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);
    bool attractTangle(const int t, const int pl, bitset &R, bitset &Z, bitset &G, const int max_prio);
    inline void attractTangles(const int pl, int v, bitset &R, bitset &Z, bitset &G, const int max_prio);

    bool search(const int player);
    void partition(bitset &R, int top, bitset &Even, bitset &Odd);
    bool sptl(bitset &V, bitset &R, const int player, int &lowest_top);
    bool extractTangles(int i, bitset &R, int *str);
};

class PDFTLSolver : public DFTLSolver
{
public:
    PDFTLSolver(Oink *oink, Game *game) : DFTLSolver(oink, game) { prepartition = true; }
    virtual ~PDFTLSolver() { }
};



}

#endif
