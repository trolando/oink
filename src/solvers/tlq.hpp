/*
 * Copyright 2024-2025 Tom van Dijk, University of Twente
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

#ifndef TLQ_HPP
#define TLQ_HPP

#include <memory>
#include <vector>

#include "oink/solver.hpp"

namespace pg {

/**
 * This is an implementation of the universal recursive algorithm run on a quasi-polynomial tree
 * The original version was based on https://www.mimuw.edu.pl/~parys/publications/2018-parity-algorithm.pdf
 * I modified the algorithm to produce a strategy as well.
 * I changed it with the optimizations from Lehtinen et al (arXiv)
 * Furthermore I added some optimizations, i.e., shortcuts in the tree
 */

class TLQSolver : public Solver
{
public:
    TLQSolver(Oink& oink, Game& game);
    virtual ~TLQSolver();

    virtual void run();

protected:
    unsigned long long iterations = 0;
    unsigned int dominions = 0;
    unsigned int tangles = 0;

    std::vector<int> tvdata; // the vertex-strategy pairs of the tangles (each ending with -1)
    std::vector<int> toutdata; // the escapes of the tangles (each ending with -1)
    std::vector<size_t> tvpos; // for each tangle, its vertex-strategy pairs in <tvdata>
    std::vector<size_t> toutpos; // for each tangle, its escapes in <toutdata>
    // The tangles that escape to a vertex, in discovery order in <tin_pool>.
    // Only pages with an escape destination have descriptors, and only vertices
    // marked in <has_tin> have an initialized descriptor.
    struct Backlinks { size_t offset; unsigned int size, capacity; };
    static constexpr unsigned int TIN_PAGE_SIZE = 256;
    std::vector<int> tin_pool;
    std::vector<std::unique_ptr<Backlinks[]>> tin_pages;
    bitset has_tin;
    std::vector<int> tpr; // priority of a tangle
    std::vector<unsigned int> tdone; // all members in Z in this attraction

    struct PearceFrame { unsigned int vertex, edge; };
    std::vector<PearceFrame> pea_state; // DFS stack
    std::vector<unsigned int> pea_S;
    std::unique_ptr<unsigned int[]> pea_vidx; // rindex
    bitset pea_root;  // root
    int pea_curidx;   // index

    std::vector<int> tangle; // stores the new tangle
    std::vector<int> tangleto;
    bitset escapes; // which escapes we considered

    uintqueue Q;
    bitset G; // the unsolved game
    bitset S; // solved vertices (in queue Q)
    std::vector<bitset*> Rs;
    std::vector<bitset*> Hs;
    std::vector<int> Ps; // priorities
    std::vector<int> heads; // record heads for the tangle learning
    bitset Processed; // already processed vertices of the current attraction
    unsigned int *escs, *esc_epoch; // lazy escape counts, scoped to one attraction
    unsigned int attract_epoch = 0;
    unsigned int *gdeg = nullptr; // number of successors in <G> of each vertex in <G>, or 0 if not known yet
                                  // (none while all successors of all vertices are in <G>)
    bool count_in_G = false; // the current attraction is in <G> (see attract)
    bool top_is_G = false; // the top level is still the whole unsolved game <G>
    std::vector<int> touched; // the vertices counted by the last attraction in <G>
    void updateDegrees(bool counted);
    void remaximizeTop(int pl, bitset &Z, bitset &R, const bitset *removed, unsigned int firstnew);
    bool secondIsWhole();
    bool attract(int pl, bitset &Z, bitset &R, const bitset *crossing = nullptr, bool resume = false);
    int *str;

    bitset W0, W1; // current approximation of winning areas

    bool attracts(const int pl, const int v, bitset &Z, bitset &R);
    void attractVertices(const int pl, const int v, bitset &Z, bitset &R);
    bool attractTangle(int t, int pl, bitset &R, bitset &Z, bitset &G);
    void attractTangles(int pl, int v, bitset &R, bitset &Z, bitset &G);
    void addBacklinks(int tidx);
    bool extractTangles(int startvertex, bitset &R);
    void solve(bitset &Subgame, int vtop, int pe, int po);
};

}

#endif
