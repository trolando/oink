/*
 * Copyright 2017-2018 Tom van Dijk, Johannes Kepler University Linz
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

#ifndef MSPM_HPP
#define MSPM_HPP

#include <vector>

#include "oink/solver.hpp"
#include "oink/uintqueue.hpp"

namespace pg {

/**
 * Small progress measures with strategy derivation for both players, as in
 * M. Gazda and T.A.C. Willemse, "Improvement in Small Progress Measures", GandALF 2015
 * (Algorithm 2, SPM-Within).
 *
 * Only the measure of one player <P> is lifted (Top means that <P> wins). When the first
 * vertex v of the current context W is lifted to Top, its guarded attractor RES (vertices
 * of priority at most that of v) is won by <P>; the context continues in REM, which is W
 * without RES and without the opponent's attractor IRR to the higher priorities. When REM
 * is done, the attractor of <P> to RES and the Top vertices of REM is removed from W.
 * The strategy of <P> follows from these steps, and that of the opponent from the measure.
 * (Priorities as in Oink: the highest priority seen infinitely often decides.)
 */
class MSPMSolver : public Solver
{
public:
    MSPMSolver(Oink& oink, Game& game);
    virtual ~MSPMSolver();

    virtual void run();
    virtual bool parseOptions(std::string& opts);

    long lift_attempt = 0;
    long lift_count = 0;

protected:
    int k = 0; // number of positions in a measure: highest priority + 1 (at least 2)
    int P = 1; // the player of the measure: pms[k*v+P] == -1 (Top) means that <P> wins v
    int force_player = -1; // the player of the measure if set ("even" or "odd"), else the smaller

    std::vector<int> pms;    // the measures, <k> per vertex (positions of parity <P>)
    std::vector<int> counts; // for each priority, the number of vertices (the bound)
    std::vector<int> tmp, best;
    std::vector<int> str;    // the strategy
    std::vector<int> depth;  // the innermost context of each vertex; -1 if removed (won by
                             // <P>), -2 if disabled
    std::vector<int> cnt;    // for attractors: successors in the context not yet attracted
    uintqueue Q;
    bitset dirty, inA, inB;

    struct Frame {
        int D;                  // depth of this context
        std::vector<int> W;     // its vertices (some may have left it)
        std::vector<int> res;   // RES of the current step
        std::vector<int> rem;   // REM of the current step (the context of the child)
        bool child = false;     // whether the child for REM has returned
    };

    bool is_top(int v) const { return pms[(size_t)k*v+P] == -1; }
    bool pm_less(const int *a, const int *b, int d) const;
    void Prog(int *dst, const int *src, int d) const;
    bool lift(int v);
    int lift_context(int D, const std::vector<int> &W);
    void attract(int pl, int D, const std::vector<int> &targets, int max_prio, bitset &in, const bitset *exclude, std::vector<int> &out, bool strategy);
};

}

#endif
