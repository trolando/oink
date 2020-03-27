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

#ifndef TL_HPP
#define TL_HPP

#include <stack>
#include <set>
#include <map>
#include <tuple>

#include "oink.hpp"
#include "solver.hpp"
#include "npp.hpp"

namespace pg {

/**
 * Implementation of the tangle learning solver.
 * [2018] Tom van Dijk, Attracting Tangles to Solve Parity Games, CAV.
 */

class TLSolver : public Solver
{
public:
    TLSolver(Oink *oink, Game *game);
    virtual ~TLSolver();

    /**
     * Run the solver.
     */
    virtual void run();

protected:
    /**
     * Statistics
     */
    int tangles;     // number of tangles (not dominions)
    int iterations;  // number of iterations (standard)
    int turns;       // number of turns (alternating)
    int dominions;   // number of dominions

    /**
     * Controls which variation of tangle learning we do...
     */
    bool alternating; // alternating tangle learning
    bool onthefly; // on-the-fly tangle learning

    int *inverse; // reverse lookup priority->vertex
    int *region; // current region of each vertex
    int *strategy; // current strategy of each vertex

    std::vector<int*> vout; // exits for each tangle
    std::vector<int> *vin; // for each normal vertex (inverse of <vout>)
    std::vector<int*> vv; // the tangle (vertex-strategy pairs, v0-s0-v1-s1-...-<-1>
    std::vector<int> vp; // priority of a tangle
    std::vector<int> vr; // current region of each tangle

    std::vector<int> *regions; // current regions
    std::vector<int> *vregions; // current regions (for tangles)

    uintqueue Q; // queue (used by attractor)
    uintqueue tarres; // helper structure for tarjan bottom SCC algorithm
    uintqueue tangleto; // helper structure to record tangle exits
    bitset bs_exits; // helper structure

    /**
     * Attract to vertex <v> as player <pl> for a region of priority <pr>.
     */
    inline void attractTo(const int pr, const int pl, int cur);

    /**
     * Run the attractor computation for current priority <prio>.
     */
    void attract(int prio); // perform the attractor computation

    /**
     * Compute the next region (starting at vertex <i>)
     * Returns <true> if the region is nonempty.
     */
    bool computeRegion(int i);

    /**
     * Extract tangles for the region whose highest vertex is vertex <i>.
     * (isHighest is set when we are the highest region in the game.)
     *
     * Returns:
     * -2 if there were no tangles
     * -1 if there was a dominion
     * <pr> the highest attracting region if there were tangles
     */
    int extractTangles(int i, bool isHighest);
};

class ATLSolver : public TLSolver
{
public:
    ATLSolver(Oink *oink, Game *game) : TLSolver(oink, game) { alternating = true; onthefly = false; }
    virtual ~ATLSolver() { }
};

class OTFATLSolver : public TLSolver
{
public:
    OTFATLSolver(Oink *oink, Game *game) : TLSolver(oink, game) { alternating = true; onthefly = true; }
    virtual ~OTFATLSolver() { }
};

class OTFTLSolver : public TLSolver
{
public:
    OTFTLSolver(Oink *oink, Game *game) : TLSolver(oink, game) { alternating = false; onthefly = true; }
    virtual ~OTFTLSolver() { }
};

}

#endif
