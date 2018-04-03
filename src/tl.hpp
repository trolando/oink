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
#include <boost/dynamic_bitset.hpp>

#include "oink.hpp"
#include "solver.hpp"
#include "npp.hpp"

namespace pg {

class TLSolver : public Solver
{
public:
    TLSolver(Oink *oink, Game *game);
    virtual ~TLSolver();

    virtual void run();

protected:
    int max_prio;

    /**
     * Statistics
     */
    int tangles;     // number of tangles (not dominions)
    int iterations;  // number of iterations (standard)
    int turns;       // number of turns (alternating)
    int dominions;   // number of dominions

    /**
     * Control whether we use a greedy approach (false) or alternating (true)
     */
    bool alternating;
    bool onthefly;

    int *inverse; // reverse lookup priority->vertex
    int *region; // current region of each vertex
    int *strategy; // current strategy of each vertex

    std::vector<int*> vout; // for each tangle
    std::vector<int> *vin; // for each normal vertex
    std::vector<int*> vv; // the tangle (vertex-strategy pairs)
    std::vector<int> vp; // priority of a tangle
    std::vector<int> vr; // region of a tangle

    std::vector<int> *regions; // current regions
    std::vector<int> *vregions; // current regions (for tangles)

    uintqueue Q;
    uintqueue tarres;
    uintqueue tangleto;
    boost::dynamic_bitset<unsigned long long> bs_exits;

    inline void attractTo(const int pr, const int pl, int cur);
    void attract(int prio); // perform the attractor computation
    bool computeRegion(int i); // compute next region, returns true if nonempty
    int extractTangles(int i, bool isHighest); // extract tangles
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
