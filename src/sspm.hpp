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

#ifndef SSPM_HPP
#define SSPM_HPP

#include "oink/solver.hpp"

namespace pg {

class SSPMSolver : public Solver
{
public:
    SSPMSolver(Oink& oink, Game& game);
    virtual ~SSPMSolver();

    virtual void run();

protected:
    /**
     * So the maximal length of bitstrings l := ceil(log2(n_vertices))
     * And (for even measures) the depth h := # of even priorities
     *
     * So we could reasonably assume 0 <= l <= 32
     *
     * Two formats are possible to encode the array, apart from the bitstring
     * - for each bit the depth: log n * log h 
     */
    int l, h;
    bitset pm_b;
    int *pm_d;

    bitset tmp_b;
    int *tmp_d;

    bitset best_b;
    int *best_d;

    bitset test_b;
    int *test_d;

    uintqueue Q;
    bitset dirty;
    bitset unstable;

    bool bounded = false;

    int *cap; // caps!
    uint64_t *lift_counters;

    // Copy pm[idx] into tmp
    void to_tmp(int idx);
    // Copy tmp into pm[idx]
    void from_tmp(int idx);
    // Copy pm[idx] into best
    void to_best(int idx);
    // Copy best into pm[idx]
    void from_best(int idx);
    // Copy tmp into best
    void tmp_to_best();
    // Copy tmp into test
    void tmp_to_test();

    // Render pm[idx] to given ostream
    void stream_pm(std::ostream &out, int idx);
    // Render tmp to given ostream
    void stream_tmp(std::ostream &out, int h);
    // Render best to given ostream
    void stream_best(std::ostream &out, int h);

    // Compare tmp to best
    int compare(int pindex);
    // Compare tmp to test
    int compare_test(int pindex);

    // Bump tmp given priority p
    void trunc_tmp(int pindex);
    void prog_tmp(int pindex, int h);
    void prog_cap_tmp(int pindex);

    // Lift node, triggered by change to target
    bool lift(int node, int target, int &str, int pl);

    inline void todo_push(int node) {
        if (dirty[node]) return;
        Q.push(node);
        dirty[node] = true;
#ifndef NDEBUG
        if (trace >= 2) logger << "push(" << node << ")" << std::endl;
#endif
    }

    inline int todo_pop() {
        int node = Q.pop();
        dirty[node] = false;
#ifndef NDEBUG
        if (trace >= 2) logger << "pop() => " << node << std::endl;
#endif
        return node;
    }

    int lift_count = 0;
    int lift_attempt = 0;

    void run(int nbits, int depth, int player);
};

class BoundedSSPMSolver : public SSPMSolver
{
public:
    BoundedSSPMSolver(Oink& oink, Game& game) : SSPMSolver(oink, game) { bounded = true; }
    virtual ~BoundedSSPMSolver() { }
};

}

#endif 
