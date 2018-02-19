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

#include <boost/dynamic_bitset.hpp>

#include "oink.hpp"
#include "solver.hpp"
#include "uintqueue.hpp"

namespace pg {

class SSPMSolver : public Solver
{
public:
    SSPMSolver(Oink *oink, Game *game);
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
    int l, h0, h1;
    boost::dynamic_bitset<unsigned long long> pm_0_b;
    int *pm_0_d;
    boost::dynamic_bitset<unsigned long long> pm_1_b;
    int *pm_1_d;

    boost::dynamic_bitset<unsigned long long> tmp_b;
    int *tmp_d;

    boost::dynamic_bitset<unsigned long long> best_b;
    int *best_d;

    boost::dynamic_bitset<unsigned long long> test_b;
    int *test_d;

    uintqueue Q;
    boost::dynamic_bitset<unsigned long long> dirty;
    boost::dynamic_bitset<unsigned long long> unstable;

    // Copy pm[idx] into tmp
    void to_tmp_0(int idx);
    void to_tmp_1(int idx);
    // Copy tmp into pm[idx]
    void from_tmp_0(int idx);
    void from_tmp_1(int idx);
    // Copy pm[idx] into best
    void to_best_0(int idx);
    void to_best_1(int idx);
    // Copy best into pm[idx]
    void from_best_0(int idx);
    void from_best_1(int idx);
    // Copy tmp into best
    void tmp_to_best();
    // Copy tmp into test
    void tmp_to_test();

    // Render pm[idx] to given ostream
    void stream_pm_0(std::ostream &out, int idx);
    void stream_pm_1(std::ostream &out, int idx);
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
    void prog_tmp(int pindex);

    // Lift node, triggered by change to target
    bool lift_0(int node, int target, int &str);
    bool lift_1(int node, int target, int &str);

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
};

}

#endif 
