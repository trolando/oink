/*
 * Copyright 2021 Paweł Parys
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
 *
 * File prepared by Paweł Parys.
 */

#ifndef ZLKPP_HPP
#define ZLKPP_HPP

#include <vector>

#include "oink/solver.hpp"

namespace pg {

#define ZLK_STANDARD 0
#define ZLK_WARSAW 1
#define ZLK_LIVERPOOL 2

class ZLKPPSolver : public Solver
{
public:
    ZLKPPSolver(Oink& oink, Game& game, int variant);

    virtual void run();

private:
    int variant, min_dominion;
    unsigned long long iterations = 0;

    // nodes currently in the game:
    bool *cur_nodes_bm; // as a bitmap
    int *cur_nodes_next, *cur_nodes_prev; // on an ordered two-way list
    int cur_first_node; // first node in the list (undefined if cur_num_nodes == 0)
    int cur_num_nodes; // number of those nodes

    //auxiliary arrays for computing attractors:
    int *num_successors;
    bool *is_in_attractor;
    
    int *strategy;
    
    bool get_attractor(int player, std::vector<int> &nodes); // returns whether the attractor is larger than the input set
    
    void remove_nodes(const std::vector<int> nodes);
    void restore_nodes(const std::vector<int> nodes);
    std::vector<int> get_cur_nodes();
    void set_cur_nodes(const std::vector<int> nodes); // it assumes that the "current game" is a subset of "nodes", and that "nodes" are sorted
    
    std::vector<int> get_nodes_of_max_priority(int max_priority);
    
    bool do_step(int max_priority, int prec_cur, int prec_opo, int &num_nodes_in_h, bool &reached_bottom_cur, bool &reached_bottom_opo);

    void solve_liverpool(int max_priority, int prec_cur, int prec_opo, bool &reached_bottom_cur, bool &reached_bottom_opo);
    std::vector<int> solve(int max_priority, int prec_cur, int prec_opo, bool &reached_bottom_cur, bool &reached_bottom_opo);
};

}

#endif
