/*
 * Copyright 2021-2022 Tom van Dijk
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

#ifndef SSI_HPP
#define SSI_HPP

#include "oink/solver.hpp"

namespace pg {

class SSISolver : public Solver
{
public:
    SSISolver(Oink& oink, Game& game);
    virtual ~SSISolver();

    virtual void run();

    uintqueue Q;         // main queue
    bitset G;            // remaining game to solve
    bitset C;            // mark vertices on a cycle
    bitset V;            // helper set
    bitset W0, W1;       // won by 0 or 1
    bitset halt0, halt1; // whether halted by player 0/1
    bitset V0, V1;       // controller by 0 or 1

    int k;               // k := 1+pr(G)
    int *str0, *str1;    // strategy (player/opponent)
    int *val;            // current value, k-tuple for each vertex

    int *first_in;       // helper for linked-list style in-edges
    int *next_in;        // helper for linked-list style in-edges

    bool si_val_less(int a, int b);    // is val[a] < val[b] from Even's perspective
    void compute_vals_ll(int side);    // compute values of finite paths
    int mark_solved(int side);         // set cycles to won by player
    int switch_opp_strategy(int side); // switch opponent strategy
    int switch_sym_strategy(int side); // switch primary strategy
};

}

#endif 
