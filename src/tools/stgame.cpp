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

#include <algorithm>
#include <iostream>
#include <random>
#include <cstring>

#include "oink/game.hpp"

using namespace std;

int
rng(int low, int high)
{
    static random_device rand_dev;
    static mt19937 generator(rand_dev());
    return uniform_int_distribution<int>(low, high)(generator);
}

int*
find(int *arr, int len, int val)
{
    for (;;) {
        if (len <= 0) return NULL;
        int mid = len/2;
        int v = arr[mid];
        if (v == val) {
            return arr + mid;
        } else if (v < val) {
            arr += (mid + 1);
            len -= (mid + 1);
        } else {
            len = mid;
        }
    }
}

int
main(int argc, char** argv)
{
    if (argc != 5 && argc != 6) {
        cout << "Syntax: " << argv[0] << " nNodes minOutDeg maxOutDeg minInDeg maxInDeg" << endl;
        return -1;
    }

    // parse arguments
    int n = stoi(argv[1]);
    int minD = stoi(argv[2]);
    int maxD = stoi(argv[3]);
    int minI = stoi(argv[4]);
    int maxI = stoi(argv[5]);

    // check bad arguments
    if (minD < 1) {
        cerr << "min degree < 1?!" << endl;
        return -1;
    }
    if (minD > n) {
        cerr << "min degree > num nodes?!" << endl;
        return -1;
    }
    if (minD > maxD) {
        cerr << "min degree > max degree?!" << endl;
        return -1;
    }
    if (minI < 1) {
        cerr << "min degree < 1?!" << endl;
        return -1;
    }
    if (minI > n) {
        cerr << "min degree > num nodes?!" << endl;
        return -1;
    }
    if (minI > maxI) {
        cerr << "min degree > max degree?!" << endl;
        return -1;
    }
 
    // allocate edge target array
    int *src = new int[n];
    int *tgt = new int[n];
    for (int j=0; j<n; j++) src[j] = j;
    for (int j=0; j<n; j++) tgt[j] = j;
    int src_count = n, tgt_count = n;
    int intodo = n, outtodo = n;

    int *outcounts = new int[n];
    int *incounts = new int[n];

    memset(outcounts, 0, sizeof(int[n]));
    memset(incounts, 0, sizeof(int[n]));

    // array for temporarily removing unwanted target nodes
    int **ptrs = new int*[maxD+1];

    // generate game
    pg::Game game(n);
    game.vec_init();

    // initialize <n> nodes
    for (int i=0; i<n; i++) game.init_vertex(i, i, rng(0, 1));

    while ((outtodo > 0 and tgt_count > 0) or (intodo > 0 and src_count > 0)) {
        // obtain source node
        int from = src[rng(0, src_count-1)];

        // (temporarily) remove target nodes and <from> from tgt array
        int removed = 0;
        ptrs[0] = find(tgt, tgt_count, from);

        int out_len = 0;
        auto edges = game.outedges() + game.firstout(from);
        while (edges[out_len] != -1) {
            ptrs[out_len+1] = find(tgt, tgt_count, edges[out_len]);
            out_len++;
        }

        std::sort(ptrs, ptrs+out_len+1, [] (int* a, int* b) { return a > b; });

        for (int i=0; i<=out_len and ptrs[i] != NULL; i++) {
            std::swap(*ptrs[i], tgt[tgt_count-1]);
            if (tgt_count <= 1) { cerr << "out of bounds\n"; return -1; }
            tgt_count--;
            removed++;
        }

        // obtain target node
        int to = tgt[rng(0, tgt_count-1)];

        // put removed nodes back in
        while (removed != 0) {
            tgt_count++;
            std::swap(*ptrs[removed-1], tgt[tgt_count-1]);
            removed--;
        }

        // add the edge
        if (game.vec_add_edge(from, to)) {
            outcounts[from]++;
            incounts[to]++;
        }

        // update source node counts
        const int from_count = outcounts[from];
        if (from_count == minD) outtodo--;
        if (from_count >= maxD) {
            // remove from src
            int *ptr = find(src, src_count, from);
            std::copy(ptr + 1, src + src_count, ptr);
            src_count--;
        }

        // update target node counts
        const int to_count = incounts[to];
        if (to_count == minI) intodo--;
        if (to_count >= maxI) {
            // remove from tgt
            int *ptr = find(tgt, tgt_count, to);
            std::copy(ptr + 1, tgt + tgt_count, ptr);
            tgt_count--;
        }
    }

    // write game
    game.vec_finish();
    game.write_pgsolver(cout);
   
    // free arrays
    delete[] src;
    delete[] tgt;
    delete[] outcounts;
    delete[] incounts;

    return 0;
}
