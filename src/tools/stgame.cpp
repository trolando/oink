#include <algorithm>
#include <iostream>
#include <random>

#include "game.hpp"

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

    // array for temporarily removing unwanted target nodes
    int **ptrs = new int*[maxD+1];

    // generate game
    pg::Game game(n);

    // initialize <n> nodes
    for (int i=0; i<n; i++) game.initNode(i, i, rng(0, 1));

    while ((outtodo > 0 and tgt_count > 0) or (intodo > 0 and src_count > 0)) {
        // obtain source node
        int from = src[rng(0, src_count-1)];

        // (temporarily) remove target nodes and <from> from tgt array
        int removed = 0;
        ptrs[0] = find(tgt, tgt_count, from);
        const int of_size = game.out[from].size();
        for (int i=0; i<of_size; i++) ptrs[i+1] = find(tgt, tgt_count, game.out[from][i]);
        std::sort(ptrs, ptrs+of_size+1, [] (int* a, int* b) { return a > b; });

        for (int i=0; i<=of_size and ptrs[i] != NULL; i++) {
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
        game.addEdge(from, to);

        // update source node counts
        const int from_count = game.out[from].size();
        if (from_count == minD) outtodo--;
        if (from_count >= maxD) {
            // remove from src
            int *ptr = find(src, src_count, from);
            std::copy(ptr + 1, src + src_count, ptr);
            src_count--;
        }

        // update target node counts
        const int to_count = game.in[to].size();
        if (to_count == minI) intodo--;
        if (to_count >= maxI) {
            // remove from tgt
            int *ptr = find(tgt, tgt_count, to);
            std::copy(ptr + 1, tgt + tgt_count, ptr);
            tgt_count--;
        }
    }

    // write game
    game.write_pgsolver(cout);
   
    // free arrays
    delete[] src;
    delete[] tgt;

    return 0;
}
