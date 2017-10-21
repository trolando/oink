#include <iostream>
#include <random>

using namespace std;

int
rng(int low, int high)
{
    static random_device rand_dev;
    static mt19937 generator(rand_dev());
    return uniform_int_distribution<int>(low, high)(generator);
}

int
main(int argc, char** argv)
{
    if (argc != 5 && argc != 6) {
        cout << "Syntax: " << argv[0] << " nNodes maxPrio minDeg maxDeg [noSelfLoops]" << endl;
        cout.flush();
        return -1;
    }

    // parse arguments
    int n = stoi(argv[1]);
    int maxP = stoi(argv[2]);
    int minD = stoi(argv[3]);
    int maxD = stoi(argv[4]);
    int SL = argc == 5 ? 1 : 0;

    // check bad arguments
    if (maxD > n) {
        cout << "max degree >= num nodes?!" << endl;
        cout.flush();
        return -1;
    }
    if (minD > maxD) {
        cout << "min degree >= max degree?!" << endl;
        cout.flush();
        return -1;
    }
    if (maxD == n && SL == 0) {
        maxD--; // just fix the number
    }

    // allocate edge target array
    int *tgt = new int[n];

    // generate game
    cout << "parity " << n << ";" << endl;
    for (int i=0; i<n; i++) {
        // generate priority and owner
        cout << i << " " << rng(0, maxP) << " " << rng(0, 1) << " ";
        // determine number of edges
        int e_count = rng(minD, maxD);
        // initialize edge target array
        for (int j=0; j<n; j++) tgt[j] = j;
        // initialize size of array
        int m = n;
        // maybe remove the self-loop
        if (SL == 0) {
            tgt[i] = n-1;
            m--;
        }
        // generate random edges
        for (int j=0; j<e_count; j++) {
            int x = rng(0, m-1);
            if (j != 0) cout << ",";
            cout << tgt[x];
            tgt[x] = tgt[m-1];
            m--;
        }
        // generate no label
        cout << ";" << endl;
    }

    // free edge target array
    delete[] tgt;

    return 0;
}
