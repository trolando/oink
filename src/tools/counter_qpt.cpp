#include <iostream>

using namespace std;

int
main(int argc, char** argv)
{
    if (argc != 2) {
        cout << "Syntax: " << argv[0] << " N" << endl;
        return -1;
    }

    int n = stoi(argv[1]);

    // generate the Nth game
    cout << "parity " << 2*n << ";" << endl;
    for (int i=1; i<=2*n; i++) {
        // node i-1 has priority i, owned by 1
        cout << i-1 << " " << i << " " << 1 << " ";
        // if even, can go to 1
        if ((i&1) == 0) {
            cout << "0";
            if (i < 2*n) cout << ",";
        }
        if (i < 2*n) {
            cout << i;
        }
        // node i has label i
        cout << " \"" << i << "\";" << endl;
    }
}
