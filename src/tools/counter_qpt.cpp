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
        cout << i-1 << " " << i << " " << 0 << " ";
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
