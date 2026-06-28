/*
 * Copyright 2024 Tom van Dijk
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

/**
 * Unit tests for the MultiStrategy class.
 */

#include <iostream>

#include "oink/multistrategy.hpp"

using namespace pg;

static int failures = 0;

static void
check(const char* name, bool ok)
{
    if (!ok) {
        std::cerr << "FAIL [" << name << "]" << std::endl;
        failures++;
    }
}

int
main()
{
    // A fresh MultiStrategy of given size has no strategy edges.
    MultiStrategy ms(10);
    check("size", ms.size() == 10);
    check("not empty", !ms.empty());
    check("initially none set", !ms.has(0) and !ms.has(5) and !ms.has(9));

    // A default-constructed one is empty.
    MultiStrategy empty;
    check("default empty", empty.empty() and empty.size() == 0);

    // add/has/remove operate per edge index.
    ms.add(2);
    ms.add(3);
    ms.add(7);
    check("add records edges", ms.has(2) and ms.has(3) and ms.has(7));
    check("others untouched", !ms.has(1) and !ms.has(4) and !ms.has(6));
    ms.remove(3);
    check("remove clears one", !ms.has(3) and ms.has(2) and ms.has(7));

    // clear_range clears a contiguous block (one vertex's edges) only.
    ms.add(3);
    ms.add(4);
    ms.add(5);  // block [2..6) now: 2,4,5,(3) set; 7 outside the block
    ms.clear_range(2, 4); // clears indices 2,3,4,5
    check("clear_range clears block", !ms.has(2) and !ms.has(3) and !ms.has(4) and !ms.has(5));
    check("clear_range keeps outside", ms.has(7));

    // bits() exposes the underlying bitset with the same contents.
    check("bits reflects edges", ms.bits().test(7) and ms.bits().count() == 1);

    // resize keeps existing bits and zeroes new ones.
    ms.resize(20);
    check("resize keeps old", ms.has(7));
    check("resize grows", ms.size() == 20 and !ms.has(15));

    // Copy is independent.
    MultiStrategy copy = ms;
    copy.add(15);
    check("copy independent", ms.has(15) == false and copy.has(15) == true);

    // reset clears everything but keeps the size.
    ms.reset();
    check("reset clears all", !ms.has(7) and ms.bits().count() == 0);
    check("reset keeps size", ms.size() == 20);

    if (failures) {
        std::cerr << failures << " multistrategy test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all multistrategy tests passed" << std::endl;
    return 0;
}
