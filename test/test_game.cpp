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
 * Unit tests for Game dynamic growth.
 */

#include <iostream>

#include "oink/game.hpp"

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
    // Start from a single-vertex game and add more vertices. This grows the
    // vertex arrays via v_sizeup; with a 1-vertex start the growth step must
    // still advance (a buggy "size += size/2" would stay at 1 and loop forever).
    Game g(1);
    g.init_vertex(0, 5, 0);
    g.init_vertex(1, 6, 1); // triggers growth from capacity 1
    g.init_vertex(2, 7, 0); // triggers growth again

    check("grew to hold all vertices", g.vertexcount() >= 3);
    check("priorities preserved", g.priority(0) == 5 and g.priority(1) == 6 and g.priority(2) == 7);
    check("owners preserved", g.owner(0) == 0 and g.owner(1) == 1 and g.owner(2) == 0);

    if (failures) {
        std::cerr << failures << " game test(s) failed" << std::endl;
        return 1;
    }
    std::cout << "all game tests passed" << std::endl;
    return 0;
}
