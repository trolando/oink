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

#include "oink/solver.hpp"
#include "oink/oink.hpp"

namespace pg {

Solver::Solver(Oink& oink, Game& game) : game(game), logger(oink.logger), trace(oink.trace), disabled(oink.disabled), oink(oink)
{
#ifndef NDEBUG
    // sanity check if the game is properly sorted
    for (int i=1; i<nodecount(); i++) assert(priority(i-1) <= priority(i));
#endif
}

}
