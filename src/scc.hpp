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

#ifndef SCC_HPP
#define SCC_HPP

#include <vector>
#include "game.hpp"

namespace pg {

/**
 * Find a bottom SCC starting from the first non-disabled node.
 * Avoids "disabled" nodes.
 * (if nonempty is set, only obtain a non-empty bottom SCC.)
 */
void getBottomSCC(Game &game, std::vector<int> &scc, bool nonempty=false);

/**
 * Find a bottom SCC starting from the given node.
 * Avoids "disabled" nodes.
 * (if nonempty is set, only obtain a non-empty bottom SCC.)
 */
void getBottomSCC(Game &game, int start, std::vector<int> &scc, bool nonempty=false);

}

#endif 
