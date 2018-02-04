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

#ifndef VERIFIER_HPP
#define VERIFIER_HPP

#include <istream>
#include "game.hpp"

namespace pg {

class Verifier
{
public:
    Verifier(Game* pg, std::ostream &lgr) : game(pg), logger(lgr) { }

    void verify(bool fullgame=true, bool even=true, bool odd=true);

    int n_strategies = 0;

protected:
    Game *game;
    std::ostream &logger;
};

}

#endif 
