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
