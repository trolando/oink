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

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include <stack>

#include "cxxopts.hpp"
#include "oink/game.hpp"
#include "oink/pgparser.hpp"

using namespace std;
using namespace pg;

int
rng(int low, int high)
{
    static random_device rand_dev;
    static mt19937 generator(rand_dev());
    return uniform_int_distribution<int>(low, high)(generator);
}

/**
 * Tarjan's SCC algorithm, modified to only compute the bottom SCC
 */
void
tarjan(Game *game, int n, std::vector<int> &res, bool nonempty)
{
    // initialize
    unsigned n_vertices = game->vertexcount();
    int *low = new int[n_vertices];
    memset(low, 0, sizeof(int[n_vertices]));
    int pre = 0;

    // search stack "st"
    std::stack<int> st;
    st.push(n);

    while (!st.empty()) { // st is never empty
        int idx = st.top();
        if (low[idx] == 0) {
            // first time we see it
            low[idx] = ++pre;
            res.push_back(idx);
        }
        int min = low[idx];
        bool pushed = false;
        for (auto curedge = game->outs(idx); *curedge != -1; curedge++) {
            int to = *curedge;
            if (low[to] == 0) {
                // not visited
                st.push(to);
                pushed = true;
                break;
            } else {
                // visited it, update min
                if (low[to] < min) min = low[to];
            }
        }
        if (pushed) continue; // we pushed...

        if (min < low[idx]) {
            low[idx] = min;
            st.pop();
            continue;
        }

        /**
         * At this point, we found a bottom SCC. Now check if it is empty.
         * A SCC is "empty" if it contains no edges, i.e., consists of 1 vertex without self-loops.
         */

        if (nonempty) {
            if (res.back() == idx and !game->has_edge(idx, idx)) {
                // it has no edges!
                res.pop_back();
                st.pop();
                if (st.empty()) break;
                continue;
            }
        }

        // found bottom SCC; remove vertices not in the bottom SCC...
        if (res.front() != idx) res.erase(res.begin(), std::find(res.begin(), res.end(), idx));
        break;
    }

    delete[] low;
}


bool
nudge(Game *game, int profile)
{
    /**
     * Select a random vertex
     */
    int n = rng(0, game->vertexcount()-1);

    /**
     * Select a random action
     * action=0: change owner of the vertex
     * action=1: change priority of the vertex
     * action=2: remove the vertex without forwarding edges
     * action=3: remove the vertex and forward edges if it has only 1 outgoing edge
     * action=4: for one predecessor replace the edge to me by all my edges
     * action=5: remove a random edge (if 2+ outgoing edges)
     * action=6: add a random edge
     */
    int action;
    if (profile == 0) {
        // only actions 0 1 2 3 5 : only change/remove, no adding edges
        action = (int[5]){0,1,2,3,5}[rng(0, 4)];
    } else if (profile == 1) {
        // only actions 0 1 5 6 : just change edges and owner/priority
        action = (int[4]){0,1,5,6}[rng(0, 3)];
    } else if (profile == 2) {
        // only actions 0 1 2 3 4 5 : only add edge when removing vertices
        action = (int[6]){0,1,2,3,4,5}[rng(0, 5)];
    } else {
        action = rng(0, 6);
    }

    /**
     * Perform the action
     */
    if (action == 0) {
        // flip owner of the vertex
        game->set_owner(n, 1-game->owner(n));
        return true;
    } else if (action == 1) {
        // change priority of the vertex
        game->set_priority(n, rng(0, game->vertexcount()*2));
        return true;
    } else if (action == 2) {
        // remove the vertex, but only if no incoming edge is a solo edge
        for (auto curedge = game->ins(n); *curedge != -1; curedge++) {
            if (game->outcount(*curedge) == 1) return false;
        }
        pg::bitset mask(game->vertexcount());
        mask.set();
        mask[n] = false;
        auto subgame = game->extract_subgame(mask);
        game->swap(*subgame);
        return true;
    } else if (action == 3) {
        // remove the vertex and forward incoming edges if it has only 1 outgoing edge
        if (game->outcount(n) == 1 and game->outs(n)[0] != n) {
            game->vec_init();
            for (auto curedge = game->ins(n); *curedge != -1; curedge++) {
                int from = *curedge;
                // add each edge "n -> to" to "from -> to"
                for (auto curedge = game->outs(n); *curedge != -1; curedge++) {
                    game->vec_add_edge(from, *curedge);
                }
            }
            game->vec_finish();
            // remove the vertex
            pg::bitset mask(game->vertexcount());
            mask.set();
            mask[n] = false;
            auto subgame = game->extract_subgame(mask);
            game->swap(*subgame);
            return true;
        }
    } else if (action == 4) {
        // for one random predecessor replace the edge to me by all my edges
        if (game->incount(n) > 0) {
            int from = game->ins(n)[rng(0, game->incount(n)-1)];
            if (from != n) {
                game->vec_init();
                game->vec_remove_edge(from, n);
                for (auto curedge = game->outs(n); *curedge != -1; curedge++) {
                    game->vec_add_edge(from, *curedge);
                }
                game->vec_finish();
                return true;
            }
        }
    } else if (action == 5) {
        // remove a random edge (if there are outgoing edges to remove)
        if (game->outcount(n) > 1) {
            int edge = rng(0, game->outcount(n)-1);
            int m = game->outs(n)[edge];
            game->vec_init();
            game->vec_remove_edge(n, m);
            game->vec_finish();
            return true;
        }
    } else if (action == 6) {
        // add a random edge
        game->vec_init();
        bool res = game->vec_add_edge(n, rng(0, game->vertexcount()-1));
        game->vec_finish();
        if (res) return true;
    }
    return false;
}

/**
 * Change the game a bit.
 */
int
main(int argc, char **argv)
{
    cxxopts::Options opts(argv[0], "Parity game solver");
    opts.add_options()
        ("input", "Input parity game", cxxopts::value<std::string>())
        ("output", "Output parity game", cxxopts::value<std::string>())
        ("help", "Print help")
        ("m,modify", "Modify graph with profile 0=only remove, 1=stable #vertices, 2=everything", cxxopts::value<int>())
        ("b,bottom-scc", "Obtain random bottom SCC before writing")
        ("i,inflate", "Inflate before writing")
        ("c,compress", "Compress before writing")
        ("r,renumber", "Renumber before writing")
        ("o,order", "Order vertices by priority before writing")
        ("u,unlabel", "Remove labels")
        ("evenodd", "Swap players")
        ("minmax", "Turn a mingame into a maxgame and vice versa")
        ;
    opts.parse_positional(std::vector<std::string>({"input", "output"}));
    auto options = opts.parse(argc, argv);

    if (options.count("help")) {
        std::cout << opts.help() << std::endl;
        return 0;
    }

    /**
     * Read a game from file or stdin.
     */
    Game *game = new Game();
    try {
        if (options.count("input")) {
            std::ifstream file(options["input"].as<std::string>());
            *game = PGParser::parse_pgsolver(file, false);
            file.close();
        } else {
            *game = PGParser::parse_pgsolver(std::cin, false);
        }
    } catch (std::runtime_error &err) {
        std::cerr << "parsing error: " << err.what() << std::endl;
        return -1;
    }

    /**
     * Check if we have to modify the game randomly.
     */
    if (options.count("modify")) {
        int profile = options["modify"].as<int>();
        while (nudge(game, profile) == false) {}
    }

    /**
     * If asked, compute a bottom SCC from a random vertex
     */
    if (options.count("b")) {
        std::vector<int> scc;
        tarjan(game, rng(0, game->vertexcount()-1), scc, true);
        auto sub = game->extract_subgame(scc);
        game->swap(*sub);
    }

    /**
     * Reindex before transformations
     */
    int *mapping = new int[game->vertexcount()];
    game->sort(mapping);

    /**
     * Perform even/odd or min/max transformations
     */
    if (options.count("evenodd")) game->evenodd();
    if (options.count("minmax")) game->minmax();

    /**
     * Inflate/compress/renumber
     */
    if (options.count("i")) game->inflate();
    if (options.count("c")) game->compress();
    if (options.count("r")) game->renumber();

    /**
     * Either reindex again (-o) or undo previous reindex
     */
    if (options.count("o")) game->sort();
    else game->permute(mapping);

    /**
     * Remove labels
     */
    if (options.count("u")) {
        for (int i=0; i<game->vertexcount(); i++) game->set_label(i, "");
    }

    /**
     * Write to output file or to stdout
     */
    if (options.count("output")) {
        std::ofstream file(options["output"].as<std::string>());
        game->write_pgsolver(file);
        file.close();
    } else {
        game->write_pgsolver(std::cout);
    }

    delete game;
    delete[] mapping;
    return 0;
}
