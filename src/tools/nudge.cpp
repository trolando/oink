#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

#include "cxxopts.hpp"
#include "game.hpp"
#include "scc.hpp"

using namespace std;
using namespace pg;

int
rng(int low, int high)
{
    static random_device rand_dev;
    static mt19937 generator(rand_dev());
    return uniform_int_distribution<int>(low, high)(generator);
}

bool
nudge(Game *game, int profile)
{
    /**
     * Select a random node
     */
    int n = rng(0, game->n_nodes-1);

    /**
     * Select a random action
     * action=0: change owner of the node
     * action=1: change priority of the node
     * action=2: remove the node without forwarding edges
     * action=3: remove the node and forward edges if it has only 1 outgoing edge
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
        // only actions 0 1 2 3 4 5 : only add edge when removing nodes
        action = (int[6]){0,1,2,3,4,5}[rng(0, 5)];
    } else {
        action = rng(0, 6);
    }

    /**
     * Perform the action
     */
    if (action == 0) {
        // change owner of the node
        game->owner[n] = 1 - game->owner[n];
        return true;
    } else if (action == 1) {
        // change priority of the node
        game->priority[n] = rng(0, game->n_nodes*2);
        return true;
    } else if (action == 2) {
        // remove the node
        std::vector<int> tokeep;
        for (int i=0; i<game->n_nodes; i++) if (n != i) tokeep.push_back(i);
        Game *sub = game->extract_subgame(tokeep, NULL);
        *game = *sub;
        delete sub;
        return true;
    } else if (action == 3) {
        // remove the node and forward edges if it has only 1 outgoing edge
        if (game->out[n].size() == 1 && game->out[n][0] != n) {
            for (auto &from : game->in[n]) {
                // add each <to> to <from>
                for (auto &to : game->out[n]) game->addEdge(from, to);
            }
            // remove the node
            std::vector<int> tokeep;
            for (int i=0; i<game->n_nodes; i++) if (n != i) tokeep.push_back(i);
            Game *sub = game->extract_subgame(tokeep, NULL);
            *game = *sub;
            delete sub;
            return true;
        }
    } else if (action == 4) {
        // for one predecessor replace the edge to me by all my edges
        if (game->in[n].size() > 0) {
            int from = game->in[n][rng(0, game->in[n].size()-1)];
            if (from != n) {
                auto &out = game->out[from];
                out.erase(std::remove(out.begin(), out.end(), n), out.end());

                // add each <to> to <from>
                for (auto &to : game->out[n]) {
                    if (std::find(game->out[from].begin(), game->out[from].end(), to) == game->out[from].end()) {
                        game->out[from].push_back(to);
                        game->in[to].push_back(from);
                    }
                }
                return true;
            }
        }
    } else if (action == 5) {
        // remove a random edge (if there are outgoing edges to remove)
        if (game->out[n].size() > 1) {
            int edge = rng(0, game->out[n].size()-1);
            int m = game->out[n][edge];
            game->out[n].erase(game->out[n].begin()+edge);
            auto &in = game->in[m];
            in.erase(std::remove(in.begin(), in.end(), n), in.end());
            return true;
        }
    } else if (action == 6) {
        // add a random edge
        if (game->addEdge(n, rng(0, game->n_nodes-1))) return true;
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
        ("m,modify", "Modify graph with profile 0=only remove, 1=stable #nodes, 2=everything", cxxopts::value<int>())
        ("b,bottom-scc", "Obtain random bottom SCC before writing")
        ("i,inflate", "Inflate before writing")
        ("c,compress", "Compress before writing")
        ("r,renumber", "Renumber before writing")
        ("o,order", "Order nodes by priority before writing")
        ("evenodd", "Swap players")
        ("minmax", "Turn a mingame into a maxgame and vice versa")
        ;
    opts.parse_positional(std::vector<std::string>({"input", "output"}));
    opts.parse(argc, argv);

    if (opts.count("help")) {
        std::cout << opts.help() << std::endl;
        return 0;
    }

    /**
     * Read a game from file or stdin.
     */
    Game *game = new Game();
    try {
        if (opts.count("input")) {
            std::ifstream file(opts["input"].as<std::string>());
            game->parse_pgsolver(file);
            file.close();
        } else {
            game->parse_pgsolver(std::cin);
        }
    } catch (const char *err) {
        std::cerr << "parsing error: " << err << std::endl;
        return -1;
    }

    /**
     * Check if we have to modify the game randomly.
     */
    if (opts.count("modify")) {
        int profile = opts["modify"].as<int>();
        while (nudge(game, profile) == false) {}
    }

    /**
     * If asked, compute a bottom SCC from a random node
     */
    if (opts.count("b")) {
        std::vector<int> scc;
        getBottomSCC(*game, rng(0, game->n_nodes-1), scc, true);
        Game *sub = game->extract_subgame(scc, NULL);
        delete game;
        game = sub;
    }

    /**
     * Reindex before transformations
     */
    int *mapping = new int[game->n_nodes];
    game->reindex(mapping);

    /**
     * Perform even/odd or min/max transformations
     */
    if (opts.count("evenodd")) game->evenodd();
    if (opts.count("minmax")) game->minmax();

    /**
     * Inflate/compress/renumber
     */
    if (opts.count("i")) game->inflate();
    if (opts.count("c")) game->compress();
    if (opts.count("r")) game->renumber();

    /**
     * Either reindex again (-o) or undo previous reindex
     */
    if (opts.count("o")) game->reindex();
    else game->permute(mapping);

    /**
     * Write to output file or to stdout
     */
    if (opts.count("output")) {
        std::ofstream file(opts["output"].as<std::string>());
        game->write_pgsolver(file);
        file.close();
    } else {
        game->write_pgsolver(std::cout);
    }

    delete game;
    delete[] mapping;
    return 0;
}
