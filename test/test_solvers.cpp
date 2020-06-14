/*
 * Copyright 2020 Tom van Dijk
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

#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/time.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include "tools/cxxopts.hpp"
#include "game.hpp"
#include "oink.hpp"
#include "solvers.hpp"
#include "verifier.hpp"

using namespace pg;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;

Solvers solvers;

bool opt_inflate = false;
bool opt_compress = false;
bool opt_single = false;
bool opt_loops = false;
bool opt_wcwc = false;
int opt_workers = 0;
int opt_trace = -1;

static double
wctime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + 1E-6 * time.tv_usec;
}


/*------------------------------------------------------------------------*/

volatile bool quit = false;

static void (*sig_int_handler)(int);

static void
catchsig(int sig)
{
    if (sig == SIGINT) {
        (void)signal(SIGINT, sig_int_handler);
        quit = true;
    }
}

static void
setsighandlers(void)
{
    sig_int_handler = signal(SIGINT, catchsig);
}

/*------------------------------------------------------------------------*/


int
test_solver(Game &game, int solverid, double &time, std::ostream &log)
{
    game.reset_solution();
    game.ensure_sorted();

    // solve a copy
    Game copy(game);
    Oink solver(copy, log);
    solver.setRenumber(); // default
    if (opt_inflate) solver.setInflate();
    if (opt_compress) solver.setCompress();
    solver.setSolveSingle(opt_single);
    solver.setRemoveLoops(opt_loops);
    solver.setRemoveWCWC(opt_wcwc);
    solver.setWorkers(opt_workers);
    solver.setSolver(solverid);
    if (opt_trace >= 0) solver.setTrace(opt_trace);
    else solver.setTrace(0);

    double begin = wctime();
    try {
        solver.run();
    } catch (pg::Error &err) {
        log << "solver error: " << err.what() << std::endl;
        return 1;
    }
    time = wctime() - begin;

    game.copy_solution(copy);

    try {
        Verifier v(&game, log);
        v.verify(true, true, true);
    } catch (const char *err) {
        log << "verification error: " << err << std::endl;
        return 2;
    }

    return 0;
}
 

int
main(int argc, char **argv)
{
    cxxopts::Options opts(argv[0], "Test parity game solvers");
    opts.custom_help("[OPTIONS...] [FILES...]");
    opts.add_options()
        ("help", "Print help")
        ;
    opts.add_options("Preprocessing")
        ("inflate", "Inflate the game before solving")
        ("compress", "Compress the game before solving")
        ("single", "Enable preprocessor \"single\" (solve single-parity games)")
        ("loops", "Enable preprocessor \"loops\" (remove/solve self-loops)")
        ("wcwc", "Enable preprocessor \"wcwc\" (solve winner-controlled winning cycles)")
        ;
    opts.add_options("Random games")
        ("count", "Number of random games", cxxopts::value<int>()->default_value("100"))
        ("seed", "Seed for the random seed generator", cxxopts::value<unsigned int>())
        ("size", "Size (number of vertices) of each random game", cxxopts::value<int>()->default_value("100"))
        ("maxp", "Maximum priority of a vertex of each random game (default: <size>)", cxxopts::value<int>())
        ("maxe", "Maximum number of edges of each random game (default: 4*<size>)", cxxopts::value<int>())
        ("gameseed", "Create just one game with the given seed and parameters size, maxp and maxe, and write this game to stdout", cxxopts::value<unsigned int>())
        ;
    opts.add_options("Solvers")
        ("all", "Run all solvers")
        ;
    for (unsigned id=0; id<solvers.count(); id++) {
        opts.add_options("Solvers")(solvers.label(id), solvers.desc(id));
    }
    opts.add_options("Solving")
        ("t,trace", "Write trace with given level (0-3) to stdout", cxxopts::value<int>())
        ("w,workers", "Number of workers for parallel algorithms, or -1 for sequential, 0 for autodetect", cxxopts::value<int>()->default_value("0"))
        ;

    /* Parse command line */
    auto options = opts.parse(argc, argv);

    if (options.count("help")) {
        std::cout << opts.help({"","Preprocessing","Random games","Solving","Solvers"}) << std::endl;
        return 0;
    }

    if (options.count("gameseed")) {
        // create random game and write to stdout
        int size = options["size"].as<int>();
        int maxP = options.count("maxp") ? options["maxp"].as<int>() : size; // default #priorities = size
        int maxE = options.count("maxe") ? options["maxe"].as<int>() : 4*size; // default avg degree [1..4]
        // limit maxE to bound [0 .. size*size]
        if (maxE < 0) maxE = 0;
        if (maxE > (size*size)) maxE = size*size;
        // create a random game
        Game g;
        g.set_random_seed(options["gameseed"].as<unsigned int>());
        g.init_random_game(size, maxP, maxE-size);
        // and write it to stdout
        g.write_pgsolver(std::cout);
        return 0;
    }

    // options
    opt_inflate = options.count("inflate") != 0;
    opt_compress = options.count("compress") != 0;
    opt_single = options.count("single") != 0;
    opt_loops = options.count("loops") != 0;
    opt_wcwc = options.count("wcwc") != 0;
    if (options.count("workers")) opt_workers = options["workers"].as<int>();
    if (options.count("trace")) opt_trace = options["trace"].as<int>();

    std::cout << "Selected solvers:";

    int solver_count = 0;
    for (unsigned id=0; id<solvers.count(); id++) {
        if (options.count("all") or options.count(solvers.label(id))) {
            solver_count++;
            std::cout << " " << solvers.label(id);
        }
    }
    if (solver_count == 0) {
        std::cout << " (none)" << std::endl << std::endl;
        std::cout << "Use --help for program options." << std::endl << std::endl;
        std::cout << "- Select one or more solvers" << std::endl;
        std::cout << "- Select either:" << std::endl;
        std::cout << "  - a list of files/directories to test" << std::endl;
        std::cout << "  - options for random games"<< std::endl;
        return 0;
    } else {
        std::cout << std::endl;
    }

    setsighandlers();

    int final_res = 0;
    std::stringstream log;
    double time;
    long total=0, good=0;

    double times[solvers.count()];
    for (unsigned i=0; i<solvers.count(); i++) times[i] = 0.0;
    int sgood[solvers.count()];
    for (unsigned i=0; i<solvers.count(); i++) sgood[i] = 0;

    if (argc > 1) {
        // obtain the list of files
        std::vector<fs::path> files;
        for (int i=1; i<argc; i++) {
            fs::path p(argv[i]);
            if (!exists(p)) {
                std::cerr << "path \"" << argv[i] << "\" not found!" << std::endl;
            } else if (is_directory(p)) {
                for (auto cp=fs::directory_iterator(p); cp != fs::directory_iterator(); cp++) {
                    if (is_regular_file(*cp)) files.push_back(*cp);
                }
            } else if (is_regular_file(p)) {
                files.push_back(p);
            }
        }
        sort(files.begin(), files.end());
        for (auto &cp : files) {
            std::string filename = cp.leaf().string();
            std::cout << filename << ": " << std::flush;
            io::filtering_istream in;
            if (boost::algorithm::ends_with(filename, ".bz2")) in.push(io::bzip2_decompressor());
            if (boost::algorithm::ends_with(filename, ".gz")) in.push(io::gzip_decompressor());
            std::ifstream inp(cp.c_str(), std::ios_base::binary);
            in.push(inp);
            try {
                Game game;
                game.parse_pgsolver(in, opt_loops);
                inp.close();
                total++;
                for (unsigned id=0; id<solvers.count(); id++) {
                    if (options.count("all") or options.count(solvers.label(id))) {
                        std::cout << std::flush;
                        log.str("");
                        int res = test_solver(game, id, time, opt_trace == -1 ? log : std::cout);
                        if (res == 0) {
                            sgood[id]++;
                            good++;
                            std::cout << "\033[38;5;82m" << solvers.label(id) << "\033[m";
                        } else {
                            final_res = res;
                            std::cout << "\033[38;5;196m" << solvers.label(id) << "\033[m";
                        }
                        std::cout << " \033[38;5;8m(" << std::fixed << std::setprecision(0) << (1000.0*time) << ")\033[m ";
                        times[id] += time;
                    }
                }
                std::cout << std::endl;
            } catch (const char *s) {
                std::cout << s << std::endl;
                std::cout << "not a parity game input?!" << std::endl;
            }
        }
    } else {
        // random
        int n = options["count"].as<int>();
        int size = options["size"].as<int>();
        int maxP = options.count("maxp") ? options["maxp"].as<int>() : size; // default #priorities = size
        int maxE = options.count("maxe") ? options["maxe"].as<int>() : 4*size; // default avg degree [1..4]

        if (maxE < 0) maxE = 0;
        if (maxE > (size*size)) maxE = size*size;

        unsigned int seriesseed;
        if (options.count("seed")) {
            seriesseed = options["seed"].as<unsigned int>();
        } else {
            std::random_device rd;
            seriesseed = rd();
        }

        std::cout << "Creating " << n << " random games: --size=" << size << " --maxp=" << maxP << " --maxe=" << maxE << " --seed=" << seriesseed << std::endl;

        std::mt19937 generator(seriesseed);
        Game g;

        for (int i=0; i<n && !quit; i++) {
            unsigned int seed;
            if (n == 1 and options.count("seed")) seed = options["seed"].as<unsigned int>();
            else seed = generator();
            g.set_random_seed(seed);
            g.init_random_game(size, maxP, maxE-size);

            std::cout << "game " << i << " (gameseed=" << seed << " size=" << g.vertexcount() << "," << g.edgecount() << "): ";
            std::cout << std::endl << std::flush;
            total++;
            for (unsigned id=0; id<solvers.count(); id++) {
                if (options.count("all") or options.count(solvers.label(id))) {
                    std::cout << std::flush;
                    log.str("");
                    int res = test_solver(g, id, time, opt_trace == -1 ? log : std::cout);
                    if (res == 0) {
                        sgood[id]++;
                        good++;
                        std::cout << "\033[38;5;82m" << solvers.label(id) << "\033[m";
                    } else {
                        final_res = res;
                        std::cout << "\033[38;5;196m" << solvers.label(id) << "\033[m";

                        std::ostringstream fn;
                        fn << "bad_" << solvers.label(id) << "_" << i << ".pg";
                        std::ofstream fout(fn.str());
                        g.write_pgsolver(fout);
                        fout.close();

                        fn.str("");
                        fn << "bad_" << solvers.label(id) << "_" << i << ".pg.log";
                        std::ofstream flog(fn.str());
                        flog << log.str();
                        flog.close();
                    }
                    std::cout << " \033[38;5;8m(" << std::fixed << std::setprecision(0) << (1000.0*time) << ")\033[m ";
                    times[id] += time;
                }
            }
            std::cout << std::endl;
        }
    }

    std::cout << "\033[38;5;226msummary\033[m: " << total << " games" << std::endl;
    std::cout << "\033[38;5;226msolvers\033[m:";
    for (unsigned id=0; id<solvers.count(); id++) {
        if (options.count("all") or options.count(solvers.label(id))) {
            if (sgood[id] == total) std::cout << " \033[38;5;82m";
            else std::cout << " \033[38;5;196m";
            std::cout << solvers.label(id) << "\033[m";
            if (sgood[id] != total) std::cout << " (" << (total-sgood[id]) << " bad)";
        }
    }
    std::cout << std::endl;
    std::cout << "\033[38;5;226mtimes\033[m:  ";
    for (unsigned id=0; id<solvers.count(); id++) {
        if (options.count("all") or options.count(solvers.label(id))) {
            std::cout << " \033[38;5;226m" << solvers.label(id) << "\033[m";
            std::cout << " (" << std::fixed << std::setprecision(0) << (1000.0*times[id]) << " ms)\033[m";
        }
    }
    std::cout << std::endl;

    return final_res;
}
