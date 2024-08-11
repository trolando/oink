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
#include <random>
#include <sstream>
#include <sys/time.h>
#include <optional>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/random/random_device.hpp>
#include <boost/process.hpp>

#include "tools/cxxopts.hpp"
#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "oink/solver.hpp"
#include "verifier.hpp"
#include "lace.h"
#include "oink/pgparser.hpp"

using namespace pg;
namespace fs = boost::filesystem;
namespace io = boost::iostreams;
namespace bp = boost::process;

bool opt_inflate = false;
bool opt_compress = false;
bool opt_single = false;
bool opt_loops = false;
bool opt_wcwc = false;
bool opt_sort = false;
int opt_workers = 0;
int opt_trace = -1;
std::optional<std::string> opt_solver_opts = {};

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

class ExternalSolver : public Solver {
public:
    ExternalSolver(Oink& oink, Game& game, std::string cmd) : Solver(oink, game), cmd(cmd) {}

    void run() {
        // TODO: handle partial files (where some nodes are disabled / only solving a subgame)

        /*
        // Extract the subgame if necessary
        if (disabled.any()) {
            // only keep the ones that are not disabled
            bitset selected(disabled);
            selected.flip();
            auto g = game->extract_subgame(selected);
        }
        */

        // Obtain temporary filenames
        auto path1 = fs::unique_path();
        auto path2 = fs::unique_path();

        // Write parity game to temporary file
        std::ofstream file(path1.native());
        game.write_pgsolver(file);
        file.close();

        // Replace %I with input filename and %O with output filename
        replaceAll(cmd, "%I", path1.native());
        replaceAll(cmd, "%O", path2.native());

        // Run the external program (no timeout)
        bp::ipstream out;
        bp::ipstream err;
        bp::system("/bin/sh", "-c", cmd, bp::std_out > out, bp::std_err > err);

        // Send any stdout/stderr contents to the logger
        std::string line;
        while (std::getline(out, line)) logger << line << std::endl;
        while (std::getline(err, line)) logger << line << std::endl;

        // Read solution from temporary file
        std::ifstream solution(path2.native());
        game.parse_solution(solution);
        solution.close();

        // Delete temporary files
        fs::remove(path1);
        fs::remove(path2);
    }

private:
    std::string cmd;

    static void replaceAll(std::string& str, const std::string& from, const std::string& to) {
        size_t start_pos = 0;
        while((start_pos = str.find(from, start_pos)) != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            start_pos += to.length(); // Move past the replaced substring
        }
    }
};


int
test_solver(Game &game, const std::string& solverid, double &time, std::ostream &log)
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
    if (opt_solver_opts.has_value()) solver.setSolverOptions(*opt_solver_opts);

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
        Verifier v(game, log);
        v.verify(true, true, true);
    } catch (std::runtime_error &err) {
        log << "verification error: " << err.what() << std::endl;
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
        ("sort", "Sort the list of files for solving (if given a list of files)")
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
    for (const auto& id : Solvers::getSolverIDs()) {
        opts.add_options("Solvers")(id, Solvers::desc(id));
    }
    opts.add_options("Solvers")("e,external", "External solver, e.g., 'python solver.py %I %O'", cxxopts::value<std::vector<std::string>>());
    opts.add_options("Solving")
        ("t,trace", "Write trace with given level (0-3) to stdout", cxxopts::value<int>())
        ("c,configure", "Additional configuration options for the solver", cxxopts::value<std::string>())
        ("w,workers", "Number of workers for parallel algorithms, or -1 for sequential, 0 for autodetect", cxxopts::value<int>()->default_value("-1"))
        ;
    opts.allow_unrecognised_options();

    /* Parse command line */
    auto options = opts.parse(argc, argv);
    auto unmatched = options.unmatched();

    if (options.count("help")) {
        std::cout << opts.help({"","Preprocessing","Random games","Solving","Solvers"}) << std::endl;
        return 0;
    }

    if (options.count("gameseed")) {
        // create random game and write to stdout
        long size = options["size"].as<int>();
        long maxP = options.count("maxp") ? options["maxp"].as<int>() : size; // default #priorities = size
        long maxE = options.count("maxe") ? options["maxe"].as<int>() : size*4L; // default avg degree [1..4]
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
    opt_sort = options.count("sort") != 0;
    if (options.count("workers")) opt_workers = options["workers"].as<int>();
    if (options.count("trace")) opt_trace = options["trace"].as<int>();

    std::cout << "Selected solvers:";

    std::vector<std::string> solvers;

    for (const auto& id : Solvers::getSolverIDs()) {
        if (options.count("all") or options.count(id)) {
            solvers.push_back(id);
            std::cout << " " << id;
        }
    }
    if (options.count("external")) {
        for (auto& e : options["external"].as<std::vector<std::string>>()) {
            Solvers::add(e, "external tool", 0, [&](Oink& oink, Game& game) {
                return std::make_unique<ExternalSolver>(oink, game, e);
            });
            solvers.push_back(e);
            std::cout << " '" << e << "'";
        }
    }
    if (solvers.size() == 0) {
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

    if (options.count("configure")) {
        if (solvers.size() > 1) {
            std::cout << std::endl << "Solver configuration options (--configure, -c) cannot be used with multiple solvers at the same time" << std::endl;
            return 0;
        }
        opt_solver_opts = { options["configure"].as<std::string>() };
    }

    setsighandlers();

    if (opt_workers >= 0) lace_start(opt_workers, 10000000UL);

    int final_res = 0;
    std::stringstream log;
    double time;
    long total=0;

    std::map<std::string, double> times;
    std::map<std::string, int> sgood;
    for (auto& id : solvers) {
        times[id] = 0.0;
        sgood[id] = 0;
    }

    if (unmatched.size() > 0) {
        // obtain the list of files
        std::vector<fs::path> files;
        for (auto item : unmatched) {
            fs::path p(item);
            if (!exists(p)) {
                std::cerr << "path \"" << item << "\" not found!" << std::endl;
            } else if (is_directory(p)) {
                for (auto cp=fs::directory_iterator(p); cp != fs::directory_iterator(); cp++) {
                    if (is_regular_file(*cp)) files.push_back(*cp);
                }
            } else if (is_regular_file(p)) {
                files.push_back(p);
            }
        }
        if (opt_sort) {
            sort(files.begin(), files.end());
        } else {
            std::mt19937 rng(std::time(nullptr));
            std::shuffle(files.begin(), files.end(), rng);
        }
        for (auto &cp : files) {
            std::string filename = cp.filename().string();
            std::cout << filename << ": " << std::flush;
            io::filtering_istream in;
            if (boost::algorithm::ends_with(filename, ".bz2")) in.push(io::bzip2_decompressor());
            if (boost::algorithm::ends_with(filename, ".gz")) in.push(io::gzip_decompressor());
            std::ifstream inp(cp.c_str(), std::ios_base::binary);
            in.push(inp);
            try {
                Game game = PGParser::parse_pgsolver_renumber(in, opt_loops);
                inp.close();
                total++;
                for (const auto& id : solvers) {
                    std::cout << std::flush;
                    log.str("");
                    int res = test_solver(game, id, time, opt_trace == -1 ? log : std::cout);
                    if (res == 0) {
                        sgood[id]++;
                        std::cout << "\033[38;5;82m" << id << "\033[m";
                    } else {
                        final_res = res;
                        std::cout << "\033[38;5;196m" << id << "\033[m";
                    }
                    std::cout << " \033[38;5;8m(" << std::fixed << std::setprecision(0) << (1000.0*time) << ")\033[m ";
                    times[id] += time;
                }
                std::cout << std::endl;
            } catch (std::runtime_error &err) {
                std::cout << err.what() << std::endl;
                std::cout << "not a parity game input?!" << std::endl;
            }
        }
    } else {
        // random
        uint64_t n = options["count"].as<int>();
        uint64_t size = options["size"].as<int>();
        uint64_t maxP = options.count("maxp") ? options["maxp"].as<int>() : size; // default #priorities = size
        uint64_t maxE = options.count("maxe") ? options["maxe"].as<int>() : 4L*size; // default avg degree [1..4]

        if (maxE > (size*size)) maxE = size*size;

        unsigned int seriesseed;
        if (options.count("seed")) {
            seriesseed = options["seed"].as<unsigned int>();
        } else {
            boost::random::random_device rd;
            seriesseed = rd();
        }

        std::cout << "Creating " << n << " random games: --size=" << size << " --maxp=" << maxP << " --maxe=" << maxE << " --seed=" << seriesseed << std::endl;

        boost::random::mt19937 generator(seriesseed);
        Game g;

        for (unsigned int i=0; i<n && !quit; i++) {
            unsigned int seed;
            if (n == 1 and options.count("seed")) seed = options["seed"].as<unsigned int>();
            else seed = generator();
            g.set_random_seed(seed);
            g.init_random_game(size, maxP, maxE-size);

            std::cout << "game " << i << " (gameseed=" << seed << " size=" << g.vertexcount() << "," << g.edgecount() << "): ";
            std::cout << std::endl << std::flush;
            total++;
            for (const auto& id : solvers) {
                std::cout << std::flush;
                log.str("");
                int res = test_solver(g, id, time, opt_trace == -1 ? log : std::cout);
                if (res == 0) {
                    sgood[id]++;
                    std::cout << "\033[38;5;82m" << id << "\033[m";
                } else {
                    final_res = res;
                    std::cout << "\033[38;5;196m" << id << "\033[m";

                    std::ostringstream fn;
                    fn << "bad_" << id << "_" << i << ".pg";
                    std::ofstream fout(fn.str());
                    g.write_pgsolver(fout);
                    fout.close();

                    fn.str("");
                    fn << "bad_" << id << "_" << i << ".pg.log";
                    std::ofstream flog(fn.str());
                    flog << log.str();
                    flog.close();
                }
                std::cout << " \033[38;5;8m(" << std::fixed << std::setprecision(0) << (1000.0*time) << ")\033[m ";
                times[id] += time;
            }
            std::cout << std::endl;
        }
    }

    if (opt_workers >= 0) lace_stop();

    std::cout << "\033[38;5;226msummary\033[m: " << total << " games" << std::endl;
    std::cout << "\033[38;5;226msolvers\033[m:";
    for (const auto& id : solvers) {
        if (sgood[id] == total) std::cout << " \033[38;5;82m";
        else std::cout << " \033[38;5;196m";
        std::cout << id << "\033[m";
        if (sgood[id] != total) std::cout << " (" << (total-sgood[id]) << " bad)";
    }
    std::cout << std::endl;
    std::cout << "\033[38;5;226mtimes\033[m:  ";
    for (const auto& id : solvers) {
        std::cout << " \033[38;5;226m" << id << "\033[m";
        std::cout << " (" << std::fixed << std::setprecision(0) << (1000.0*times[id]) << " ms)\033[m";
    }
    std::cout << std::endl;

    return final_res;
}
