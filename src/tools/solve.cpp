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

#include <csignal>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sys/time.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "cxxopts.hpp" 
#include "game.hpp"
#include "oink.hpp"
#include "solvers.hpp"
#include "verifier.hpp"
#include "tools/getrss.h"

using namespace pg;

/*------------------------------------------------------------------------*/

static double
wctime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec + 1E-6 * time.tv_usec;
}

static double t_start;

/*------------------------------------------------------------------------*/

// timestamp_filter adds a timestamp at the beginning of every line.
namespace io = boost::iostreams;
class timestamp_filter : public io::output_filter
{
public:
    timestamp_filter() {}
    struct category : io::output_filter::category, io::flushable_tag { };

    template<typename Sink> bool put(Sink& snk, char c);
    template<typename Device> void close(Device&);
    template<typename Sink> bool flush(Sink& snk);

private:
    bool is_start = true;
    char sz[16];
    const char* pos = NULL;
    const char* end = NULL;
};

template<typename Sink>
bool timestamp_filter::put(Sink& dest, char c)
{
    if (is_start) {
        if (c == '\n') return true;  // ignore consecutive endl
        is_start = false;
        pos = sz;
        end = sz + snprintf(sz, 16, "[% 8.2f] ", wctime() - t_start);
    }

    while (pos != end) {
        if (!io::put(dest, *pos)) return false;
        pos++;
    }

    if (!io::put(dest, c)) return false;
    if (c == '\n') is_start = true;

    return true;
}

template<typename Sink>
bool timestamp_filter::flush(Sink& dest)
{
    while (pos != end) {
        if (!io::put(dest, *pos)) return false;
        pos++;
    }

    return io::flush(dest);
}

template<typename Device>
void timestamp_filter::close(Device&)
{
    is_start = true;
    pos = end = NULL;
}

// global variable so signal handlers can work with it
io::filtering_ostream out;

/*------------------------------------------------------------------------*/

static void (*sig_int_handler)(int);
static void (*sig_segv_handler)(int);
static void (*sig_abrt_handler)(int);
static void (*sig_term_handler)(int);
static void (*sig_alrm_handler)(int);

static void
resetsighandlers(void)
{
    (void)signal(SIGINT, sig_int_handler);
    (void)signal(SIGSEGV, sig_segv_handler);
    (void)signal(SIGABRT, sig_abrt_handler);
    (void)signal(SIGTERM, sig_term_handler);
    (void)signal(SIGALRM, sig_alrm_handler);
}

static void
catchsig(int sig)
{
    // note: this can actually deadlock because we are writing to stdout...

    if (sig == SIGALRM) {
        resetsighandlers();
        out << std::endl << "terminated due to timeout" << std::endl;
        out.flush();
        exit(-1);
    } else if (sig == SIGINT) {
        // CTRL-C
        resetsighandlers();
        out << std::endl << "received INT signal" << std::endl;
        out.flush();
        exit(-SIGINT);
    } else if (sig == SIGABRT) {
        // really bad
        resetsighandlers();
        out << std::endl << "terminated due to ABORT signal" << std::endl;
        out.flush();
        exit(-SIGABRT);
        raise(sig);
    } else {
        resetsighandlers();
        out << std::endl << "terminated due to signal " << sig << std::endl;
        out.flush();
        exit(-sig);
    }
}

static void
setsighandlers(void)
{
    sig_int_handler = signal(SIGINT, catchsig);
    sig_segv_handler = signal(SIGSEGV, catchsig);
    sig_abrt_handler = signal(SIGABRT, catchsig);
    sig_term_handler = signal(SIGTERM, catchsig);
    sig_alrm_handler = signal(SIGALRM, catchsig);
}

/*------------------------------------------------------------------------*/

static char*
to_h(double size, char *buf)
{
    const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    int i = 0;
    for (;size>1024;size/=1024) i++;
    sprintf(buf, "%.*f %s", i, size, units[i]);
    return buf;
}

/*------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    t_start = wctime();

    setsighandlers();

    cxxopts::Options opts(argv[0], "Parity game solver");
    opts.add_options()
        ("help", "Print help")
        ("t,trace", "Generate trace (with increasing verbosity)")
        ("v,verify", "Verify solution")
        ("p,print", "Print solution to stdout")
        ("i,input", "Input parity game", cxxopts::value<std::string>())
        ("sol", "Input (partial) solution", cxxopts::value<std::string>())
        ("o,output", "Output game or solution", cxxopts::value<std::string>())
        ("dot", "Write .dot file (before preprocessing)", cxxopts::value<std::string>())
        /* Preprocessing */
        ("inflate", "Inflate game")
        ("compress", "Compress game")
        ("no-single", "Do not solve single-parity games during preprocessing")
        ("no-loops", "Do not remove self-loops during preprocessing (default behavior)")
        ("no-wcwc", "Do not solve winner-controlled winning cycles during preprocessing")
        ("no", "Do not touch the game at all")
        /* Solving */
        ("scc", "Iteratively solve bottom SCCs")
        ("s,solver", "Use given solver (--solvers for info)", cxxopts::value<std::string>())
        ("solvers", "List available solvers")
        ("w,workers", "Number of workers for parallel code", cxxopts::value<int>())
        ("z,timeout", "Number of seconds for timeout", cxxopts::value<int>())
        ;

    /* Add solvers */
    Solvers solvers;
    for (unsigned id=0; id<solvers.count(); id++) {
        opts.add_options()(solvers.label(id), solvers.desc(id));
    }

    /* Parse command line */
    opts.parse_positional(std::vector<std::string>({"input", "output"}));
    auto options = opts.parse(argc, argv);

    if (options.count("help")) {
        std::cout << opts.help() << std::endl;
        return 0;
    }

    if (options.count("solvers")) {
        solvers.list(std::cout);
        return 0;
    }

    /* Setup timestamp filter */

    out.push(timestamp_filter());
    out.push(std::cout);

    /**
     * STEP 1
     * Read the game that must be solved.
     * (Supports bz2 and gz compression.)
     */

    Game pg;

    try {
        if (options.count("input")) {
            std::string filename = options["input"].as<std::string>();
            io::filtering_istream in;
            if (boost::algorithm::ends_with(filename, ".bz2")) in.push(io::bzip2_decompressor());
            if (boost::algorithm::ends_with(filename, ".gz")) in.push(io::gzip_decompressor());
            std::ifstream file(filename, std::ios_base::binary);
            in.push(file);
            pg.parse_pgsolver(in, options.count("no-loops") == 0 and options.count("no") == 0);
            file.close();
        } else {
            pg.parse_pgsolver(std::cin, options.count("no-loops") == 0);
        }
        out << "parity game with " << pg.nodecount() << " nodes and " << pg.edgecount() << " edges." << std::endl;
    } catch (const char *err) {
        out << "parsing error: " << err << std::endl;
        return -1;
    }

    /**
     * STEP 2
     * Parse the (partial) solution.
     */

    try {
        if (options.count("sol")) {
            std::ifstream file(options["sol"].as<std::string>());
            pg.parse_solution(file);
            file.close();
            out << "solution parsed." << std::endl;
        }
    } catch (const char *err) {
        out << "parsing error: " << err << std::endl;
        return -1;
    }

    /**
     * STEP 3
     * If requested, write .dot file
     */

    if (options.count("dot")) {
        std::ofstream file(options["dot"].as<std::string>());
        pg.write_dot(file);
        file.close();
        out << "dot file written." << std::endl;
    }

    /**
     * STEP 4
     * Reindex the game so all nodes are in order of priority.
     * (Remember the mapping to reverse the reindex later.)
     */

    int *mapping = new int[pg.nodecount()];
    pg.sort(mapping);
    out << "parity game reindexed" << std::endl;

    /**
     * STEP 5
     * Configure the solver.
     */

    Oink en(pg, out);
    en.setTrace(options.count("t"));

    // preprocessing options
    bool no = options.count("no");
    if (options.count("inflate")) en.setInflate();
    else if (options.count("compress")) en.setCompress();
    else if (!no) en.setRenumber();
    if (no or options.count("no-single")) en.setSolveSingle(false);
    if (no or options.count("no-loops")) en.setRemoveLoops(false);
    if (no or options.count("no-wcwc")) en.setRemoveWCWC(false);

    // solver
    if (options.count("solver")) {
        en.setSolver(solvers.id(options["solver"].as<std::string>()));
    } else {
        en.setSolver("atl"); // default solver :)
        for (unsigned id=0; id<solvers.count(); id++) {
            if (options.count(solvers.label(id))) en.setSolver(id);
        }
    }

    // solving options
    if (options.count("scc")) en.setBottomSCC(true);
    if (options.count("workers")) en.setWorkers(options["workers"].as<int>());

    /**
     * STEP 6
     * Run the solver and report the time.
     */

    if (options.count("timeout")) alarm(options["timeout"].as<int>());

    try {
        double begin = wctime();
        en.run();
        double end = wctime();
        out << "total solving time: " << std::fixed << (end-begin) << " sec." << std::endl;
    } catch (pg::Error &err) {
        out << "solving error: " << err.what() << std::endl;
        return -1;
    }

    /**
     * STEP 7
     * Verify the solution.
     */

    if (options.count("v")) {
        try {
            out << "verifying solution..." << std::endl;
            Verifier v(&pg, out);
            double vbegin = wctime();
            v.verify(true, true, true);
            double vend = wctime();
            out << "solution verified (" << v.n_strategies << " strategies)." << std::endl;
            out << "verification took " << std::fixed << (vend - vbegin) << " sec." << std::endl;
        } catch (const char *err) {
            out << "verification error: " << err << std::endl;
            return -1;
        }
    }

    char buf[32];
    to_h(getCurrentRSS(), buf);
    out << "current memory usage: " << buf << std::endl;
    to_h(getPeakRSS(), buf);
    out << "peak memory usage: " << buf << std::endl;

    /**
     * STEP 8
     * Revert reindex if we need to output.
     */

    if (options.count("output") or options.count("p")) pg.permute(mapping);

    if (options.count("output")) {
        // write solution to file
        if (options.count("output")) {
            std::ofstream file(options["output"].as<std::string>());
            pg.write_sol(file);
        }
    }

    if (options.count("p")) {
        // print winning nodes
        bool banner = false;
        for (int i=0; i<pg.nodecount(); i++) {
            if (pg.solved[i] and pg.winner[i] == 0) {
                if (!banner) out << "won by even:";
                banner = true;
                // out << " " << i; // << "(" << pg.priority(i) << ")";
                out << " " << pg.label_vertex(i);
            }
        }
        if (banner) out << std::endl;
        banner = false;
        for (int i=0; i<pg.nodecount(); i++) {
            if (pg.solved[i] and pg.winner[i] == 1) {
                if (!banner) out << "won by odd:";
                banner = true;
                out << " " << pg.label_vertex(i);
            }
        }
        if (banner) out << std::endl;
    }

    delete[] mapping;

    resetsighandlers();
    return 0;
}
