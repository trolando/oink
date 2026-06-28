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
#include <unistd.h>

#include "cxxopts.hpp"
#include "oink/oink.hpp"
#include "oink/solvers.hpp"
#include "oink/pgparser.hpp"
#include "verifier.hpp"
#include "oink/io.hpp"
#include "oink/version.hpp"
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

// timestamp_streambuf prefixes a timestamp to the beginning of every line
// before forwarding to an underlying ostream.
class timestamp_streambuf : public std::streambuf
{
public:
    explicit timestamp_streambuf(std::ostream& sink) : sink(sink) {}

protected:
    int_type overflow(int_type ch) override
    {
        if (traits_type::eq_int_type(ch, traits_type::eof())) return ch;
        const char c = traits_type::to_char_type(ch);
        if (is_start) {
            if (c == '\n') return ch;  // ignore consecutive endl
            is_start = false;
            char sz[16];
            int n = snprintf(sz, sizeof(sz), "[% 8.2f] ", wctime() - t_start);
            sink.write(sz, n);
        }
        sink.put(c);
        if (c == '\n') is_start = true;
        return ch;
    }

    int sync() override { sink.flush(); return sink ? 0 : -1; }

private:
    std::ostream& sink;
    bool is_start = true;
};

// An ostream that prefixes every line with a timestamp.
class timestamp_ostream : public std::ostream
{
public:
    explicit timestamp_ostream(std::ostream& sink)
        : std::ostream(nullptr), buf(sink) { rdbuf(&buf); }

private:
    timestamp_streambuf buf;
};

// global variable so signal handlers can work with it
timestamp_ostream out(std::cout);

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

std::string to_h(double size) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
    int i = 0;
    while (size > 1024 && i < 8) {
        size /= 1024;
        i++;
    }
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(i) << size << " " << units[i];
    return oss.str();
}

/*------------------------------------------------------------------------*/

int main(int argc, char **argv)
{
    t_start = wctime();

    setsighandlers();

    cxxopts::Options opts(argv[0], "Parity game solver");
    opts.add_options()
        ("help", "Print help")
        ("version", "Print version and exit")
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
        ("c,configure", "Additional configuration options for the solver", cxxopts::value<std::string>())
        ("w,workers", "Number of workers for parallel code", cxxopts::value<int>())
        ("z,timeout", "Number of seconds for timeout", cxxopts::value<int>())
        ;

    /* Add solvers */
    for (const auto & id : Solvers::getSolverIDs()) {
        opts.add_options()(id, Solvers::desc(id));
    }

    /* Parse command line */
    opts.parse_positional(std::vector<std::string>({"input", "output"}));
    auto options = opts.parse(argc, argv);

    if (options.count("help")) {
        std::cout << opts.help() << std::endl;
        return 0;
    }

    if (options.count("version")) {
        std::cout << "oink " << pg::version() << std::endl;
        std::cout << "compressed input: " << supported_input_formats() << std::endl;
        return 0;
    }

    if (options.count("solvers")) {
        Solvers::list(std::cout);
        return 0;
    }

    /* Setup timestamp filter */


    /**
     * STEP 1
     * Read the game that must be solved.
     * (Transparently decompresses .gz/.bz2/.xz when built with the
     *  corresponding compression library.)
     */

    Game pg;

    try {
        if (options.count("input")) {
            std::string filename = options["input"].as<std::string>();
            // time opening (and, for compressed input, decompressing) the file
            auto open_begin = wctime();
            auto in = open_input(filename);
            auto open_end = wctime();
            out << "reading input took " << std::fixed << (open_end-open_begin) << " sec." << std::endl;
            // time parsing
            auto begin = wctime();
            pg = PGParser::parse_pgsolver_renumber(*in, options.count("no-loops") == 0 and options.count("no") == 0);
            auto end = wctime();
            out << "parsing took " << std::fixed << (end-begin) << " sec." << std::endl;
        } else {
            pg = PGParser::parse_pgsolver_renumber(std::cin, options.count("no-loops") == 0 and options.count("no") == 0);
        }
        out << "parity game with " << pg.nodecount() << " nodes and " << pg.edgecount() << " edges." << std::endl;
    } catch (std::runtime_error &err) {
        out << "parsing error: " << err.what() << std::endl;
        return -1;
    }

    /**
     * STEP 2
     * Parse the (partial) solution.
     */

    // A partial solution (--sol) is parsed into Oink's solution after Oink is
    // constructed (see below), since the solution now lives with the solver.

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
        en.setSolver(options["solver"].as<std::string>());
    } else {
        en.setSolver("tl"); // default solver
        for (const auto& id : Solvers::getSolverIDs()) {
            if (options.count(id)) en.setSolver(id);
        }
    }

    // solving options
    if (options.count("scc")) en.setBottomSCC(true);
    if (options.count("workers")) en.setWorkers(options["workers"].as<int>());

    if (options.count("configure")) {
        en.setSolverOptions(options["configure"].as<std::string>());
    }

    // Parse a partial solution (--sol) into Oink's solution (in the now-sorted
    // game's numbering). Oink picks up pre-solved vertices when run() is called.
    if (options.count("sol")) {
        try {
            std::ifstream file(options["sol"].as<std::string>());
            parse_solution(pg, en.solution(), file);
            out << "solution parsed." << std::endl;
        } catch (std::runtime_error &err) {
            out << "parsing error: " << err.what() << std::endl;
            return -1;
        }
    }

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
            pg.ensure_sorted(); // Verifier requires a game sorted by priority
            Verifier v(pg, en.solution(), out);
            double vbegin = wctime();
            v.verify(true, true, true);
            double vend = wctime();
            out << "solution verified (" << v.numberOfStrategies() << " strategies)." << std::endl;
            out << "verification took " << std::fixed << (vend - vbegin) << " sec." << std::endl;
        } catch (std::runtime_error &err) {
            out << "verification error: " << err.what() << std::endl;
            return -1;
        }
    }

    out << "current memory usage: " << to_h(getCurrentRSS()) << std::endl;
    out << "peak memory usage: " << to_h(getPeakRSS()) << std::endl;

    /**
     * STEP 8
     * Revert reindex if we need to output.
     */

    // Take a copy of the (Oink-owned) solution; revert the reindex on both the
    // game and the solution together so output is in the input numbering.
    pg::Solution sol = en.solution();
    if (options.count("output") or options.count("p")) {
        sol.permute(mapping); // does not consume mapping
        pg.permute(mapping);  // consumes mapping
    }

    if (options.count("output")) {
        // write solution to file
        if (options.count("output")) {
            std::ofstream file(options["output"].as<std::string>());
            write_solution(pg, sol, file);
        }
    }

    if (options.count("p")) {
        // print winning nodes
        bool banner = false;
        for (int i=0; i<pg.nodecount(); i++) {
            if (sol.is_solved(i) and sol.winner(i) == 0) {
                if (!banner) out << "won by even:";
                banner = true;
                // out << " " << i; // << "(" << pg.priority(i) << ")";
                out << " " << pg.label_vertex(i);
            }
        }
        if (banner) out << std::endl;
        banner = false;
        for (int i=0; i<pg.nodecount(); i++) {
            if (sol.is_solved(i) and sol.winner(i) == 1) {
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
