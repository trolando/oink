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
#include "verifier.hpp"

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

static void
resetsighandlers(void)
{
    (void)signal(SIGINT, sig_int_handler);
    (void)signal(SIGSEGV, sig_segv_handler);
    (void)signal(SIGABRT, sig_abrt_handler);
    (void)signal(SIGTERM, sig_term_handler);
}

static void
catchsig(int sig)
{
    resetsighandlers();
    out << "terminated due to signal " << sig << std::endl;
    out.flush();
    exit(-1);
    // raise(sig);
}

static void
setsighandlers(void)
{
    sig_int_handler = signal(SIGINT, catchsig);
    sig_segv_handler = signal(SIGSEGV, catchsig);
    sig_abrt_handler = signal(SIGABRT, catchsig);
    sig_term_handler = signal(SIGTERM, catchsig);
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
        ("sol", "Input solution", cxxopts::value<std::string>())
        ("o,output", "Output game or solution", cxxopts::value<std::string>())
        ("dot", "Write .dot file (before preprocessing)", cxxopts::value<std::string>())
        /* Preprocessing */
        ("inflate", "Inflate game")
        ("compress", "Compress game")
        ("no-single", "Do not solve single-parity games during preprocessing")
        ("no-loops", "Do not process self-loops during preprocessing")
        ("no-wcwc", "Do not solve winner-controlled winning cycles during preprocessing")
        /* Solving */
        ("scc", "Iteratively solve bottom SCCs")
        ("pp", "Use PP solver")
        ("ppp", "Use PP+ solver")
        ("rr", "Use RR solver")
        ("dp", "Use DP solver")
        ("rrdp", "Use RRDP solver")
        ("qpt", "Use QPT PM solver")
        ("psi", "Use Parallel SI solver")
        ("zlk", "Use improved Zielonka solver")
        ("spm", "Use SPM solver")
        ("mspm", "Use Maciej' Modified SPM solver")
        ("tspm", "Use Traditional SPM solver")
        ("s,solver", "Use given solver (--solvers for info)", cxxopts::value<std::string>())
        ("solvers", "List available solvers")
        ("w,workers", "Number of workers for parallel code", cxxopts::value<int>())
        ;
    opts.parse_positional(std::vector<std::string>({"input", "output"}));
    opts.parse(argc, argv);

    if (opts.count("help")) {
        std::cout << opts.help() << std::endl;
        return 0;
    }

    if (opts.count("solvers")) {
        listSolvers(std::cout);
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
        size_t edgecount;
        if (opts.count("input")) {
            std::string filename = opts["input"].as<std::string>();
            std::ifstream file(filename, std::ios_base::binary);
            io::filtering_istream in;
            if (boost::algorithm::ends_with(filename, ".bz2")) in.push(io::bzip2_decompressor());
            if (boost::algorithm::ends_with(filename, ".gz")) in.push(io::gzip_decompressor());
            in.push(file);
            edgecount = pg.parse_pgsolver(in);
            file.close();
        } else {
            edgecount = pg.parse_pgsolver(std::cin);
        }
        out << "parity game with " << pg.n_nodes << " nodes and " << edgecount << " edges." << std::endl;
    } catch (const char *err) {
        out << "parsing error: " << err << std::endl;
        return -1;
    }

    /**
     * STEP 2
     * Parse the (partial) solution.
     */

    try {
        if (opts.count("sol")) {
            std::ifstream file(opts["sol"].as<std::string>());
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

    if (opts.count("dot")) {
        std::ofstream file(opts["dot"].as<std::string>());
        pg.write_dot(file);
        file.close();
        out << "dot file written." << std::endl;
    }

    /**
     * STEP 4
     * Reindex the game so all nodes are in order of priority.
     * (Remember the mapping to reverse the reindex later.)
     */

    int *mapping = new int[pg.n_nodes];
    pg.reindex(mapping);
    out << "parity game reindexed" << std::endl;

    /**
     * STEP 5
     * Configure the solver.
     */

    Oink en(pg, out);
    en.setTrace(opts.count("t"));

    // preprocessing options
    if (opts.count("inflate")) en.setInflate();
    else if (opts.count("compress")) en.setCompress();
    else en.setRenumber();
    if (opts.count("no-single")) en.setSolveSingle(false);
    if (opts.count("no-loops")) en.setRemoveLoops(false);
    if (opts.count("no-wcwc")) en.setRemoveWCWC(false);

    // solver
    if (opts.count("solver")) en.setSolver(solverToId(opts["solver"].as<std::string>()));
    else en.setSolver(Solvers::TL);

    if (opts.count("zlk")) en.setSolver(Solvers::ZLK);
    if (opts.count("pp")) en.setSolver(Solvers::PP);
    if (opts.count("ppp")) en.setSolver(Solvers::PPP);
    if (opts.count("rr")) en.setSolver(Solvers::RR);
    if (opts.count("dp")) en.setSolver(Solvers::DP);
    if (opts.count("rrdp")) en.setSolver(Solvers::RRDP);
    if (opts.count("psi")) en.setSolver(Solvers::PSI);
    if (opts.count("spm")) en.setSolver(Solvers::SPM);
    if (opts.count("tspm")) en.setSolver(Solvers::TSPM);
    if (opts.count("mspm")) en.setSolver(Solvers::MSPM);
    if (opts.count("qpt")) en.setSolver(Solvers::QPT);
    if (opts.count("tl")) en.setSolver(Solvers::TL);

    // solving options
    if (opts.count("scc")) en.setBottomSCC(true);
    if (opts.count("workers")) en.setWorkers(opts["workers"].as<int>());

    /**
     * STEP 6
     * Run the solver and report the time.
     */

    double begin = wctime();
    en.run();
    double end = wctime();
    out << "solving took " << std::fixed << (end-begin) << " sec." << std::endl;

    /**
     * STEP 7
     * Verify the solution.
     */

    if (opts.count("v")) {
        try {
            out << "verifying solution..." << std::endl;
            Verifier v(&pg, out);
            double vbegin = wctime();
            v.verify(true, opts.count("qpt")==0, true);
            double vend = wctime();
            out << "solution verified (" << v.n_strategies << " strategies)." << std::endl;
            out << "verification took " << std::fixed << (vend - vbegin) << " sec." << std::endl;
        } catch (const char *err) {
            out << "verification error: " << err << std::endl;
            return -1;
        }
    }

    /**
     * STEP 8
     * Revert reindex if we need to output.
     */

    if (opts.count("output") or opts.count("p")) pg.permute(mapping);

    if (opts.count("output")) {
        // write solution to file
        if (opts.count("output")) {
            std::ofstream file(opts["output"].as<std::string>());
            pg.write_sol(file);
        }
    }

    if (opts.count("p")) {
        // print winning nodes
        bool banner = false;
        for (int i=0; i<pg.n_nodes; i++) {
            if (pg.dominion[i] == 0) {
                if (!banner) out << "won by even:";
                banner = true;
                out << " " << i; // << "(" << pg.priority[i] << ")";
            }
        }
        if (banner) out << std::endl;
        banner = false;
        for (int i=0; i<pg.n_nodes; i++) {
            if (pg.dominion[i] == 1) {
                if (!banner) out << "won by odd:";
                banner = true;
                out << " " << i; // << "(" << pg.priority[i] << ")";
            }
        }
        if (banner) out << std::endl;
    }

    delete[] mapping;

    resetsighandlers();
    return 0;
}
