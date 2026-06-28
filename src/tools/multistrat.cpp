/*
 * Copyright 2024 Tom van Dijk
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

/**
 * multistrat: inspect and benchmark the multi-strategy solvers (fpim, fpjm).
 *
 * Usage:
 *   multistrat --demo                 print small worked examples
 *   multistrat <file-or-dir> ...      table comparing fpi/fpim and fpj/fpjm over
 *                                     the given .pg files (and directories of them)
 *
 * For each game the table reports, beyond the single strategy of fpi/fpj, how
 * many additional strategy edges the multi-strategy solvers find (im_extra /
 * jm_extra), how many vertices have a genuine choice (im_ch / jm_ch), and the
 * solving times. Each solution is verified (winners must agree and every
 * recorded strategy edge must verify).
 */

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "oink/oink.hpp"
#include "oink/game.hpp"
#include "oink/game_builder.hpp"
#include "oink/pgparser.hpp"
#include "verifier.hpp"

using namespace pg;
namespace fs = std::filesystem;

// Solve a copy of <g> with <solver> (sequential, no preprocessing so the whole
// game is handled by the solver and the multi-strategy is complete). <ms> is set
// to the milliseconds spent in run().
static Game
solve_with(const Game& g, const std::string& solver, double& ms)
{
    Game copy(g);
    std::stringstream log;
    Oink ok(copy, log);
    ok.setSolver(solver);
    ok.setRenumber();
    ok.setSolveSingle(false);
    ok.setRemoveLoops(false);
    ok.setRemoveWCWC(false);
    ok.setWorkers(-1);
    auto t0 = std::chrono::high_resolution_clock::now();
    ok.run();
    auto t1 = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return copy;
}

struct Stats {
    long strat_vertices = 0; // winner-owned solved vertices (>=1 strategy edge)
    long total_edges = 0;    // total strategy edges across those vertices
    long choice_vertices = 0;// vertices with >=2 strategy edges
    int max_set = 0;         // largest strategy set at any vertex
    long even_won = 0;       // vertices won by Even (winner sanity check)
};

static Stats
collect(const Game& g, bool multi)
{
    Stats s;
    const int* base = g.outedges();
    for (int v = 0; v < g.nodecount(); v++) {
        if (!g.isSolved(v)) continue;
        if (g.getWinner(v) == 0) s.even_won++;
        if (g.getWinner(v) != g.owner(v)) continue; // only winner-owned carry a strategy
        int cnt = 0;
        if (multi) {
            for (auto e = g.outs(v); *e != -1; e++) if (g.isStrategyEdgeIndex(e - base)) cnt++;
        } else {
            cnt = (g.getStrategy(v) != -1) ? 1 : 0;
        }
        if (cnt > 0) {
            s.strat_vertices++;
            s.total_edges += cnt;
            if (cnt > 1) s.choice_vertices++;
            if (cnt > s.max_set) s.max_set = cnt;
        }
    }
    return s;
}

static bool
verify_ok(Game& g)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier v(g, log);
    try { v.verify(true, true, true); } catch (std::runtime_error&) { return false; }
    return true;
}

static void
print_strategy(const Game& g, const char* tag)
{
    const int* base = g.outedges();
    std::cout << "  [" << tag << "]\n";
    for (int v = 0; v < g.nodecount(); v++) {
        std::cout << "    v" << v << " (prio " << g.priority(v) << ", "
                  << (g.owner(v) ? "Odd" : "Even") << ") won by "
                  << (g.getWinner(v) ? "Odd" : "Even");
        if (g.getWinner(v) == g.owner(v)) {
            std::cout << "  strategy {";
            bool first = true;
            for (auto e = g.outs(v); *e != -1; e++) {
                if (g.isStrategyEdgeIndex(e - base)) { std::cout << (first ? "" : ",") << *e; first = false; }
            }
            std::cout << "}";
        }
        std::cout << "\n";
    }
}

static void
demo()
{
    double ms;

    // K4 of even-priority Even vertices: every move stays in the (entirely Even)
    // dominion, so every vertex reports all of its outgoing edges as strategies.
    {
        GameBuilder b(4);
        for (int v = 0; v < 4; v++) {
            b.set_priority(v, 2);
            b.set_owner(v, Player::Even);
            for (int w = 0; w < 4; w++) if (w != v) b.add_edge(v, w);
        }
        Game g = b.build();
        std::cout << "Example 1: K4 of even vertices (each has 3 winning moves)\n";
        Game s = solve_with(g, "fpim", ms); print_strategy(s, "fpim");
        Game s2 = solve_with(g, "fpjm", ms); print_strategy(s2, "fpjm");
        std::cout << "\n";
    }

    // The self-loop trap: v0 (odd priority, Even-owned) has a self-loop AND an
    // escape into the even region. The self-loop is a LOSING move (odd cycle); the
    // multi-strategy must keep only the escape, even though the self-loop "stays
    // in" the dominion.
    {
        GameBuilder b(2);
        b.set_priority(0, 1); b.set_owner(0, Player::Even); b.add_edge(0, 0); b.add_edge(0, 1);
        b.set_priority(1, 2); b.set_owner(1, Player::Even); b.add_edge(1, 1);
        Game g = b.build();
        std::cout << "Example 2: self-loop trap (v0 odd self-loop excluded, only escape kept)\n";
        Game s = solve_with(g, "fpim", ms); print_strategy(s, "fpim");
        Game s2 = solve_with(g, "fpjm", ms); print_strategy(s2, "fpjm");
        std::cout << "\n";
    }

    // A mix: some vertices have a choice, some are forced.
    {
        GameBuilder b(4);
        b.set_priority(2, 2); b.set_owner(2, Player::Even); b.add_edge(2, 3); b.add_edge(2, 2);
        b.set_priority(3, 2); b.set_owner(3, Player::Even); b.add_edge(3, 2); b.add_edge(3, 3);
        b.set_priority(0, 0); b.set_owner(0, Player::Even); b.add_edge(0, 2); b.add_edge(0, 3);
        b.set_priority(1, 0); b.set_owner(1, Player::Even); b.add_edge(1, 2);
        Game g = b.build();
        std::cout << "Example 3: mixed (v0 has choice {2,3}, v1 forced {2})\n";
        Game s = solve_with(g, "fpim", ms); print_strategy(s, "fpim");
        std::cout << "\n";
    }
}

static void
benchmark(const std::vector<fs::path>& files)
{
    std::cout << std::left << std::setw(40) << "benchmark"
              << std::right << std::setw(8) << "n" << std::setw(9) << "m"
              << std::setw(8) << "strat" << std::setw(9) << "im_extra" << std::setw(8) << "im_ch"
              << std::setw(9) << "jm_extra" << std::setw(8) << "jm_ch"
              << std::setw(8) << "t_fpi" << std::setw(8) << "t_fpim"
              << std::setw(8) << "t_fpj" << std::setw(8) << "t_fpjm" << "  ok\n";

    long tot_strat=0, tot_im_extra=0, tot_im_ch=0, tot_jm_extra=0, tot_jm_ch=0;
    double tt_fpi=0, tt_fpim=0, tt_fpj=0, tt_fpjm=0;
    int n_files=0, n_with_choice=0, n_bad=0, im_max=0, jm_max=0;

    for (const auto& p : files) {
        Game g;
        try {
            std::ifstream f(p.string(), std::ios_base::binary);
            g = PGParser::parse_pgsolver(f, false);
        } catch (...) { continue; }

        double t_fpi, t_fpim, t_fpj, t_fpjm;
        Game gi  = solve_with(g, "fpi",  t_fpi);
        Game gim = solve_with(g, "fpim", t_fpim);
        Game gj  = solve_with(g, "fpj",  t_fpj);
        Game gjm = solve_with(g, "fpjm", t_fpjm);

        Stats si = collect(gi, false), sim = collect(gim, true);
        Stats sj = collect(gj, false), sjm = collect(gjm, true);

        bool ok = (si.even_won == sim.even_won) && (si.even_won == sj.even_won)
               && (si.even_won == sjm.even_won)
               && verify_ok(gi) && verify_ok(gj) && verify_ok(gim) && verify_ok(gjm);
        if (!ok) n_bad++;

        long im_extra = sim.total_edges - sim.strat_vertices;
        long jm_extra = sjm.total_edges - sjm.strat_vertices;

        std::string name = p.filename().string();
        if (name.size() > 39) name = name.substr(0, 36) + "...";
        std::cout << std::left << std::setw(40) << name
                  << std::right << std::setw(8) << g.nodecount() << std::setw(9) << g.edgecount()
                  << std::setw(8) << sim.strat_vertices << std::setw(9) << im_extra << std::setw(8) << sim.choice_vertices
                  << std::setw(9) << jm_extra << std::setw(8) << sjm.choice_vertices
                  << std::setw(8) << std::fixed << std::setprecision(1) << t_fpi << std::setw(8) << t_fpim
                  << std::setw(8) << t_fpj << std::setw(8) << t_fpjm
                  << "  " << (ok ? "y" : "N") << "\n";

        tot_strat += sim.strat_vertices; tot_im_extra += im_extra; tot_im_ch += sim.choice_vertices;
        tot_jm_extra += jm_extra; tot_jm_ch += sjm.choice_vertices;
        tt_fpi += t_fpi; tt_fpim += t_fpim; tt_fpj += t_fpj; tt_fpjm += t_fpjm;
        if (sim.choice_vertices > 0 || sjm.choice_vertices > 0) n_with_choice++;
        if (sim.max_set > im_max) im_max = sim.max_set;
        if (sjm.max_set > jm_max) jm_max = sjm.max_set;
        n_files++;
    }

    std::cout << "\n=== TOTALS over " << n_files << " benchmarks ("
              << n_with_choice << " have multi-strategy choice, " << n_bad << " bad) ===\n";
    std::cout << "strategy-vertices=" << tot_strat << "\n";
    std::cout << "fpim: extra-edges=" << tot_im_extra << " choice-vertices=" << tot_im_ch << " max-set=" << im_max << "\n";
    std::cout << "fpjm: extra-edges=" << tot_jm_extra << " choice-vertices=" << tot_jm_ch << " max-set=" << jm_max << "\n";
    std::cout << std::fixed << std::setprecision(1)
              << "time(ms): fpi=" << tt_fpi << " fpim=" << tt_fpim << " fpj=" << tt_fpj << " fpjm=" << tt_fpjm << "\n";
}

int
main(int argc, char** argv)
{
    std::vector<std::string> args(argv + 1, argv + argc);
    if (args.empty()) {
        std::cout << "Usage: " << argv[0] << " --demo | <file-or-dir>...\n";
        return 0;
    }
    if (args[0] == "--demo") { demo(); if (args.size() == 1) return 0; }

    std::vector<fs::path> files;
    for (const auto& a : args) {
        if (a == "--demo") continue;
        fs::path p(a);
        if (fs::is_directory(p)) {
            for (auto& e : fs::directory_iterator(p))
                if (e.is_regular_file() && e.path().extension() == ".pg") files.push_back(e.path());
        } else if (fs::is_regular_file(p)) {
            files.push_back(p);
        }
    }
    std::sort(files.begin(), files.end());
    if (!files.empty()) benchmark(files);
    return 0;
}
