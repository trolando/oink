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
// game is handled by the solver and the multi-strategy is complete). The solved
// game is written to the caller-owned <out> (so it outlives the local Oink), and
// the solution is returned (its multi-strategy references <out>). <ms> is set to
// the milliseconds spent in run().
static Solution
solve_with(const Game& g, const std::string& solver, Game& out, double& ms)
{
    out = g;
    std::stringstream log;
    Oink ok(out, log);
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
    return ok.solution(); // value type; <out> holds the matching (solved, sorted) game
}

struct Stats {
    long strat_vertices = 0; // winner-owned solved vertices (>=1 strategy edge)
    long total_edges = 0;    // total strategy edges across those vertices
    long choice_vertices = 0;// vertices with >=2 strategy edges
    int max_set = 0;         // largest strategy set at any vertex
    long even_won = 0;       // vertices won by Even (winner sanity check)
};

static Stats
collect(const Game& g, const Solution& sol)
{
    Stats s;
    std::vector<int> moves;
    for (int v = 0; v < g.nodecount(); v++) {
        if (!sol.is_solved(v)) continue;
        if (sol.winner(v) == 0) s.even_won++;
        if (sol.winner(v) != g.owner(v)) continue; // only winner-owned carry a strategy
        // strategy_targets dispatches single vs multi uniformly (the consumer API)
        moves.clear();
        strategy_targets(g, sol, v, moves);
        const int cnt = (int)moves.size();
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
verify_ok(Game& g, const Solution& sol)
{
    std::stringstream log;
    g.ensure_sorted();
    Verifier v(g, sol, log);
    try { v.verify(true, true, true); } catch (std::runtime_error&) { return false; }
    return true;
}

static void
print_strategy(const Game& g, const Solution& sol, const char* tag)
{
    std::cout << "  [" << tag << "]\n";
    std::vector<int> moves;
    for (int v = 0; v < g.nodecount(); v++) {
        std::cout << "    v" << v << " (prio " << g.priority(v) << ", "
                  << (g.owner(v) ? "Odd" : "Even") << ") won by "
                  << (sol.winner(v) ? "Odd" : "Even");
        if (sol.winner(v) == g.owner(v)) {
            moves.clear();
            strategy_targets(g, sol, v, moves);
            std::cout << "  strategy {";
            for (size_t i = 0; i < moves.size(); i++) std::cout << (i ? "," : "") << moves[i];
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
        Game g_s; Solution s = solve_with(g, "fpim", g_s, ms); print_strategy(g_s, s, "fpim");
        Game g_s2; Solution s2 = solve_with(g, "fpjm", g_s2, ms); print_strategy(g_s2, s2, "fpjm");
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
        Game g_s; Solution s = solve_with(g, "fpim", g_s, ms); print_strategy(g_s, s, "fpim");
        Game g_s2; Solution s2 = solve_with(g, "fpjm", g_s2, ms); print_strategy(g_s2, s2, "fpjm");
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
        Game g_s; Solution s = solve_with(g, "fpim", g_s, ms); print_strategy(g_s, s, "fpim");
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
        Game gi, gim, gj, gjm;
        Solution si_sol  = solve_with(g, "fpi",  gi,  t_fpi);
        Solution sim_sol = solve_with(g, "fpim", gim, t_fpim);
        Solution sj_sol  = solve_with(g, "fpj",  gj,  t_fpj);
        Solution sjm_sol = solve_with(g, "fpjm", gjm, t_fpjm);

        Stats si = collect(gi, si_sol), sim = collect(gim, sim_sol);
        Stats sj = collect(gj, sj_sol), sjm = collect(gjm, sjm_sol);

        bool ok = (si.even_won == sim.even_won) && (si.even_won == sj.even_won)
               && (si.even_won == sjm.even_won)
               && verify_ok(gi, si_sol) && verify_ok(gj, sj_sol)
               && verify_ok(gim, sim_sol) && verify_ok(gjm, sjm_sol);
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
