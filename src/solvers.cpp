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

#include "oink/solvers.hpp"
#include <algorithm>
#include <initializer_list>

#include "solvers/zlk.hpp"
#include "solvers/zlkj.hpp"
#include "solvers/pp.hpp"
#include "solvers/ppj.hpp"
#include "solvers/ppp.hpp"
#include "solvers/rr.hpp"
#include "solvers/dp.hpp"
#include "solvers/rrdp.hpp"
#include "solvers/fpi.hpp"
#include "solvers/fpim.hpp"
#include "solvers/fpj.hpp"
#include "solvers/fpjm.hpp"
#include "solvers/psi.hpp"
#include "solvers/ssi.hpp"
#include "solvers/spm.hpp"
#include "solvers/tspm.hpp"
#include "solvers/mspm.hpp"
#include "solvers/qpt.hpp"
#include "solvers/tl.hpp"
#include "solvers/tlpar.hpp"
#include "solvers/tlj.hpp"
#include "solvers/rtl.hpp"
#include "solvers/npp.hpp"
#include "solvers/sspm.hpp"
#include "solvers/zlkpp.hpp"
#include "solvers/zlkq.hpp"
#include "solvers/zlkw.hpp"
#include "solvers/ppq.hpp"
#include "solvers/dtl.hpp"
#include "solvers/tlq.hpp"
#include "solvers/gpm.hpp"
#include "solvers/pmtl.hpp"

namespace pg {

Solvers::Solvers()
{
    _add("zlkq", "qpt Zielonka", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKQSolver>(oink, game); });
    _add("zlkw", "Warsaw quasi-polynomial Zielonka", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKWSolver>(oink, game); });
    _add("zlk", "parallel Zielonka", 1, [] (Oink& oink, Game& game) { return std::make_unique<ZLKSolver>(oink, game); });
    _add("zlkj", "Zielonka with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKJSolver>(oink, game); });
    _add("uzlk", "unoptimized Zielonka", 0, [] (Oink& oink, Game& game) { return std::make_unique<UnoptimizedZLKSolver>(oink, game); }, true);
    _add("zlkpp-std", "Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_STANDARD); }, true);
    _add("zlkpp-waw", "Warsaw quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_WARSAW); }, true);
    _add("zlkpp-liv", "Liverpool quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_LIVERPOOL); }, true);
    _add("npp", "priority promotion NPP", 0, [] (Oink& oink, Game& game) { return std::make_unique<NPPSolver>(oink, game); });
    _add("pp", "priority promotion PP", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPSolver>(oink, game); });
    _add("ppj", "priority promotion with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPJSolver>(oink, game); });
    _add("ppp", "priority promotion PP+", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPPSolver>(oink, game); }, true);
    _add("rr", "priority promotion RR", 0, [] (Oink& oink, Game& game) { return std::make_unique<RRSolver>(oink, game); }, true);
    _add("dp", "priority promotion PP+ with DP strategy", 0, [] (Oink& oink, Game& game) { return std::make_unique<DPSolver>(oink, game); }, true);
    _add("rrdp", "priority promotion RR with DP strategy", 0, [] (Oink& oink, Game& game) { return std::make_unique<RRDPSolver>(oink, game); }, true);
    _add("ppq", "qpt Zielonka accelerated by priority promotion", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPQSolver>(oink, game); });
    _add("fpi", "fixpoint iteration", 1, [] (Oink& oink, Game& game) { return std::make_unique<FPISolver>(oink, game); });
    _add("fpim", "fixpoint iteration (multi-strategy)", 0, [] (Oink& oink, Game& game) { return std::make_unique<FPIMSolver>(oink, game); });
    _add("fpj", "fixpoint iteration with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<FPJSolver>(oink, game); });
    _add("fpjm", "fixpoint iteration with justifications (multi-strategy)", 0, [] (Oink& oink, Game& game) { return std::make_unique<FPJMSolver>(oink, game); });
    _add("fpjg", "greedy fixpoint iteration with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<FPJGSolver>(oink, game); });
    _add("psi", "parallel strategy improvement", 1, [] (Oink& oink, Game& game) { return std::make_unique<PSISolver>(oink, game); });
    _add("ssi", "symmetric strategy improvement", 0, [] (Oink& oink, Game& game) { return std::make_unique<SSISolver>(oink, game); }, true);
    _add("spm", "accelerated small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<SPMSolver>(oink, game); }, true);
    _add("tspm", "traditional small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<TSPMSolver>(oink, game); }, true);
    _add("mspm", "small progress measures with both strategies from one measure (Gazda, Willemse)", 0, [] (Oink& oink, Game& game) { return std::make_unique<MSPMSolver>(oink, game); }, true);
    _add("sspm", "succinct small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<SSPMSolver>(oink, game); }, true);
    _add("bsspm", "bounded succinct small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<BoundedSSPMSolver>(oink, game); }, true);
    _add("qpt", "quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<QPTSolver>(oink, game); }, true);
    _add("bqpt", "bounded quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<BoundedQPTSolver>(oink, game); }, true);
    _add("dtl", "distance tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<DTLSolver>(oink, game); }, true);
    _add("idtl", "interleaved distance tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<IDTLSolver>(oink, game); }, true);
    _add("rtl", "recursive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<RTLSolver>(oink, game); });
    _add("ortl", "one-sided recursive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<ORTLSolver>(oink, game); }, true);
    _add("tl", "tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<TLSolver>(oink, game); });
    _add("tlpar", "parallel tangle learning", 1, [] (Oink& oink, Game& game) { return std::make_unique<TLParSolver>(oink, game); });
    _add("tlj", "tangle learning with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<TLJSolver>(oink, game); });
    _add("tlq", "qpt recursive with tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<TLQSolver>(oink, game); });
    _add("pmtl", "progress measures and tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game); }, true);
    _add("pmtl-small", "progress measures and tangle learning, small measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Small); });
    _add("pmtl-succinct", "progress measures and tangle learning, succinct measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Succinct); });
    _add("pmtl-ordered", "progress measures and tangle learning, ordered measures (as pmtl)", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Ordered); });
    _add("pmtl-strahler", "progress measures and tangle learning, Strahler measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Strahler); });
    _add("gpm", "generic progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<GPMSolver>(oink, game); }, true);
    _add("gpm-small", "generic progress measures, small measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<GPMSolver>(oink, game, MeasureKind::Small); });
    _add("gpm-succinct", "generic progress measures, succinct measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<GPMSolver>(oink, game, MeasureKind::Succinct); });
    _add("gpm-ordered", "generic progress measures, ordered measures (as gpm)", 0, [] (Oink& oink, Game& game) { return std::make_unique<GPMSolver>(oink, game, MeasureKind::Ordered); });
    _add("gpm-strahler", "generic progress measures, Strahler measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<GPMSolver>(oink, game, MeasureKind::Strahler); });

    _add("pma", "progress measures with attractors", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Ordered, false); }, true);
    _add("pma-small", "progress measures with attractors, small measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Small, false); }, true);
    _add("pma-succinct", "progress measures with attractors, succinct measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Succinct, false); }, true);
    _add("pma-ordered", "progress measures with attractors, ordered measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Ordered, false); }, true);
    _add("pma-strahler", "progress measures with attractors, Strahler measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<PMTLSolver>(oink, game, MeasureKind::Strahler, false); }, true);

    auto group = [&](Family family, std::initializer_list<const char*> ids) {
        int order = 0;
        for (const auto* id : ids) {
            auto& info = solvers.at(id);
            info.family = family;
            info.familyOrder = order++;
        }
    };
    group(Family::TangleLearning, {"tl", "tlpar", "tlj", "tlq", "rtl", "ortl", "dtl", "idtl",
          "pmtl", "pmtl-small", "pmtl-succinct", "pmtl-ordered", "pmtl-strahler"});
    group(Family::ProgressMeasures, {"spm", "tspm", "mspm", "sspm", "bsspm", "qpt", "bqpt",
          "gpm", "gpm-small", "gpm-succinct", "gpm-ordered", "gpm-strahler",
          "pma", "pma-small", "pma-succinct", "pma-ordered", "pma-strahler"});
    group(Family::FixpointIteration, {"fpi", "fpim", "fpj", "fpjg", "fpjm"});
    group(Family::PriorityPromotion, {"npp", "pp", "ppj", "ppp", "rr", "dp", "rrdp", "ppq"});
    group(Family::StrategyImprovement, {"psi", "ssi"});
    group(Family::Zielonka, {"zlk", "zlkj", "uzlk", "zlkq", "zlkw", "zlkpp-std", "zlkpp-waw", "zlkpp-liv"});

    solvers["mspm"].optionsHelp = "even | odd (default: smaller measure domain)";
    solvers["psi"].optionsHelp = "parmin=N,chunk=N (positive integers; defaults: 2048,256)";
    solvers["tlpar"].optionsHelp = "parmin=N,chunk=N (positive integers; defaults: 2048,256)";
    for (const char* id : {"gpm", "gpm-small", "gpm-succinct", "gpm-ordered", "gpm-strahler"})
        solvers[id].optionsHelp = "small, succinct, ordered, strahler, no-shortcuts (comma separated)";
    for (const char* id : {"pmtl", "pmtl-small", "pmtl-succinct", "pmtl-ordered", "pmtl-strahler"})
        solvers[id].optionsHelp = "small, succinct, ordered, strahler, up/noup, tangles/notangles, no-shortcuts (comma separated)";

    for (const char* id : {"pma", "pma-small", "pma-succinct", "pma-ordered", "pma-strahler"})
        solvers[id].optionsHelp = "small, succinct, ordered, strahler, up/noup, tangles/notangles (default: notangles), no-shortcuts (comma separated)";

    // No incoming edges required up front; sequential FPI may build them lazily.
    for (const char* id : {"fpi", "fpim", "psi"}) solvers[id].needsInEdges = false;

    solvers["tlpar"].sequentialFallback = "tl";
}

std::string
Solvers::familyName(Family family)
{
    switch (family) {
        case Family::TangleLearning: return "Tangle learning";
        case Family::ProgressMeasures: return "Progress measures";
        case Family::FixpointIteration: return "Fixpoint iteration";
        case Family::PriorityPromotion: return "Priority promotion";
        case Family::StrategyImprovement: return "Strategy improvement";
        case Family::Zielonka: return "Zielonka recursion";
        case Family::Other: return "Other / external";
    }
    return "Other / external";
}

std::vector<std::string>
Solvers::getSolverIDsByFamily()
{
    std::vector<std::string> ids;
    for (const auto& entry : instance().solvers) ids.push_back(entry.first);
    const auto& registry = instance().solvers;
    std::sort(ids.begin(), ids.end(), [&](const std::string& a, const std::string& b) {
        const auto& x = registry.at(a);
        const auto& y = registry.at(b);
        if (x.family != y.family) return x.family < y.family;
        if (x.familyOrder != y.familyOrder) return x.familyOrder < y.familyOrder;
        return a < b;
    });
    return ids;
}

void
Solvers::list(std::ostream &out)
{
    out << "List of solvers:" << std::endl;
    for (const auto& entry : instance().solvers) {
        const auto& label = entry.first;
        const auto& info = entry.second;
        out << "* " << label << ":\t" << info.description << std::endl;
        if (!info.optionsHelp.empty()) out << "    options: " << info.optionsHelp << std::endl;
    }
}

/**
 * Construct solver with the given parameters
 */
std::unique_ptr<Solver>
Solvers::construct(const std::string& id, Oink& oink, Game& game) 
{
    return instance().solvers[id].constructor(oink, game); 
}

bool
Solvers::validateOptions(const std::string& id, const std::string& options, std::ostream& errors)
{
    if (options.empty()) return true;
    if (!hasOptions(id)) {
        errors << id << " does not accept configuration options" << std::endl;
        return false;
    }
    Game game(1);
    game.vec_init();
    game.init_vertex(0, 0, 0);
    game.vec_add_edge(0, 0);
    game.vec_finish();
    game.ensure_sorted();
    Oink controller(game, errors, needsInEdges(id));
    auto solver = construct(id, controller, game);
    std::string text = options;
    return solver->parseOptions(text);
}

void
Solvers::add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor)
{
    add(id, description, isParallel, constructor, "");
}

void
Solvers::add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor, const std::string& optionsHelp, bool experimental)
{
    instance().solvers[id] = {description, isParallel, constructor, experimental};
    instance().solvers[id].optionsHelp = optionsHelp;
}

}
