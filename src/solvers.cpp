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

#include "solvers/zlk.hpp"
#include "solvers/pp.hpp"
#include "solvers/ppp.hpp"
#include "solvers/rr.hpp"
#include "solvers/dp.hpp"
#include "solvers/rrdp.hpp"
#include "solvers/fpi.hpp"
#include "solvers/fpj.hpp"
#include "solvers/psi.hpp"
#include "solvers/ssi.hpp"
#include "solvers/spm.hpp"
#include "solvers/tspm.hpp"
#include "solvers/mspm.hpp"
#include "solvers/qpt.hpp"
#include "solvers/tl.hpp"
#include "solvers/rtl.hpp"
#include "solvers/npp.hpp"
#include "solvers/sspm.hpp"
#include "solvers/zlkpp.hpp"
#include "solvers/zlkq.hpp"
#include "solvers/ppq.hpp"
#include "solvers/ptl.hpp"
#include "solvers/dtl.hpp"
#include "solvers/tlq.hpp"

namespace pg {

Solvers::Solvers()
{
    _add("zlkq", "qpt Zielonka", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKQSolver>(oink, game); });
    _add("zlk", "parallel Zielonka", 1, [] (Oink& oink, Game& game) { return std::make_unique<ZLKSolver>(oink, game); });
    _add("uzlk", "unoptimized Zielonka", 0, [] (Oink& oink, Game& game) { return std::make_unique<UnoptimizedZLKSolver>(oink, game); });
    _add("zlkpp-std", "Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_STANDARD); });
    _add("zlkpp-waw", "Warsaw quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_WARSAW); });
    _add("zlkpp-liv", "Liverpool quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return std::make_unique<ZLKPPSolver>(oink, game, ZLK_LIVERPOOL); });
    _add("npp", "priority promotion NPP", 0, [] (Oink& oink, Game& game) { return std::make_unique<NPPSolver>(oink, game); });
    _add("pp", "priority promotion PP", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPSolver>(oink, game); });
    _add("ppp", "priority promotion PP+", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPPSolver>(oink, game); });
    _add("rr", "priority promotion RR", 0, [] (Oink& oink, Game& game) { return std::make_unique<RRSolver>(oink, game); });
    _add("dp", "priority promotion PP+ with DP strategy", 0, [] (Oink& oink, Game& game) { return std::make_unique<DPSolver>(oink, game); });
    _add("rrdp", "priority promotion RR with DP strategy", 0, [] (Oink& oink, Game& game) { return std::make_unique<RRDPSolver>(oink, game); });
    _add("ppq", "qpt Zielonka accelerated by priority promotion", 0, [] (Oink& oink, Game& game) { return std::make_unique<PPQSolver>(oink, game); });
    _add("fpi", "fixpoint iteration", 1, [] (Oink& oink, Game& game) { return std::make_unique<FPISolver>(oink, game); });
    _add("fpj", "fixpoint iteration with justifications", 0, [] (Oink& oink, Game& game) { return std::make_unique<FPJSolver>(oink, game); });
    _add("fpjg", "greedy fixpoint iteration with justifications", 1, [] (Oink& oink, Game& game) { return std::make_unique<FPJGSolver>(oink, game); });
    _add("psi", "parallel strategy improvement", 1, [] (Oink& oink, Game& game) { return std::make_unique<PSISolver>(oink, game); });
    _add("ssi", "symmetric strategy improvement", 0, [] (Oink& oink, Game& game) { return std::make_unique<SSISolver>(oink, game); });
    _add("spm", "accelerated small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<SPMSolver>(oink, game); });
    _add("tspm", "traditional small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<TSPMSolver>(oink, game); });
    _add("mspm", "Maciej' modified small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<MSPMSolver>(oink, game); });
    _add("sspm", "succinct small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<SSPMSolver>(oink, game); });
    _add("bsspm", "bounded succinct small progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<BoundedSSPMSolver>(oink, game); });
    _add("qpt", "quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<QPTSolver>(oink, game); });
    _add("bqpt", "bounded quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return std::make_unique<BoundedQPTSolver>(oink, game); });
    _add("ptl", "progressive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<PTLSolver>(oink, game); });
    _add("spptl", "single-player progressive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<SPPTLSolver>(oink, game); });
    _add("dtl", "distance tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<DTLSolver>(oink, game); });
    _add("idtl", "interleaved distance tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<IDTLSolver>(oink, game); });
    _add("rtl", "recursive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<RTLSolver>(oink, game); });
    _add("ortl", "one-sided recursive tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<ORTLSolver>(oink, game); });
    _add("tl", "tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<TLSolver>(oink, game); });
    _add("tlq", "qpt recursive with tangle learning", 0, [] (Oink& oink, Game& game) { return std::make_unique<TLQSolver>(oink, game); });
}       

void
Solvers::list(std::ostream &out)
{
    out << "List of solvers:" << std::endl;
    for (const auto& entry : instance().solvers) {
        const auto& label = entry.first;
        const auto& info = entry.second;
        out << "* " << label << ":\t" << info.description << std::endl;
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

void
Solvers::add(const std::string& id, const std::string& description, bool isParallel, const SolverConstructor& constructor) 
{
    instance().solvers[id] = {description, isParallel, constructor};
}

}
