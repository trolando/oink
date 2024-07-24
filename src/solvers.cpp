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

#include "zlk.hpp"
#include "pp.hpp"
#include "ppp.hpp"
#include "rr.hpp"
#include "dp.hpp"
#include "rrdp.hpp"
#include "fpi.hpp"
#include "fpj.hpp"
#include "psi.hpp"
#include "ssi.hpp"
#include "spm.hpp"
#include "tspm.hpp"
#include "mspm.hpp"
#include "qpt.hpp"
#include "tl.hpp"
#include "rtl.hpp"
#include "npp.hpp"
#include "sspm.hpp"
#include "zlkpp.hpp"
#include "zlkq.hpp"
#include "ppq.hpp"
#include "ptl.hpp"
#include "dtl.hpp"

namespace pg {

Solvers::Solvers()
{
    _add("zlkq", "qpt Zielonka", 0, [] (Oink& oink, Game& game) { return new ZLKQSolver(&oink, &game); });
    _add("zlk", "parallel Zielonka", 1, [] (Oink& oink, Game& game) { return new ZLKSolver(&oink, &game); });
    _add("uzlk", "unoptimized Zielonka", 0, [] (Oink& oink, Game& game) { return new UnoptimizedZLKSolver(&oink, &game); });
    _add("zlkpp-std", "Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return new ZLKPPSolver(&oink, &game, ZLK_STANDARD); });
    _add("zlkpp-waw", "Warsaw quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return new ZLKPPSolver(&oink, &game, ZLK_WARSAW); });
    _add("zlkpp-liv", "Liverpool quasipolynomial Zielonka (implementation by Paweł Parys)", 0, [] (Oink& oink, Game& game) { return new ZLKPPSolver(&oink, &game, ZLK_LIVERPOOL); });
    _add("npp", "priority promotion NPP", 0, [] (Oink& oink, Game& game) { return new NPPSolver(&oink, &game); });
    _add("pp", "priority promotion PP", 0, [] (Oink& oink, Game& game) { return new PPSolver(&oink, &game); });
    _add("ppp", "priority promotion PP+", 0, [] (Oink& oink, Game& game) { return new PPPSolver(&oink, &game); });
    _add("rr", "priority promotion RR", 0, [] (Oink& oink, Game& game) { return new RRSolver(&oink, &game); });
    _add("dp", "priority promotion PP+ with DP strategy", 0, [] (Oink& oink, Game& game) { return new DPSolver(&oink, &game); });
    _add("rrdp", "priority promotion RR with DP strategy", 0, [] (Oink& oink, Game& game) { return new RRDPSolver(&oink, &game); });
    _add("ppq", "qpt Zielonka accelerated by priority promotion", 0, [] (Oink& oink, Game& game) { return new PPQSolver(&oink, &game); });
    _add("fpi", "fixpoint iteration", 1, [] (Oink& oink, Game& game) { return new FPISolver(&oink, &game); });
    _add("fpj", "fixpoint iteration with justifications", 0, [] (Oink& oink, Game& game) { return new FPJSolver(&oink, &game); });
    _add("fpjg", "greedy fixpoint iteration with justifications", 1, [] (Oink& oink, Game& game) { return new FPJGSolver(&oink, &game); });
    _add("psi", "parallel strategy improvement", 1, [] (Oink& oink, Game& game) { return new PSISolver(&oink, &game); });
    _add("ssi", "symmetric strategy improvement", 0, [] (Oink& oink, Game& game) { return new SSISolver(&oink, &game); });
    _add("spm", "accelerated small progress measures", 0, [] (Oink& oink, Game& game) { return new SPMSolver(&oink, &game); });
    _add("tspm", "traditional small progress measures", 0, [] (Oink& oink, Game& game) { return new TSPMSolver(&oink, &game); });
    _add("mspm", "Maciej' modified small progress measures", 0, [] (Oink& oink, Game& game) { return new MSPMSolver(&oink, &game); });
    _add("sspm", "succinct small progress measures", 0, [] (Oink& oink, Game& game) { return new SSPMSolver(&oink, &game); });
    _add("bsspm", "bounded succinct small progress measures", 0, [] (Oink& oink, Game& game) { return new BoundedSSPMSolver(&oink, &game); });
    _add("qpt", "quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return new QPTSolver(&oink, &game); });
    _add("bqpt", "bounded quasi-polynomial time progress measures", 0, [] (Oink& oink, Game& game) { return new BoundedQPTSolver(&oink, &game); });
    _add("ptl", "progressive tangle learning", 0, [] (Oink& oink, Game& game) { return new PTLSolver(&oink, &game); });
    _add("spptl", "single-player progressive tangle learning", 0, [] (Oink& oink, Game& game) { return new SPPTLSolver(&oink, &game); });
    _add("dtl", "distance tangle learning", 0, [] (Oink& oink, Game& game) { return new DTLSolver(&oink, &game); });
    _add("idtl", "interleaved distance tangle learning", 0, [] (Oink& oink, Game& game) { return new IDTLSolver(&oink, &game); });
    _add("rtl", "recursive tangle learning", 0, [] (Oink& oink, Game& game) { return new RTLSolver(&oink, &game); });
    _add("ortl", "one-sided recursive tangle learning", 0, [] (Oink& oink, Game& game) { return new ORTLSolver(&oink, &game); });
    _add("tl", "tangle learning", 0, [] (Oink& oink, Game& game) { return new TLSolver(&oink, &game); });
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

}
