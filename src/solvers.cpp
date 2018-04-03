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

#include "solvers.hpp"

#include "zlk.hpp"
#include "pp.hpp"
#include "ppp.hpp"
#include "rr.hpp"
#include "dp.hpp"
#include "rrdp.hpp"
#include "apt.hpp"
#include "psi.hpp"
#include "spm.hpp"
#include "tspm.hpp"
#include "mspm.hpp"
#include "qpt.hpp"
#include "tl.hpp"
#include "npp.hpp"
#include "sspm.hpp"

namespace pg {

Solvers::Solvers()
{
    add("zlk", "parallel Zielonka", 1, [] (Oink* oink, Game* game) { return new ZLKSolver(oink, game); });
    add("uzlk", "unoptimized Zielonka", 1, [] (Oink* oink, Game* game) { return new UnoptimizedZLKSolver(oink, game); });
    add("npp", "priority promotion NPP", 0, [] (Oink* oink, Game* game) { return new NPPSolver(oink, game); });
    add("pp", "priority promotion PP", 0, [] (Oink* oink, Game* game) { return new PPSolver(oink, game); });
    add("ppp", "priority promotion PP+", 0, [] (Oink* oink, Game* game) { return new PPPSolver(oink, game); });
    add("rr", "priority promotion RR", 0, [] (Oink* oink, Game* game) { return new RRSolver(oink, game); });
    add("dp", "priority promotion PP+ with DP strategy", 0, [] (Oink* oink, Game* game) { return new DPSolver(oink, game); });
    add("rrdp", "priority promotion RR with DP strategy", 0, [] (Oink* oink, Game* game) { return new RRDPSolver(oink, game); });
    add("apt", "APT (no strategy)", 0, [] (Oink* oink, Game* game) { return new APTSolver(oink, game); });
    add("psi", "parallel strategy improvement", 1, [] (Oink* oink, Game* game) { return new PSISolver(oink, game); });
    add("spm", "accelerated small progress measures", 0, [] (Oink* oink, Game* game) { return new SPMSolver(oink, game); });
    add("tspm", "traditional small progress measures", 0, [] (Oink* oink, Game* game) { return new TSPMSolver(oink, game); });
    add("mspm", "Maciej' modified small progress measures", 0, [] (Oink* oink, Game* game) { return new MSPMSolver(oink, game); });
    add("sspm", "succinct small progress measures", 0, [] (Oink* oink, Game* game) { return new SSPMSolver(oink, game); });
    add("qpt", "quasi-polynomial time progress measures", 0, [] (Oink* oink, Game* game) { return new QPTSolver(oink, game); });
    add("tl", "tangle learning", 0, [] (Oink* oink, Game* game) { return new TLSolver(oink, game); });
    add("atl", "alternating tangle learning", 0, [] (Oink* oink, Game* game) { return new ATLSolver(oink, game); });
    add("otftl", "on-the-fly tangle learning", 0, [] (Oink* oink, Game* game) { return new OTFTLSolver(oink, game); });
    add("otfatl", "on-the-fly alternating tangle learning", 0, [] (Oink* oink, Game* game) { return new OTFATLSolver(oink, game); });
}       

void
Solvers::add(std::string the_label, std::string the_desc, int the_ispar, std::function<Solver*(Oink*, Game*)> the_cons)
{
    labels.push_back(the_label);
    descriptions.push_back(the_desc);
    ispar.push_back(the_ispar);
    constructors.push_back(the_cons);
}

int
Solvers::id(std::string lbl)
{
    int id = 0;
    for (auto s : labels) {
        if (s == lbl) return id;
        id++;
    }
    return -1;
}

void
Solvers::list(std::ostream &out)
{
    out << "List of solvers:" << std::endl;
    for (unsigned i=0; i<count(); i++) {
        out << "* " << label(i) << ":\t" << desc(i) << std::endl;
    }
}

}
