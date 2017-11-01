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

namespace pg {

Solvers::Solvers()
{
    add("zlk", "parallel Zielonka", 1, [] (Oink* oink, Game* game, std::ostream& lgr) { return new ZLKSolver(oink, game, lgr); });
    add("pp", "priority promotion PP", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new PPSolver(oink, game, lgr); });
    add("ppp", "priority promotion PP+", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new PPPSolver(oink, game, lgr); });
    add("rr", "priority promotion RR", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new RRSolver(oink, game, lgr); });
    add("dp", "priority promotion PP+ with DP strategy", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new DPSolver(oink, game, lgr); });
    add("rrdp", "priority promotion RR with DP strategy", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new RRDPSolver(oink, game, lgr); });
    add("apt", "APT (no strategy)", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new APTSolver(oink, game, lgr); });
    add("psi", "parallel strategy improvement", 1, [] (Oink* oink, Game* game, std::ostream& lgr) { return new PSISolver(oink, game, lgr); });
    add("spm", "accelerated small progress measures", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new SPMSolver(oink, game, lgr); });
    add("tspm", "traditional small progress measures", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new TSPMSolver(oink, game, lgr); });
    add("mspm", "Maciej' modified small progress measures", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new MSPMSolver(oink, game, lgr); });
    add("qpt", "quasi-polynomial time progress measures", 0, [] (Oink* oink, Game* game, std::ostream& lgr) { return new QPTSolver(oink, game, lgr); });
}       

void
Solvers::add(std::string the_label, std::string the_desc, int the_ispar, std::function<Solver*(Oink*, Game*, std::ostream&)> the_cons)
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
