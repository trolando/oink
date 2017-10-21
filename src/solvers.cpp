#include "solvers.hpp"

#include "pp.hpp"
#include "ppp.hpp"
#include "rr.hpp"
#include "dp.hpp"
#include "rrdp.hpp"
#include "qpt.hpp"
#include "psi.hpp"
#include "zlk.hpp"
#include "spm.hpp"
#include "mspm.hpp"
#include "tspm.hpp"

namespace pg {

Solvers solverToId(std::string s)
{
         if (s == "zlk")  return ZLK;
    else if (s == "pp")   return PP;
    else if (s == "ppp")  return PPP;
    else if (s == "rr")   return RR;
    else if (s == "dp")   return DP;
    else if (s == "rrdp") return RRDP;
    else if (s == "psi")  return PSI;
    else if (s == "spm")  return SPM;
    else if (s == "tspm") return TSPM;
    else if (s == "mspm") return MSPM;
    else if (s == "qpt")  return QPT;
    else                  return NONE;
}

std::string solverToString(Solvers s)
{

    switch (s) {
    case ZLK:  return "parallel Zielonka";
    case PP:   return "PP";
    case PPP:  return "PP+";
    case RR:   return "RR";
    case DP:   return "DP";
    case RRDP: return "RR/DP";
    case PSI:  return "parallel SI";
    case SPM:  return "small progress measures";
    case TSPM: return "traditional small progress measures";
    case MSPM: return "Maciej' modified small progress measures";
    case QPT:  return "QPT progress measures";
    default:   return "(none)";
    }
}

bool solverIsParallel(Solvers s)
{
    return s == ZLK or s == PSI;
}

void
listSolvers(std::ostream &out)
{
    out << "List of solvers:" << std::endl;
    out << "* zlk:  parallel Zielonka" << std::endl;
    out << "* pp:   priority promotion (basic)" << std::endl;
    out << "* ppp:  priority promotion PP+" << std::endl;
    out << "* rr:   priority promotion RR" << std::endl;
    out << "* dp:   priority promotion PP+ with DP strategy" << std::endl;
    out << "* rrdp: priority promotion RR with DP strategy" << std::endl;
    out << "* psi:  parallel strategy improvement" << std::endl;
    out << "* spm:  accelerated small progress measures" << std::endl;
    out << "* tspm: traditional small progress measures" << std::endl;
    out << "* mspm: Maciej' modified small progress measures" << std::endl;
    out << "* qpt:  quasi-polynomial time progress measures" << std::endl;
}       

Solver *
constructSolver(Solvers s, Oink *oink, Game *game, std::ostream &lgr)
{
    switch (s) {
    case ZLK:  return new ZLKSolver(oink, game, lgr);
    case PP:   return new PPSolver(oink, game, lgr);
    case PPP:  return new PPPSolver(oink, game, lgr);
    case RR:   return new RRSolver(oink, game, lgr);
    case DP:   return new DPSolver(oink, game, lgr);
    case RRDP: return new RRDPSolver(oink, game, lgr);
    case PSI:  return new PSISolver(oink, game, lgr);
    case SPM:  return new SPMSolver(oink, game, lgr);
    case TSPM: return new TSPMSolver(oink, game, lgr);
    case MSPM: return new MSPMSolver(oink, game, lgr);
    case QPT:  return new QPTSolver(oink, game, lgr);
    default:   return NULL;
    }
}

}
