#ifndef PSI_HPSI
#define PSI_HPSI

#include <queue>

#include "solver.hpp"
#include "lace.h"

namespace pg {

class PSISolver : public Solver
{
public:
    PSISolver(Oink *oink, Game *game, std::ostream &lgr);
    virtual ~PSISolver();

    virtual void run();
    void setWorkers(int count) { workers = count; }

    int minor = 0, major = 0;

    // variables are public to avoid warnings due to C/C++ issues
    int workers = 0; // number of workers for Lace

    bool si_val_less(int a, int b);
    void print_debug();

    void compute_vals_seq(void);
    int mark_solved_seq(void);
    int switch_strategy_seq(int pl);

    friend void compute_val_WORK(_WorkerP*, _Task*, int, PSISolver*);
    friend void compute_all_val_WORK(_WorkerP*, _Task*, PSISolver*);
    friend int mark_solved_rec_WORK(_WorkerP*, _Task*, PSISolver*, int, int);
    friend int switch_strategy_WORK(_WorkerP*, _Task*, PSISolver*, int, int, int);
};

}

#endif 
