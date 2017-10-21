Oink
====
Oink is an implementation of modern parity game solvers written in C++.
Oink aims to provide high-performance implementations of state-of-the-art
algorithms representing different approaches to solving parity games.

Oink is licensed with the Apache 2.0 license and is aimed at both researchers
and practitioners. For convenience, Oink compiles into a library that can
be used by other projects.

Oink is developed (&copy; 2017) by Tom van Dijk and the
[Formal Methods and Verification](http://fmv.jku.at/)
group at the Johannes Kepler University Linz as part of the RiSE project.

The main author of Oink is Tom van Dijk who can be reached via <tom@tvandijk.nl>.
Please let us know if you use Oink in your projects.

The main repository of Oink is https://github.com/trolando/oink.

Implemented algorithms
----------------------

Oink implements various modern algorithms.

1. Various algorithms from the ''priority promotion'' family
  * PP (standard priority promotion)
  * PPP (priority promotion plus, improved reset heuristic)
  * RR (PP with an even better improved reset heuristic)
  * DP (PPP with the delayed promotion strategy)
  * RRDP (RR with the delayed promotion strategy)
2. Zielonka's recursive algorithm (recursion-free)
3. QPT progress measures (theoretically quasi-polynomial)
4. Parallel strategy improvement
5. Three variations of small progress measures
  * Traditional SPM, extended with updates to priority caps and regular cross-parity updates
  * Maciej' Modified SPM algorithm
  * Fast SPM, which is traditional SPM but instead of slowly increasing winning measures,
    it directly computes the best escape.

The priority promotion family of algorithms have been proposed in 2016.
The Zielonka implementation is inspired by work in 2017.
The parallel strategy improvement implementation is inspired by
work in 2017 but uses a different approach with work-stealing.
The QPT progress measures algorithm was published in 2016.
The fast SPM approach is a new approach.

The parallel Zielonka algorithm and parallel strategy improvement
are parallelized using the work-stealing framework Lace.

The solver can further be tuned using several pre-processors.

1. Removing all self-loops (recommended)
2. Removing winner-controlled winning cycles (recommended)
3. Inflating or compressing the game before solving it
4. SCC decomposition to solve the parity game one SCC at a time.  
   May either improve or deteriorate performance

Tools
-----

Oink comes with several simple tools that are built around the library
liboink. 

oink

verify

For generating games:

rngame
Similar to the random game generator that comes with PGSOLVER but faster.

stgame
Similar to the steady game generator that comes with PGSOLVER but faster.

counter\_rr

counter\_dp

counter\_m

counter\_qpt

counter\_core

For transforming games:

nudge
