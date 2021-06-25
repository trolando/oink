Oink
====
Oink is an implementation of modern parity game solvers written in C++.
Oink aims to provide high-performance implementations of state-of-the-art
algorithms representing different approaches to solving parity games.

Oink is licensed with the Apache 2.0 license and is aimed at both researchers
and practitioners. For convenience, Oink compiles into a library that can
be used by other projects.

Oink was initially developed (&copy; 2017-2021) by Tom van Dijk and the
[Formal Models and Verification](http://fmv.jku.at/)
group at the Johannes Kepler University Linz as part of the RiSE project,
and now in the [Formal Methods and Tools](https://fmt.ewi.utwente.nl) group
at the University of Twente.

The main author of Oink is Tom van Dijk who can be reached via <tom@tvandijk.nl>.
Please let us know if you use Oink in your projects.
If you are a researcher, please cite the following publication.

Tom van Dijk (2018) [Oink: An Implementation and Evaluation of Modern Parity Game Solvers](https://doi.org/10.1007/978-3-319-89960-2_16). In: TACAS 2018.

The main repository of Oink is https://github.com/trolando/oink.

Implemented algorithms
----------------------

Oink implements various modern algorithms.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
TL              | Tangle learning (standard variant)
RTL             | Recursive Tangle Learning (research variant)
PTL             | Progressive Tangle Learning (research variant)
SPPTL           | Single-player Progressive Tangle Learning (research variant)
DTL             | Distraction-free Tangle Learning (research variant)
IDTL            | Interleaved Distraction-free Tangle Learning (research variant)
FPI             | (parallel) Distraction Fixpoint Iteration (similar to APT algorithm and to Bruse/Falk/Lange)
FPJ             | Fixpoint Iteration using Justifications
FPJG            | Fixpoint Iteration using Justifications, greedy variant
NPP             | Priority promotion (implementation by authors Benerecetti et al)
PP              | Priority promotion (basic algorithm)
PPP             | Priority promotion PP+ (better reset heuristic)
RR              | Priority promotion RR (even better improved reset heuristic)
DP              | Priority promotion PP+ with the delayed promotion strategy
RRDP            | Priority promotion RR with the delayed promotion strategy
ZLK             | (parallel) Zielonka's recursive algorithm
UZLK            | an unoptimized version of Zielonka for research purposes
ZLKQ            | Quasi-polynomial time recursive algorithm (optimized implementation by TvD)
ZLKPP-STD       | Zielonka's recursive algorithm (standard version, by Pawel Parys)
ZLKPP-WAW       | Zielonka's recursive algorithm (Warsaw version, by Pawel Parys)
ZLKPP-LIV       | Zielonka's recursive algorithm (Liverpool version, by Pawel Parys)
PSI             | (parallel) strategy improvement
TSPM            | Traditional small progress measures (standard variant without cycle acceleration)
SPM             | Accelerated version of small progress measures (cycle acceleration variant by TvD)
MSPM            | Maciej Gazda's modified small progress measures
SSPM            | Quasi-polynomial time succinct progress measures
BSSPM           | Quasi-polynomial time succinct progress measures (bounding variant by TvD and MJ, often faster)
QPT             | Quasi-polynomial time ordered progress measures
BQPT            | Quasi-polynomial time ordered progress measures (bounding variant by TvD and MJ, often faster)

* The Zielonka implementation is inspired by work in 2017 with additional improvements by TvD.
* The priority promotion family of algorithms has been proposed in 2016.
* The parallel strategy improvement implementation is inspired by the work of Fearnley in CAV 2017 but uses a different approach with work-stealing.
* The accelerated SPM approach is a novel approach developed by Tom van Dijk.  
  The idea is to let progress measures increase until a bound, then analyse the halted result. The measures are then immediately increased to a higher value, skipping many small increases.
* The QPT ordered progress measures algorithm was published by Fearnley et al in SPIN 2016.
* The QPT succinct progress measures algorithm was published by Jurdzinski et al in LICS 2017.
* The bounded variants `BSSPM` and `BQPT` are an idea by Tom van Dijk and Marcin Jurdzinski.
* The ZLKQ algorithm is an implementation by TvD based on work by Parys, Lehtinen et al. https://arxiv.org/abs/2104.09717 with additional optimizations, including a shortcut in the tree to the next priority in the subgame and skipping recursions if the remaining game is smaller than the precision parameter.
* The ZLKPP algorithms are the implementations by Pawel Parys accompanying the paper "Karoliina Lehtinen, Pawel Parys, Sven Schewe, Dominik Wojtczak: A Recursive Approach to Solving Parity Games in Quasipolynomial Time." https://arxiv.org/abs/2104.09717

The parallel algorithms use the work-stealing framework Lace.

The solver can further be tuned using several pre-processors:

1. Removing all self-loops (recommended)
2. Removing winner-controlled winning cycles (recommended)
3. Inflating or compressing the priorities before solving it (compression is useful for many solvers)
4. SCC decomposition to solve the parity game one SCC at a time.  
   May either improve or deteriorate performance

Tools
-----

Oink comes with several simple tools that are built around the library `liboink`.

Main tools:

Tool          | Description
:------------ | :-------------
oink          | The main tool for solving parity games
verify        | Small tool that just verifies a solution (can be done by Oink too)
nudge         | Swiss army knife for transforming parity games
dotty         | Small tool that just generates a .dot graph of a parity game
test\_solvers | Main testing tool for parity game solvers and benchmarking

Tools to generate games:

Tool           | Description
:------------- | :----------
rngame         | Faster version of the random game generator of PGSolver
stgame         | Faster version of the steady game generator of PGSolver
counter\_rr    | Counterexample to the RR solver
counter\_dp    | Counterexample to the DP solver
counter\_m     | Counterexample of Maciej Gazda, PhD Thesis, Sec. 3.5
counter\_qpt   | Counterexample of Fearnley et al, An ordered approach to solving parity games in quasi polynomial time and quasi linear space, SPIN 2017
counter\_core  | Counterexample of Benerecetti et al, Robust Exponential Worst Cases for Divide-et-Impera Algorithms for Parity Games, GandALF 2017
counter\_rob   | SCC version of counter\_core.
counter\_dtl   | Counterexample to the DTL solver
tc             | Two binary counters generator (game family that is an exponential lower bound for many of the recursive algorithms). See https://doi.org/10.4204/EPTCS.305.8 for details.
tc+            | TC modified to defeat the RTL solver

Usage
-----------

Oink is compiled using CMake.
Optionally, use `ccmake` to set options.
By default, Oink does not compile the extra tools, only the library `liboink` and the main tools `oink` and `test_solvers`.
Oink requires several Boost libraries.
```
mkdir build && cd build
cmake .. && make
ctest
```

Oink provides usage instructions via `oink --help`. Typically, Oink is provided a parity game either
via stdin (default) or from a file. The file may be zipped using the gzip or bzip2 format, which is detected if the
filename ends with `.gz` or `.bz2`.

What you want?                          | But how then?
:-------------------------------------- | :---------------------------------
To quickly solve a gzipped parity game: | `oink -v game.pg.gz game.sol`
To verify some solution:                | `oink -v game.pg.gz --sol game.sol`

A typical call to Oink is: `oink [options] [solver] <filename> [solutionfile]`. This reads a parity game from `filename`, solves it with the chosen solver (default: `--npp`), then writes the solution to `<solutionfile>` (default: don't write).
Typical options are:
- `-v` verifies the solution after solving the game.
- `-w <workers>` sets the number of worker threads for parallel solvers. By default, these solvers run their sequential version. Use `-w 0` to automatically determine the maximum number of worker threads.
- `--inflate` and `--compress` inflate/compress the game before solving it.
- `--scc` repeatedly solves a bottom SCC of the parity game.
- `--no-wcwc`, `--no-loops` and `--no-single` disable preprocessors that eliminate winner-controlled winning cycles, self-loops and single-parity games. Use `--no` to disable all preprocessors.
- `-z <seconds>` kills the solver after the given time.
- `--sol <filename>` loads a partial or full solution.
- `--dot <dotfile>` writes a .dot file of the game as loaded.
- `-p` writes the vertices won by even/odd to stdout.
- `-t` (once or multiple times) increases verbosity level.
