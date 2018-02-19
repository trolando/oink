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

Algorithm       | Description
:-------------- | :----------
NPP  | Priority promotion (bitvector implementation by Benerecetti et al)
PP   | Priority promotion (basic algorithm)
PPP  | Priority promotion PP+ (better reset heuristic)
RR   | Priority promotion RR (even better improved reset heuristic)
DP   | Priority promotion PP+ with the delayed promotion strategy
RRDP | Priority promotion RR with the delayed promotion strategy
ZLK  | (parallel) Zielonka's recursive algorithm
APT  | APT algorithm by Kupferman and Vardi
PSI  | (parallel) strategy improvement
TSPM | Traditional small progress measures
SPM  | Accelerated version of small progress measures
MSPM | Maciej Gazda's modified small progress measures
SSPM | Succinct small progress measures
QPT  | Quasi-polynomial time progress measures

* The priority promotion family of algorithms has been proposed in 2016.
* The Zielonka implementation is inspired by work in 2017.
* The APT algorithm was published in 1998 and again in 2016.
* The parallel strategy improvement implementation is inspired by work in 2017 but uses a different approach with work-stealing.
* The accelerated SPM approach is a novel approach.
* The succinct small progress measures algorithm was published in 2017. (Jurdzinski and Lazic, at LICS 2017)
* The QPT progress measures algorithm was published in 2017. (Fearnley et al, at SPIN 2017.)

The parallel algorithms use the work-stealing framework Lace.

The solver can further be tuned using several pre-processors:

1. Removing all self-loops (recommended)
2. Removing winner-controlled winning cycles (recommended)
3. Inflating or compressing the game before solving it
4. SCC decomposition to solve the parity game one SCC at a time.  
   May either improve or deteriorate performance

Tools
-----

Oink comes with several simple tools that are built around the library
liboink.

Main tools:

Tool    | Description
:------ | :-------------
oink    | The main tool for solving parity games
verify  | Small tool that just verifies a solution (can be done by Oink too)
nudge   | Swiss army knife for transforming parity games
dotty   | Small tool that just generates a .dot graph of a parity game

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

Instructions
-----------

Oink is compiled using CMake.
Optionally, use `ccmake` to set options.
By default, Oink does not compile the extra tools, only the library `liboink` and the tool `oink`.
Oink requires the Boost libraries, in particular `boost_iostreams`.
```
mkdir build && cd build
cmake .. && make && make install
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
- `--no-wcwc`, `--no-loops` and `--no-single` disable preprocessors that eliminate winner-controlled winning cycles, self-loops and single-parity games.
- `-z <seconds>` kills the solver after the given time.
- `--sol <filename>` loads a partial or full solution.
- `--dot <dotfile>` writes a .dot file of the game as loaded.
- `-p` writes the vertices won by even/odd to stdout.
- `-t` (once or multiple times) increases verbosity level.
