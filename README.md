# Oink
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![CI testing](https://github.com/trolando/oink/actions/workflows/ci-build.yml/badge.svg)](https://github.com/trolando/oink/actions/workflows/ci-build.yml)

Oink is a modern implementation of parity game solvers written in C++.
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

## Implemented algorithms

Oink implements various modern algorithms.
Several of the algorithms can be run with multiple threads.
These parallel algorithms use the work-stealing framework Lace.

### Tangle learning family

The tangle learning algorithms are a focus of research for finding a polynomial-time algorithm to solve parity games.

See also:  
* Tom van Dijk (2018) [Attracting Tangles to Solve Parity Games](https://doi.org/10.1007/978-3-319-96142-2/_14). In: CAV 2018.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
TL              | Tangle learning (standard variant)
RTL             | Recursive Tangle Learning (research variant)
ORTL            | One-sided Recursive Tangle Learning (research variant)
PTL             | Progressive Tangle Learning (research variant)
SPPTL           | Single-player Progressive Tangle Learning (research variant)
DTL             | Distance Tangle Learning (research variant based on 'good path length')
IDTL            | Interleaved Distance Tangle Learning (interleaved variant of DTL)
TLQ             | Tangle Learning applied to the quasi-polynomial time recursive algorithm

### Fixpoint algorithms

The fixpoint algorithms are based on a translation of parity games to a formula of the modal mu-calculus, which can be solved
using a naive fixpoint algorithm. Despite its atrocious behavior on constructed hard games, for practical parity games
they are among the fastest algorithms.
The FPI algorithm has a scalable parallel implementation.

See also:  
* Tom van Dijk and Bob Rubbens (2019) [Simple Fixpoint Iteration To Solve Parity Games](https://doi.org/10.4204/EPTCS.305.9). In: GandALF 2019.  
* Ruben Lapauw, Maurice Bruynooghe, Marc Denecker (2020) [Improving Parity Game Solvers with Justifications](https://doi.org/10.1007/978-3-030-39322-9_21). In: VMCAI 2020.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
FPI             | (parallel) Distraction Fixpoint Iteration (similar to APT algorithm and to Bruse/Falk/Lange)
FPJ             | Fixpoint Iteration using Justifications
FPJG            | Fixpoint Iteration using Justifications (greedy variant by TvD)

### Priority promotion family

The priority promotion algorithms are efficient algorithms described in several papers by Benerecetti, Dell'Erba and Mogavero.
The implementations PP, PPP, RR, DP and RRDP are by TvD; the authors provided their own implementation NPP.

See also:
* Massimo Benerecetti, Daniele Dell'Erba, Fabio Mogavero (2016) [Solving Parity Games via Priority Promotion](https://doi.org/10.1007/978-3-319-41540-6_15). In: CAV 2016.
* Massimo Benerecetti, Daniele Dell'Erba, Fabio Mogavero (2016) [A Delayed Promotion Policy for Parity Games](https://doi.org/10.4204/EPTCS.226.3). In: GandALF 2016.
* Massimo Benerecetti, Daniele Dell'Erba, Fabio Mogavero (2016) [Improving Priority Promotion for Parity Games](https://doi.org/10.1007/978-3-319-49052-6_8). In: HVC 2016.
* Massimo Benerecetti, Daniele Dell'Erba, Fabio Mogavero, Sven Schewe, Dominik Wojtczak [Priority Promotion with Parysian Flair](https://arxiv.org/abs/2105.01738). In: arXiv.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
NPP             | Priority promotion (implementation by authors BDM)
PP              | Priority promotion (basic algorithm)
PPP             | Priority promotion PP+ (better reset heuristic)
RR              | Priority promotion RR (even better reset heuristic)
DP              | Priority promotion PP+ with the delayed promotion strategy
RRDP            | Priority promotion RR with the delayed promotion strategy
PPQ             | Pawel Parys's QP recursive algorithm augmented with priority promotion

### Zielonka's recursive algorithm

One of the oldest algorithms for solving parity games is the recursive algorithm due to Zielonka, based on work by McNaughton.
The standard implementation ZLK is described in the TACAS 2018 paper about Oink. It is an optimized version based on earlier optimizations by various authors.
The UZLK variant removes some of the optimizations for research purposes.

The ZLKQ algorithm is an implementation by TvD based on work by Parys and by Lehtinen et al., with additional optimizations,
including a shortcut in the tree to the next priority in the subgame and skipping recursions if the remaining game is smaller than the precision parameter.

The ZLKPP algorithms are the implementations by Pawel Parys accompanying the paper with Lehtinen et al.

See also:
* Pawel Parys (2019) [Parity Games: Zielonka's Algorithm in Quasi-Polynomial Time](https://doi.org/10.4230/LIPIcs.MFCS.2019.10). In: MFCS 2019.
* Karoliina Lehtinen, Pawel Parys, Sven Schewe, Dominik Wojtczak (2021) [A Recursive Approach to Solving Parity Games in Quasipolynomial Time](https://arxiv.org/abs/2104.09717).

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
ZLK             | (parallel) Zielonka's recursive algorithm
UZLK            | an unoptimized version of Zielonka for research purposes
ZLKQ            | Quasi-polynomial time recursive algorithm (optimized implementation by TvD)
ZLKPP-STD       | Zielonka's recursive algorithm (standard version, by Pawel Parys)
ZLKPP-WAW       | Quasi-polynomial time recursive algorithm (Warsaw version, by Pawel Parys)
ZLKPP-LIV       | Quasi-polynomial time recursive algorithm (Liverpool version, by Pawel Parys)

### Strategy improvement

The strategy improvement algorithm is one of the older algorithms, receiving a lot of attention from the academic community in the past.
The parallel strategy improvement implementation is inspired by the work of Fearnley in CAV 2017 but uses a different approach with work-stealing.
This is described in the TACAS 2018 paper about Oink.
The symmetric strategy improvement algorithm is my own take on the 2015 ICALP paper by Schewe et al. It applies the principle of only updating strategies that are also the best response to the SI variation of PSI, that is, considering only finite paths simplifying the valuation.

See also:
* John Fearnley (2017) [Efficient Parallel Strategy Improvement for Parity Games](https://doi.org/10.1007/978-3-319-63390-9_8). In: CAV 2017.
* Sven Schewe, Ashutosh Trivedi, Thomas Varghese (2015) [Symmetric Strategy Improvement](https://doi.org/10.1007/978-3-662-47666-6_31). In: ICALP 2015.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
PSI             | (parallel) strategy improvement
SSI             | symmetric strategy improvement

### Progress measures

Progress measures algorithms are based on a monotonically increasing value attached to every vertex of the parity game.
Several of the quasi-polynomial solutions to parity games are modifications of the small progress measures algorithm.

The original algorithm by Jurdzinski in 2000 is implemented as the TSPM algorithm.
The accelerated SPM approach is a novel approach developed by TvD.
The idea is to let progress measures increase until a bound, then analyse the halted result.
The measures are then immediately increased to a higher value, skipping many small increases.

The QPT ordered progress measures algorithm was published by Fearnley et al in SPIN 2016.  The QPT succinct progress measures algorithm was published by Jurdzinski et al in LICS 2017.

The bounded variants `BSSPM` and `BQPT` are an idea by Tom van Dijk and Marcin Jurdzinski.

See also:
* Marcin Jurdzinski (2000) [Small Progress Measures for Solving Parity Games](https://doi.org/10.1007/3-540-46541-3_24). In: STACS 2000.
* Maciej Gazda, Tim Willemse (2015) [Improvement in Small Progress Measures](https://doi.org/10.4204/EPTCS.193.12). In: GandALF 2015.
* John Fearnley, Sanjay Jain, Sven Schewe, Frank Stephan, Dominik Wojtczak (2017) [An ordered approach to solving parity games in quasi polynomial time and quasi linear space](https://doi.org/10.1145/3092282.3092286). In: SPIN 2017.
* Marcin Jurdzinski, Ranko Lazic (2017) [Succinct progress measures for solving parity games](https://doi.org/10.1109/LICS.2017.8005092). In: LICS 2017.

Algorithm       | Description
:-------------- | :-------------------------------------------------------------------------------------------------
TSPM            | Traditional small progress measures (standard variant without cycle acceleration)
SPM             | Accelerated version of small progress measures (cycle acceleration variant by TvD)
MSPM            | Maciej Gazda's modified small progress measures
SSPM            | Quasi-polynomial time succinct progress measures
BSSPM           | Quasi-polynomial time succinct progress measures (bounding variant by TvD and MJ, often faster)
QPT             | Quasi-polynomial time ordered progress measures
BQPT            | Quasi-polynomial time ordered progress measures (bounding variant by TvD and MJ, often faster)

### Preprocessing

Oink can also apply the following preprocessors before solving the game:

* Removing all self-loops (recommended)
* Removing winner-controlled winning cycles (recommended)
* Inflating or compressing the priorities before solving the parity game.  
  Compression is useful for several solvers. By default, priorities are simply renumbered (to remove gaps). Several algorithms actually perform on-the-fly compression and are not affected by renumbering or compression.
* SCC decomposition to solve the parity game one SCC at a time.  
  May either improve or deteriorate performance.

## Tools

Oink comes with several simple tools that are built around the library `liboink`.

### Main tools

Tool          | Description
:------------ | :-------------
oink          | The main tool for solving parity games
verify        | Small tool that just verifies a solution (can be done by Oink too)
nudge         | Swiss army knife for transforming parity games
dotty         | Small tool that just generates a .dot graph of a parity game
test\_solvers | Main testing tool for parity game solvers and benchmarking

### Game generators

Tool           | Description
:------------- | :----------
rngame         | Faster version of the random game generator of PGSolver
stgame         | Faster version of the steady game generator of PGSolver
counter\_rr    | Counterexample to the RR solver
counter\_dp    | Counterexample to the DP solver
counter\_m     | Counterexample of Maciej Gazda, PhD Thesis, Sec. 3.5
counter\_qpt   | Counterexample of Fearnley et al, [An ordered approach to solving parity games in quasi polynomial time and quasi linear space](https://doi.org/10.1145/3092282.3092286). In: SPIN 2017.
counter\_core  | Counterexample of Benerecetti, Dell'Erba, Mogavero, [Robust Exponential Worst Cases for Divide-et-Impera Algorithms for Parity Games](https://doi.org/10.4204/EPTCS.256.9). In: GandALF 2017.
counter\_rob   | SCC version of counter\_core.
counter\_dtl   | Counterexample to the DTL solver
counter\_ortl  | Counterexample to the ORTL solver
counter\_symsi | Counterexample of Matthew Maat to (standard) symmetric strategy improvement
tc             | Two binary counters generator (game family that is an exponential lower bound for many algorithms). See also Tom van Dijk (2019) [A Parity Game Tale of Two Counters](https://doi.org/10.4204/EPTCS.305.8). In: GandALF 2019.
tc+            | TC modified to defeat the RTL solver

### Two binary counters

The two binary counters game family is an exponential lower bound for many algorithms:
* The tangle learning algorithm TL
* All fixpoint algorithms FPI, FPJ, FPJG
* All priority promotion algorithms NPP, PP, PPP, DP, RRDP
* All normal Zielonka algorithms ZLK, UZLK, ZLKPP-STD
* The progress measures algorithms TSPM, SPM, MSPM

The two binary counters game family appears to be a quasi-polynomial lower bound for the QP algorithms:
* The Zielonka variations ZLKQ, ZLKPP-WAW, ZLKPP-LIV
* The progress measures variations SSPM, QPT

Some algorithms can solve the two binary counters games in polynomial time:
* The tangle learning algorithms RTL, ORTL, PTL, SPPTL, DTL, IDTL
* The strategy improvement algorithm PSI
* Unsure: the BSSPM and BQPT algorithms

## Usage

### Prerequisites

- A C++17 compiler (GCC 9 or newer, Clang, or MSVC) and CMake 3.14 or newer.
- [Lace](https://github.com/trolando/lace), the work-stealing framework used by the parallel solvers. By default CMake downloads it automatically with `git` during configuration, so the first build needs network access. To build offline, install Lace yourself and CMake will pick it up via `find_package(lace)`.
- Optional: the `zlib`, `libbz2` and `liblzma` development packages enable transparent decompression of `.gz`, `.bz2` and `.xz` parity games. They are detected automatically; a missing library simply disables that input format.

Oink does **not** require Boost.

### Building

Oink uses an out-of-tree CMake build. By default it builds the library `liboink`, the main tool `oink`, the test runner `test_solvers` and the example program (use `ccmake` or `-D` options such as `OINK_BUILD_EXTRA_TOOLS=ON` to change this).

```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build
```

### Getting started

Solve one of the bundled example games and verify the solution:

```
./build/oink -v examples/abcg_arbiter.tlsf.ehoa.pg
```

Add `-p` to print which vertices each player wins.

Oink provides usage instructions via `oink --help`. Typically, Oink is provided a parity game either
via stdin (default) or from a file. The file may be zipped using the gzip or bzip2 format, which is detected if the
filename ends with `.gz` or `.bz2`.

What you want?                          | But how then?
:-------------------------------------- | :---------------------------------
To quickly solve a gzipped parity game: | `oink -v game.pg.gz game.sol`
To verify some solution:                | `oink -v game.pg.gz --sol game.sol`

A typical call to Oink is: `oink [options] [solver] <filename> [solutionfile]`. This reads a parity game from `filename`, solves it with the chosen solver (default: `--tl`), then writes the solution to `<solutionfile>` (default: don't write).
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

## Using Oink as a library

After `cmake --install build`, Oink ships a CMake package and a pkg-config file, so a consuming project can use it directly:

```cmake
find_package(oink REQUIRED)
target_link_libraries(my_app PRIVATE oink::oink)
```

A minimal end-to-end example — parse a parity game, solve it, and print the winner and strategy of each vertex — is in [`examples/simple.cpp`](examples/simple.cpp), which the build compiles to `oink-example-simple`:

```
./build/oink-example-simple examples/abcg_arbiter.tlsf.ehoa.pg
```

For complete standalone consumer projects (CMake and pkg-config), see `test/install/consumer-cmake` and `test/install/consumer-pkgconfig`.

## Testing an external solver

The `test_solvers` tool can run and **verify** an external solver — not only Oink's built-in ones — against a directory (or a list) of parity games. This is the easiest way to check your own solver for correctness against the bundled corpus.

Use `--external name:command` (or `-e`), where `command` is how your solver is invoked, with two placeholders:

- `%I` — **input**: the parity game (in PGSolver format) that `test_solvers` writes for each game.
- `%O` — **output**: the file your solver must write its solution to.

For every game, `test_solvers` writes the game to `%I`, runs your command, reads the solution back from `%O`, and verifies the winning regions and strategies. For example:

```
./build/test_solvers tests --external "mysolver:mysolver %I %O"
./build/test_solvers tests -e "py:python3 solve.py %I %O"
```

`name` is just a label used in the report (the first `:` separates it from the command). You can list several `--external` solvers, and mix in built-in ones (e.g. add `--tl`) to compare in the same run.

### Solution format

Your solver reads a parity game in PGSolver format (the same format as the bundled `.pg` games) from `%I`, and writes a solution to `%O`. The solution has one line per vertex:

```
paritysol <count>;            # optional header line, ignored on input
<vertex> <winner>;            # winner is 0 (Even) or 1 (Odd)
<vertex> <winner> <strategy>; # include <strategy> only when the winner owns the vertex
```

`<strategy>` is the successor the winning player moves to; it is required exactly when `<winner>` equals the vertex's owner, and omitted otherwise. Every vertex must appear (a full solution is expected).
