# Changelog

All notable changes to Oink will be documented in this file.

## [Unreleased]

### Added

- Solver ZLKJ, Zielonka with justifications, after Algorithm 5 of
  Lapauw, Bruynooghe and Denecker (VMCAI 2020).
- Solver PPJ, priority promotion with justifications, after Algorithm 4.4
  of the thesis of Ruben Lapauw.
- Solver PMA, progress measures with attractors, formerly only available
  as the `notangles` option of PMTL. PMA is experimental and offers small,
  succinct, ordered and Strahler measures.
- The tool `bench_solvers` benchmarks the solvers with verified, repeated
  and interleaved trials, in worker processes with timeouts, optional
  once-per-game preprocessing, a live table with a review, JSON Lines and
  CSV reports, and named experiments that can be resumed. Its summary
  includes the geometric mean relative to the fastest solver per game.
  See `test/bench/README.md`.
- `test_solvers -i` chooses the solvers and their configurations in an
  interactive screen.
- Solution files can hold multi-strategies: all winning moves of a vertex,
  separated by commas. Oink writes the multi-strategies of FPIM and FPJM
  this way (before, one of their moves), and reads and verifies them.
- The solver registry reports solver families, experimental status and
  configuration syntax. `Oink::effectiveSolver()` and
  `Oink::effectiveSolverOptions()` report the solver and options used.

### Changed

- Reading a PGSolver game now rejects repeated edges, which were previously
  accepted and could lead to wrong solver results. Library callers get a
  `std::runtime_error` from either PGSolver parser, `GameBuilder::build()`,
  `Game::e_add()`, or the adjacency-list `Game` constructor for a repeated
  edge. The error names both endpoints.
- Random games are now generated the same way on every platform. Games
  selected with `--gameseed` may differ from those in earlier releases.
- Sorting a game by priority now keeps vertices with equal priority in
  their original order.
- ZLK and ZLKQ do less work in recursive calls. Priority promotion and
  tangle learning scan fewer edges or revisit fewer regions. The fixpoint
  solvers reuse evaluations and process updates in batches. GPM and PMTL
  do fewer measure operations.
- FPJM keeps the justifications of the new winner within a batch. Its
  multi-strategies can differ from before, and usually have more moves.
- `oink --solvers` shows the available solver options. Oink validates
  options before loading the game. Library callers can catch invalid
  solver options as exceptions.
- Per-step solver traces with `-t` are available in Debug builds. Release
  builds still print the solver summaries.

### Removed

- The solvers PTL and SPPTL (progressive tangle learning), unpublished
  research variants; SPPTL could stop with vertices unsolved.

### Fixed

- A successor beyond the last vertex was accepted when the header of a
  PGSolver game gives the number of vertices; it is now an error.
- Parallel FPI could race while reading and updating its winner bitsets
  or writing Debug traces.
- Oink rejects a partial solution if a strategy move or an opponent move
  leaves its winning region. With `--sol`, errors give the vertex numbers
  from the input file.
- Oink now excludes vertices from an installed partial solution during
  preprocessing, so they are only attracted to once.
- Small measures for Odd could read past the end of the priority counts
  when the highest priority was even.
- `test_solvers` reports a failing external solver's exit status.
- `oink`, `nudge` and `test_solvers` report option parsing errors and exit
  with an error status.
- Fixed a Clang build failure with libstdc++ on Ubuntu 26.04.
- Clang can build the library and extra tools with warnings treated as
  errors.
- PMA and PMTL could leave vertices unsolved: the attractor shortcut
  (optimisation 2) ignored the regions whose lifts were postponed to the
  next round. This bug dates from the first PMTL.
- PMTL could skip lifts in a region that started with a known tangle,
  because it added the head twice.
- Reading a solution file ignored the rest of a malformed line (so
  `1,2x,0` was read as `1`) and the moves of a vertex won by its opponent;
  both are now errors.
- SPM could give a vertex to both players or end with a losing strategy:
  its analysis of measures beyond the counts made two unsound
  conclusions. It now keeps every measure at or below the least fixed
  point.
- The verifier was slow on games with vertices with many moves: its search
  scanned the moves of a vertex from the start each time it returned to
  that vertex. On the game of `counter_ortl 3000` it took over 10 s, now
  0.1 s. It could also accept a wrong solution after more than 2^31 steps
  of the search, as it stored a 64-bit number in an `int`.

## [1.1.1] - 2026-11-14

### Fixed

- TLQ could give wrong solutions: after removing a dominion, it did not
  always revisit the older regions that attracted to the removed heads.
  This bug dates from the first TLQ.

## [1.1.0] - 2026-11-09

This release makes tangle learning and PMTL much faster, adds parallel
tangle learning, tangle learning with justifications and multi-strategies,
and removes the dependency on Boost. The library API changes: the solution
moved from `Game` to the new class `Solution`.

### Added

- Solver TLPAR, tangle learning that attracts large regions in parallel with
  Lace. On small games it runs TL. About 3.5x faster than TL with 64 workers
  on the largest games.
- Solver TLJ, tangle learning with justifications, after Algorithm 7 of
  Lapauw, Bruynooghe and Denecker (VMCAI 2020).
- Multi-strategies: a solution can contain all winning moves of a vertex.
  The solvers FPIM and FPJM compute them, the verifier checks them, and the
  tool `multistrat` shows them.
- Strahler measures for GPM and PMTL (`--gpm-strahler`, `--pmtl-strahler`),
  the measures of Daviaud, Jurdziński and Thejaswini (ICALP 2020), adapted
  from the STRPM solver by Vanessa Flügel and Guillermo A. Pérez. The other
  measures have shorthands as well, such as `--gpm-small` and
  `--pmtl-succinct`.
- Input compressed with xz, next to gzip and bzip2. Each format is enabled
  when its library (zlib, libbz2, liblzma) is found.
- `oink --version` and the header `oink/version.hpp`.
- A pkg-config file `oink.pc`, next to the CMake package.
- In the library: `GameBuilder` to build games, `Solution`, `SolverConfig`,
  `SubgameView` and the enum `Player`.
- The build options `OINK_ENABLE_SANITIZERS`, `OINK_ENABLE_TSAN`,
  `OINK_WARNINGS_AS_ERRORS`, `OINK_BUILD_TOOLS` and `OINK_BUILD_NATIVE`.
- `tc --tangle-loop`, a variant of the two counters games.
- CI on Linux, macOS, FreeBSD and Windows (MSYS2).

### Changed

- The solution (winners and strategies) moved from `Game` to the class
  `Solution`, which `Oink` owns; the game is only input. Games are built with
  `GameBuilder`. Programs that use the library need to be updated.
- TL is about twice as fast over the benchmark sets: model checking 3.2x,
  random games 2.0x, equivalence checking 1.8x. It is now the fastest
  sequential solver of Oink on each set.
- PMTL is 3-5x faster on the games that both versions solve, and solves more
  of the large games within the time limit: one region for all tops with the
  same measure and priority, tangles in partially closed regions, and closed
  regions solved as dominions. GPM and PMTL with succinct and ordered
  measures are faster too.
- FPI and FPJ are faster, and PSI computes its values in rounds instead of
  recursively.
- MSPM implements the algorithm of Gazda and Willemse (GandALF 2015).
- Parsing is faster: the parser reads its input in blocks, and large bzip2
  files are decompressed in parallel. Loading large games is about 4x faster.
- Parallel solvers solve small games sequentially, unless `--force-parallel`
  is given.
- `test_solvers` runs the solvers sequentially unless `-w` is given, reports
  the time to load and parse each game, and accepts a configuration per
  solver, such as `--tlpar=parmin=1,chunk=4`.
- The default build type is Release, and `-march=native` is off by default
  (`OINK_BUILD_NATIVE`), so builds are portable.
- Incoming edges are built only for solvers that use them.
- Lace 2.6.0.

### Removed

- The dependency on Boost.

### Fixed

- PMTL could leave vertices unsolved: the shortcuts of optimisation 2 missed
  vertices that are lifted by attracting a tangle, and opponent vertices
  whose escape was solved in a dominion. Both bugs date from the first PMTL.
- TSPM could stop with a logic error, when a measure computed for a larger
  bound overflowed after the bound was lowered.
- MSPM gave wrong solutions when a successor was solved for the opponent, and
  read past its measures on games with only priority 0.
- RRDP could loop forever.
- A data race in the parallel attractor of Zielonka, and Zielonka in parallel
  on large games, where the Lace task deque overflowed.
- `--scc` attracted only within the SCC after solving a bottom SCC; it now
  attracts over the whole game.
- SPM read out of bounds and could leave the highest priority of a player
  uncomputed; SSI hung on platforms where `char` is unsigned (Linux on
  aarch64).
- Memory leaks in PPQ, TLQ and PMTL, uninitialised vertex arrays in `Game`,
  and `Game::sizeup` from size 1.
