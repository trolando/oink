# Changelog

All notable changes to Oink will be documented in this file.

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
