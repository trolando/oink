# bench_solvers

`bench_solvers` benchmarks Oink's solvers on parity game files or on random
games. Every trial solves a fresh copy of the game and is verified; the
metric is the wall time of `Oink::run`. It replaces the benchmark mode of the
earlier `test_solvers --benchmark`; `test_solvers` itself tests solvers.

```
./build/bench_solvers examples                         # the default solvers on a directory
./build/bench_solvers --tl --tlq --ppq --reference=tl --jsonl=run.jsonl --csv=run.csv examples
./build/bench_solvers -i examples                      # setup screen, live table, review
./build/bench_solvers --live --review --heatmap examples
./build/bench_solvers --experiment=experiments/synthesis-01 --bench-seed=42 examples
./build/bench_solvers --resume=experiments/synthesis-01
./build/bench_solvers --tl --zlk --count=50 --size=200 # random games
```

Without solver flags it benchmarks tl, tlj, tlq, npp, fpi, fpj and zlk, in
that order. Solvers are selected as in `test_solvers`: `--tl`,
`--tlpar=parmin=1,chunk=8` for a configured variant (several variants of one
solver may be given), `-c` to configure the only selected solver, `--all`,
and `-e name:command` for an external solver (`%I` is the game, `%O` the
solution). `bench_solvers --help` lists all options.

## What is measured

The games are read one at a time. Files are canonicalised, deduplicated,
sorted and then shuffled with the schedule seed (`--sort` keeps the sorted
order). In directories, files with the extensions .pg, .gm, .pgsolver, .gz,
.bz2 and .xz count, and files without an extension that start with `parity`;
`--recursive` searches subdirectories. Explicit files can have any name. A
file that cannot be read is a failure of every solver on that game.

For every game:

1. It is parsed (priorities above INT_MAX are renumbered, self-loops are
   kept) and sorted, outside any timing.
2. With `--preprocess` it is preprocessed once, outside the timing: Oink
   renumbers the priorities and solves single-parity games, self-loops and
   winner-controlled winning cycles. `--single`, `--loops`, `--wcwc`,
   `--inflate` and `--compress` choose another recipe. The partial solution is
   verified, and the solvers solve the residual game with all preprocessing of
   Oink disabled. A game that preprocessing solves completely has no trials.
   Preprocessing is off by default.
3. The rounds run: `--warmup` warm-up rounds (default 1), then three pilot
   trials when calibrating, then the measured repetitions. Every round runs
   the solvers in a shuffled order. With `--min-solve-time` (default 0.001
   seconds; 0 disables it) the median of the pilots fixes the sample count of
   each solver on the game, max(`--repeat`, ceil(target/pilot)), at most
   `--max-repeat` (default 10000), before the measured trials start.
   `--repeat` defaults to 7. `--quick` means one repetition, no warm-up and no
   calibration.
4. Every trial copies the residual game and constructs Oink (with the
   eager in-edges) outside the solve interval. The solve interval runs from
   immediately before calling `Oink::run` to its return. It includes Oink's
   solver selection, construction and bookkeeping, `--scc`, and any
   incoming edges built inside the call, including lazy construction. The
   same boundary applies to every solver. With preprocessing, the solution
   and all moves of a multi-strategy are mapped
   back and verified against the original game. A failed trial (solver
   error, invalid solution, timeout, crash) fails the solver on that game;
   its later trials are reported as skipped, and the game and log are saved
   in `--artifacts` (default `benchmark-artifacts`).

On POSIX systems the trials of a game run in a forked worker process in its
own process group, which serves all trials of the game. `--timeout` limits
each solve to 30 wall seconds by default. Use 0 to disable the solve limit.
`--safety-timeout` gives setup and verification 300 wall seconds each by
default. The first setup watchdog includes worker and Lace startup. The
reported setup time measures the per-trial copy and configuration. The
verification phase includes collecting solver details and checking the
solution. Input loading and once-per-game preprocessing precede the
worker's deadlines.
Either timeout stops the worker and its subprocesses. After a failure the
next trial gets a new worker. `--in-process` runs the trials in the tool itself,
without crash isolation or timeouts; this is the only mode on Windows.

Trial and result records give the failure phase (`setup`, `solve` or
`verification`). A solve deadline gives `timeout`, and a safety deadline
gives `safety_timeout`. Allocation failures with explicit evidence give
`memory_failure`; unexplained worker exits and signals give `crash`.
Trial records retain the timings already reported by the worker and use
null for unavailable timings. CSV includes the failed cell's phase.

The schedule seed (`--bench-seed`, random by default and always reported)
reproduces the file order, the solver order of every round and the solver
RNG seeds; the series of random games is derived from it unless `--seed` is
given. The solver seed is fixed per game unless `--vary-solver-seed`. With
`-w`, Lace runs in the worker, once per game.

## Output

As each game finishes, stderr gets a line with the median and interquartile
range of every solver. At the end, stdout gets the summary: per solver the
coverage (verified games out of the games that need a solver; `!N` counts
failures), the sum of the per-game medians, `geom. VBS`, the geometric mean
of its median over the fastest median of each game (the virtual best
solver; a solve timeout counts as the time limit), and `VBS W/T`, the games where
it is uniquely fastest or within 1% of the fastest median. With `--reference`
(a selected name, its id if it has one variant, or `first`), the summary adds
the geometric mean speedup against the reference over the games both
verified, and wins/ties/losses with a 1% tie band. Solvers with full coverage
come first, by their sum; on a terminal the sums and `geom. VBS` are coloured
green (fastest) to red. These are point estimates, not significance tests.

`--jsonl=FILE` writes JSON Lines, flushed as they happen: the metadata (build,
machine, settings, solvers, a clock probe), every game, preparation, warm-up,
pilot and measured trial (status, setup, solve, verification and CPU seconds,
the solver that actually ran), the calibration decisions, the result per game
and solver (median, quartiles, MAD), pauses, extractions, the summary and an
end marker. `--csv=FILE` writes one row per game and solver with the same
statistics and the residual size. Existing reports are only replaced with
`--overwrite`. The exit status is 0 on success, 1 if a trial failed, 130 when
interrupted, and 3 for a bad command line.

## Live table and review

`--live` shows a table of the medians that updates as the trials finish, with
the sums, the coverage, a progress bar and the current trial below it. A cyan
`>` marks the selected game, the fastest complete median of a row is green,
failures are red. `h` or `--heatmap` colours every timing by its rank in the
row instead (green for the fastest 5%, red for the slowest 5%).

The table is drawn and the keys are read between trials only; keys pressed
during a trial wait for it to finish.

| Key | Action |
| :-- | :----- |
| Up/Down, k/j | select a game (the newest is followed until you move up) |
| Left/Right, [ ] | page the solver columns |
| h | heatmap on or off |
| e | extract the selected game |
| Space | resume (browsing, h and e pause the measurements) |
| q, Ctrl-C | stop the benchmark; the results so far are kept |

With `--review` the table stays open after the benchmark. Then also:

| Key or mouse | Action |
| :----------- | :----- |
| click a solver header, 1-9 | sort by that solver (fastest first); again reverses; failures stay last |
| click GAME, g | sort by name |
| o | the original order |
| click a row, wheel | select a game |
| q | leave the review and print the summary |

`e` writes the game to `--exports` (default `benchmark-extracts`) as a `.pg`
file with a `.jsonl` sidecar: the metadata, the game, its results and the
measured solve times. With preprocessing, the residual game, the partial
solution and the vertex map are written too. Pauses are reported in the JSON
Lines; measure for publication without browsing.

Without a terminal on stdin and stderr, with `--plain`, `TERM=dumb` or
`NO_COLOR`, the output is plain text without escape sequences.

## Setup screen

`-i` starts with the solver chooser (the same as `test_solvers -i`), then
runs the benchmark with the live table and the review. AVAILABLE lists the
solvers by family with their descriptions, SELECTED is the roster in its
order.

| Key | Action |
| :-- | :----- |
| Up/Down, PageUp/PageDown, Home/End | move in the active pane |
| Tab, Left/Right | switch panes |
| Enter | add the solver; in SELECTED, edit its configuration |
| Space | add or remove the solver; in SELECTED, remove the row |
| Delete, Backspace, - | remove the selected row |
| c / d | configure / add a configured variant (solvers with options) |
| / | filter by id, description or family |
| x | show or hide experimental solvers (hidden by default) |
| R | mark the selected row as the reference, or clear it |
| r w t b p | repetitions, warm-ups, timeout, schedule seed, sampling target |
| P | preprocessing on or off |
| s | start |
| q, Escape | cancel |

Configurations are checked by the solver's own parser; the editor shows the
accepted syntax and refuses a duplicate variant. In an editor, Enter saves,
Escape cancels and Ctrl-U clears the text.

## Named experiments

`--experiment=DIR` runs a benchmark that can be resumed, with a progress bar
instead of the table (`--plain` prints a line per game). `-i` shows the setup
screen first. DIR holds:

| File | Contents |
| :--- | :------- |
| plan.txt | settings, solvers, games in their order with file hashes, program and host |
| done.txt | one line per finished game, with all its observations |
| games/N.jsonl | the events of finished game N |
| sessions/N.jsonl | everything session N did, also trials of an interrupted game |
| results.jsonl, results.csv | the reports of all finished games, rebuilt after every session |
| status.json | how many games are finished |
| failures/ | the games and logs of failed trials |

A game is committed when it finishes. `--resume=DIR` (or `--experiment=DIR`
again, after a prompt on a terminal) continues with the saved plan; the only
option it takes besides is `--plain`. Finished games are kept, an interrupted
game is measured again from its warm-ups, and failures are not retried. The
resume is refused if the program was rebuilt, the host differs or an input
file changed. Only one benchmark at a time can work on an experiment.

A resumed experiment is not an uninterrupted one: the sessions are reported
apart. For measurements to publish, use an idle machine, run independent
experiments for replication, and report the coverage and the failures with
every total.
