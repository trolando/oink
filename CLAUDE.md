# Oink Agent Instructions

## Project

Oink is a high-performance C++ parity-game solver library and command-line tool. It implements many solver families, and some parallel solvers use Lace.

The goal of current work is incremental modernization, not a rewrite.

## Current architecture notes

Important files and directories:

- `include/oink/game.hpp`: main parity-game representation and public game API.
- `src/game.cpp`: implementation of `Game`.
- `include/oink/solver.hpp`: base solver interface.
- `src/solver.cpp`, `src/oink.cpp`, `src/solvers.cpp`: solver orchestration.
- `src/solvers/*.cpp`: individual algorithms.
- `include/oink/bitset.hpp`: custom dynamic bitset; treat as performance-sensitive.
- `include/oink/uintqueue.hpp`: custom queue; treat as performance-sensitive.
- `src/pgparser.cpp`: pgsolver parser.
- `src/verifier.cpp`: solution verification.
- `CMakeLists.txt`: build, tools, tests, and dependencies.
- `tests/`: parity-game test corpus.
- `test/`: test drivers.

## Modernization goal

Use C++ where it improves ownership, invariants, APIs, and maintainability.

Do not modernize for style alone. Preserve the performance-oriented character of the codebase.

Prefer:

- RAII ownership.
- `std::vector` for contiguous owned storage.
- Explicit immutable/mutable separation.
- Small helper types where they clarify invariants.
- Simple non-owning views or raw pointers in hot loops.
- Small, local refactors with tests.

Avoid:

- Wholesale rewrites.
- Rust or C migration.
- Template-heavy abstractions.
- `shared_ptr` ownership graphs.
- Exception-based solver control flow.
- `std::list`, `std::map`, or ranges-heavy code in hot solver paths unless justified.
- Replacing the custom `bitset` with `std::vector<bool>`.

## Refactoring rules

- Make the smallest coherent change.
- Do not refactor unrelated files opportunistically.
- Do not change solver algorithms unless the task explicitly asks for it.
- Keep public API compatibility where practical.
- Preserve existing raw pointer edge iteration until solvers are deliberately migrated.
- Build and test after each coherent change.
- If tests fail, stop and explain the failure before continuing.
- If a performance-sensitive representation changes, note the possible performance impact.

## Build commands

Use out-of-tree builds.

Debug configure:

```bash
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DOINK_BUILD_TESTS=ON
```

Debug build:

```bash
cmake --build build-debug -j
```

Run tests:

```bash
ctest --test-dir build-debug --output-on-failure
```

Release configure:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DOINK_BUILD_TESTS=ON
```

Release build:

```bash
cmake --build build-release -j
```

Release tests:

```bash
ctest --test-dir build-release --output-on-failure
```

## Git workflow

Use small commits.

Before editing:

```bash
git status --short
```

Never mix unrelated changes in one commit.

Commit messages must be short and scoped:

```text
scope: imperative summary
```

Examples:

```text
docs: add refactoring contract
cmake: add sanitizer option
game: initialize priorities
game: simplify label ownership
solver: add config struct
bitset: isolate intrinsics
tests: document disabled solvers
```

Rules:

- Keep subject line under 60 characters.
- Use lowercase scope.
- Use imperative mood.
- No trailing period.
- Do not include generated/co-authored footers unless asked.
- Do not commit unless explicitly told to commit.
- After each patch, show `git diff --stat` and propose the commit message.

## Task discipline

When given a roadmap phase, implement only that phase.

For each task:

1. Inspect relevant files.
2. State a concise implementation plan.
3. Make the smallest patch.
4. Build.
5. Run relevant tests.
6. Show the diff summary.
7. Propose one short commit message.

If a task is too broad, split it into numbered smaller commits and ask which one to do first.
