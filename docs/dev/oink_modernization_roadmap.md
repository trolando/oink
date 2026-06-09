# Oink Modernization Roadmap for Agentic Refactoring

## Objective

Modernize Oink by using C++ where it improves correctness, maintainability, ownership, and API clarity, while preserving the current performance-oriented character of the project.

The goal is **not** to rewrite Oink in “fancy modern C++.” The goal is to make the existing C++ codebase safer, cleaner, easier to extend, and easier to package, without disrupting solver performance or research usability.

## Repository

Repository:

```text
https://github.com/trolando/oink
```

Assume the current codebase is a high-performance C++ parity-game solver library/tool with many solver implementations. The code currently uses a mixture of C++ abstractions and C-style storage, especially in core graph and solver data structures.

## High-Level Principles

Follow these principles throughout the work:

1. Preserve solver behavior.
2. Preserve performance unless a regression is explicitly justified.
3. Prefer small, reviewable PRs/commits.
4. Modernize ownership before changing algorithms.
5. Keep hot loops simple and benchmarkable.
6. Do not introduce template-heavy abstractions in solver internals.
7. Do not rewrite all solvers at once.
8. Do not migrate to Rust or C.
9. Do not bump the C++ standard beyond C++17 unless explicitly approved.
10. Treat Oink as a research-grade high-performance library, not as a toy project.

## Non-Goals

Do **not** do the following unless explicitly instructed:

- Do not perform a wholesale rewrite.
- Do not replace all pointer-based iteration in hot solver loops immediately.
- Do not replace the custom bitset with `std::vector<bool>`.
- Do not introduce `shared_ptr`-based ownership.
- Do not use exceptions for normal solver control flow.
- Do not make solver code “clever” with heavy templates, ranges, or metaprogramming.
- Do not refactor all solver families in one pass.
- Do not remove Lace support.
- Do not remove Boost.Iostreams compression support without checking its usage.

## Target End State

The desired end state is:

```text
- Game structure is owned through RAII containers.
- Game structure and solver output are separated.
- Solvers receive immutable game data and explicit mutable solution/scratch state.
- Raw owning pointers are removed from core public types.
- Raw non-owning pointers/views remain allowed in performance-critical loops.
- Parsing/building is separated from immutable game representation.
- CI catches memory errors, undefined behavior, and basic regressions.
- Packaging and CMake targets are cleaner.
```

## Working Method

For each phase:

1. Inspect the relevant files first.
2. Make the smallest useful change.
3. Build.
4. Run tests.
5. Run sanitizer tests where applicable.
6. Compare behavior against the existing implementation.
7. Avoid large unrelated cleanups.
8. Document any intentional behavior change.

Prefer commit messages like:

```text
game: replace label ownership with vector storage
solution: introduce explicit Solution class
cmake: add sanitizer build options
bitset: move compiler intrinsics behind portability helpers
```

---

# Phase 0 — Establish Refactoring Contract

## Goal

Create a clear contract for future refactors.

## Tasks

Create:

```text
docs/dev/refactoring-contract.md
```

Suggested content:

```markdown
# Oink Refactoring Contract

Refactors must preserve solver behavior and performance unless explicitly stated.

## Correctness

- Existing tests must pass.
- Solver output must verify successfully.
- Any changed solver behavior must be documented.
- Disabled tests must be explicitly linked to a known issue or TODO.

## Performance

- Hot solver loops may use raw pointers or simple index loops.
- Owning raw pointers should be removed from public/core types.
- Performance regressions above 3–5% on representative benchmarks require explanation.
- Micro-optimizing before correctness cleanup is discouraged.

## C++ Usage

Use C++ for:

- RAII ownership
- `std::vector` contiguous storage
- move semantics
- explicit immutable/mutable separation
- `enum class` for small closed domains
- local view/span abstractions
- safer construction APIs

Avoid C++ for:

- template-heavy solver abstractions
- `std::list`/`std::map` in hot paths unless measured
- `shared_ptr` ownership graphs
- exception-driven solver logic
- ranges-heavy inner loops
```

## Acceptance Criteria

- Document exists.
- It is referenced from the main developer documentation or README if appropriate.
- No code behavior changes in this phase.

---

# Phase 1 — Safety and CI Baseline

## Goal

Make future refactors safe by improving diagnostics and test coverage.

## Tasks

### 1. Add sanitizer options

Add CMake options:

```cmake
option(OINK_ENABLE_SANITIZERS "Enable address and undefined behavior sanitizers" OFF)
option(OINK_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
```

For GCC/Clang, when sanitizers are enabled:

```text
-fsanitize=address,undefined
-fno-omit-frame-pointer
-g
```

Do not force sanitizers on by default.

### 2. Add CI jobs

Add or update CI jobs for:

```text
- Debug build
- Release build
- Debug + ASan/UBSan
- GCC
- Clang
```

Do not add MSVC yet unless the code is already close to supporting it.

### 3. Fix known memory initialization problems

Look for manually allocated arrays that are not initialized before reads.

Pay special attention to fields like:

```text
_priority
_owner
_strategy
_solved
_winner
_outedges
_inedges
_label
```

Prefer RAII/zero-initializing containers in later phases, but fix known initialization bugs immediately.

### 4. Classify disabled tests

Search CMake and test files for disabled/commented tests.

For each disabled test:

```text
- Re-enable if fixed.
- Otherwise add a clear TODO with reason.
- Link to an issue if possible.
```

### 5. Add malformed-input parser tests

Add small parser tests for:

```text
- empty file
- missing header
- invalid vertex id
- invalid owner
- invalid priority
- malformed edge list
- duplicate or out-of-range vertex
```

## Acceptance Criteria

- Existing tests pass.
- Sanitizer build completes.
- Known uninitialized read issues are fixed or documented.
- Disabled tests are no longer silently ignored.

---

# Phase 2 — Separate Game from Solution

## Goal

Separate immutable parity-game data from mutable solver output.

Currently, `Game` likely contains both game structure and solution-related state such as solved vertices, winners, and strategies. This makes it harder to reuse the same game across solvers and blurs invariants.

## New Types

Introduce a `Solution` class.

Possible header:

```text
include/oink/solution.hpp
```

Suggested shape:

```cpp
#pragma once

#include <vector>
#include "oink/bitset.hpp"

namespace pg {

using Vertex = int;

enum class Player : unsigned char {
    Even = 0,
    Odd = 1
};

class Solution {
public:
    explicit Solution(int vertex_count = 0);

    void resize(int vertex_count);
    void reset();

    int vertex_count() const noexcept;

    bool is_solved(Vertex v) const noexcept;
    Player winner(Vertex v) const noexcept;
    Vertex strategy(Vertex v) const noexcept;

    void solve(Vertex v, Player winner, Vertex strategy);
    void set_winner(Vertex v, Player winner);
    void set_strategy(Vertex v, Vertex strategy);

private:
    bitset solved_;
    bitset winner_;
    std::vector<Vertex> strategy_;
};

} // namespace pg
```

Adjust names/namespaces to match the existing project.

## Migration Strategy

Do **not** remove old `Game` solution methods immediately.

Instead:

1. Add `Solution`.
2. Let `Game` temporarily forward old solution methods to an internal `Solution`.
3. Gradually update solvers to take/use `Solution`.
4. Once solvers are migrated, remove solution state from `Game`.

## Solver Interface Target

Eventually move toward:

```cpp
class Solver {
public:
    Solver(Oink& oink, const Game& game, Solution& solution);
    virtual ~Solver() = default;

    virtual void run() = 0;

protected:
    Oink& oink;
    const Game& game;
    Solution& solution;
};
```

## Acceptance Criteria

- `Solution` exists.
- Existing solver behavior is unchanged.
- Old APIs still compile for now.
- At least one nontrivial solver uses `Solution` directly.

---

# Phase 3 — Replace Owning Raw Arrays in Game

## Goal

Move `Game` ownership to RAII containers while preserving memory layout and hot-loop access patterns.

## Current Problem

The current `Game` type appears to own many arrays manually, such as priorities, edges, labels, counts, solved/winner/strategy state, and adjacency offsets.

Owning raw pointers are error-prone:

```cpp
int* _priority;
int* _outedges;
int* _firstouts;
std::string** _label;
```

The target is not to eliminate raw pointers entirely. The target is to eliminate **owning raw pointers**.

## Proposed Internal Storage

Use something like:

```cpp
class Game {
public:
    int vertex_count() const noexcept;
    int edge_count() const noexcept;

    int priority(Vertex v) const noexcept;
    Player owner(Vertex v) const noexcept;

    const int* outs_raw(Vertex v) const noexcept;
    const int* ins_raw(Vertex v) const noexcept;

    int out_count(Vertex v) const noexcept;
    int in_count(Vertex v) const noexcept;

private:
    std::vector<int> priority_;
    bitset owner_;

    std::vector<int> first_out_;
    std::vector<int> out_count_;
    std::vector<int> out_edges_;

    std::vector<int> first_in_;
    std::vector<int> in_count_;
    std::vector<int> in_edges_;

    std::vector<std::string> labels_;
};
```

## Important Compatibility Rule

Initially preserve existing sentinel-based edge iteration if the solvers depend on it.

For example:

```cpp
const int* outs_raw(Vertex v) const noexcept {
    return out_edges_.data() + first_out_[v];
}
```

Existing loops of the form:

```cpp
for (const int* to = game.outs(v); *to != -1; ++to) {
    ...
}
```

may continue to work during the transition.

## Label Storage

Replace nullable heap-allocated labels with one of:

### Preferred simple version

```cpp
std::vector<std::string> labels_;
```

Interpret empty string as no label.

### More explicit version

```cpp
std::vector<std::optional<std::string>> labels_;
```

Use this only if empty labels are meaningful.

## Migration Strategy

1. Convert labels first.
2. Convert priorities/owners.
3. Convert edge offsets/counts.
4. Convert edge arrays.
5. Remove obsolete destructor logic.
6. Remove `malloc`, `calloc`, `realloc`, and `free` from `Game`.

## Acceptance Criteria

- `Game` no longer manually frees its primary storage.
- Existing solver loops still compile.
- Existing tests pass.
- Sanitizer build passes.
- No meaningful performance regression on small benchmark corpus.

---

# Phase 4 — Introduce GameBuilder

## Goal

Separate mutable construction/parsing from immutable game usage.

## Problem

`Game` currently appears to mix:

```text
- storage
- mutation
- parsing
- writing
- edge construction
- sorting
- priority compression
- subgame extraction
- solution state
```

A builder clarifies the lifecycle.

## New Type

Add:

```text
include/oink/game_builder.hpp
src/game_builder.cpp
```

Suggested shape:

```cpp
class GameBuilder {
public:
    explicit GameBuilder(int vertex_count);

    void set_priority(Vertex v, int priority);
    void set_owner(Vertex v, Player owner);
    void set_label(Vertex v, std::string label);

    void add_edge(Vertex from, Vertex to);

    Game build();

private:
    int vertex_count_;
    std::vector<int> priority_;
    bitset owner_;
    std::vector<std::string> labels_;
    std::vector<std::vector<Vertex>> outgoing_;
};
```

This builder does not have to be maximally optimized. It is used during construction/parsing, not in solver hot loops.

## Parser Migration

Route parser code through `GameBuilder`:

```cpp
Game parse_pgsolver(std::istream& in) {
    GameBuilder builder(n);
    ...
    return builder.build();
}
```

## Preserve Old APIs Temporarily

If existing code uses methods like:

```cpp
e_start(...)
e_add(...)
e_finish(...)
```

keep them temporarily and implement them internally using the builder or a compatibility path.

Mark as deprecated only after the new builder path is stable.

## Acceptance Criteria

- Parser can construct a `Game` through `GameBuilder`.
- Existing parser behavior is unchanged.
- Construction-time validation improves where easy.
- `Game` becomes closer to immutable after construction.

---

# Phase 5 — Add Lightweight Views for Edge Access

## Goal

Provide safer non-owning access without hurting hot loops.

## C++17 Constraint

Since Oink currently targets C++17, do not use `std::span` unless the project moves to C++20.

Options:

1. Add a tiny local span-like type.
2. Use a pair of pointer + size.
3. Keep raw pointer API and add count accessors.

## Suggested Minimal Type

```cpp
template <typename T>
class Span {
public:
    Span(T* data, std::size_t size) noexcept
        : data_(data), size_(size) {}

    T* data() const noexcept { return data_; }
    std::size_t size() const noexcept { return size_; }

    T* begin() const noexcept { return data_; }
    T* end() const noexcept { return data_ + size_; }

    T& operator[](std::size_t i) const noexcept { return data_[i]; }

private:
    T* data_;
    std::size_t size_;
};
```

Then expose:

```cpp
Span<const Vertex> out_edges(Vertex v) const noexcept;
Span<const Vertex> in_edges(Vertex v) const noexcept;
```

Keep:

```cpp
const Vertex* outs_raw(Vertex v) const noexcept;
const Vertex* ins_raw(Vertex v) const noexcept;
```

## Acceptance Criteria

- New view API exists.
- Old raw API still exists.
- At least one non-hot or easy solver/helper uses the view API.
- Hot solvers do not need to be rewritten yet.

---

# Phase 6 — Introduce Domain Types Carefully

## Goal

Use C++ types to make invalid states harder to express.

## First Step: Player

Introduce:

```cpp
enum class Player : unsigned char {
    Even = 0,
    Odd = 1
};
```

Provide helpers:

```cpp
inline Player opponent(Player p) noexcept {
    return p == Player::Even ? Player::Odd : Player::Even;
}

inline int player_index(Player p) noexcept {
    return static_cast<int>(p);
}
```

Use this at API boundaries first:

```cpp
Player owner(Vertex v) const noexcept;
void set_owner(Vertex v, Player owner);
void solve(Vertex v, Player winner, Vertex strategy);
```

Do not attempt to replace all internal `int` player variables at once.

## Second Step: Aliases

Use aliases initially:

```cpp
using Vertex = int;
using Priority = int;
```

Do not introduce a strong `Vertex` wrapper yet. That will touch too much code.

## Acceptance Criteria

- `Player` exists.
- Public APIs begin to use `Player`.
- Internal hot loops may still use integers where conversion would be noisy.
- No solver behavior changes.

---

# Phase 7 — Modernize bitset Internals

## Goal

Keep the custom bitset, but make its ownership and portability cleaner.

## Do Not

Do not blindly replace the custom bitset with:

```cpp
std::vector<bool>
boost::dynamic_bitset
std::set
std::unordered_set
```

The custom bitset likely exists for performance reasons.

## Refactor Direction

Replace manual ownership with:

```cpp
std::vector<uint64_t> bits_;
std::size_t size_;
```

Expose low-level access where needed:

```cpp
uint64_t* data() noexcept;
const uint64_t* data() const noexcept;
std::size_t block_count() const noexcept;
```

## Intrinsics Isolation

Move compiler-specific operations into a helper header/source:

```text
include/oink/intrinsics.hpp
src/intrinsics.cpp
```

Provide:

```cpp
namespace pg::intrinsics {

int popcount64(uint64_t x) noexcept;
int countr_zero64(uint64_t x) noexcept;
int countl_zero64(uint64_t x) noexcept;

}
```

Use compiler builtins internally:

```text
GCC/Clang: __builtin_popcountll, __builtin_ctzll, __builtin_clzll
MSVC later: _popcnt64, _tzcnt_u64, __lzcnt64 or safe fallbacks
```

Always handle zero arguments safely where required.

## Acceptance Criteria

- `bitset` owns memory through RAII.
- Existing bitset API mostly remains.
- Compiler-specific code is isolated.
- Existing tests pass.
- Bitset-heavy benchmarks do not regress significantly.

---

# Phase 8 — Clarify Solver Context and Scratch State

## Goal

Make solver dependencies explicit.

## Target Structure

Solvers should eventually have access to:

```text
- const Game& game
- Solution& solution
- SolverConfig config
- SolverScratch or algorithm-local state
- Oink/logging/orchestration only where needed
```

Introduce:

```cpp
struct SolverConfig {
    int workers = 0;
    bool trace = false;
    bool preprocess = true;
    std::optional<uint64_t> seed;
};
```

Avoid making `Oink` a catch-all dependency.

## Solver Refactor Order

Do not start with the most complicated parallel/tangle solvers.

Suggested order:

```text
1. verifier / simple helper code
2. one simple sequential solver
3. Zielonka
4. FPI
5. progress-measure solvers
6. priority-promotion solvers
7. tangle-learning solvers
8. parallel solvers
```

## Acceptance Criteria

- `SolverConfig` exists.
- At least one solver uses explicit config/state cleanly.
- No solver-family-wide rewrite yet.
- Existing CLI behavior remains unchanged.

---

# Phase 9 — Add SubgameView

## Goal

Avoid unnecessary game copying and clarify active-subgame operations.

## New Type

Add:

```cpp
class SubgameView {
public:
    SubgameView(const Game& game, const bitset& active);

    const Game& game() const noexcept;
    bool contains(Vertex v) const noexcept;
    const bitset& active_vertices() const noexcept;

private:
    const Game* game_;
    const bitset* active_;
};
```

Use it first in non-invasive places such as preprocessing, attractor computation, or SCC-related helper code.

## Do Not

Do not remove existing `extract_subgame` immediately.

Some algorithms may legitimately benefit from extracted/reindexed games.

## Acceptance Criteria

- `SubgameView` exists.
- One helper/preprocessor uses it.
- Existing extraction code remains available.
- No broad solver rewrite.

---

# Phase 10 — Reduce Boost Dependencies Where Easy

## Goal

Reduce dependency friction without wasting time.

## Investigate Uses

Search for:

```text
boost::filesystem
boost::regex
boost::random
boost::system
boost::iostreams
```

## Likely Replacements

Prefer:

```text
Boost.Filesystem -> std::filesystem
Boost.Regex      -> std::regex or simpler parsing
Boost.Random     -> <random>
Boost.System     -> may disappear when filesystem disappears
Boost.Iostreams  -> keep for gzip/bzip2 unless replacing deliberately
```

## Rules

- Do not remove Boost.Iostreams casually.
- Do not break compressed input support.
- Do not introduce a larger dependency to remove a smaller one.
- Keep changes separate from core solver refactors.

## Acceptance Criteria

- At least one obsolete Boost dependency is removed if straightforward.
- CMake dependency list is simplified.
- Existing compressed input behavior is preserved.

---

# Phase 11 — Improve CMake Packaging

## Goal

Make Oink easier to consume as a library.

## Tasks

1. Ensure install targets are correct.
2. Export CMake package targets.
3. Support downstream usage:

```cmake
find_package(oink CONFIG REQUIRED)
target_link_libraries(mytool PRIVATE oink::oink)
```

4. Make CPU-specific optimization opt-in for packaged builds.

If currently defaulting to native CPU optimization, change to:

```cmake
option(OINK_ENABLE_NATIVE "Build with native CPU optimizations" OFF)
```

For local benchmarking, users can enable:

```text
-DOINK_ENABLE_NATIVE=ON
```

## Suggested Options

```cmake
option(OINK_BUILD_TOOLS "Build Oink command-line tools" ON)
option(OINK_BUILD_TESTS "Build Oink tests" ON)
option(OINK_BUILD_BENCHMARKS "Build Oink benchmarks" OFF)
option(OINK_ENABLE_NATIVE "Build with native CPU optimizations" OFF)
option(OINK_ENABLE_LTO "Enable link-time optimization" OFF)
option(OINK_ENABLE_SANITIZERS "Enable sanitizers" OFF)
```

## Acceptance Criteria

- Project builds as before.
- Install/export works.
- Native CPU optimization is not forced on downstream packagers.
- Tools/tests/benchmarks can be toggled independently.

---

# Phase 12 — Naming and Public API Cleanup

## Goal

Make the public API easier to understand without unnecessary churn.

## Rules

- Do this after structural refactors.
- Keep compatibility wrappers where reasonable.
- Avoid renaming internals purely for style.
- Rename public concepts only when the new name is clearly better.

## Suggested API Direction

Prefer:

```cpp
vertex_count()
edge_count()
priority(v)
owner(v)
out_edges(v)
in_edges(v)
out_count(v)
in_count(v)
```

over ambiguous or historical names.

If old names exist, keep wrappers temporarily:

```cpp
int nodecount() const noexcept {
    return vertex_count();
}
```

Later mark as deprecated:

```cpp
[[deprecated("Use vertex_count()")]]
int nodecount() const noexcept;
```

## Acceptance Criteria

- Public API has clearer naming.
- Old names still work during transition.
- Migration path is documented.

---

# Suggested PR Sequence

## Milestone 1 — Safety Baseline

```text
PR 1: Add refactoring contract document.
PR 2: Add sanitizer CMake option and CI job.
PR 3: Fix known initialization issues.
PR 4: Classify disabled tests and add malformed parser tests.
```

## Milestone 2 — Ownership Cleanup

```text
PR 5: Replace label pointer storage with vector storage.
PR 6: Introduce Solution class.
PR 7: Move Game solution state behind Solution forwarding.
PR 8: Replace Game-owned raw arrays with vector-backed storage.
```

## Milestone 3 — Construction and Access APIs

```text
PR 9: Introduce GameBuilder.
PR 10: Route parser through GameBuilder.
PR 11: Add edge view/span API while preserving raw pointer API.
PR 12: Introduce Player enum class at API boundaries.
```

## Milestone 4 — Core Utility Cleanup

```text
PR 13: Modernize bitset ownership.
PR 14: Move intrinsics into portability helpers.
PR 15: Add SolverConfig.
PR 16: Introduce SubgameView.
```

## Milestone 5 — Library Quality

```text
PR 17: Reduce easy Boost dependencies.
PR 18: Improve CMake install/export.
PR 19: Make native CPU optimization opt-in.
PR 20: Add release/versioning hygiene.
```

## Milestone 6 — Solver-by-Solver Refactors

```text
PR 21+: Refactor one solver family at a time.

Order:
1. verifier/helper code
2. simple sequential solver
3. Zielonka
4. FPI
5. progress-measure solvers
6. priority-promotion solvers
7. tangle-learning solvers
8. parallel solvers
```

---

# Validation Commands

Adjust these to the actual repository layout.

## Configure Debug

```bash
cmake -S . -B build-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DOINK_BUILD_TESTS=ON
```

## Build Debug

```bash
cmake --build build-debug -j
```

## Run Tests

```bash
ctest --test-dir build-debug --output-on-failure
```

## Configure Sanitizers

```bash
cmake -S . -B build-asan \
  -DCMAKE_BUILD_TYPE=Debug \
  -DOINK_BUILD_TESTS=ON \
  -DOINK_ENABLE_SANITIZERS=ON
```

## Build Sanitizers

```bash
cmake --build build-asan -j
```

## Run Sanitizer Tests

```bash
ctest --test-dir build-asan --output-on-failure
```

## Configure Release

```bash
cmake -S . -B build-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DOINK_BUILD_TESTS=ON
```

## Build Release

```bash
cmake --build build-release -j
```

## Run Release Tests

```bash
ctest --test-dir build-release --output-on-failure
```

---

# Benchmarking Guidance

Before and after major representation changes, compare performance on a fixed corpus.

Create or identify a corpus such as:

```text
benchmarks/small/
benchmarks/medium/
benchmarks/hard/
```

For each benchmark, record:

```text
- input file
- solver
- wall-clock time
- user CPU time
- memory usage if easy
- output verification result
```

Use a simple script such as:

```bash
scripts/benchmark-smoke.sh
```

This does not need to be a full scientific benchmarking framework. It only needs to catch obvious regressions.

Do not block all refactors on tiny noise-level changes. Treat regressions above roughly 3–5% as requiring investigation.

---

# Files to Inspect First

Start by inspecting these likely core files:

```text
include/oink/game.hpp
src/game.cpp

include/oink/solver.hpp
src/solver.cpp

include/oink/bitset.hpp
include/oink/uintqueue.hpp

src/oink.cpp
include/oink/oink.hpp

CMakeLists.txt
.github/workflows/*
```

Then inspect representative solvers:

```text
src/solvers/*
include/oink/solvers/*
```

Choose one simple sequential solver as the first solver migration target.

---

# Final Desired Architecture

A reasonable final architecture is:

```text
Game
  Immutable parity-game structure:
  - priorities
  - owners
  - outgoing edges
  - incoming edges
  - labels

GameBuilder
  Mutable construction API:
  - parser uses this
  - tests can use this
  - validates construction

Solution
  Mutable solver result:
  - solved vertices
  - winners
  - strategies

Solver
  Abstract solver interface:
  - const Game&
  - Solution&
  - SolverConfig
  - algorithm-local scratch

Oink
  Orchestration:
  - solver selection
  - logging
  - preprocessing coordination
  - CLI/library integration

SubgameView
  Non-owning active-subgame view:
  - const Game&
  - active bitset

bitset
  Custom high-performance bitset:
  - RAII-owned storage
  - compiler intrinsics isolated
```

The result should still be recognizably Oink: compact, fast, and suitable for parity-game research. The main difference should be that ownership, construction, and solver state are explicit instead of implicit in raw pointers and lifecycle conventions.
