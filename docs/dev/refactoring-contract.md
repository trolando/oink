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
