# Tests Overview

This directory contains a small Fortran-based test harness and several
standalone test programs. There is no external test framework.

## Current Workflow (as-is)

Build the harness:

```
make -C tests
```

Run all tests (builds and executes the standalone programs too):

```
make -C tests run
```

Clean build artifacts and dump files:

```
make -C tests clean
```

Notes:
- The harness reads `tests/mendel.in`.
- The harness currently calls only `test_restart` from `test_main.f90`.
- The `run` target executes all `test_*.f90` programs.

Standalone examples (build manually):

```
mpif90 -O3 -I../src -I/usr/local/include -J../src -c ../src/random_pkg.f90
mpif90 -O3 -I../src -I/usr/local/include -J../src -o test_poisson test_poisson.f90 random_pkg.o
```

## What Tests Exist

- `test_main.f90`: harness; currently runs only `test_restart`.
- `test_restart.f90`: standalone program for restart I/O.
- `test_competition.f90`: standalone program for tribal competition logic.
- `test_inputs.f90`: reads and writes parameters.
- `test_poisson.f90`: checks Poisson RNG behavior.
- `test_polygenic.f90`: quick polygenic string/encoding scratchpad.

## Small Gains (Low Effort)

1) Make the harness run more tests without major refactors.
   - Add calls in `test_main.f90` for `test_inputs`-style checks or port logic
     into subroutines in `unittest.f90`.
2) Add Makefile convenience targets for standalone tests (done).
3) Add a `make run` target to build and execute everything (done).
4) Add a basic "smoke test" check to `test_restart` so it validates a small
   invariant beyond "no crash".

## Longer-Term Plan (Optional)

1) Consolidate tests into the harness:
   - Convert standalone programs into `subroutine test_*` functions and call
     them from `test_main.f90`.
2) Add a minimal runner:
   - Track pass/fail counts in `unittest.f90` instead of early `stop`.
3) Add deterministic tests:
   - Fix RNG seeds per test to make runs reproducible.
4) Add a CI-friendly "quick" config:
   - Smaller `mendel.in` with reduced populations and generations.
