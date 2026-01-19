# Repository Guidelines

## Project Structure & Module Organization
- Fortran source lives in `src/` (e.g., `src/mendel.f90`, `src/inputs.f90`, `src/genome.f90`). Shared headers live in `include/` (`include/common.h`).
- Tests and test harness sources live in `tests/` (e.g., `test_main.f90`, `unittest.f90`, `test_*.f90`).
- Input templates and runtime files for the SPC UI are in `app/` (e.g., `app/mendel.in`, `app/mendel.j2`, `app/spc.json`, `app/about.html`).
- Build and container tooling remain at the repo root: `Makefile`, `tests/Makefile`, `Dockerfile`.

## Build, Test, and Development Commands
- `make` or `make release`: build the default parallel executable `mendel` using `mpif90`.
- `make debug`: build with debug flags and checks.
- `make preserial` then `make serial`: build `mendel_serial` without MPI code sections.
- `make test`: build a test executable named `test` from core sources.
- `make -C tests`: build the dedicated test harness `tests/test_main`.
- `make dist`: build SPC deployment packages into `dist/` (requires binaries in `bin/*/mendel`).
- `docker build . -t mendel:2.7.3` and `docker run -p 8580:8580 mendel:2.7.3`: run the legacy SPC UI.

## Coding Style & Naming Conventions
- Fortran 90 free-form style; follow existing conventions (lowercase keywords, 3-space indentation inside blocks).
- Module and file names are descriptive and lowercase (e.g., `polygenic.f90`, `selection.f90`).
- Keep interfaces in `common.h` and module declarations in corresponding `.f90` files.

## Testing Guidelines
- Tests are Fortran-based with a custom `unittest.f90` harness; no external framework.
- Test files are named `test_*.f90` or live in `tests/` with `test_main.f90` as the entry point.
- Run `make -C tests` to compile the harness; run `./tests/test_main` manually.

## Commit & Pull Request Guidelines
- Git history uses short, imperative, lowercase summaries (e.g., “fix issue…”, “output gen…”). Keep messages under ~70 chars.
- PRs should describe the change, reference any issues, and note how tests were run.
- Include input/output examples or relevant `mendel.in` snippets when behavior changes.

## Configuration Tips
- Runtime configuration is primarily via `app/mendel.in`. Start from `app/mendel.j2` or `tests/mendel.in` when adding new parameters.
- MPI support is optional; ensure serial builds still work when editing MPI sections.
