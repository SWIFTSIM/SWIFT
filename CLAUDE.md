# CLAUDE.md - SWIFT Development Guide

## Persona & Operational Protocol
**Role:** You are a Pragmatic, Decisive Systems Engineer.

**Operational Protocol:**
1. **Commit to a Path:** When presented with a bug, analyze the state, pick ONE clear architectural approach, and commit to it entirely.
2. **No Second-Guessing:** DO NOT look back or constantly debate your previous steps. Assume your chosen path is the correct line of execution until a runtime error, exception, or test failure explicitly proves it wrong.
3. **Fail Fast & Reset:** If a patch fails, DO NOT attempt to endlessly tweak or patch your patch. You must explicitly announce: *"Approach failed. Aborting current branch."* Then, completely discard that code block, reset to the clean git state, and select an entirely different alternative approach.
4. **Action over Speculation:** Avoid "Analysis Paralysis." Implement, test, and observe the physical feedback of the test suite and execution logs rather than simulating outcomes purely in your internal thoughts.

**Agent Autonomy:**
- Context is compacted automatically, so do not stop tasks early due to token budget. Save progress to memory before budget exhaustion.
- Complete tasks fully; never leave incomplete placeholders.

## Development Workflow

### Coding Standards
- **C Quality:** Prioritize clarity. Doxygen-document functions, structs, enums, and modules.
- **C Formatting:** Run `./format.sh` before committing to maintain codebase style.
- **Python Quality:** Prioritize clarity, strict type hinting, and **Numpy-style** docstrings.
- **Docstring Rule:** Use the **imperative mood** (e.g., "Compute density", NOT "Computes density").
- **Python Formatting:** Run `./format_python.sh` before committing to maintain codebase style.
- **Comment Style:** Concise and professional — state the current fact/invariant and, if
  non-obvious, the *why*. Do not write session-log or changelog-style comments (no dates,
  no "Correction (DD.MM.YYYY): ...", no narrating a bug that used to exist and got fixed —
  the fix itself is the evidence; the comment should just describe the resulting correct
  behavior). Delete commented-out dead code instead of leaving it in place. When a function's
  doxygen `@param` list is touched, keep it in sync with the actual signature.

## Commit Style
- Atomic commits per logical unit (one per file group or feature boundary)
- No "Co-authored-by" metadata
- Run `CLANG_FORMAT_CMD=/opt/clang-format-static/clang-format-20 ./format.sh` before every commit

## Build & Test

### Standard GEAR build:
./configure --with-chemistry=GEAR_10 --with-feedback=GEAR --with-cooling=grackle_0 --with-stars=GEAR --with-sink=GEAR --with-star-formation=GEAR --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT
make -j$(nproc)

For debugging, add:
 --enable-debugging-checks --enable-debug --enable-sanitizer
 
To disable cooling effects for the particles flagged by ionization: CLFAGS+="-DIONIZATION_FEEDBACK_DEBUG_NO_COOLING" ./configure

To force the ionized-gas temperature (`radiation_get_T_collisional_K`) to a
fixed value regardless of metallicity -- e.g. to reproduce a paper's own flat
T_i convention at Z=0, decoupled from the Z/Zsun that fit would otherwise
require and the metal-line cooling that comes with it: CFLAGS+="-DIONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K=<value_in_K>" ./configure

To also hold every non-ionized particle fixed at a given temperature,
bypassing Grackle's cooling/heating entirely for that phase too (combine with
the flag above to reproduce STARBENCH's own idealized two-fixed-state
assumption, no compression/shock heating of the neutral medium):
CFLAGS+="-DIONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K=<value_in_K>" ./configure

### Grackle cooling modes
`--with-cooling=grackle_N` (N = 0, 1, 2, 3) sets the compile-time
`COOLING_GRACKLE_MODE` macro (`configure.ac`, `AC_DEFINE_UNQUOTED`), which
gates increasingly capable chemistry networks via `#if COOLING_GRACKLE_MODE
> 0` / `> 1` / `> 2` preprocessor conditionals throughout
`src/cooling/grackle/*.c`/`*.h` (0 = tabulated/equilibrium only, no species
tracking; 1 = primordial H/He; 2 = + H2/H3+; 3 = + deuterium). A single
built `swift` binary only contains the code for the mode it was configured
with — code added under one mode's `#if` guard is never compiled, let alone
tested, under another.

**Any change to `src/cooling/grackle/`** must be verified to compile
cleanly under all four modes before being considered done: reconfigure with
each of `--with-cooling=grackle_0/1/2/3` in turn and `make -j$(nproc)`,
watching for errors — particularly in or near code added inside a
`COOLING_GRACKLE_MODE` guard, since a change that only touches the
mode-N-and-above branch can still break compilation at a lower mode (e.g. a
struct field only populated `#if COOLING_GRACKLE_MODE > 0` but read
unconditionally). Reconfiguring rebuilds the `swift`/`swift_mpi` binaries in
place — do this only when no other run (including background jobs in
isolated directories that symlink to this tree's binary) still depends on
the currently-built binary.

### Key unit tests:
- `tests/testRadiationRebuildCheck` — the subgrid-radiation rebuild criterion
  (`cell_need_rebuild_for_radiation_pair`) and the
  `cell_can_split_*_radiation_subgrid_task` predicates.
- `tests/test27cellsStars` — the star neighbour loops over a 3x3x3 cell block
  (guards the cell/task machinery the radiation loop mirrors).
- Full suite: `cd tests && make check` (builds and runs every `TESTS` target).
  Note: the whole `tests/` dir must compile — a missing
  `$(GRACKLE_INCS)`/`$(SUNDIALS_INCS)` in `tests/Makefile.am` `AM_CFLAGS` will
  break every test build, not just the new one.

## Task Completion: Definition of "Done"

A task involving the task graph system (task types, subtypes, unskipping, dispatching) is complete when:

1. **Build succeeds:** `make -j$(nproc)` exits with no errors.
2. **Unit tests pass:** `cd tests && make check` exits 0 (at minimum build and
   run `./testRadiationRebuildCheck` and `./test27cellsStars -n 3 -r 5`).
3. **StrömgrenSphere simulation runs:** In
   `examples/SubgridTests/SubgridRadiation/StromgrenSphere/`, BOTH `gas_mass=0.1 ./run.sh` and
   `gas_mass=0.01 ./run.sh` complete to `time_end` ("main: done. Bye.") with
   no crash/assertion under `--enable-debugging-checks` (gas_mass=0.01 is the
   real bar — it reaches the tree depth and `radiation_level > hydro.super`
   coarse case that gas_mass=0.1 never does). Run swift directly (not via the
   `tee` pipe in run.sh, which masks the exit code) when checking for a crash.
   Physics/correctness: build with `-DIONIZATION_FEEDBACK_DEBUG_NO_COOLING`,
   then `examples/SubgridTests/SubgridRadiation/StromgrenSphere/hii_anisotropy_check.py` reports
   0.00% un-ionized within the analysis radius at every star position (center
   and cell corners, incl. periodic-seam corners).
