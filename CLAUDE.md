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
- **Example READMEs:** Written for an external computational astrophysicist who wants to
  configure and run the example — not for us, and not a lab notebook. Concise, self-contained,
  no "bible length" documents. Must cover: what it reproduces (a one-line pointer to the
  paper is fine), how to configure, how to run. Must NOT cover: validation status, what was
  investigated or bisected, comparisons against deprecated/old results, or anything that
  assumes the reader watched this project's own debugging happen. That content belongs in
  `theory/GEAR/Radiation/` or the agent's own memory, never in the shipped README. Before
  writing or editing one, re-read a sibling README (e.g. `StromgrenSphereClump/README`) for
  length and tone calibration. This has been a repeat problem — take it seriously.

## Commit Style
- Atomic commits per logical unit (one per file group or feature boundary)
- No "Co-authored-by" metadata
- Run `CLANG_FORMAT_CMD=/opt/clang-format-static/clang-format-20 ./format.sh` before every commit

## Build & Test

### Standard GEAR build:
./configure --with-chemistry=GEAR_10 --with-feedback=GEAR --with-cooling=grackle_0 --with-stars=GEAR --with-sink=GEAR --with-star-formation=GEAR --with-tracers=GEAR --with-kernel=wendland-C2 --with-grackle=$GRACKLE_ROOT
make -j$(nproc)

`--with-tracers=GEAR` is required, not optional: the HII ionization tag lives
in `tracers_xpart_data`, so `src/feedback/GEAR/radiation.c` does not compile
without it.

For debugging, add `--enable-debugging-checks --enable-debug --enable-sanitizer`.

GEAR radiation debug flags, all via `CFLAGS+="..." ./configure`:
`-DIONIZATION_FEEDBACK_DEBUG_NO_COOLING` (tagged gas never cools, so tags never
expire), `-DIONIZATION_FEEDBACK_DEBUG_FIXED_IONIZED_TEMPERATURE_K=<K>` (flat
T_i regardless of Z), `-DIONIZATION_FEEDBACK_DEBUG_FIXED_NEUTRAL_TEMPERATURE_K=<K>`
(also pins neutral gas; with the previous flag this is STARBENCH's idealized
two-state setup).

### Grackle cooling modes
`--with-cooling=grackle_N` sets `COOLING_GRACKLE_MODE` (0 = tabulated, no
species; 1 = H/He; 2 = + H2; 3 = + D), which gates code via `#if
COOLING_GRACKLE_MODE > N` throughout `src/cooling/grackle/`. A binary contains
only its own mode's code. **Any change under `src/cooling/grackle/` must be
compile-checked under all four modes** — a mode-gated line is otherwise never
compiled, let alone tested. Reconfiguring rebuilds `swift` in place, so don't
do it while a run depends on the current binary.

### Key unit tests
- `tests/testRadiationRebuildCheck` — `cell_need_rebuild_for_radiation_pair`
  and the `cell_can_split_*_radiation_subgrid_task` predicates.
- `tests/test27cellsStars` — star neighbour loops over a 3x3x3 cell block.
- Full suite: `cd tests && make check`.

## Task Completion: Definition of "Done"

A task touching the task graph (types, subtypes, unskipping, dispatch) is done when:

1. **Builds:** `make -j$(nproc)` clean.
2. **Unit tests:** `./testRadiationRebuildCheck` and
   `./test27cellsStars -n 3 -N 3 -r 5` both exit 0 (`-N` is required; omitting
   it exits 1 on usage, not on a real failure). `make check` as a whole cannot
   complete when configure picks `-march=core2 -mavx2` — `testKernel.c` needs
   FMA, which `core2` lacks. Pre-existing and unrelated; run the two directly
   and say so.
3. **StrömgrenSphere runs:** in
   `examples/SubgridTests/SubgridRadiation/StromgrenSphere/`, both
   `gas_mass=0.1 ./run.sh` and `gas_mass=0.01 ./run.sh` reach `time_end` under
   `--enable-debugging-checks`. `gas_mass=0.01` is the real bar (it alone hits
   `radiation_level > hydro.super`). Invoke swift directly — run.sh's `tee`
   pipe masks the exit code. Then build with
   `-DIONIZATION_FEEDBACK_DEBUG_NO_COOLING` and run
   `hii_anisotropy_check.py -s 'snap/snapshot_*.hdf5'` (a glob, not one
   snapshot: the star dies at the end and `h_hii` is legitimately 0 there).
   Pass bar, per resolution:
   - `gas_mass=0.1`: 0.00% un-ionized in every octant (the region saturates
     the half-box well before star death, so full ionization is achievable).
   - `gas_mass=0.01`: an isotropic rim residue is EXPECTED, not a bug — the
     region is still growing at star death and the check measures against
     the peak radius, so the never-reached outer shell stays neutral (finite
     photon budget). Pass = no octant/axis outlier (all octants within a few
     tenths of a percent of each other, mean-cos-to-diagonal ~ 0) AND overall
     residue at or below the 2.5% archived reference (run_A2_gm001,
     `dod_6c96aaeeb_phaseC`; re-measured 2026-08-12). A directional
     signature or a residue well above the reference indicates a real
     coverage/task-graph problem.
   (History: the criterion previously said "0.00% in every octant" at both
   resolutions — that number was only ever real at `gas_mass=0.1`; the
   phaseC archive's 0.01 "pass" cited a mislabeled 0.1-resolution log.
   Amended by Darwin's decision, 2026-08-12.)

Criteria 1-3 all run under `IONIZATION_FEEDBACK_DEBUG_NO_COOLING`, which
compiles out every path that depends on tags expiring. Work touching the
**photon budget, tag lifetimes, or region maintenance** is therefore not
covered by them and also needs:

4. **A cooling-enabled run actually exercises maintenance.** Build without
   `NO_COOLING`, with `CFLAGS+="-DSWIFT_DEBUG_CHECKS_VERBOSE"` *and*
   `--enable-debugging-checks` (the verbose macro has no configure flag and
   does not build alone — its blocks use `cell->cellID`). In the log, `held`
   and `recomb photons` must both go non-zero, and the star must not pin at
   `min_star_timestep`. `held` stuck at 0 all run = the mechanism silently
   never ran. If `claimed` collapses to 0 while the `has NOT exhausted` line
   still shows photons left, the gather is finding no eligible gas — check
   `feedback_part_can_be_ionized`'s density/temperature gates against the
   run's real gas state before blaming the budget.
5. **Cadence-independence**, for any photon-budget or time-integration change:
   `husmith2017_grackle_coupling` at `HII_rebuild_time_Myr` ∈ {0.5, 0.1, 0.03,
   0.01, 0.003, 0.001}, then
   `hu_smith_analytic_check.py`. The six must agree — cadence is a numerical
   parameter. Report the discarded-budget fraction too: at coarse cadence the
   `max_ngbs` ceiling binds and caps how far agreement can extend.
6. **Restart works**, for any change to the radiation sub-struct of
   `feedback_spart_data` or to the radiation tables. Both have caused
   restart-only segfaults.
7. **All four Grackle modes compile**, for any change under
   `src/cooling/grackle/`.
