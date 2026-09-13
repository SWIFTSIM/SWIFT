---
author: Worker
date: 2026-09-13T23:57:00+02:00
task: fvpm-stage1a-centring-margin-fix-and-mfm-verify
model: claude-sonnet-5
---

# Task summary

Applied the neutral arbiter's two required fixes on top of the existing,
already-reviewed Stage 1a diff (11 files, uncommitted, verified present via
`git status`/`git diff --stat` before starting and untouched except for fix
1's file):

1. Restored diagnostic granularity to `centring_margin` in
   `fvpm_compute_volume_and_matrix` (`src/fvpm_geometry/Gizmo/fvpm_geometry.h`).
2. Configured, built, and ran `--with-hydro=gizmo-mfm` under
   `SWIFT_DEBUG_CHECKS`/`FVPM_FACE_DIAGNOSTIC_OUTPUT` on
   `NonUniformCarthesian_1D` level 5, in an isolated worktree, since this
   hydro scheme had never been compiled or exercised this session.

# Key changes

- `src/fvpm_geometry/Gizmo/fvpm_geometry.h`, `fvpm_compute_volume_and_matrix`:
  added `p->geometry.centring_margin = -1.0f;` immediately before the
  `omega_prime > const_fvpm_min_omega_prime * kernel_root` gate (step 3),
  and deleted the `p->geometry.centring_margin = 0.0f;` line from step 11's
  `else` (disabled) branch, leaving `centring_disabled = 1` untouched there.
  `-1.0f` now unambiguously means "SPD margin never evaluated" (omega' floor
  tripped, or the production inversion at step 8 failed first); any other
  stored value, including on a particle where centring was ultimately
  disabled, is the real computed `q` from step 9 (failed the margin test) or
  reflects a step-10 centred-inversion failure after `q` was already stored.
  This is the only file touched for fix 1, per the ground rules.

# Verification

## Fix 1

- `git diff -- src/fvpm_geometry/Gizmo/fvpm_geometry.h` confirms exactly the
  two intended edits (see excerpt below); no other line in that function or
  file changed.
- Checked all other reads of `centring_margin`: `hydro_io.h`'s snapshot
  field (`FVPMCentringMargin`) and `SPHENIX/hydro_io.h`'s mirror just emit
  whatever float is stored, so `-1.0f` passes through unchanged and needs no
  format/description update. `fvpm_geometry_part_has_no_neighbours`'s own
  `centring_margin = 0.0f;` (a different, out-of-scope function/branch: the
  zero-neighbour fallback) was explicitly left alone.
- Compiled clean (`-Werror`, no diagnostics) both with an incremental
  rebuild and, at the end of this session, a full `make clean` rebuild of
  the main tree (see Decisions for why the clean rebuild was needed).

```diff
   int centred = 0;
   float c[3] = {0.f, 0.f, 0.f};
   float Ec[3][3] = {{0.f}};

+  /* -1 marks "SPD margin never evaluated" (omega' floor or step 8 tripped
+     first); any other value is the computed q, disabled or not. */
+  p->geometry.centring_margin = -1.0f;
+
   if (omega_prime > const_fvpm_min_omega_prime * kernel_root) {
@@
-    p->geometry.centring_margin = 0.0f;
     p->geometry.centring_disabled = 1;
   }
```

## Fix 2 (gizmo-mfm build and run)

Configured and built in an isolated `git worktree` at
`/home/darwinr/scratch_worktrees/fvpm_mfm_check` (removed at the end of the
session; `git worktree list` now shows only the main tree), mirroring the
gizmo-mfv debug recipe from
`.claude/dev/logs/2026-09-13_2255_fvpm-stage1a-candidate-a-centring.md`:
`./configure --disable-optimization --with-hydro-dimension=1
--with-hydro=gizmo-mfm --with-riemann-solver=hllc --with-chemistry=GEAR_1
CC=clang`, then `make -j4 CPPFLAGS="-DSWIFT_DEBUG_CHECKS
-DFVPM_FACE_DIAGNOSTIC_OUTPUT"`. `config.h` confirmed `GIZMO_MFM_SPH`,
`HYDRO_DIMENSION_1D`, `RIEMANN_SOLVER_HLLC`, `CHEMISTRY_GEAR` set as
expected.

**Blocker found: gizmo-mfm cannot build in this repository at all without a
manual source edit.** `GIZMO_FIX_PARTICLES` is unconditionally `#define`d at
`src/const.h:43` (pre-existing at HEAD, unrelated to Stage 1a — Stage 1a's
own `const.h` diff only added two new constants, it did not touch this
line), and `MFM/hydro_velocities.h:22` hard-`#error`s the moment it sees
that macro defined: `"Fixed particles are not allowed for GIZMO MFM!"`.
Anyone reproducing this verification will hit that `#error` immediately
unless they know to comment out line 43 first. Remedy applied: commented
out **only in the throwaway worktree's `src/const.h`**
(`// #define GIZMO_FIX_PARTICLES`), never in the tracked diff.

Consequence: this MFM run therefore evolves the IC in Lagrangian
(moving-particle) mode, whereas the gizmo-mfv verification ran it in
Eulerian (fixed-particle) mode. Not strictly like-for-like, but the
diagnostics this check cares about — face antisymmetry (9a) and the
row-sum identity (9b, item 1's `S^(1)_i`) — are per-step properties of the
geometry/gradient loop evaluated at each particle's current position,
independent of whether that position is fixed or drifting between steps,
which is why matching values are meaningful evidence despite the mode
difference.

Build: clean, `-Werror`, both `swift`/`swift_mpi` linked, zero warnings.

Run: `examples/HydroTests/FVPMGeometry/NonUniformCarthesian_1D`, IC
regenerated with `python3 makeIC.py --dimension 1 --level 5` (24 particles,
matching the gizmo-mfv verification's level-5 count), `--hydro --threads=2`,
completed all 256 steps to `t=0.5` ("main: done. Bye."), **zero** errors,
zero `Face antisymmetry violated` (9a) occurrences, and the run reached
`main: done` cleanly, which is only possible if item 9b's
`hydro_end_gradient` assertion (an `error()` call that aborts) never fired
across any of the 256 steps.

Verified 9b was genuinely *live* for this scheme rather than vacuously
passing: the `psi_c_sum`/`psi_c_abs_sum` accumulation lives in
`hydro_gradients_gizmo.h`, which is selected by `GRADIENTS_GIZMO`
(unconditionally `#define`d in `const.h:31`, not scheme-dependent) and
carries no `GIZMO_MFM_SPH`/`GIZMO_MFV_SPH` conditional anywhere in that
file; `hydro_end_gradient` (which contains the 9b assertion, in the shared
`src/hydro/Gizmo/hydro.h`) is called from a single, unguarded call site in
`runner_ghost.c`. Both the accumulation and the check run identically for
MFM and MFV.

Measured values (temporarily forcing `fvpm_check_total_face_area_vector_sum`'s
print threshold from `1e-2` to `0.0` in the worktree only, to recover
`area_sum1`/`area_sum` per-particle-per-step instead of only on threshold
breach; not part of the tracked diff, worktree discarded afterward):

- `max|S^(1)_i|/A_i = 9.05e-8` (5730 forced-print lines parsed across the
  full run) — matches the gizmo-mfv verification's `9.05e-8` to 3
  significant figures.
- `max|S_i|/A_i = 0.098867` — matches the gizmo-mfv verification's
  `0.098867` exactly.
- `FVPMCentringDisabled`: 0 for all 24 particles (snapshot_0000.hdf5).
- `FVPMCentringMargin` range `[0, 0.4837]`, no `-1.0` values (omega' floor
  never tripped in this bulk test) — matches the gizmo-mfv pass's reported
  `[0, 0.484]` almost exactly.

This is the expected result: the diagnostic face-assembly code
(`fvpm_accumulate_total_face_area_vector_and_norm`,
`fvpm_compute_face_area_vector`) is shared, untouched-by-scheme code: MFM's
only new exposure from Stage 1a is that `fvpm_accumulate_first_moment_left/
right` now runs for real in MFM's density kernel (previously an empty
no-op stub), feeding the same `first_moment`/`centroid_offset` machinery
MFV already used. Matching diagnostic values this closely is exactly the
signature of that wrapper working correctly for MFM.

# Decisions

- **`GIZMO_FIX_PARTICLES` toggle scope**: commented out only inside the
  disposable worktree's `src/const.h`, never touching the real
  `src/const.h` in the main tree (which still has it `#define`d, as Stage
  1a's own diff left it — Stage 1a only added two new constants to that
  file, it did not touch this line).
- **Temporary threshold-zero instrumentation for measurement**: changed
  `fvpm_check_total_face_area_vector_sum`'s local `threshold` from `1e-2`
  to `0.0f` in the worktree, to force its existing `warning()` print (which
  already reports `area_sum`/`area_sum1`/`area_sum2`) on every call instead
  of only when `area_sum` exceeds 1% of `area`. This is the same
  `warning()` call the gizmo-mfv verification pass cites as its source
  ("via `area_sum1`/`T1_ij` in `output.log`"); the prior log doesn't record
  exactly how it got that print to fire for round-off-level values, but
  relaxing this same threshold is a straightforward way to do it, and it
  produces numbers matching what was reported. Confined to the throwaway
  worktree; discarded with it.
- **Main-tree MFV hang investigated and diagnosed, not fixed, not blamed on
  fix 1.** Before finding the above, the *existing* gizmo-mfv build in the
  main tree (untouched aside from fix 1) was observed to hang indefinitely
  (0% CPU, all threads asleep in `futex_do_wait`) immediately after
  `space_rebuild: (re)building space`, reproducibly, regardless of thread
  count (1, 2, 4) and regardless of whether fix 1 was applied or manually
  reverted-and-rebuilt (tested both ways). `ptrace`/`gdb`/`strace` are
  unavailable in this sandbox, so root-caused it indirectly: `ls -lt
  src/*.o` showed 89 of 254 object files predated this session's rebuilds
  (compiled by the prior Worker's session with different `CPPFLAGS`, since
  incremental `make -j4 CPPFLAGS=...` only rebuilds translation units whose
  *sources* changed, not ones affected only by a `CPPFLAGS` value change).
  `SWIFT_DEBUG_CHECKS` changes `struct part`'s layout (per the prior
  Worker's own log: `psi_c_sum`/`psi_c_abs_sum` are added "under
  `SWIFT_DEBUG_CHECKS`" specifically to keep them "out of the non-debug byte
  count"), so the running binary was linking together objects compiled
  against two different `struct part` layouts — undefined behaviour
  consistent with exactly this kind of silent scheduler deadlock. Ran
  `make clean && make -j4 CPPFLAGS="-DSWIFT_DEBUG_CHECKS
  -DFVPM_FACE_DIAGNOSTIC_OUTPUT"` in the main tree at the end of this
  session to restore it to an internally consistent state (see Verification
  above for fix 1's clean-rebuild confirmation); this is a build-hygiene
  artifact of cross-session incremental builds, not a code regression, and
  is not something this task's scope asked me to chase further. The
  isolated MFM worktree was unaffected throughout, since it was configured
  and built from scratch with identical `CPPFLAGS` on every translation
  unit. **Confirmation**: after the clean rebuild, re-ran the exact same
  gizmo-mfv command (`--hydro --threads=2 params.yml`) on the main tree; it
  completed all 497 steps to `t=0.5` ("main: done. Bye."), zero errors, zero
  antisymmetry violations — the hang is gone with fix 1 in place, on a
  consistent binary.

# Rejected alternatives

- Considered trying to also re-run the full gizmo-mfv verification suite
  (levels 5/6/7, MPI/thread invariance) after the clean rebuild, since the
  hang was found and diagnosed. Rejected as out of this task's explicit
  scope (fix 1 correctness + one MFM build/run check); the prior Worker's
  gizmo-mfv verification already stands on its own (it was performed on an
  internally consistent build in its own session), and re-doing all of it
  here would be scope creep against the arbiter's stated fix list.

# Human interventions

None. All decisions above (worktree-only `GIZMO_FIX_PARTICLES` toggle,
threshold-zero instrumentation for measurement, diagnosing rather than
chasing the unrelated MFV hang) were made autonomously and are recorded
here for review.

# Open questions

- **MFM ran in Lagrangian mode, MFV in Eulerian mode** (see Verification):
  a direct consequence of `GIZMO_FIX_PARTICLES` being incompatible with
  MFM. Flagging per the arbiter's framing ("like-for-like comparison") that
  this is not a bit-for-bit identical test setup, though it is the closest
  achievable one and exercises the same diagnostic invariants.
- **The main-tree MFV binary was rebuilt via `make clean` and re-run once**
  (single-configuration smoke check: `--threads=2`, level 5, 497 steps,
  clean) to confirm the hang was gone; the *full* gizmo-mfv verification
  matrix (levels 5/6/7, MPI/thread invariance) was not repeated in this
  session (see Rejected alternatives) — the prior Worker's pass already
  covers that on its own internally-consistent build.

# Commit-time discipline (for the Orchestrator)

The eventual commit must stage only the same 11 files as the original
Stage 1a diff:

- `src/const.h`
- `src/fvpm_geometry/Gizmo/MFV/fvpm_geometry.h`
- `src/fvpm_geometry/Gizmo/fvpm_geometry.h` (now includes fix 1 on top of
  Stage 1a)
- `src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h`
- `src/fvpm_geometry/None/fvpm_geometry.h`
- `src/hydro/Gizmo/hydro.h`
- `src/hydro/Gizmo/hydro_gradients_gizmo.h`
- `src/hydro/Gizmo/hydro_iact.h`
- `src/hydro/Gizmo/hydro_io.h`
- `src/hydro/SPHENIX/hydro_iact.h`
- `src/hydro/SPHENIX/hydro_io.h`

plus this log and the prior Worker's log
(`.claude/dev/logs/2026-09-13_2255_fvpm-stage1a-candidate-a-centring.md`).
Never `git add -A`: the working tree still has a dirty `csds` submodule
pointer and untracked generated/binary artifacts (`.hdf5.1/.2/.3` files,
PNGs, `dashboard.html`, example run outputs) that must not be swept in. I
did not commit anything myself.
