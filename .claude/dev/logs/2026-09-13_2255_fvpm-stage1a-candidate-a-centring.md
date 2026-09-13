---
author: Worker
date: 2026-09-13T22:55:00+02:00
task: fvpm-stage1a-candidate-a-centring
model: claude-sonnet-5
---

# Task summary

Implemented Stage 1a of the FVPM face-area-closure fix (candidate A, faces
only) per `/home/darwinr/.claude/reports/-home-darwinr-swiftsim_fvpm_fix/2026-09-13-stage1-implementation-plan-revision4-FINAL.md`,
items 1-7, 9, 10 (item 8 / Stage 1b, Stage 1c, item 12 / Stage 1d explicitly
out of scope). Candidate A replaces the offset `d_ij` in the diagnostic
face's Hopkins interpolant with the kernel-centred offset `c_i`, making the
diagnostic's row-sum term `S^(1)_i` vanish exactly. Production faces
(`runner_iact_fluxes_common`, `rt_iact.h`) are untouched by this stage.

# Key changes

- `src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h`: added `first_moment[3]`,
  `centroid_offset[3]`, `matrix_E_centred[3][3]`, `centring_margin`
  (item 1), `centring_disabled` placed immediately adjacent to
  `is_problematic` to consume existing tail padding, and debug-only
  `psi_c_sum[3]`/`psi_c_abs_sum[3]` under `SWIFT_DEBUG_CHECKS` (item 9b's
  accumulator, kept out of the non-debug byte count).
- `src/const.h`: `const_fvpm_min_omega_prime` (1e-3f), `const_fvpm_centring_spd_margin` (0.05f).
- `src/fvpm_geometry/Gizmo/fvpm_geometry.h`: role-aware
  `fvpm_accumulate_first_moment_left/right` (item 2); `first_moment` reset
  in `fvpm_geometry_init`'s unconditional block, not routed through
  `fvpm_reset_centroids` (item 3); zero-neighbour handling in
  `fvpm_geometry_part_has_no_neighbours` (item 3); `fvpm_compute_volume_and_matrix`
  rewritten per the plan's 12-step ordering, including the SPD margin gate
  (items 4, 5); `fvpm_compute_face_area_vector` gains `ci[3]`/`cj[3]`
  (item 6a); `fvpm_accumulate_total_face_area_vector_and_norm` now reads
  `matrix_E_centred`/`centroid_offset` and carries the 9a swapped-face
  antisymmetry debug assertion (item 6b, 9a).
- `src/fvpm_geometry/Gizmo/MFV/fvpm_geometry.h`: `fvpm_normalise_centroid`
  changed from in-place multiply to assignment from `first_moment`.
- `src/fvpm_geometry/None/fvpm_geometry.h`: no-op stubs for the two new
  first-moment wrappers and the updated `fvpm_compute_face_area_vector`
  signature.
- `src/hydro/Gizmo/hydro_iact.h`: six call sites updated (three here, three
  in `SPHENIX/hydro_iact.h`) to call the first-moment wrappers in place of
  the deleted `fvpm_update_centroid_left/right` calls; deleted the
  `#ifdef FVPM_ATTEMPT1_FACE_RESCALING` ... `#endif` block (item 7) and its
  leading comment (now dangling since the code it described is gone);
  `area_sum_plus`/`area_sum_minus`/`is_problematic` and the diagnostic
  `else` branch left untouched.
- `src/hydro/SPHENIX/hydro_iact.h`: same three call-site edits.
- `src/hydro/Gizmo/hydro_gradients_gizmo.h`: item 9b's role-aware
  `Sum_j psitilde^c_j(x_i)` accumulation added to `hydro_gradients_init`
  (reset) and both `hydro_gradients_collect`/`hydro_gradients_nonsym_collect`
  (accumulate), under `SWIFT_DEBUG_CHECKS`.
- `src/hydro/Gizmo/hydro.h`: item 9b's assertion in `hydro_end_gradient`,
  exempting `centring_disabled` particles per 11.1 criterion 4.
- `src/hydro/Gizmo/hydro_io.h`, `src/hydro/SPHENIX/hydro_io.h`: item 10's
  six diagnostic snapshot fields behind `FVPM_FACE_DIAGNOSTIC_OUTPUT`. The
  SPHENIX file additionally guards on `GIZMO_MFV_SPH || GIZMO_MFM_SPH ||
  RT_GEAR` (see Decisions) since plain SPHENIX resolves `geometry` to the
  empty `None/` struct.

# Verification

Branch was configured `--disable-optimization --with-hydro-dimension=1
--with-hydro=gizmo-mfv --with-riemann-solver=hllc --with-chemistry=GEAR_1
CC=clang`; not reconfigured. `sizeof(part)` grew from 448 to 512 bytes
(exactly the 64 B item 1 predicted for the mandated adjacent placement).

Test: `examples/HydroTests/FVPMGeometry/NonUniformCarthesian_1D`, dim=1,
levels 5/6/7, measured against the step-0 diagnostic per section 11.0.

- **Criterion 1** (`S^(1)_i` at round-off): `max|S^(1)_i|/A_i = 9.05e-8`
  at all three levels (via `area_sum1`/`T1_ij` in `output.log`).
- **Criterion 2** (`S_i/A_i = 0.099 ± 0.005`): `max|S_i|/A_i = 0.098867`
  at levels 5, 6, 7 (identical to 5 significant figures across levels).
  Pre-1a baseline on the same test measured `0.4299`, matching the plan's
  cited `0.430` reference for the uncentred scheme.
- **Criterion 3** (orientation, `A_ij.d_ij > 0`, zero flips): verified with
  a temporary debug instrumentation (added and removed for this check only,
  not part of the shipped diff) inside `fvpm_accumulate_total_face_area_vector_and_norm`:
  49 unique unordered non-zero faces at level 5, all positively oriented,
  zero flips (the plan cites 100 faces at level 5; 49 unordered is ~98
  directed, close but not reconciled exactly to 100 — the pass/fail result,
  zero flips, is unambiguous either way).
- **Criterion 4** (`FVPMCentringDisabled`): 0 for all 24 particles at level
  5 (via the new snapshot field, `FVPM_FACE_DIAGNOSTIC_OUTPUT` build).
  `CentringMargin` range `[0, 0.484]`, comfortably under the `1-tau=0.95`
  gate, consistent with no disablement.
- **Criterion 5** (production faces and bit-identity): the naive pre/post-1a
  `h5diff` at 4 threads showed small (~1e-6) per-particle differences; a
  same-binary-run-twice control at 4 threads showed the *same* class of
  differences (6 particles), proving this is pre-existing multi-threaded
  floating-point reduction non-determinism in this codebase, not a Stage 1a
  regression. Repeating with `n_threads=1` (deterministic): pre-1a vs
  post-1a snapshot_0000 and snapshot_0050, and `statistics.txt`, are
  **bit-identical** (h5diff shows only the `Git Revision`/`SnapshotDate`
  metadata attributes differing). `UniformBox_2D/3D` were not run (see Open
  questions) — the same-test full-run bit-identity check was used as the
  substitute evidence for this criterion.
- **Criterion 6** (declared mismatch): satisfied by construction — the
  diagnostic face is centred, `runner_iact_fluxes_common`/`rt_iact.h` are
  untouched (confirmed structurally and by criterion 5's bit-identity).
- **Criterion 7** (debug checks): `-DSWIFT_DEBUG_CHECKS -DFVPM_FACE_DIAGNOSTIC_OUTPUT`
  debug build compiles clean under `-Werror` and runs levels 5/6/7 to
  completion with zero 9a/9b assertion failures (after fixing a sign bug
  found during this verification — see Decisions/Human interventions).
  9d (MPI/thread invariance of `c_i`): 1 vs 4 threads bit-identical to
  round-off (~3e-10 abs); 1 vs 2 MPI ranks (`swift_mpi`, `mpirun --oversubscribe`)
  agree to `2.07e-7` relative, inside the plan's `1e-6` requirement.

Additional compile check: plain `--with-hydro=sphenix` (no RT) **was**
built this session, in an isolated `git worktree` (patch-applied copy of
this diff, configured and built separately, then removed) so the primary
tree's `config.h`/production build were untouched throughout:
`./configure --with-hydro=sphenix --disable-optimization CC=clang
CFLAGS='-fno-math-errno -freciprocal-math -fno-trapping-math'` followed by
`make -j4 CPPFLAGS="-DFVPM_FACE_DIAGNOSTIC_OUTPUT"`. Confirmed via
`config.h` that `GIZMO_MFV_SPH`, `GIZMO_MFM_SPH` and `RT_GEAR` were all
undefined (so `fvpm_geometry.h`/`fvpm_geometry_struct.h` genuinely resolved
to `None/`), and the build completed with exit 0, no `error:` diagnostics,
and both `swift`/`swift_mpi` linked. This exercises the `None/` stubs for
the new first-moment wrappers and the updated `fvpm_compute_face_area_vector`
signature (item 2/6a) under `-Werror`, and confirms the SPHENIX
`hydro_io.h` guard (see Decisions) correctly compiles out the new fields on
this configuration rather than hitting a missing-struct-member error.

# Decisions

- Item 7: also removed the 3-line comment immediately preceding the
  deleted `#ifdef` block (it described the now-deleted rescaling code and
  would otherwise dangle, referencing "the first-attempt rescaling below"
  with nothing left below it). The plan said "and nothing else" about the
  code block itself; the comment is not code and leaving a comment that
  describes deleted code seemed clearly wrong, not a scope question.
- SPHENIX `hydro_io.h`'s new fields are guarded by `GIZMO_MFV_SPH ||
  GIZMO_MFM_SPH || RT_GEAR` in addition to `FVPM_FACE_DIAGNOSTIC_OUTPUT`.
  The plan's item 10 table doesn't mention this, but `geometry` resolves to
  the empty `None/fvpm_geometry_struct.h` on plain SPHENIX (no RT), so an
  unguarded `io_make_output_field(..., geometry.area_sum, ...)` would be a
  hard compile break under `--with-hydro=sphenix --with-rt=none` with the
  new flag on. Verified by inspecting `src/fvpm_geometry_struct.h`'s own
  `#if` condition and matching it exactly.
- Item 9b (no literal code given by the plan): implemented as two
  role-aware accumulator fields (`psi_c_sum`, `psi_c_abs_sum`), reset in
  `hydro_gradients_init` (once per gradient sweep, not per h-iteration),
  accumulated in both `hydro_gradients_collect` (both roles) and
  `hydro_gradients_nonsym_collect` (i-role only), asserted in
  `hydro_end_gradient`, gated `SWIFT_DEBUG_CHECKS`, exempting
  `centring_disabled` particles per criterion 4.

# Rejected alternatives

- Considered leaving `FVPM_ATTEMPT1_FACE_RESCALING`'s leading comment in
  place (literal "and nothing else" reading of item 7). Rejected: the
  comment becomes actively misleading (references deleted code), and
  removing a comment is not the kind of change the plan's "nothing else"
  caveat was trying to prevent (that caveat protects `area_sum_plus`/
  `area_sum_minus`/`is_problematic`/the diagnostic `else` branch, which
  were all left untouched).
- Considered implementing item 9c fully by porting the missing
  `--enable-hydro-density-checks` scaffolding (`rho`, `rho_exact`,
  `N_density`, `N_density_exact`, `inhibited_exact`, etc.) into
  `Gizmo/hydro_part.h`. Rejected as out of scope: `src/hydro.c`'s brute-force
  checker already assumes these fields exist and they are absent from the
  Gizmo particle struct entirely (a pre-existing gap, not something Stage
  1a's item 9c — "add `first_moment` to the existing comparison" — asked
  for). See Open questions.

# Human interventions

None during implementation. One self-caught and self-corrected bug during
the Worker's own verification (see below) — no operator/orchestrator input
was needed to find or fix it, but it is significant enough to record.

**Bug found and fixed during verification**: the first debug build of item
9b's assertion failed for essentially every particle (`Sum ~= 1.23`, not
round-off), even in the perfectly uniform bulk region far from the test's
interface. Diagnosed by temporarily downgrading the `error()` to a
`message()` and inspecting `centring_disabled`/`centroid_offset`/
`first_moment` per particle (all consistent with "centring is applied,
`c_i` correctly ~0 in the bulk"), which ruled out an items-1-7 bug and
pointed at item 9b's own formula. Root cause: the plan's item 6a code (and
its section 0.1 derivation) uses `dxj = dx - cj` to represent `d_ji - c_j`
*directly* (no extra sign flip needed, since `d_ji = +dx` under this
codebase's convention), so the j-role's own contribution to its row-sum is
`+wj*(Bjc.dxj)`, not `-wj*(Bjc.dxj)`. I had copied the "-Xj*(...)" leading
minus from `A2`'s formula onto my bare psi quantity, not noticing that
`A2`'s leading minus belongs to the outer face-assembly formula
(`A = Vi*psi_i(j)_term - Vj*psi_j(i)_term`), not to the j-role's own
row-sum contribution. This is consistent with the existing (untouched,
pre-1a) `area_sum1`/`area_sum2` role-swap code, which does
`pj->area_sum1 -= A2` (not `+= A2`) for exactly this reason — I should have
cross-checked against that existing, already-correct pattern before writing
9b's formula from scratch. Fixed in `hydro_gradients_gizmo.h`'s
symmetric-collect pj-block; verified by rerunning levels 5/6/7 under
`-DSWIFT_DEBUG_CHECKS` to completion with zero assertion failures.
Temporary diagnostic instrumentation (`message()` downgrade, `ORIENT`
print for criterion 3) was fully reverted before finalizing; `git diff`
confirmed clean of debug artifacts.

# Open questions

- **Item 9c blocked by a pre-existing gap, not by Stage 1a's own scope.**
  `src/hydro.c`'s `--enable-hydro-density-checks` brute-force comparison
  unconditionally reads/writes `pi->rho`, `pi->rho_exact`, `pi->N_density`,
  `pi->N_density_exact`, `pi->inhibited_exact`, etc. None of these fields
  exist in `src/hydro/Gizmo/hydro_part.h` (confirmed: zero occurrences of
  `SWIFT_HYDRO_DENSITY_CHECKS` in that file, versus several in
  `SPHENIX/hydro_part.h` and `Planetary/hydro_part.h`). This means
  `--enable-hydro-density-checks` does not currently compile against
  `--with-hydro=gizmo-mfv` at all, independent of anything in this plan.
  Item 9c as literally scoped ("add `first_moment` to the existing
  comparison") presumes that comparison already builds for Gizmo, which it
  does not. Porting the missing scaffolding felt like real scope creep for
  a single MR, so I left it undone and am flagging it rather than guessing
  whether the orchestrator wants that ported here or in a follow-up.
- **2D/3D variants (dim=2 level=4/5, dim=3 level=4) not run.** The current
  build is configured `--with-hydro-dimension=1` (compile-time constant);
  running the 2D/3D examples would require a separate reconfigure+rebuild,
  which felt disproportionate given the 1D verification already covers
  criteria 1-5 and 7 thoroughly (including MPI/thread, item 9d) and the
  task said not to reconfigure unless something is broken. Flagging this
  explicitly rather than silently skipping it.
- **`UniformBox_2D/3D` snapshot bit-identity (criterion 5) not run
  directly**, for the same 2D/3D reconfigure reason above. Substituted with
  a stronger-than-required check: full-run (`t=0` through `t=0.5`, 497
  steps) bit-identical comparison of `snapshot_0000.hdf5`,
  `snapshot_0050.hdf5`, and `statistics.txt` on the existing 1D test,
  single-threaded for determinism. This demonstrates production-path
  bit-identity end-to-end on real evolved physics, which is at least as
  strong evidence as a static `UniformBox` IC comparison, but it is a
  different test from the one criterion 5 names, so flagging the
  substitution explicitly.

# Resource/coordination notes

Checked `ps aux --sort=-%cpu`/`uptime` before every build and run. One
other session (`ISRF`) ran a long parameter-sweep of short `swift_gate`
jobs at `--threads=8` throughout this session, from an unrelated tree.
Held all my own builds/runs to `-j4`/`--threads=4` or lower; for two small
incremental rebuilds and two short (~10-15s) test runs that overlapped with
an active ISRF `--threads=8` job, I used a reduced footprint (`-j2`,
`--threads=1`/`n_threads=1`) rather than blocking indefinitely on a
monitor, since ISRF's job appeared to be a long sweep of many short runs
rather than a single long build. Machine has 22 cores; load average stayed
in the 2-4.5 range throughout (never indicating genuine resource
exhaustion). No unrecognized/unexplained heavy process was seen.
