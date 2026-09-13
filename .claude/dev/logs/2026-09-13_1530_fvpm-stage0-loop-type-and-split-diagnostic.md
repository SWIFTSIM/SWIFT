---
author: Worker
date: 2026-09-13T15:31:02+02:00
task: fvpm-stage0-loop-type-and-split-diagnostic
model: claude-sonnet-5
---

# Task summary

Implemented Stage 0 items 3 and 4 of the FVPM face-closure design plan
(`design_plan.tex`, `\subsection{Staged implementation plan}`):

- Item 3: switch Gizmo's gradient loop from type-1 (gather) to type-2
  (union) via `EXTRA_HYDRO_LOOP_TYPE2` in `src/part.h`.
- Item 4: split the per-particle face-area closure defect `Sum_j A_ij`
  into its two additive terms `T1` (row sum, particle i's own data) and
  `T2` (column sum, neighbour j's data), wired into real code (not the
  scratch patch).

# Key changes

- `src/part.h`: added `#define EXTRA_HYDRO_LOOP_TYPE2` to the
  `GIZMO_MFV_SPH || GIZMO_MFM_SPH` block, same line the previous session
  reported a permission refusal on. The edit went through cleanly this
  session with no refusal.
- `src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h`: added
  `float area_sum1[3]` and `float area_sum2[3]` with one-line doc
  comments matching the style of `area`/`area_sum`.
- `src/fvpm_geometry/Gizmo/fvpm_geometry.h`:
  - Reset `area_sum1`/`area_sum2` in `fvpm_geometry_init` (inside the
    existing `time_bin == 0` block, alongside `area_sum`) and in
    `fvpm_geometry_part_has_no_neighbours` (unconditional, matching
    `area_sum` there).
  - `fvpm_compute_face_area_vector`: added `A1[3]`, `A2[3]` output
    parameters. Both branches (well-behaved Hopkins face and the
    ill-conditioned SPH-style fallback) now compute the two additive
    terms directly instead of only their sum, so no kernel evaluation is
    duplicated -- `wi`, `wj`, `hi_inv_dim`, `hj_inv_dim`, `Bi`, `Bj` are
    reused as already computed.
  - `fvpm_accumulate_total_face_area_vector_and_norm`: accumulates
    `area_sum1 += A1`, `area_sum2 += A2` for `pi`; for `pj` in symmetric
    mode, `area_sum1 -= A2`, `area_sum2 -= A1` (terms swap roles and flip
    sign for `A_ji`; independently re-derived from `A_ji = -A_ij` term
    by term, matches the scratch patch's logic).
  - `fvpm_check_total_face_area_vector_sum`: extended the existing
    `warning(...)` call (only fires when a particle is already flagged)
    to also print `Sum_j T1_ij` and `Sum_j T2_ij`. No second unconditional
    `message()` dump was added.
- `src/fvpm_geometry/None/fvpm_geometry.h`: updated the
  `fvpm_compute_face_area_vector` stub signature to match (empty body,
  no struct fields touched -- `area`/`area_sum`/`area_sum1`/`area_sum2`
  are Gizmo-only, confirmed by grep before editing).

`git diff --stat`:
```
 src/fvpm_geometry/Gizmo/fvpm_geometry.h        | 70 +++++++++++++++++++++-----
 src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h |  6 +++
 src/fvpm_geometry/None/fvpm_geometry.h         |  4 +-
 src/part.h                                     |  1 +
 4 files changed, 67 insertions(+), 14 deletions(-)
```

# Verification

Build: `make -j4` from the existing configure (`--disable-optimization
--with-hydro-dimension=1 --with-hydro=gizmo-mfv --with-riemann-solver=hllc
--with-chemistry=GEAR_1`, per `config.log`; not reconfigured). Clean
build, no warnings/errors from the touched files; second `make -j4` run
reported nothing to rebuild.

Acceptance test: `examples/HydroTests/FVPMGeometry/NonUniformCarthesian_1D`,
`n_threads=2 uniform=0 dim=1 level=5 ./run.sh`.

- Particle id 16 (x = 0.46875, the "x = 0.469" particle in the design
  doc): `Sum_j A_ij` went from `+2.523341e-01` (baseline, pre-existing
  `1D_nonuniform_eulerian_particles_l5/output.log` from before this
  session's changes, gather-set loop) to `+3.485932e-01` in the new run
  (union-set loop). Design doc acceptance: `+0.252` -> `+0.349`. Matches
  to 3 significant figures.
- Particle id 17 (x = 0.5, the interface particle in `fourth_attempt.tex`
  Table `tab:split`): new run gives `A_tot = 1.947409e+00`,
  `Sum_j A_ij = -8.370214e-01`, `Sum_j T1_ij = -7.522899e-01`,
  `Sum_j T2_ij = -8.473156e-02`. Normalised: `S/A = -0.430`,
  `S1/A = -0.386`, `S2/A = -0.044`. Matches Table `tab:split`
  (`-0.430`, `-0.386`, `-0.044`) essentially exactly, and `|S1| >> |S2|`
  as expected at the interface.

Regression sanity check: same test, `uniform=1` (Cartesian, level 5).
Ran to completion (`main: done. Bye.`, exit 0), zero `WARNING` lines in
`output.log` -- a uniform grid still closes near round-off, as expected.

# Decisions

- Modified `fvpm_compute_face_area_vector`'s signature (added `A1`, `A2`
  output params) rather than recomputing kernels in the caller, since
  grep confirmed it has exactly one call site
  (`fvpm_accumulate_total_face_area_vector_and_norm`); the design_plan.tex
  table itself notes "None of the solvers calls
  `fvpm_compute_face_area_vector`" -- it is diagnostic-only, so widening
  its signature carries no risk to the real hydro force/gradient loops
  (those are separate copies in `hydro_iact.h`/`hydro_gradients_gizmo.h`
  per `sec:copies`).
- Extended the ill-conditioned SPH-style fallback branch of
  `fvpm_compute_face_area_vector` to also emit an analogous `A1`/`A2`
  split (i-only term vs j-only term of the SPH-style face), even though
  the design doc's `T1`/`T2` definitions are stated for the well-behaved
  Hopkins-face branch. This keeps `area_sum1 + area_sum2 == area_sum`
  exactly in all cases, including if a particle's geometry matrix is ever
  flagged ill-conditioned during the test.
- Verified the `pj` sign/role-swap convention independently against
  `A_ji = -A_ij` term by term (see reasoning in this session) before
  trusting the scratch patch's logic; it checked out and matches the
  patch's `pj->area_sum1 -= T2; pj->area_sum2 -= T1;` construction
  exactly.
- Left the `stage0_verify_nonuniform`/`stage0_verify_uniform` run
  directories in place under `examples/.../NonUniformCarthesian_1D/`
  (untracked, same convention as the pre-existing `1D_nonuniform_*`
  directories already in the working tree) rather than deleting them, so
  the reviewer can inspect the raw logs if needed.

# Human interventions

None. The item-3 edit to `src/part.h` that a previous session reported
as refused by the permission system went through without any refusal or
prompt in this session.

# Open questions

None outstanding for Stage 0 items 3/4. Stage 1 (candidate A, faces only,
all solvers) is out of scope for this task per the assignment.
