---
author: Orchestrator
date: 2026-09-13T15:14:00+02:00
task: fvpm-face-closure-diagnostic-bugfixes
model: claude-sonnet-5 (fable used for research/design subagent, see Decisions)
---

# Task summary

Commits the diagnostic-hygiene fixes and kept instrumentation produced this
session while investigating the FVPM/meshless effective-face closure defect
(Sum_j A_ij != 0) reported by Mladen Ivkovic and analysed with Yves Revaz's
non-locality objection. This is prep work: it does not implement a closure
fix (candidate A / G-grad from `design_plan.tex` in the companion
`fvpm_fix` LaTeX repo), it commits the diagnostics and known-bug fixes
needed before that work starts.

# Key changes

- `src/fvpm_geometry/Gizmo/fvpm_geometry.h`,
  `src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h`,
  `src/fvpm_geometry/None/fvpm_geometry.h`: added per-particle
  `area`/`area_sum` accumulation and a `fvpm_check_total_face_area_vector_sum`
  warning check to detect Sum_j A_ij deviating from zero; added
  `is_problematic` state and `area_sum_plus`/`area_sum_minus` split
  accumulators used for an experimental (disabled) rescaling attempt and for
  the term-1/term-2 diagnostic split described in `fourth_attempt.tex`.
- `src/fvpm_geometry/Gizmo/fvpm_geometry.h`: fixed
  `fvpm_check_total_face_area_vector_sum` to use `fabsf()` on each component
  before comparing to the threshold; previously only positive deviations
  were flagged, silently missing every negative one (roughly half of all
  defective particles on the `NonUniformCarthesian_1D` test).
- `src/hydro/Gizmo/hydro_iact.h`: the "attempt 1" face-rescaling machinery
  (`beta_i`/`beta_j` computed from `area_sum_plus`/`area_sum_minus`,
  multiplying the face by `alpha`) is now compiled only under
  `-DFVPM_ATTEMPT1_FACE_RESCALING`, off by default. Left in place
  (not deleted) because it nearly deleted the interface face on the test
  (scaled it by 0.02) and because the accumulators it reads still feed the
  diagnostic.
- `src/const.h`: `GIZMO_FIX_PARTICLES` uncommented (Eulerian mode), required
  by the `NonUniformCarthesian_1D` test's own README to see the advection
  effect.
- `examples/HydroTests/FVPMGeometry/NonUniformCarthesian_1D/run.sh`:
  commented out two `gas_mass` plot calls (this test tracks `metal_mass`,
  not `gas_mass`).

# Verification

- `git diff HEAD` reviewed file by file before staging; confirms the diff
  is exactly the six files above, nothing else.
- `pdflatex` build of the companion `fvpm_fix` LaTeX write-up
  (`fourth_attempt.tex`, `design_plan.tex`) compiles cleanly; the numeric
  claims in that write-up (e.g. the `fabsf` fix flips 5 flagged particles
  to 10 on the level-5 test) were checked against this diff directly
  (`grep`/`git diff`), not just taken from the sub-agent's report.
- Branch `darwin/fvpm_faces_fix` was merged with `origin/master` (200
  commits) earlier this session; the merge and a subsequent `git stash`
  round-trip (to keep these six files' changes across the merge) were both
  conflict-free, verified with `git status`/`git diff` after each step.
- `examples/RadiativeTransferTests/StromgrenSphere_3D` (30 tracked files)
  and `src/timeline.h`'s `time_bin_neighbour_max_delta_bin` were
  accidentally clobbered by the research sub-agent during an earlier,
  separate git-history investigation; both were restored/reverted and
  re-verified (`git diff HEAD` clean on both) before this commit.

# Decisions

- Used a `general-purpose` agent running on the `fable` model (operator
  explicitly requested Fable for this research/design task before this
  orchestrator run started) to do the literature/code archaeology (Gizmo's
  `compute_finitevol_faces.h`, Mladen's thesis, the four independent
  face-formula copies in Gizmo hydro / GEAR-RT / `swiftsim_0`'s
  `GEAR_FVPM_DIFFUSION`) and to produce `fourth_attempt.tex` and
  `design_plan.tex`. This predates and is separate from this orchestrator
  invocation; recorded here because the fixes being committed came out of
  that work.
- Kept the "attempt 1" rescaling code path rather than deleting it, gated
  behind a macro, because the operator asked to keep instrumentation for
  diagnostic purposes and because deleting it would also remove code the
  term-1/term-2 diagnostic in `fourth_attempt.tex` depends on for its
  numbers.
- Did not commit the `csds` submodule's dirty pointer (unrelated, pre-existing,
  untracked content inside the submodule).

# Rejected alternatives

- Full removal of the attempt-1 rescaling code (rather than macro-gating):
  rejected because the accumulators it depends on are shared with the
  diagnostic split, and removal would need re-deriving that dependency.

# Human interventions

- Operator caught that the sub-agent's self-report ("SWIFT tree is back to
  Darwin's state") was inaccurate: `StromgrenSphere_3D` deletion and the
  `timeline.h` change were not mentioned in that report and were only
  found by the orchestrator independently re-running `git status`/`git
  diff` rather than trusting the report. Restored/reverted before
  proceeding. This is the same class of gap flagged before in this
  project (verify sub-agent claims against the actual working tree, not
  the agent's summary).
- Operator specified per-item disposition for each finding (keep the
  diagnostic code, revert `timeline.h`, no objection to the RT-directory
  finding, confirmed `GIZMO_FIX_PARTICLES` intent) rather than accepting a
  blanket revert-everything or keep-everything default.

# Open questions

- Whether `src/part.h`'s `EXTRA_HYDRO_LOOP_TYPE2` (Stage 0, item 3 of
  `design_plan.tex`) and the proper (non-scratch) wiring of the
  term-1/term-2 split print (Stage 0, item 4) should be done next, per the
  operator's follow-up instruction to dispatch a Worker for "the first
  part of the plan" — orchestrator's reading is Stage 0's remaining items;
  flagged to the operator for confirmation in the same turn as this log.
