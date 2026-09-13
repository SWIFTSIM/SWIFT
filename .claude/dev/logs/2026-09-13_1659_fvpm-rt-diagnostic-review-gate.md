---
author: Orchestrator
date: 2026-09-13T16:59:00+02:00
task: fvpm-force-loop-accumulation-rt-diffusion
model: claude-sonnet-5 (red-team-reviewer and neutral-arbiter subagents for the gate)
---

# Task summary

Records the review gate (red-team + neutral-arbiter) run on the Worker
deliverable logged in
`2026-09-13_1625_fvpm-force-loop-accumulation-rt-diffusion.md`, before
committing it. Covers the RT-side face-closure diagnostic (accumulated in
RT's flux/transport loop, not its gradient loop, to avoid needing a new
type-2 RT-gradient dispatcher or touching SPHENIX's shared gradient loop)
and the design-doc write-up of the equivalent (unimplemented) plan for the
diffusion solver in `swiftsim_0`.

# Key changes

None beyond what the Worker already implemented (see the referenced log).
This entry documents the gate outcome that authorizes committing that work.

# Verification

- Orchestrator independently confirmed before dispatching the gate: exact
  3-file diff scope (`fvpm_geometry_struct.h`, `rt.h`, `rt_iact.h`), no
  leftover git worktrees or stray branches from the Worker's disposable
  build-verification worktrees (`git worktree list`/`git branch` matched
  the pre-existing 29-branch baseline exactly), `swiftsim_0` untouched
  (status hash unchanged), `FVPM_RT_FACE_CLOSURE_DIAGNOSTIC` confirmed
  undefined anywhere in the build (off by default), companion LaTeX
  (`design_plan.tex`) recompiles cleanly (33 pages).
- Red-team review: no blocking defects. Independently re-derived and
  confirmed the antisymmetric role-swap for the RT split under RT's own
  gating condition (`mode==1 || pj->rt_data.flux_dt<0`, not a copy-paste
  of Gizmo's different `interaction_mode` convention); confirmed the
  `RT_GEAR`-on-`gizmo-mfv` struct-selector guard in `fvpm_geometry.h` is
  consistent (no field-access mismatch); re-verified in `cell_unskip.c`
  that `task_subtype_rt_transport` already gets correct MPI cross-rank
  activation (`MPI_SYMMETRIC_FORCE_INTERACTION_RT`) — the exact treatment
  found *missing* for Gizmo's gradient loop in the prior gate, present
  here; confirmed the reset/check lifecycle across RT sub-cycles has no
  early-return gap; confirmed the new `rt.h` -> `rt_iact.h` include is
  load-bearing and non-circular. One cosmetic nit (the two new lifecycle
  helpers live in `rt_iact.h`, a leaf interaction header, but are only
  called from `rt.h`, inverting the codebase's normal layering) — not
  blocking. One open item: never compiled with the diagnostic flag on in
  *this* sandbox (`config.h` here is `RT_NONE`+`GIZMO_MFV_SPH`, no
  grackle available).
- Neutral-arbiter arbitration: **APPROVED-WITH-CAVEATS**. Reconciled the
  Worker's build claim ("four full builds in disposable worktrees") against
  red-team's "untested by any build in this repo": not a contradiction —
  the Worker's builds (two non-MPI `RT_GEAR`+`sphenix`+grackle configs, one
  MPI build, one targeted recompile of `runner_doiact_hydro.c` with the
  diagnostic flag, independently confirmed by the arbiter via `grep` to be
  the exact translation unit instantiating `TASK_LOOP_RT_TRANSPORT`) did
  compile-verify the new code; red-team's statement is only about this
  specific sandbox's currently-active config, not a rebuttal of the
  Worker's claim. Both agree: no *runtime* smoke test occurred; that gap
  is real, undisputed, and must be recorded rather than silently dropped.
  Also corrected one overstatement in the red-team report (claimed the
  ill-conditioned-branch `A1+A2` was "exact, not approximate"; the diff's
  own comment correctly says otherwise — a report inaccuracy, not a diff
  defect, and far below the diagnostic's 1e-2 threshold either way).
  Independently found and confirmed a point neither review had stated
  explicitly: the new `rt_area*` fields are write-only into a `warning()`
  call and read nowhere else in the codebase, so even a latent bug in this
  diagnostic cannot corrupt simulation state, only misleading log output,
  and only when manually enabled.

# Decisions

- Commit only the three reviewed files by explicit path, plus this and
  the Worker's log entry — never `git add -A`/`git add .` — per the
  arbiter's explicit instruction, given the same class of untracked HDF5
  binaries sitting in the working tree as the previous gate.
- Record the diagnostic flag's build-verification status (four
  configurations) and the absence of a runtime smoke test in the commit
  message itself, per the arbiter's required caveat, mirroring how the
  prior commit recorded its MPI limitation in the message body rather
  than only in this log.
- Left the include-layering cosmetic nit (moving the two helpers from
  `rt_iact.h` into `rt.h`) unaddressed in this commit; non-blocking per
  both reviews, noted as a fast-follow.

# Rejected alternatives

- None beyond what the Worker's own log already records (switching RT's
  gradient task to type-2 from scratch, rejected in favour of
  force-loop accumulation; conflating RT's diagnostic fields with
  Gizmo's, rejected because `RT_GEAR`+`gizmo-mfv` is a real, buildable
  combination that would otherwise corrupt one diagnostic or the other).

# Human interventions

- Operator directed this whole line of work ("you mentioned the
  possibility to accumulate in the force loop to not alter gradients of
  other schemes. So, dispatch an agent to work it out"), redirecting from
  the gradient-loop-type approach used for Gizmo to a force-loop
  accumulation approach for RT/diffusion specifically, after the
  Orchestrator explained the SPHENIX/RT/diffusion generalization problem
  in conversation.

# Open questions

- No runtime smoke test of `FVPM_RT_FACE_CLOSURE_DIAGNOSTIC` has been
  performed; running `examples/RadiativeTransferTests/Advection_1D` with
  the flag enabled is a named follow-up (already written up in
  `design_plan.tex` and the Worker's dev log).
- The diffusion-solver equivalent (in `swiftsim_0`) remains a design
  only, not implemented; `swiftsim_0` was correctly left untouched
  throughout (read-only), per its own pre-existing uncommitted WIP that
  has not been cleared for editing.
- The include-layering cosmetic nit (`rt_iact.h` helper placement) is
  unaddressed; low priority.
