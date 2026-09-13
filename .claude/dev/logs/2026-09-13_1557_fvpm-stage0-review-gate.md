---
author: Orchestrator
date: 2026-09-13T15:57:00+02:00
task: fvpm-stage0-loop-type-and-split-diagnostic
model: claude-sonnet-5 (red-team-reviewer and neutral-arbiter subagents for the gate)
---

# Task summary

Records the review gate (red-team + neutral-arbiter) run on the Worker
implementation logged in
`2026-09-13_1530_fvpm-stage0-loop-type-and-split-diagnostic.md`, before
committing it. Covers Stage 0 items 3 and 4 of `design_plan.tex` (in the
companion `fvpm_fix` LaTeX repo): `EXTRA_HYDRO_LOOP_TYPE2` for Gizmo, and
wiring the term-1/term-2 closure-defect split into real (non-scratch) code.

# Key changes

None beyond what the Worker already implemented (see the referenced log).
This entry documents the gate outcome that authorizes committing that work.

# Verification

- Independent orchestrator re-derivation of the antisymmetric role-swap
  (`A1_ji = -A2_ij`, `A2_ji = -A1_ij`) confirmed correct before the gate
  was even run.
- Red-team review (`red-team-reviewer` subagent) of the diff: found the
  `A = A1 + A2` split exact in the well-conditioned branch, mathematically
  (not bitwise) exact in the ill-conditioned SPH-fallback branch (ULP-scale
  only, diagnostic-only impact, `FVPM_ATTEMPT1_FACE_RESCALING` — the only
  code that reads `area_sum`/`is_problematic` in the actual solver path —
  confirmed undefined in this build), struct growth safe (no hardcoded
  `sizeof` assumptions found), sym/nonsym split safe (no double-counting).
  One critical finding: `src/cell_unskip.c`'s MPI cross-rank task-unskip
  logic (~lines 1963-1990) only walks `c->hydro.force`, never
  `c->hydro.gradient`; combining `EXTRA_HYDRO_LOOP_TYPE2` with the
  pre-existing `MPI_SYMMETRIC_FORCE_INTERACTION` in the Gizmo block (a
  combination `REMIX_SPH` does not exercise) means a boundary pair with an
  inactive local cell and an active remote cell satisfying the union
  condition only via the remote particle's `h` will silently drop that
  gradient contribution under MPI. Reproducible today: `WITH_MPI` is
  unconditionally on for the `swift_mpi` target in this checkout and a
  built `swift_mpi` binary exists.
- Neutral-arbiter arbitration: **APPROVED-WITH-CAVEATS**. Confirmed the
  critical finding firsthand (not just re-quoting the report), and
  corrected the red-team's proposed fix as insufficient on its own: the
  four `MPI_SYMMETRIC_FORCE_INTERACTION` inactive-local branches
  (`cell_unskip.c` ~1746-1761, 1795-1808, 1845-1854, 1889-1902) currently
  only exchange the finished `task_subtype_gradient` result under
  `EXTRA_HYDRO_LOOP` (never `task_subtype_rho`, by an explicit NOTE at
  ~line 1748), which was correct only because gradient was gather-type
  before this change. Simply unskipping the gradient task cross-rank
  without also exchanging `task_subtype_rho` there would run the
  newly-unskipped remote gradient task against particles whose density
  output was never received: a worse bug than the one being fixed. Also
  confirmed: no `SWIFT_DEBUG_CHECKS`/`DO_DRIFT_DEBUG_CHECKS` assertion
  would catch this (silent wrongness, not a crash); the >80-column
  comment-style finding is real but non-blocking (file already had
  pre-existing violations, no CI enforcement); the untracked
  HDF5/`csds`-submodule/generated-output files in the working tree are
  unrelated to the diff and must not be swept into this commit.
- `pdflatex` build of the companion LaTeX write-up unaffected by this gate
  (no LaTeX changes in this round).

# Decisions

- Commit only the four reviewed files by explicit path (`git add
  src/part.h src/fvpm_geometry/Gizmo/fvpm_geometry.h
  src/fvpm_geometry/Gizmo/fvpm_geometry_struct.h
  src/fvpm_geometry/None/fvpm_geometry.h`, plus the two `.claude/dev/logs`
  entries), never `git add -A`/`git commit -a`, per the arbiter's explicit
  instruction, given the untracked HDF5 binaries and dirty `csds`
  submodule pointer sitting in the same working tree.
- Record the MPI limitation and the corrected two-part fast-follow fix
  (gradient unskip AND `task_subtype_rho` exchange) in the commit message
  itself, verbatim per the arbiter's required text, rather than only in
  this log, so the limitation travels with `git log`/`git blame` even if
  this log file is later moved or pruned.
- Did not apply the optional `clang-format` cleanup of the four new
  >80-column comment lines in this commit; non-blocking per the arbiter,
  left as a follow-up before this branch is pushed for shared review.
- Did not implement the `cell_unskip.c` fast-follow fix in this task: out
  of scope for Stage 0 (single-rank diagnostic work), and the operator
  explicitly said "we will take care of MPI at a later point."

# Rejected alternatives

- Applying the red-team's literally-quoted one-line fix ("add a
  `c->hydro.gradient` walk") was considered and rejected by the
  arbiter as insufficient by itself; a correct fix needs the paired
  `task_subtype_rho` exchange fix too. Not attempted in this task; left
  as a fully-specified fast-follow for whoever picks up the MPI fix.

# Human interventions

- Operator explicitly deferred the MPI fix ("We will take care of MPI at
  a later point") rather than blocking on it, after seeing the red-team
  finding relayed in conversation. This authorizes committing Stage 0's
  Gizmo change with the limitation documented, rather than reverting or
  gating further on an MPI fix.

# Open questions

- The `cell_unskip.c` fast-follow (gradient-loop cross-rank unskip +
  `task_subtype_rho` exchange fix) remains unimplemented; needs its own
  task before any multi-rank Gizmo production use of this branch.
