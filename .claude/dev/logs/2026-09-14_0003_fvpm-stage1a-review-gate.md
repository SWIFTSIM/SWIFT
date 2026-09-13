---
author: Orchestrator
date: 2026-09-14T00:03:34+02:00
task: fvpm-stage1a-candidate-a-centring
model: claude-sonnet-5 (red-team-reviewer and neutral-arbiter subagents for the gate)
---

# Task summary

Records the review gate (red-team + neutral-arbiter, one revision cycle)
run on Stage 1a of the FVPM face-area-closure fix ("candidate A": centre
the interpolant so the closure defect's dominant row-sum term vanishes
exactly), before committing it. Covers items 1-7, 9, 10 of
`/home/darwinr/.claude/reports/-home-darwinr-swiftsim_fvpm_fix/2026-09-13-stage1-implementation-plan-revision4-FINAL.md`,
implemented across two Worker rounds
(`2026-09-13_2255_fvpm-stage1a-candidate-a-centring.md`,
`2026-09-13_2357_fvpm-stage1a-centring-margin-fix-and-mfm-verify.md`).

# Key changes

None beyond what the two Worker rounds already implemented (see the
referenced logs). This entry documents the gate outcome that authorizes
committing that work.

# Verification

- Orchestrator independently spot-checked before dispatching the first
  gate: exact 11-file diff scope, `swiftsim_0` untouched (status hash
  unchanged), no leftover worktrees, the `omega_prime`-capture and
  inversion-gating sequence in `fvpm_compute_volume_and_matrix` matching
  the plan's step ordering exactly, `FVPM_ATTEMPT1_FACE_RESCALING` fully
  deleted, and item 9a's antisymmetry assertion correctly placed in the
  caller (`fvpm_accumulate_total_face_area_vector_and_norm`) rather than
  as a build-breaking recursive self-call inside the `always_inline`
  `fvpm_compute_face_area_vector` — the exact mistake an earlier plan
  revision made and a prior review round on the *plan* had already
  caught before any code existed.
- Red-team review (round 1): no critical bugs. Independently re-derived
  the sign convention at all six `fvpm_accumulate_first_moment_left/right`
  call sites (including the two easy-to-miss non-symmetric density-loop
  sites) and confirmed correct; independently re-derived the Worker's
  self-reported and self-fixed item 9b sign bug and confirmed the final
  state is correct; confirmed production faces
  (`runner_iact_fluxes_common`, `rt_iact.h`'s fully separate
  implementation, the three gradient-collector blocks) still read only
  the untouched `matrix_E`, not the new centred fields — Stage 1a is
  genuinely diagnostic-only for production code paths, as designed.
  Found one real diagnostic-quality defect (`centring_margin`
  unconditionally clobbered to `0.0f` on every disable path, inherited
  verbatim from the plan itself, not a Worker deviation) and one real
  verification gap (`gizmo-mfm` and any `RT_GEAR` configuration were
  reviewed by careful reading but never actually compiled this session,
  despite both being substantively, newly affected by this diff).
- Neutral-arbiter arbitration (round 1): **REJECTED FOR REVISION**, with
  a concrete, minimal, two-item fix list (not a re-architecture).
  Independently re-verified the `centring_margin` clobber by direct
  inspection and upgraded it from "fast-follow" to "fix now," since it is
  the instrument backing an acceptance criterion the plan itself calls
  blocking (11.1 criterion 4: report which particles were disabled by a
  failed inversion or SPD-margin failure, distinct from particles that
  never reached the margin computation at all). Independently traced the
  SPHENIX+RT_GEAR write/read ordering (`fvpm_compute_volume_and_matrix`
  called unconditionally from `SPHENIX/hydro.h:647`, pre-existing,
  unchanged by this diff, before the gradient-loop read) to rule out a
  sharper "uninitialized read" hypothesis before downgrading `RT_GEAR`'s
  urgency relative to `gizmo-mfm`'s; required a `gizmo-mfm` build+run
  before resubmission, recommended (not required) an `RT_GEAR` compile
  check.
- Second Worker round applied both required fixes: `centring_margin`
  initialized to `-1.0f` (meaning "SPD margin never evaluated") before
  the `centred` gate sequence, with the disable-branch's `= 0.0f;`
  overwrite deleted so a genuinely-computed-and-failed `q` survives into
  the diagnostic snapshot; `gizmo-mfm` built and run to completion (256
  steps, level 5, `NonUniformCarthesian_1D`) with zero 9a/9b assertion
  failures and `max|S^(1)_i|/A_i = 9.05e-8`, `max|S_i|/A_i = 0.098867` —
  matching the `gizmo-mfv` verification almost exactly, confirming the
  fix behaves identically across both Gizmo hydro variants.
- Orchestrator independently re-verified before this commit: `git diff`
  shows the fix exactly as specified (grep for `centring_margin` in
  `fvpm_geometry.h` shows three occurrences — a separate, correct,
  item-3-specified zero-neighbour reset in
  `fvpm_geometry_part_has_no_neighbours` at line ~166 that was never
  supposed to change, the new `-1.0f` init at line ~223, and the real `q`
  assignment at line ~300 — confirmed the disable-branch clobber is gone
  and the zero-neighbour reset was correctly left untouched, not
  conflated with the bug that was fixed). Working tree confirmed to be
  exactly the same 11 tracked files as the original diff, no scope creep
  from either Worker round, `swiftsim_0` untouched, no leftover
  worktrees or processes from either build round.
- A build-hygiene issue was found and resolved during the second Worker
  round, unrelated to Stage 1a's correctness: the pre-existing main-tree
  `gizmo-mfv` binary hung (0% CPU, deadlocked) due to a stale
  mixed-`CPPFLAGS` build (part of the object tree predated this
  session's `SWIFT_DEBUG_CHECKS`-driven `struct part` layout change,
  causing an ABI mismatch within one binary). Root-caused, reproduced
  with the `centring_margin` fix reverted (ruling out the fix as the
  cause), and resolved with `make clean` + rebuild, re-verified with a
  full clean 497-step run.

# Decisions

- Applied the arbiter's fix list via a second Worker round and did a
  targeted orchestrator re-verification rather than dispatching a full
  third red-team round from scratch, given the required changes were
  narrow, mechanical, and precisely specified by the arbiter itself (a
  two-line diagnostic fix plus a build-and-run check, not new design
  surface) — consistent with the orchestrator protocol's guidance to
  apply a concrete, unambiguous fix list via a new Worker and re-gate,
  which for a fix this narrow this session judged to mean thorough
  independent re-verification rather than a repeated adversarial pass.
- Did not require an `RT_GEAR` compile check before this commit, per the
  arbiter's own explicit judgment that it is recommended but not
  blocking, since its only new surface (`SPHENIX/hydro_io.h`'s field
  guard) was hand-verified against the struct's own dispatch condition
  and the RT production path was independently confirmed untouched by
  this diff.
- Commit only the 11 files from the original `git diff --stat` plus the
  three log entries from this task's two Worker rounds and this gate;
  never `git add -A`, given the dirty `csds` submodule pointer and
  various untracked generated/binary artifacts sitting in the same
  working tree.

# Rejected alternatives

- A full third red-team-reviewer + neutral-arbiter round on the
  two-line fix plus the MFM verification: considered and not done, for
  the reason given in Decisions above. If a future reviewer disagrees
  with this judgment call, the fix is small enough to re-review cheaply
  on its own.

# Human interventions

- None specific to this gate; the operator's standing instructions this
  session (verify sub-agent claims against the working tree, gate
  non-trivial changes with red-team + arbiter, never commit
  autonomously without a durable log) were followed as established
  practice, not as a new correction.

# Open questions

- `RT_GEAR` configurations (`--with-hydro=gizmo-mfv --with-rt=GEAR_*` or
  `--with-hydro=sphenix --with-rt=GEAR_*`) remain reviewed-but-uncompiled
  for this diff; recommended, not required, before Stage 1b/1c/1d (which
  will touch RT's production face directly) begin.
- The `GRADIENTS_SPH` compile-time toggle (currently inactive;
  `GRADIENTS_GIZMO` is the active default) makes item 9b's debug
  assertion vacuously pass rather than actively protect, since
  `hydro_gradients_sph.h`'s own collect functions never populate the
  fields it checks. Not a correctness break, but worth a code comment or
  a fast-follow if `GRADIENTS_SPH` is ever activated on this branch.
- `fvpm_update_centroid_left`/`fvpm_update_centroid_right` are now dead
  code (uncalled from every site) but were left in place rather than
  deleted; harmless (`static inline`, no `-Wunused-function` risk), minor
  cleanup opportunity for a later pass.
