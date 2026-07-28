/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include <config.h>

/* Some standard headers. */
#include <signal.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"

/* Not declared in cell.h (matches cell_set_super_hydro/_gravity, which are
 * likewise only called from within cell.c) -- declared here so the
 * orphaned-link regression test below can exercise it directly. */
void cell_set_super_radiation_subgrid(struct cell *c,
                                      struct cell *super_radiation);

/**
 * @brief Build a bare pair of cells with the fields relevant to the
 * radiation rebuild/split predicates set, everything else zeroed.
 */
static void setup_pair(struct cell *ci, struct cell *cj, float dmin,
                       float stars_h_max, float stars_h_hii_max,
                       float hydro_h_max, float stars_dx_max_part,
                       float hydro_dx_max_part) {
  bzero(ci, sizeof(struct cell));
  bzero(cj, sizeof(struct cell));

  ci->dmin = dmin;
  cj->dmin = dmin;
  ci->stars.h_max = stars_h_max;
  ci->stars.h_hii_max = stars_h_hii_max;
  cj->hydro.h_max = hydro_h_max;
  ci->stars.dx_max_part = stars_dx_max_part;
  cj->hydro.dx_max_part = hydro_dx_max_part;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  struct cell ci, cj;

  /* Case 1: everything comfortably smaller than the cell -- no rebuild
   * needed for either the plain stars check or the radiation-aware one. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/0.01f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "Unexpected rebuild flagged by cell_need_rebuild_for_stars_pair "
        "for a small, static configuration.");
  if (cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "Unexpected rebuild flagged by cell_need_rebuild_for_radiation_pair "
        "for a small, static configuration.");

  /* Case 2: h_hii_max alone grown large enough to exceed the cell size,
   * while ordinary h_max/dx_max stay small. This is the scenario the two
   * checks must disagree on: the plain stars check has no h_hii term and
   * must stay silent, while the radiation-aware check must flag a
   * rebuild. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/2.0f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_stars_pair must not react to h_hii_max "
        "growth (it has no h_hii term by design).");
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after h_hii_max grew past the cell size.");

  /* Case 3: ordinary h_max grown large (not h_hii) -- both checks must
   * agree and flag a rebuild, since this is the term they share. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/2.0f,
             /*stars_h_hii_max=*/0.0f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
  if (!cell_need_rebuild_for_stars_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_stars_pair failed to flag a rebuild after "
        "stars.h_max grew past the cell size.");
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after stars.h_max grew past the cell size.");

  /* Case 4: drift alone (dx_max_part) pushes the pair over the edge, with
   * h_max/h_hii_max both otherwise small. */
  setup_pair(&ci, &cj, /*dmin=*/1.0f, /*stars_h_max=*/0.01f,
             /*stars_h_hii_max=*/0.01f, /*hydro_h_max=*/0.01f,
             /*stars_dx_max_part=*/0.6f, /*hydro_dx_max_part=*/0.6f);
  if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
    error(
        "cell_need_rebuild_for_radiation_pair failed to flag a rebuild "
        "after dx_max_part drift pushed the pair past the cell size.");

  /* Case 5: h_hii_max sits in the 1.0x-radiation_search_radius_factor(x)
   * band -- inside the search radius the runner actually uses
   * (radiation_search_radius_factor * kernel_gamma * h_hii_max), but
   * outside the unfactored kernel_gamma * h_hii_max the check used before
   * the fold-in fix. Midpoint between 1.0 and the factor, so the case sits
   * comfortably inside the band regardless of the factor's exact value. */
  {
    const float mid_factor = 0.5f * (1.0f + radiation_search_radius_factor);
    const float dmin = 1.0f;
    const float h_hii_band = dmin / (mid_factor * kernel_gamma);
    setup_pair(&ci, &cj, dmin, /*stars_h_max=*/0.01f,
               /*stars_h_hii_max=*/h_hii_band, /*hydro_h_max=*/0.01f,
               /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
    if (kernel_gamma * h_hii_band >= dmin)
      error(
          "Test setup error: h_hii_band must sit BELOW the unfactored "
          "threshold (kernel_gamma * h_hii_max < dmin) for this case to "
          "discriminate the fold-in fix.");
    if (kernel_gamma * radiation_search_radius_factor * h_hii_band <= dmin)
      error(
          "Test setup error: h_hii_band must sit ABOVE the factored "
          "threshold (kernel_gamma * factor * h_hii_max > dmin) for this "
          "case to discriminate the fold-in fix.");
    if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
      error(
          "cell_need_rebuild_for_radiation_pair missed a rebuild for "
          "h_hii_max in the 1.0x-factor(x) search-radius band: the check "
          "must bound the SAME radius the runner searches with "
          "(radiation_search_radius_factor * kernel_gamma * h_hii_max), "
          "not the unfactored kernel_gamma * h_hii_max.");
  }

  /* Now check cell_can_split_pair/self_radiation_subgrid_task(): radiation
   * must NEVER split once the cell is at or below hydro.super (non-NULL),
   * for any combination of the other fields (including small h_max/
   * h_hii_max that would otherwise geometrically allow splitting). This is
   * the fix for the duplicate-unlock crash at gas_mass=0.01: once a cell
   * IS (or lies below) a hydro.super, cell_set_super_hydro() has already
   * stamped that SAME hydro.super pointer onto every one of its progeny,
   * so splitting further would only produce sibling radiation tasks that
   * all wire dependencies to the identical hydro.super. See
   * cell_can_split_pair_radiation_subgrid_task()'s docstring in cell.h for
   * the full history. */
  struct cell c;
  const float dmins[] = {0.5f, 1.0f, 2.0f};
  const float h_values[] = {0.0f, 0.01f, 0.3f, 1.0f, 5.0f};
  struct cell dummy_super;
  bzero(&dummy_super, sizeof(struct cell));
  for (int is = 0; is < 2; ++is) {
    for (int id = 0; id < 3; ++id) {
      for (int ih = 0; ih < 5; ++ih) {
        bzero(&c, sizeof(struct cell));
        c.split = is;
        c.dmin = dmins[id];
        c.hydro.h_max = h_values[ih];
        c.stars.h_max = h_values[ih];
        c.stars.h_hii_max = h_values[(ih + 1) % 5];
        c.sinks.h_max = h_values[ih];
        c.black_holes.h_max = h_values[ih];

        /* With hydro.super already set (at-or-below hydro.super), must
         * always refuse to split regardless of geometry. */
        c.hydro.super = &dummy_super;
        if (cell_can_split_pair_radiation_subgrid_task(&c))
          error(
              "cell_can_split_pair_radiation_subgrid_task must never "
              "split once hydro.super is set "
              "(split=%d dmin=%e h=%e h_hii=%e).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max);

        if (cell_can_split_self_radiation_subgrid_task(&c))
          error(
              "cell_can_split_self_radiation_subgrid_task must never "
              "split once hydro.super is set "
              "(split=%d dmin=%e h=%e h_hii=%e).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max);

        /* With hydro.super still NULL (strictly above any hydro.super),
         * splitting must be governed purely by the geometric criterion --
         * identical to what the (now-removed) NULL-gated call would have
         * returned before this test replaced it. */
        c.hydro.super = NULL;
        const int expect_geometric =
            c.split &&
            (space_stretch * kernel_gamma * c.hydro.h_max < 0.5f * c.dmin) &&
            (space_stretch * radiation_search_radius_factor * kernel_gamma *
                 c.stars.h_max <
             0.5f * c.dmin) &&
            (space_stretch * radiation_search_radius_factor * kernel_gamma *
                 c.stars.h_hii_max <
             0.5f * c.dmin) &&
            (space_stretch * kernel_gamma * c.sinks.h_max < 0.5f * c.dmin) &&
            (space_stretch * kernel_gamma * c.black_holes.h_max <
             0.5f * c.dmin);

        if (cell_can_split_pair_radiation_subgrid_task(&c) != expect_geometric)
          error(
              "cell_can_split_pair_radiation_subgrid_task disagreed with "
              "the geometric criterion while hydro.super == NULL "
              "(split=%d dmin=%e h=%e h_hii=%e expected=%d).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max,
              expect_geometric);

        if (cell_can_split_self_radiation_subgrid_task(&c) != expect_geometric)
          error(
              "cell_can_split_self_radiation_subgrid_task disagreed with "
              "the geometric criterion while hydro.super == NULL "
              "(split=%d dmin=%e h=%e h_hii=%e expected=%d).",
              c.split, c.dmin, c.hydro.h_max, c.stars.h_hii_max,
              expect_geometric);
      }
    }
  }

  /* Confirm hydro's own split criteria stay decoupled from h_hii_max:
   * growing h_hii_max alone must NOT change cell_can_split_pair_hydro_task
   * -- hydro splitting is a pure hydro concept (see cell.h). */
  bzero(&c, sizeof(struct cell));
  c.split = 1;
  c.dmin = 1.0f;
  c.hydro.h_max = 0.01f;
  c.stars.h_max = 0.01f;
  c.stars.h_hii_max = 0.01f;
  c.sinks.h_max = 0.01f;
  c.black_holes.h_max = 0.01f;
  if (!cell_can_split_pair_hydro_task(&c))
    error(
        "Expected cell_can_split_pair_hydro_task to allow splitting for "
        "a small configuration.");
  c.stars.h_hii_max = 2.0f;
  if (!cell_can_split_pair_hydro_task(&c))
    error(
        "cell_can_split_pair_hydro_task must stay decoupled from "
        "h_hii_max -- it should still allow splitting after h_hii_max "
        "alone grew past the cell size.");

  /* Regression test for cell_set_super_radiation_subgrid()'s "orphaned
   * radiation_in link" SWIFT_DEBUG_CHECKS assertion (cell.c). Invariant
   * discovered 2026-07-12: the flat radiation_in walk in
   * runner_dosub_stars_hii_ionization_feedback only ever reads
   * radiation_level's OWN link list -- it never descends into
   * radiation_level's children to pick up links attached deeper. So every
   * cell holding a radiation_in link must BE its own radiation_level.
   * This can be violated when a pair task with one neighbour stops
   * splitting early (its link lands on a coarse ancestor, e.g. because
   * that neighbour's hydro.super sits at a shallower depth) while a
   * DIFFERENT, matching-depth neighbour's pair (or the cell's own
   * self-splitting) attaches its link to a deeper descendant: the
   * "topmost link found first wins, propagated to all descendants" rule
   * in cell_set_super_radiation_subgrid then makes the deeper link
   * permanently invisible to the flat walk -- a silent, directional
   * search-completeness bug, not a crash. */
  {
    struct cell top, child;
    struct link top_link, child_link;
    bzero(&top, sizeof(struct cell));
    bzero(&child, sizeof(struct cell));
    bzero(&top_link, sizeof(struct link));
    bzero(&child_link, sizeof(struct link));

    /* Case A: no orphaning -- only the top cell holds a link (the common,
     * correct case: a pair/self task that never split below the top
     * cell). radiation_level must land on top for both cells; no error. */
    top.split = 1;
    top.progeny[0] = &child;
    top.stars.radiation_in = &top_link;
    child.stars.radiation_in = NULL;
    cell_set_super_radiation_subgrid(&top, NULL);
    if (top.stars.radiation_level != &top)
      error(
          "cell_set_super_radiation_subgrid: expected radiation_level == "
          "top for the unsplit case.");
    if (child.stars.radiation_level != &top)
      error(
          "cell_set_super_radiation_subgrid: expected the child to "
          "inherit radiation_level == top.");

    /* Case B: orphaning -- BOTH top and child hold their own link (top's
     * simulates an unsplit, mismatched-neighbour pair; child's simulates a
     * properly-split sibling/neighbour pair or self-task). radiation_level
     * is pulled to top (topmost wins), silently orphaning child's link.
     * Must be caught by the SWIFT_DEBUG_CHECKS assertion -- run in a
     * forked child so the expected abort doesn't kill the test. If this
     * test binary was built without SWIFT_DEBUG_CHECKS, the assertion
     * doesn't exist and this case is skipped (nothing to catch it with). */
#ifdef SWIFT_DEBUG_CHECKS
    child.stars.radiation_in = &child_link;

    const pid_t pid = fork();
    if (pid == 0) {
      /* Child process: silence stderr (the error() message is expected
       * output, not a test failure), then trigger the assertion. */
      if (freopen("/dev/null", "w", stderr) == NULL) _exit(43);
      cell_set_super_radiation_subgrid(&top, NULL);
      /* Reached only if the assertion did NOT fire -- signal failure with
       * a distinguishable exit code (real error() exits with status 1). */
      _exit(42);
    } else if (pid > 0) {
      int status;
      waitpid(pid, &status, 0);
      /* error() aborts via swift_abort(), which is exit(1) normally but
       * abort() (SIGABRT) under SWIFT_DEVELOP_MODE -- accept either. */
      const int exited_with_error =
          WIFEXITED(status) && WEXITSTATUS(status) == 1;
      const int aborted = WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT;
      if (!exited_with_error && !aborted)
        error(
            "cell_set_super_radiation_subgrid failed to detect an "
            "orphaned radiation_in link (WIFEXITED=%d WEXITSTATUS=%d "
            "WIFSIGNALED=%d WTERMSIG=%d) -- the SWIFT_DEBUG_CHECKS "
            "assertion did not fire as expected.",
            WIFEXITED(status), WIFEXITED(status) ? WEXITSTATUS(status) : -1,
            WIFSIGNALED(status), WIFSIGNALED(status) ? WTERMSIG(status) : -1);
    } else {
      error("fork() failed in the orphaned-link regression test.");
    }
#endif
  }

  /* Regression tests for the periodic cdim==3 rebuild-floor fix: at that
   * grid, the 27-cell radiation stencil wires every top-level pair, so
   * cell_unskip.c's call site skips cell_need_rebuild_for_radiation_pair
   * for a pair of UNSPLIT top-level cells once
   * space_radiation_top_stencil_covers_box() is true (see its docstring in
   * space.h and the guard in cell_unskip_radiation_tasks()). That predicate
   * is called directly here (not threaded through
   * cell_need_rebuild_for_radiation_pair, which stays cell.h-only and
   * space/engine-agnostic), so it is exercised as the real function, not a
   * re-implementation that could silently drift from space_regrid.c. */

  /* Case 6: space_radiation_top_stencil_covers_box() itself -- true only
   * for a periodic box at exactly cdim == {3,3,3}. */
  {
    struct space test_space;
    const struct {
      int periodic;
      int cdim[3];
      int expected;
    } flag_cases[] = {
        {1, {3, 3, 3}, 1}, {1, {4, 3, 3}, 0}, {1, {3, 4, 3}, 0},
        {1, {3, 3, 4}, 0}, {1, {2, 3, 3}, 0}, {0, {3, 3, 3}, 0},
    };
    for (size_t i = 0; i < sizeof(flag_cases) / sizeof(flag_cases[0]); ++i) {
      bzero(&test_space, sizeof(struct space));
      test_space.periodic = flag_cases[i].periodic;
      test_space.cdim[0] = flag_cases[i].cdim[0];
      test_space.cdim[1] = flag_cases[i].cdim[1];
      test_space.cdim[2] = flag_cases[i].cdim[2];
      const int computed = space_radiation_top_stencil_covers_box(&test_space);
      if (computed != flag_cases[i].expected)
        error(
            "space_radiation_top_stencil_covers_box mismatch for case %zu "
            "(periodic=%d cdim=[%d %d %d]): got %d, expected %d.",
            i, flag_cases[i].periodic, flag_cases[i].cdim[0],
            flag_cases[i].cdim[1], flag_cases[i].cdim[2], computed,
            flag_cases[i].expected);
    }
  }

  /* Case 7: the cell_unskip.c guard's full condition, mirrored here since
   * exercising the real call site needs a fully-wired engine/scheduler/task
   * graph. A pair whose h_hii_max violates the rebuild criterion must have
   * its rebuild demand suppressed ONLY when BOTH the box-covers flag is set
   * AND both cells are unsplit top-level cells (ci->top == ci); a split
   * sub-cell pair (ci->top != ci) must still demand a rebuild even with the
   * flag set, since its own smaller dmin came from a geometric split
   * decision the top-level stencil argument says nothing about. */
  {
    const float dmin = 1.0f;
    setup_pair(&ci, &cj, dmin, /*stars_h_max=*/0.01f,
               /*stars_h_hii_max=*/2.0f, /*hydro_h_max=*/0.01f,
               /*stars_dx_max_part=*/0.0f, /*hydro_dx_max_part=*/0.0f);
    if (!cell_need_rebuild_for_radiation_pair(&ci, &cj))
      error(
          "Test setup error: this pair must violate the radiation rebuild "
          "criterion on its own for the cases below to mean anything.");

    struct space test_space;
    bzero(&test_space, sizeof(struct space));
    test_space.periodic = 1;
    test_space.cdim[0] = test_space.cdim[1] = test_space.cdim[2] = 3;

    for (int box_covers = 0; box_covers <= 1; ++box_covers) {
      for (int top_level = 0; top_level <= 1; ++top_level) {
        /* box_covers toggles between the case-6-verified true/false space;
         * top_level toggles whether ci/cj present as unsplit top cells. */
        struct space *s = box_covers ? &test_space : NULL;
        ci.top = top_level ? &ci : NULL;
        cj.top = top_level ? &cj : NULL;

        const int guard_skips = s != NULL &&
                                space_radiation_top_stencil_covers_box(s) &&
                                ci.top == &ci && cj.top == &cj;
        const int would_rebuild =
            !guard_skips && (cell_need_rebuild_for_radiation_pair(&ci, &cj) ||
                             cell_need_rebuild_for_radiation_pair(&cj, &ci));

        const int expect_skip = box_covers && top_level;
        if (expect_skip && would_rebuild)
          error(
              "cell_unskip_radiation_tasks's guard must not demand a "
              "rebuild for a criterion-violating pair of unsplit top-level "
              "cells once the box-covers flag is set.");
        if (!expect_skip && !would_rebuild)
          error(
              "cell_unskip_radiation_tasks's guard must still demand a "
              "rebuild for a criterion-violating pair when box_covers=%d "
              "top_level=%d (must only skip when BOTH hold).",
              box_covers, top_level);
      }
    }
  }

  return 0;
}
