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
#include <string.h>

/* Local headers. */
#include "swift.h"

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
            (space_stretch * kernel_gamma * c.stars.h_max < 0.5f * c.dmin) &&
            (space_stretch * kernel_gamma * c.stars.h_hii_max <
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

  return 0;
}
