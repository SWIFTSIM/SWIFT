/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* T3: task-set completeness for the fixed-aperture gas-gas sink formation
 * preparation loop (Fix A: grid sizing, Fix B: shared split predicate).
 *
 * Building the full engine_maketasks/scheduler wiring in a unit test is too
 * heavy (see sink_formation_gas_fix_plan.org, T3), so this asserts the two
 * pure decisions Fix A and Fix B are actually responsible for, exercising
 * the real production code paths rather than re-deriving the formulas:
 *
 *  (a) Fix A -- space_init() floors s->cell_min to (space_stretch * r_cut)
 *      when the fixed-aperture loop is active, independent of whether any
 *      sink particle exists in the run (with_sink=0 here). Run with
 *      dry_run=1 so the (heavy, IC-dependent) rest of space_init/space_regrid
 *      is skipped; s->cell_min is set unconditionally before that point.
 *
 *  (b) Fix B -- cell_can_split_{pair,self}_hydro_task() stop splitting once
 *      the aperture no longer fits within half a sub-cell, exactly like the
 *      pre-existing h_max terms, so the shared decomposition (and hence
 *      hydro.super) respects r_cut and the piggybacked formation_gas leaves
 *      stay complete.
 */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <string.h>

/* Local headers. */
#include "cell.h"
#include "cosmology.h"
#include "error.h"
#include "hydro_properties.h"
#include "parser.h"
#include "sink_properties.h"
#include "space.h"

/* ============================================================
 * (a) Fix A -- grid floor via space_init(dry_run=1)
 * ============================================================ */

static void test_grid_floor(void) {

  struct swift_params params;
  parser_init("testSinkFormationGasCompleteness-in-memory", &params);

  struct hydro_props hp;
  hydro_props_init_no_hydro(&hp);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  double dim[3] = {10., 10., 10.};

  /* --- Gate ON: fixed aperture active, r_cut larger than the naturally
   * implied grid spacing (tol * dmax / maxtcells ~= 0.99 * 10 / 12 ~= 0.82).
   * Expect s->cell_min to be exactly floored to space_stretch * r_cut. --- */
  {
    struct sink_props sink_properties;
    memset(&sink_properties, 0, sizeof(struct sink_props));
    sink_properties.use_fixed_r_cut = 1;
    sink_properties.cut_off_radius = 3.0f;

    struct space s;
    space_init(&s, &params, &cosmo, dim, &hp, &sink_properties,
               /*parts=*/NULL, /*gparts=*/NULL, /*sinks=*/NULL,
               /*sparts=*/NULL, /*bparts=*/NULL, /*Npart=*/0, /*Ngpart=*/0,
               /*Nsink=*/0, /*Nspart=*/0, /*Nbpart=*/0, /*Nnupart=*/0,
               /*periodic=*/1, /*replicate=*/1, /*remap_ids=*/0,
               /*generate_gas_in_ics=*/0, /*hydro=*/1, /*self_gravity=*/0,
               /*star_formation=*/0, /*with_sink=*/0, /*with_DM=*/0,
               /*with_DM_background=*/0, /*neutrinos=*/0, /*verbose=*/0,
               /*dry_run=*/1, /*nr_nodes=*/1);

    const double expected =
        (double)(space_stretch * sink_properties.cut_off_radius);
    if (fabs(s.cell_min - expected) > 1e-6)
      error(
          "Fix A: gate ON, cell_min=%.6f does not match the expected r_cut "
          "floor %.6f.",
          s.cell_min, expected);
    if (s.cell_min < sink_properties.cut_off_radius)
      error(
          "Fix A: gate ON, cell_min=%.6f is smaller than r_cut=%.6f -- the "
          "immediate-neighbour task stencil would be incomplete.",
          s.cell_min, sink_properties.cut_off_radius);

    message("Fix A (gate ON) passed: cell_min=%.6f >= r_cut=%.3f.", s.cell_min,
            sink_properties.cut_off_radius);
  }

  /* --- Gate OFF: fixed aperture not active (use_fixed_r_cut=0). Expect the
   * natural, r_cut-independent cell_min (no flooring applied). --- */
  {
    struct sink_props sink_properties;
    memset(&sink_properties, 0, sizeof(struct sink_props));
    sink_properties.use_fixed_r_cut = 0;
    sink_properties.cut_off_radius = 3.0f;

    struct space s;
    space_init(&s, &params, &cosmo, dim, &hp, &sink_properties,
               /*parts=*/NULL, /*gparts=*/NULL, /*sinks=*/NULL,
               /*sparts=*/NULL, /*bparts=*/NULL, /*Npart=*/0, /*Ngpart=*/0,
               /*Nsink=*/0, /*Nspart=*/0, /*Nbpart=*/0, /*Nnupart=*/0,
               /*periodic=*/1, /*replicate=*/1, /*remap_ids=*/0,
               /*generate_gas_in_ics=*/0, /*hydro=*/1, /*self_gravity=*/0,
               /*star_formation=*/0, /*with_sink=*/0, /*with_DM=*/0,
               /*with_DM_background=*/0, /*neutrinos=*/0, /*verbose=*/0,
               /*dry_run=*/1, /*nr_nodes=*/1);

    if (s.cell_min >= sink_properties.cut_off_radius)
      error(
          "Fix A: gate OFF, cell_min=%.6f was floored to (or above) "
          "r_cut=%.6f even though the fixed-aperture loop is inactive.",
          s.cell_min, sink_properties.cut_off_radius);

    message(
        "Fix A (gate OFF) passed: cell_min=%.6f is NOT floored by r_cut=%.3f.",
        s.cell_min, sink_properties.cut_off_radius);
  }
}

/* ============================================================
 * (b) Fix B -- cell_can_split_{pair,self}_hydro_task() gate behaviour
 * ============================================================ */

static struct cell make_bare_split_cell(int split, double dmin) {
  struct cell c;
  memset(&c, 0, sizeof(struct cell));
  c.split = split;
  c.dmin = dmin;
  return c;
}

static void test_split_predicate(void) {

  /* Boundary: space_stretch * r_cut == 0.5 * dmin, i.e.
   * dmin_boundary = 2 * space_stretch * r_cut. */
  const float r_cut = 0.1f;
  const double dmin_boundary = 2.0 * space_stretch * r_cut;

  const double dmin_above = dmin_boundary + 0.05; /* can still split */
  const double dmin_below = dmin_boundary - 0.05; /* cannot split anymore */

  /* --- Split cell, r_cut well within half a sub-cell: can split. --- */
  {
    struct cell c = make_bare_split_cell(/*split=*/1, dmin_above);
    if (!cell_can_split_pair_hydro_task(&c, r_cut))
      error(
          "Fix B: above the r_cut boundary, PAIR predicate incorrectly "
          "returned 'cannot split'.");
    if (!cell_can_split_self_hydro_task(&c, r_cut))
      error(
          "Fix B: above the r_cut boundary, SELF predicate incorrectly "
          "returned 'cannot split'.");
  }

  /* --- Split cell, exactly at the boundary: strict '<', so cannot split. ---
   */
  {
    struct cell c = make_bare_split_cell(/*split=*/1, dmin_boundary);
    if (cell_can_split_pair_hydro_task(&c, r_cut))
      error(
          "Fix B: at the r_cut boundary, PAIR predicate incorrectly "
          "returned 'can split' (should be a strict inequality).");
    if (cell_can_split_self_hydro_task(&c, r_cut))
      error(
          "Fix B: at the r_cut boundary, SELF predicate incorrectly "
          "returned 'can split' (should be a strict inequality).");
  }

  /* --- Split cell, r_cut no longer fits within half a sub-cell: cannot
   * split -- this is the actual Fix B property: the shared decomposition
   * (and hence hydro.super) stops at the r_cut level. --- */
  {
    struct cell c = make_bare_split_cell(/*split=*/1, dmin_below);
    if (cell_can_split_pair_hydro_task(&c, r_cut))
      error(
          "Fix B: below the r_cut boundary, PAIR predicate incorrectly "
          "returned 'can split'.");
    if (cell_can_split_self_hydro_task(&c, r_cut))
      error(
          "Fix B: below the r_cut boundary, SELF predicate incorrectly "
          "returned 'can split'.");
  }

  /* --- Non-split (leaf) cell: must never report 'can split', regardless of
   * r_cut. --- */
  {
    struct cell c = make_bare_split_cell(/*split=*/0, dmin_above);
    if (cell_can_split_pair_hydro_task(&c, r_cut))
      error(
          "Fix B: a non-split cell's PAIR predicate incorrectly returned "
          "'can split'.");
    if (cell_can_split_self_hydro_task(&c, r_cut))
      error(
          "Fix B: a non-split cell's SELF predicate incorrectly returned "
          "'can split'.");
  }

  /* --- Gate inactive (r_cut = 0, the documented "no fixed aperture"
   * sentinel): behaves exactly as the pre-Fix-B predicate (h_max terms
   * only; all zero here), so a split cell can always split regardless of
   * dmin. --- */
  {
    struct cell c = make_bare_split_cell(/*split=*/1, dmin_below);
    if (!cell_can_split_pair_hydro_task(&c, 0.f))
      error(
          "Fix B: with the gate inactive (r_cut=0), PAIR predicate "
          "incorrectly returned 'cannot split'.");
    if (!cell_can_split_self_hydro_task(&c, 0.f))
      error(
          "Fix B: with the gate inactive (r_cut=0), SELF predicate "
          "incorrectly returned 'cannot split'.");
  }

  message("Fix B (cell_can_split_{pair,self}_hydro_task r_cut gate) passed.");
}

int main(int argc, char *argv[]) {

  test_grid_floor();
  test_split_predicate();

  message("All T3 (task-set completeness) checks passed.");

  return 0;
}
