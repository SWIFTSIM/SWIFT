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
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @brief Set up a cell with a singular particle.
 *
 * @param c The #cell pointer to set up.
 * @param loc The location of the cell.
 * @param width The width of the cell.
 * @param x The position of the particle in the cell.
 */
static void setup_cell(struct cell *c, const double loc[3],
                       const double width[3], const double x[3]) {
  bzero(c, sizeof(struct cell));

  c->loc[0] = loc[0];
  c->loc[1] = loc[1];
  c->loc[2] = loc[2];
  c->width[0] = width[0];
  c->width[1] = width[1];
  c->width[2] = width[2];

  c->grav.count = 1;

  c->grav.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  if (c->grav.multipole == NULL) error("Failed to allocate multipole");
  bzero(c->grav.multipole, sizeof(struct gravity_tensors));

  if (posix_memalign((void **)&c->grav.parts, gpart_align,
                     c->grav.count * sizeof(struct gpart)) != 0)
    error("Error allocating gparts");
  bzero(c->grav.parts, c->grav.count * sizeof(struct gpart));

  c->grav.parts[0].x[0] = x[0];
  c->grav.parts[0].x[1] = x[1];
  c->grav.parts[0].x[2] = x[2];
  c->grav.parts[0].mass = 1.0;
  c->grav.parts[0].time_bin = 1;
  c->grav.parts[0].type = swift_type_dark_matter;
}

/**
 * @brief Free the memory allocated for a cell.
 *
 * @param c The #cell pointer to free.
 */
static void free_cell(struct cell *c) {
  free(c->grav.multipole);
  free(c->grav.parts);
}

/**
 * @brief Test the mesh rebuild criterion checks.
 *
 * This test builds two leaf cells with one particle each, builds their
 * multipoles, and checks that the mesh criterion is satisfied at rebuild and
 * every value is as expected. It then drifts the multipoles to accumulate
 * to modify the dx_max variables used in the mesh criterion and checks that
 * the criterion is correctly flagged after the drift. Finally, it rebuilds the
 * multipoles again to check that the dx_max variables are cleared at rebuild
 * for the next criterion check.
 */
int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  /* Set up an engine with a periodic space, a PM mesh, and grav properties. */
  struct engine e;
  struct space s;
  struct pm_mesh mesh;
  struct gravity_props props;
  bzero(&e, sizeof(struct engine));
  bzero(&s, sizeof(struct space));
  bzero(&mesh, sizeof(struct pm_mesh));
  bzero(&props, sizeof(struct gravity_props));

  s.periodic = 1;
  s.dim[0] = 10.0;
  s.dim[1] = 10.0;
  s.dim[2] = 10.0;

  mesh.periodic = 1;
  mesh.dim[0] = s.dim[0];
  mesh.dim[1] = s.dim[1];
  mesh.dim[2] = s.dim[2];
  mesh.r_cut_min = 0.0;
  mesh.r_cut_max = 2.5;
  mesh.r_s = 1.0;
  mesh.r_s_inv = 1.0;

  props.epsilon_DM_cur = 0.1f;
  props.epsilon_baryon_cur = 0.1f;

  e.s = &s;
  e.mesh = &mesh;
  e.gravity_properties = &props;

  /* The cell pair to check. */
  struct cell ci;
  struct cell cj;

  const double width[3] = {1.0, 1.0, 1.0};
  const double loc_i[3] = {0.0, 0.0, 0.0};
  const double loc_j[3] = {3.4, 0.0, 0.0};
  const double x_i[3] = {0.5, 0.5, 0.5};
  const double x_j[3] = {3.5, 0.5, 0.5};

  /* Build two leaf cells with one particle each for mesh criterion checks. */
  setup_cell(&ci, loc_i, width, x_i);
  setup_cell(&cj, loc_j, width, x_j);

  /* Seed dx_max with junk to ensure rebuild zeroes it. */
  ci.grav.multipole->dx_max[0] = 1.0f;
  ci.grav.multipole->dx_max[1] = 2.0f;
  ci.grav.multipole->dx_max[2] = 3.0f;
  cj.grav.multipole->dx_max[0] = 4.0f;
  cj.grav.multipole->dx_max[1] = 5.0f;
  cj.grav.multipole->dx_max[2] = 6.0f;

  /* Rebuild multipoles as done during space_rebuild (but before mesh criterion
   * checks.) */
  cell_make_multipoles(&ci, /*ti_current=*/0, &props);
  cell_make_multipoles(&cj, /*ti_current=*/0, &props);

  /* Confirm the multipoles were built. */
  if (ci.grav.multipole->m_pole.M_000 <= 0.) error("ci multipole mass not set");
  if (cj.grav.multipole->m_pole.M_000 <= 0.) error("cj multipole mass not set");

  /* Check that the dx_max variables were cleared at rebuild. */
  for (int k = 0; k < 3; ++k) {
    if (ci.grav.multipole->dx_max[k] != 0.f)
      error("dx_max not zeroed for ci: %e", ci.grav.multipole->dx_max[k]);
    if (cj.grav.multipole->dx_max[k] != 0.f)
      error("dx_max not zeroed for cj: %e", cj.grav.multipole->dx_max[k]);
  }

  /* At rebuild, the minimal distance with dx_max should be below the mesh
   * cutoff so the mesh criterion is satisfied. */
  const double max_distance2 = mesh.r_cut_max * mesh.r_cut_max;
  const double min_radius2_rebuild =
      cell_min_dist2_with_max_dx(&ci, &cj, s.periodic, s.dim);
  if (min_radius2_rebuild >= max_distance2)
    error(
        "Cells should be able to use the mesh at rebuild: min_radius2=%e "
        "max_distance2=%e",
        min_radius2_rebuild, max_distance2);

  /* Drift the multipoles to accumulate dx_max and re-evaluate the criterion. */
  ci.grav.multipole->m_pole.vel[0] = 0.1f;
  ci.grav.multipole->m_pole.max_delta_vel[0] = 0.2f;
  ci.grav.multipole->m_pole.min_delta_vel[0] = -0.2f;
  cj.grav.multipole->m_pole.vel[1] = -0.1f;
  cj.grav.multipole->m_pole.max_delta_vel[1] = 0.3f;
  cj.grav.multipole->m_pole.min_delta_vel[1] = -0.3f;

  /* Drift the multipoles to accumulate dx_max. */
  gravity_drift(ci.grav.multipole, 1.0);
  gravity_drift(cj.grav.multipole, 2.0);

  /* Check that the dx_max variables were updated by the drift. */
  if (ci.grav.multipole->dx_max[0] <= 0.f)
    error("dx_max did not increase after drift for ci");
  if (cj.grav.multipole->dx_max[1] <= 0.f)
    error("dx_max did not increase after drift for cj");

  /* After drift, the minimal distance with dx_max should exceed the mesh
   * cutoff so the mesh criterion is no longer satisfied. */
  const double min_radius2_drift =
      cell_min_dist2_with_max_dx(&ci, &cj, s.periodic, s.dim);
  if (min_radius2_drift <= max_distance2)
    error(
        "dx_max large enough to need pair task but mesh rebuild not flagged: "
        "min_radius2=%e max_distance2=%e",
        min_radius2_drift, max_distance2);

  /* Rebuild again and ensure the dx_max variables are cleared. */
  cell_make_multipoles(&ci, /*ti_current=*/1, &props);
  cell_make_multipoles(&cj, /*ti_current=*/1, &props);
  for (int k = 0; k < 3; ++k) {
    if (ci.grav.multipole->dx_max[k] != 0.f)
      error("dx_max not zeroed after rebuild for ci: %e",
            ci.grav.multipole->dx_max[k]);
    if (cj.grav.multipole->dx_max[k] != 0.f)
      error("dx_max not zeroed after rebuild for cj: %e",
            cj.grav.multipole->dx_max[k]);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Sanity-check the rebuild invariants. */
  cell_check_multipole(&ci, &props);
  cell_check_multipole(&cj, &props);
#endif

  free_cell(&ci);
  free_cell(&cj);

  return 0;
}
