/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026 Will J. Roper (w.roper@sussex.ac.uk).
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

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "swift.h"

/* Standard headers. */
#include <math.h>
#include <stdlib.h>
#include <string.h>

/**
 * @brief Set up a gravity cell with one particle and a multipole.
 *
 * @param c The #cell pointer to set up.
 * @param loc The lower corner of the cell.
 * @param width The width of the cell.
 * @param x The position of the particle in the cell.
 */
static void setup_cell(struct cell *c, const double loc[3],
                       const double width[3], const double x[3]) {

  bzero(c, sizeof(struct cell));

  for (int i = 0; i < 3; ++i) {
    c->loc[i] = loc[i];
    c->width[i] = width[i];
  }

  c->grav.count = 1;
  c->grav.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  if (c->grav.multipole == NULL) error("Failed to allocate multipole");
  bzero(c->grav.multipole, sizeof(struct gravity_tensors));

  if (posix_memalign((void **)&c->grav.parts, gpart_align,
                     c->grav.count * sizeof(struct gpart)) != 0)
    error("Error allocating gparts");
  bzero(c->grav.parts, c->grav.count * sizeof(struct gpart));

  for (int i = 0; i < 3; ++i) c->grav.parts[0].x[i] = x[i];
  c->grav.parts[0].mass = 1.0;
  c->grav.parts[0].time_bin = 1;
  c->grav.parts[0].type = swift_type_dark_matter;
}

/**
 * @brief Free the memory allocated for a test cell.
 *
 * @param c The #cell pointer to free.
 */
static void free_cell(struct cell *c) {

  free(c->grav.multipole);
  free(c->grav.parts);
}

/**
 * @brief Compute the analytic zoom mesh correction kernel at a separation.
 *
 * This matches the real-space correction Green's function used by the Hockney
 * zero-padded zoom mesh solve.
 *
 * @param r The separation.
 * @param r_s_zoom The zoom mesh smoothing scale.
 * @param r_s_global The regular PM mesh smoothing scale.
 * @return The correction potential for unit mass and G=1.
 */
static double zoom_mesh_test_kernel(const double r, const double r_s_zoom,
                                    const double r_s_global) {

  if (r == 0.)
    return -M_2_SQRTPI * 0.5 / r_s_zoom +
           M_2_SQRTPI * 0.5 / r_s_global;

  return -erf(0.5 * r / r_s_zoom) / r +
         erf(0.5 * r / r_s_global) / r;
}

/**
 * @brief Compute the expected mesh acceleration from the analytic kernel.
 *
 * This uses the same five-point stencil as the zoom mesh interpolation routine,
 * evaluated for a source and probe particle placed exactly on mesh nodes.
 *
 * @param cell_fac The mesh cell conversion factor.
 * @param source_i The source mesh index along x.
 * @param probe_i The probe mesh index along x.
 * @param r_s_zoom The zoom mesh smoothing scale.
 * @param r_s_global The regular PM mesh smoothing scale.
 * @return The expected x acceleration for unit mass and G=1.
 */
static double zoom_mesh_expected_acceleration_x(const double cell_fac,
                                                const int source_i,
                                                const int probe_i,
                                                const double r_s_zoom,
                                                const double r_s_global) {

  double phi[5];
  for (int n = -2; n <= 2; ++n) {
    const double dx = (probe_i + n - source_i) / cell_fac;
    phi[n + 2] = zoom_mesh_test_kernel(fabs(dx), r_s_zoom, r_s_global);
  }

  double a = 0.;
  a += (1. / 12.) * phi[4];
  a -= (2. / 3.) * phi[3];
  a += (2. / 3.) * phi[1];
  a -= (1. / 12.) * phi[0];

  return cell_fac * a;
}

/**
 * @brief Test the zoom mesh pair-pruning predicates.
 *
 * This constructs a mock periodic engine with both a regular PM mesh and a
 * zoom mesh. The regular PM cutoff is deliberately set too large to prune the
 * pair so that #cell_can_use_mesh must rely on the zoom mesh fallback. The test
 * then checks inside/outside bounds handling and the between-rebuild dx_max
 * safety path.
 */
int main(int argc, char *argv[]) {

  (void)argc;
  (void)argv;

  struct engine e;
  struct space s;
  struct pm_mesh mesh;
  struct zoom_pm_mesh zoom_mesh;
  struct gravity_props props;
  struct phys_const phys_const;
  bzero(&e, sizeof(struct engine));
  bzero(&s, sizeof(struct space));
  bzero(&mesh, sizeof(struct pm_mesh));
  bzero(&zoom_mesh, sizeof(struct zoom_pm_mesh));
  bzero(&props, sizeof(struct gravity_props));
  bzero(&phys_const, sizeof(struct phys_const));

  s.periodic = 1;
  s.dim[0] = 100.;
  s.dim[1] = 100.;
  s.dim[2] = 100.;

  mesh.periodic = 1;
  mesh.dim[0] = s.dim[0];
  mesh.dim[1] = s.dim[1];
  mesh.dim[2] = s.dim[2];
  mesh.r_cut_max = 50.;
  mesh.r_cut_min = 0.;
  mesh.r_s = 10.;
  mesh.r_s_inv = 0.1;

  zoom_mesh.enabled = 1;
  zoom_mesh.N[0] = 16;
  zoom_mesh.N[1] = 16;
  zoom_mesh.N[2] = 16;
  zoom_mesh.loc[0] = 10.;
  zoom_mesh.loc[1] = 10.;
  zoom_mesh.loc[2] = 10.;
  zoom_mesh.dim[0] = 30.;
  zoom_mesh.dim[1] = 30.;
  zoom_mesh.dim[2] = 30.;
  zoom_mesh.cell_fac[0] = zoom_mesh.N[0] / zoom_mesh.dim[0];
  zoom_mesh.cell_fac[1] = zoom_mesh.N[1] / zoom_mesh.dim[1];
  zoom_mesh.cell_fac[2] = zoom_mesh.N[2] / zoom_mesh.dim[2];
  zoom_mesh.r_cut_max = 2.;
  zoom_mesh.r_cut_min = 0.;
  zoom_mesh.r_s = 1.;
  zoom_mesh.r_s_inv = 1.;

  props.epsilon_DM_cur = 0.1f;
  props.epsilon_baryon_cur = 0.1f;
  phys_const.const_newton_G = 1.;

  e.s = &s;
  e.mesh = &mesh;
  e.zoom_mesh = &zoom_mesh;
  e.gravity_properties = &props;
  e.physical_constants = &phys_const;
  s.e = &e;

  const double width[3] = {1., 1., 1.};
  const double loc_i[3] = {16., 16., 16.};
  const double loc_j[3] = {21., 16., 16.};
  const double x_i[3] = {16.5, 16.5, 16.5};
  const double x_j[3] = {21.5, 16.5, 16.5};

  struct cell ci;
  struct cell cj;
  setup_cell(&ci, loc_i, width, x_i);
  setup_cell(&cj, loc_j, width, x_j);

  cell_make_multipoles(&ci, /*ti_current=*/0, &props);
  cell_make_multipoles(&cj, /*ti_current=*/0, &props);

  if (!zoom_mesh_cell_is_covered(&zoom_mesh, &ci, /*use_max_dx=*/0))
    error("ci should be covered by the zoom mesh");
  if (!zoom_mesh_cell_is_covered(&zoom_mesh, &cj, /*use_max_dx=*/0))
    error("cj should be covered by the zoom mesh");

  if (!zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj))
    error("Cells should be able to use the zoom mesh");

  if (!cell_can_use_mesh(&e, &ci, &cj))
    error("cell_can_use_mesh should fall back to the zoom mesh");

  /* Move cj outside the zoom mesh and check the predicate rejects it. */
  cj.loc[0] = 39.5;
  cell_make_multipoles(&cj, /*ti_current=*/1, &props);
  if (zoom_mesh_cell_is_covered(&zoom_mesh, &cj, /*use_max_dx=*/0))
    error("cj should not be covered by the zoom mesh");
  if (zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj))
    error("Cells should not use the zoom mesh when one cell is outside");

  /* Restore cj inside the mesh and use dx_max to make it unsafe between
   * rebuilds. */
  cj.loc[0] = loc_j[0];
  cell_make_multipoles(&cj, /*ti_current=*/2, &props);
  ci.grav.multipole->dx_max[0] = 2.5f;
  cj.grav.multipole->dx_max[0] = 2.5f;

  if (!zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj))
    error("Cells should use the zoom mesh at rebuild before dx_max inflation");
  if (zoom_mesh_can_use_mesh_between_rebuilds(&zoom_mesh, &s, &ci, &cj))
    error("dx_max should make the zoom mesh unsafe between rebuilds");

  /* Check the Hockney solve against the analytic correction kernel sampled
   * with the same finite-difference stencil. */
  const int source_i = 6;
  const int probe_i = 9;
  const double y_node = zoom_mesh.loc[1] + source_i / zoom_mesh.cell_fac[1];
  const double z_node = zoom_mesh.loc[2] + source_i / zoom_mesh.cell_fac[2];

  struct cell mesh_cells[2];
  int local_cells_with_particles_top[2] = {0, 1};
  const double source_x[3] = {zoom_mesh.loc[0] + source_i / zoom_mesh.cell_fac[0],
                              y_node, z_node};
  const double probe_x[3] = {zoom_mesh.loc[0] + probe_i / zoom_mesh.cell_fac[0],
                             y_node, z_node};
  const double source_loc[3] = {source_x[0] - 0.5, source_x[1] - 0.5,
                                source_x[2] - 0.5};
  const double probe_loc[3] = {probe_x[0] - 0.5, probe_x[1] - 0.5,
                               probe_x[2] - 0.5};

  setup_cell(&mesh_cells[0], source_loc, width, source_x);
  setup_cell(&mesh_cells[1], probe_loc, width, probe_x);
  mesh_cells[1].grav.parts[0].mass = 0.;

  s.cells_top = mesh_cells;
  s.nr_local_cells_with_particles = 2;
  s.local_cells_with_particles_top = local_cells_with_particles_top;

  threadpool_init(&e.threadpool, /*num_threads=*/1);
  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);

  const double expected_ax = zoom_mesh_expected_acceleration_x(
      zoom_mesh.cell_fac[0], source_i, probe_i, zoom_mesh.r_s, mesh.r_s);
  const double measured_ax = mesh_cells[1].grav.parts[0].a_grav_mesh[0];
  const double tolerance = 1e-7 * max(1., fabs(expected_ax));

  if (fabs(measured_ax - expected_ax) > tolerance)
    error(
        "Zoom mesh Hockney force mismatch: measured=%e expected=%e diff=%e "
        "tol=%e",
        measured_ax, expected_ax, measured_ax - expected_ax, tolerance);
  if (fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[1]) > tolerance ||
      fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[2]) > tolerance)
    error("Zoom mesh force should be aligned with the x-axis");

  zoom_mesh_clean(&zoom_mesh);
  threadpool_clean(&e.threadpool);

  free_cell(&mesh_cells[0]);
  free_cell(&mesh_cells[1]);

  free_cell(&ci);
  free_cell(&cj);

  return 0;
}
