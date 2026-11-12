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
  c->top = c;

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
                     c->grav.count * sizeof(struct gpart)) != 0) {
    error("Error allocating gparts");
  }
  bzero(c->grav.parts, c->grav.count * sizeof(struct gpart));

  for (int i = 0; i < 3; ++i) {
    c->grav.parts[0].x[i] = x[i];
  }
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

#ifdef HAVE_FFTW
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

  if (r == 0.) {
    return -M_2_SQRTPI * 0.5 / r_s_zoom + M_2_SQRTPI * 0.5 / r_s_global;
  }

  return -erf(0.5 * r / r_s_zoom) / r + erf(0.5 * r / r_s_global) / r;
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
#endif

/**
 * @brief Test zoom mesh pair pruning and, with FFTW, its force correction.
 *
 * This constructs a mock periodic engine with both a regular PM mesh and a
 * zoom mesh. The regular PM cutoff is deliberately set too large to prune the
 * pair so that #cell_can_use_mesh must rely on the zoom mesh fallback. The test
 * then checks inside/outside bounds handling and the between-rebuild dx_max
 * safety path. With FFTW it also tests the isolated solve,
 * deposition methods and kernel invalidation after regridding.
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
  struct neutrino_props neutrino_props;
  struct cosmology cosmology;
  struct zoom_region_properties zoom_props;
  bzero(&e, sizeof(struct engine));
  bzero(&s, sizeof(struct space));
  bzero(&mesh, sizeof(struct pm_mesh));
  bzero(&zoom_mesh, sizeof(struct zoom_pm_mesh));
  bzero(&props, sizeof(struct gravity_props));
  bzero(&phys_const, sizeof(struct phys_const));
  bzero(&neutrino_props, sizeof(struct neutrino_props));
  bzero(&cosmology, sizeof(struct cosmology));
  bzero(&zoom_props, sizeof(struct zoom_region_properties));

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
  zoom_mesh.use_local_patches = 1;
  zoom_mesh.N = 16;
  zoom_mesh.loc[0] = 10.;
  zoom_mesh.loc[1] = 10.;
  zoom_mesh.loc[2] = 10.;
  zoom_mesh.dim[0] = 30.;
  zoom_mesh.dim[1] = 30.;
  zoom_mesh.dim[2] = 30.;
  zoom_mesh.cell_fac = zoom_mesh.N / zoom_mesh.dim[0];
  for (int i = 0; i < 3; ++i) {
    const double cell_width = 1. / zoom_mesh.cell_fac;
    const double tol = 1e-6 * zoom_mesh.dim[i];
    zoom_mesh.active_min[i] = zoom_mesh.loc[i] + 2. * cell_width - tol;
    zoom_mesh.active_max[i] =
        zoom_mesh.loc[i] + zoom_mesh.dim[i] - 3. * cell_width + tol;
  }
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
  e.neutrino_properties = &neutrino_props;
  e.cosmology = &cosmology;
  s.e = &e;
  s.with_zoom_region = 1;
  s.zoom_props = &zoom_props;

#ifdef HAVE_FFTW
  const int test_mode[3] = {3, -2, 5};
  double expected_deconvolution = 1.;
  for (int axis = 0; axis < 3; ++axis) {
    const double x = M_PI * test_mode[axis] / 32.;
    expected_deconvolution *= pow(x / sin(x), 4.);
  }
  const double measured_deconvolution =
      zoom_mesh_CIC_deconvolution(test_mode[0], test_mode[1], test_mode[2], 32);
  if (fabs(measured_deconvolution - expected_deconvolution) > 1e-12) {
    error("Shared CIC deconvolution differs from the global PM convention");
  }
#endif

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

  if (!zoom_mesh_cell_is_covered(&zoom_mesh, &ci, /*use_max_dx=*/0)) {
    error("ci should be covered by the zoom mesh");
  }
  if (!zoom_mesh_cell_is_covered(&zoom_mesh, &cj, /*use_max_dx=*/0)) {
    error("cj should be covered by the zoom mesh");
  }

  /* A leaf inside the fine mesh cannot use it when its top cell straddles
   * the mesh boundary: long-range tree interactions use that top cell. */
  struct cell straddling_top = {0};
  for (int axis = 0; axis < 3; ++axis) {
    straddling_top.loc[axis] = 12.;
    straddling_top.width[axis] = 6.;
  }
  ci.top = &straddling_top;
  if (zoom_mesh_cell_is_covered(&zoom_mesh, &ci, /*use_max_dx=*/0)) {
    error("Leaf in a partially covered top cell used the zoom mesh");
  }
  ci.top = &ci;

  if (!zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj)) {
    error("Cells should be able to use the zoom mesh");
  }

  if (!cell_can_use_mesh(&e, &ci, &cj)) {
    error("cell_can_use_mesh should fall back to the zoom mesh");
  }

  /* Move cj outside the zoom mesh and check the predicate rejects it. */
  cj.loc[0] = 39.5;
  cell_make_multipoles(&cj, /*ti_current=*/1, &props);
  if (zoom_mesh_cell_is_covered(&zoom_mesh, &cj, /*use_max_dx=*/0)) {
    error("cj should not be covered by the zoom mesh");
  }
  if (zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj)) {
    error("Cells should not use the zoom mesh when one cell is outside");
  }

  /* Restore cj inside the mesh and use dx_max to make it unsafe between
   * rebuilds. */
  cj.loc[0] = loc_j[0];
  cell_make_multipoles(&cj, /*ti_current=*/2, &props);
  ci.grav.multipole->dx_max[0] = 2.5f;
  cj.grav.multipole->dx_max[0] = 2.5f;

  if (!zoom_mesh_can_use_mesh(&zoom_mesh, &s, &ci, &cj)) {
    error("Cells should use the zoom mesh at rebuild before dx_max inflation");
  }
  if (zoom_mesh_can_use_mesh_between_rebuilds(&zoom_mesh, &s, &ci, &cj)) {
    error("dx_max should make the zoom mesh unsafe between rebuilds");
  }

#ifdef HAVE_FFTW
  /* Check the Hockney solve against the analytic correction kernel sampled
   * with the same finite-difference stencil. CIC deconvolution modifies this
   * sampled reference slightly at this deliberately small mesh size. */
  const int source_i = 6;
  const int probe_i = 9;
  const double y_node = zoom_mesh.loc[1] + source_i / zoom_mesh.cell_fac;
  const double z_node = zoom_mesh.loc[2] + source_i / zoom_mesh.cell_fac;

  struct cell mesh_cells[2];
  int local_cells_with_particles_top[2] = {0, 1};
  const double source_x[3] = {zoom_mesh.loc[0] + source_i / zoom_mesh.cell_fac,
                              y_node, z_node};
  const double probe_x[3] = {zoom_mesh.loc[0] + probe_i / zoom_mesh.cell_fac,
                             y_node, z_node};
  const double source_loc[3] = {source_x[0] - 0.5, source_x[1] - 0.5,
                                source_x[2] - 0.5};
  const double probe_loc[3] = {probe_x[0] - 0.5, probe_x[1] - 0.5,
                               probe_x[2] - 0.5};

  setup_cell(&mesh_cells[0], source_loc, width, source_x);
  setup_cell(&mesh_cells[1], probe_loc, width, probe_x);
  mesh_cells[1].grav.parts[0].mass = 0.;

  s.cells_top = mesh_cells;
  s.local_cells_with_particles_top = local_cells_with_particles_top;
  zoom_props.nr_local_zoom_cells = 1;
  zoom_props.local_zoom_cells_top = &local_cells_with_particles_top[0];
  zoom_props.nr_local_bkg_cells = 1;
  zoom_props.local_bkg_cells_top = &local_cells_with_particles_top[1];

  threadpool_init(&e.threadpool, /*num_threads=*/4);

  zoom_mesh_allocate(&zoom_mesh);

  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);

  const double expected_ax = zoom_mesh_expected_acceleration_x(
      zoom_mesh.cell_fac, source_i, probe_i, zoom_mesh.r_s, mesh.r_s);
  const double measured_ax = mesh_cells[1].grav.parts[0].a_grav_mesh[0];
  const double analytic_tolerance = 4e-3 * max(1., fabs(expected_ax));

  if (fabs(measured_ax - expected_ax) > analytic_tolerance) {
    error(
        "Zoom mesh Hockney force mismatch: measured=%e expected=%e diff=%e "
        "tol=%e",
        measured_ax, expected_ax, measured_ax - expected_ax,
        analytic_tolerance);
  }
  if (fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[1]) > analytic_tolerance ||
      fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[2]) > analytic_tolerance) {
    error("Zoom mesh force should be aligned with the x-axis");
  }

  /* Exercise the cached-kernel path and require the same force. */
  for (int axis = 0; axis < 3; ++axis) {
    mesh_cells[1].grav.parts[0].a_grav_mesh[axis] = 0.f;
  }
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  mesh_cells[1].grav.parts[0].potential_mesh = 0.f;
#endif
  zoom_mesh.use_local_patches = 0;
  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);
  if (fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[0] - measured_ax) > 1e-7) {
    error("Atomic deposition or cached kernel changed the computed force");
  }
  zoom_mesh.use_local_patches = 1;

  /* Shift both particles off the nodes to exercise fractional CIC assignment
   * and interpolation with the cached, deconvolved kernel. */
  const double shift = 0.37 / zoom_mesh.cell_fac;
  for (int cell = 0; cell < 2; ++cell) {
    mesh_cells[cell].loc[0] += shift;
    mesh_cells[cell].grav.parts[0].x[0] += shift;
  }
  for (int axis = 0; axis < 3; ++axis) {
    mesh_cells[1].grav.parts[0].a_grav_mesh[axis] = 0.f;
  }
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  mesh_cells[1].grav.parts[0].potential_mesh = 0.f;
#endif
  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);
  const double expected_fractional_ax = -3.946721553802e-2;
  const double measured_fractional_ax =
      mesh_cells[1].grav.parts[0].a_grav_mesh[0];
  if (fabs(measured_fractional_ax - expected_fractional_ax) > 1e-7) {
    error("Off-grid zoom mesh force changed: measured=%e expected=%e",
          measured_fractional_ax, expected_fractional_ax);
  }

  /* Delta-f deposition is linear in the mesh-only particle weight. */
  double neutrino_mass = 1.;
  double neutrino_degeneracy = 1.;
  cosmology.a = 1.;
  cosmology.N_nu = 1;
  cosmology.M_nu_eV = &neutrino_mass;
  cosmology.deg_nu = &neutrino_degeneracy;
  cosmology.T_nu_0_eV = 1.;
  phys_const.const_speed_light_c = 1.;
  e.neutrino_mass_conversion_factor = 1.;
  neutrino_props.use_delta_f_mesh_only = 1;
  mesh_cells[0].grav.parts[0].type = swift_type_neutrino;
  mesh_cells[0].grav.parts[0].id_or_neg_offset = 1;
  struct neutrino_model nu_model;
  gather_neutrino_consts(&s, &nu_model);
  double delta_f_weight = 1.;
  gpart_neutrino_weight_mesh_only(&mesh_cells[0].grav.parts[0], &nu_model,
                                  &delta_f_weight);
  for (int axis = 0; axis < 3; ++axis) {
    mesh_cells[1].grav.parts[0].a_grav_mesh[axis] = 0.f;
  }
#ifndef SWIFT_GRAVITY_NO_POTENTIAL
  mesh_cells[1].grav.parts[0].potential_mesh = 0.f;
#endif
  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);
  if (fabs(mesh_cells[1].grav.parts[0].a_grav_mesh[0] -
           delta_f_weight * measured_fractional_ax) > 1e-7) {
    error("Zoom mesh did not apply delta-f deposition weighting");
  }
  neutrino_props.use_delta_f_mesh_only = 0;
  mesh_cells[0].grav.parts[0].type = swift_type_dark_matter;

  /* A zoom regrid changes the mesh origin and spacing and invalidates the
   * cached Fourier kernel before the next force calculation. */
  struct cell bkg_cell = {0};
  for (int axis = 0; axis < 3; ++axis) {
    bkg_cell.width[axis] = 10.;
    zoom_props.dim[axis] = 12.;
    zoom_props.region_lower_bounds[axis] = 21.;
  }
  zoom_props.bkg_cells_top = &bkg_cell;
  props.a_smooth = zoom_mesh.r_s * zoom_mesh.cell_fac;
  props.r_s = mesh.r_s;
  props.r_cut_max_ratio = mesh.r_cut_max / mesh.r_s;
  zoom_mesh.buffer_cells = 1;
  if (zoom_mesh.kernel_fft == NULL) {
    error("Zoom mesh kernel should have been cached before regrid");
  }
  zoom_mesh_regrid(&zoom_mesh, &s);
  if (zoom_mesh.loc[0] != 11. || zoom_mesh.dim[0] != 32. ||
      fabs(zoom_mesh.r_s - props.a_smooth * 2.) > 1e-12 ||
      zoom_mesh.kernel_fft != NULL) {
    error("Zoom mesh regrid did not update geometry and discard cached kernel");
  }
  for (int axis = 0; axis < 3; ++axis) {
    mesh_cells[1].grav.parts[0].a_grav_mesh[axis] = 0.f;
  }
  zoom_mesh_compute_potential(&zoom_mesh, &s, &e.threadpool, /*verbose=*/0);
  if (zoom_mesh.kernel_fft == NULL ||
      !isfinite(mesh_cells[1].grav.parts[0].a_grav_mesh[0]) ||
      mesh_cells[1].grav.parts[0].a_grav_mesh[0] == 0.f) {
    error("Regridded zoom mesh failed to reconstruct its kernel and force");
  }

  zoom_mesh_clean(&zoom_mesh);
  threadpool_clean(&e.threadpool);

  free_cell(&mesh_cells[0]);
  free_cell(&mesh_cells[1]);
#endif

  free_cell(&ci);
  free_cell(&cj);

  return 0;
}
