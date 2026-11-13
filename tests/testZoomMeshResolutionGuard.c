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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/**
 * @brief Build the minimal engine-adjacent structures #zoom_mesh_init needs.
 *
 * The smoothing scale of the global mesh is fixed. We vary the side length of
 * the zoom mesh to make it higher or lower resolution than the global mesh.
 *
 * @param params The #swift_params to fill in.
 * @param props The #gravity_props to fill in.
 * @param s The #space to fill in.
 * @param zoom_props The #zoom_region_properties to fill in.
 * @param bkg_cell A mock background top-level #cell.
 * @param zoom_cell A mock zoom top-level #cell.
 * @param zoom_mesh_side_length The value of Gravity:zoom_mesh_side_length.
 */
static void setup_mock_engine(struct swift_params *params,
                              struct gravity_props *props, struct space *s,
                              struct zoom_region_properties *zoom_props,
                              struct cell *bkg_cell, struct cell *zoom_cell,
                              const int zoom_mesh_side_length) {

  /* Start from nothing */
  bzero(params, sizeof(struct swift_params));
  bzero(props, sizeof(struct gravity_props));
  bzero(s, sizeof(struct space));
  bzero(zoom_props, sizeof(struct zoom_region_properties));
  bzero(bkg_cell, sizeof(struct cell));
  bzero(zoom_cell, sizeof(struct cell));

  /* The zoom mesh parameters */
  parser_init("test_zoom_mesh_resolution_guard", params);
  parser_set_param(params, "Gravity:zoom_mesh:1");
  {
    char desc[128];
    snprintf(desc, sizeof(desc), "Gravity:zoom_mesh_side_length:%d",
             zoom_mesh_side_length);
    parser_set_param(params, desc);
  }
  parser_set_param(params, "Gravity:zoom_r_cut_max:4.5");

  /* The global mesh: 128 cells over a box of 100 with a_smooth = 1.25, so
   * r_s = 1.25 * 100 / 128 = 0.9765625. */
  props->a_smooth = 1.25;
  props->r_s = 0.9765625;
  props->r_s_inv = 1. / props->r_s;
  props->r_cut_max_ratio = 4.5;

  /* A periodic box with a zoom region */
  s->periodic = 1;
  s->with_zoom_region = 1;
  s->dim[0] = 100.;
  s->dim[1] = 100.;
  s->dim[2] = 100.;

  /* A background top-level grid of 8 cells per axis over the 100-wide box. */
  for (int i = 0; i < 3; ++i) {
    bkg_cell->width[i] = 100. / 8.;
  }

  /* A zoom region 10 units across, centred in the box. */
  for (int i = 0; i < 3; ++i) {
    zoom_props->dim[i] = 10.;
    zoom_props->region_lower_bounds[i] = 45.;
    zoom_cell->width[i] = zoom_props->dim[i] / 4.;
  }
  zoom_props->nr_bkg_cells = 1;
  zoom_props->bkg_cells_top = bkg_cell;
  zoom_props->nr_zoom_cells = 1;
  zoom_props->zoom_cells_top = zoom_cell;

  s->zoom_props = zoom_props;
}

/**
 * @brief Check mesh setup and optional covered-particle force validation.
 *
 * These paths call #error() when they reject the configuration, so we fork to
 * observe that without killing the test runner.
 *
 * @param zoom_mesh_side_length The value of Gravity:zoom_mesh_side_length.
 * @param softening The current internal DM gravitational softening radius.
 * @param regrid_bkg_width New background-cell width, or 0 for no regrid.
 * @param check_covered Whether to test the force guard with a covered DM
 * particle.
 * @return 1 if setup completed normally, 0 if it errored.
 */
static int zoom_mesh_init_survives(const int zoom_mesh_side_length,
                                   const float softening,
                                   const double regrid_bkg_width,
                                   const int check_covered) {

  /* Set up the mock engine */
  struct swift_params params;
  struct gravity_props props;
  struct space s;
  struct zoom_region_properties zoom_props;
  struct cell bkg_cell, zoom_cell;
  setup_mock_engine(&params, &props, &s, &zoom_props, &bkg_cell, &zoom_cell,
                    zoom_mesh_side_length);
  props.epsilon_DM_cur = softening;
  props.epsilon_baryon_cur = 0.6f;
  struct engine e = {0};
  e.gravity_properties = &props;
  s.e = &e;

  const pid_t pid = fork();
  if (pid < 0) error("fork() failed");

  if (pid == 0) {
    /* In the child we either get through (accepted) or hit an error
     * (rejected). The error message printed in that case is expected. */
    struct zoom_pm_mesh zoom_mesh;
    bzero(&zoom_mesh, sizeof(struct zoom_pm_mesh));
    zoom_mesh_init(&zoom_mesh, &params, &props, &s, /*nr_threads=*/1,
                   /*verbose=*/0);

    /* Do we want to regrid with different background cells? */
    if (regrid_bkg_width > 0.) {
      for (int i = 0; i < 3; ++i) {
        bkg_cell.width[i] = regrid_bkg_width;
      }
      zoom_mesh_regrid(&zoom_mesh, &s);
    }

    /* Do we want to compute forces for a particle covered by the mesh? */
    if (check_covered) {
      struct gpart gp = {0};
      gp.type = swift_type_dark_matter;
      gp.time_bin = 1;
#ifdef MULTI_SOFTENING_GRAVITY
      gp.epsilon = softening;
#endif
      for (int i = 0; i < 3; ++i) {
        zoom_cell.loc[i] = 48.;
        gp.x[i] = 49.;
      }
      zoom_cell.top = &zoom_cell;
      zoom_cell.grav.parts = &gp;
      zoom_cell.grav.count = 1;
      struct gravity_tensors multipole = {0};
      multipole.m_pole.max_softening = softening;
      zoom_cell.grav.multipole = &multipole;
      int cell_index = 0;
      s.cells_top = &zoom_cell;
      zoom_props.local_zoom_cells_top = &cell_index;
      zoom_props.nr_local_zoom_cells = 1;
      if (!zoom_mesh_cell_is_covered(&zoom_mesh, &zoom_cell, 0)) {
        error("Mock DM particle must be covered by the zoom mesh");
      }
#ifdef WITH_MPI
      MPI_Init(NULL, NULL);
#endif
      zoom_mesh_compute_potential(&zoom_mesh, &s, NULL, 0);
      _exit(0);
    }
    zoom_mesh_clean(&zoom_mesh);
    _exit(0);
  }

  /* Back in the parent, how did the child get on? */
  int status;
  if (waitpid(pid, &status, 0) < 0) error("waitpid() failed");
  if (!WIFEXITED(status)) {
    error("Zoom mesh guard test child terminated with signal %d",
          WTERMSIG(status));
  }
  if (check_covered && WEXITSTATUS(status) != 1) {
    error("Covered particle guard exited with status %d instead of 1",
          WEXITSTATUS(status));
  }

  return WEXITSTATUS(status) == 0;
}

int main(int argc, char *argv[]) {

  (void)argc;
  (void)argv;

  /* 128 cells over the zoom mesh (35 across: the region of 10 plus a buffer
   * cell of 12.5 either side) gives r_s = 1.25 * 35 / 128 ~ 0.342, well below
   * the global mesh. This is fine, even with the unused baryon softening of
   * 0.6. */
  if (!zoom_mesh_init_survives(128, 0.3f, 0., 0)) {
    error("A genuinely higher-resolution zoom mesh should not be rejected.");
  }

  /* 8 cells over the same mesh gives r_s = 1.25 * 35 / 8 ~ 5.47, which is
   * lower resolution than the global mesh. We can't have that. */
  if (zoom_mesh_init_survives(8, 0.3f, 0., 0)) {
    error(
        "A zoom mesh lower resolution than the background mesh should have "
        "been rejected by zoom_mesh_init.");
  }

  /* A covered DM particle with a softening above r_s can't be allowed. */
  if (zoom_mesh_init_survives(128, 0.4f, 0., 1)) {
    error("Covered DM softening larger than the smoothing scale was accepted.");
  }

  /* Regridding the mesh from 35 to 24 across takes r_s from ~0.342 to ~0.234,
   * which is now below the softening of 0.3. */
  if (zoom_mesh_init_survives(128, 0.3f, 7., 1)) {
    error("Covered DM softening above the regridded scale was accepted.");
  }

  message("Zoom mesh resolution guard behaves as expected.");
  return 0;
}
