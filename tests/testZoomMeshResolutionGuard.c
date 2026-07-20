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
 * The background mesh's smoothing scale (@p props->r_s) is fixed by the
 * caller; @p zoom_mesh_side_length is varied to make the zoom mesh either
 * higher or lower resolution than that background mesh.
 *
 * @param params The #swift_params to fill in.
 * @param props The #gravity_props to fill in (background mesh already
 * "initialised").
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

  bzero(params, sizeof(struct swift_params));
  bzero(props, sizeof(struct gravity_props));
  bzero(s, sizeof(struct space));
  bzero(zoom_props, sizeof(struct zoom_region_properties));
  bzero(bkg_cell, sizeof(struct cell));
  bzero(zoom_cell, sizeof(struct cell));

  parser_init("test_zoom_mesh_resolution_guard", params);
  parser_set_param(params, "Gravity:zoom_mesh:1");
  {
    char desc[128];
    snprintf(desc, sizeof(desc), "Gravity:zoom_mesh_side_length:%d",
             zoom_mesh_side_length);
    parser_set_param(params, desc);
  }
  parser_set_param(params, "Gravity:zoom_r_cut_max:4.5");

  /* A pre-"initialised" background (global) periodic mesh: a_smooth=1.25,
   * box=100, mesh_side_length=128 -> r_s = 1.25 * 100 / 128 = 0.9765625. */
  props->a_smooth = 1.25;
  props->r_s = 0.9765625;
  props->r_s_inv = 1. / props->r_s;

  s->periodic = 1;
  s->with_zoom_region = 1;
  s->dim[0] = 100.;
  s->dim[1] = 100.;
  s->dim[2] = 100.;

  /* A background top-level grid of 8 cells per axis over the 100-wide box. */
  for (int i = 0; i < 3; ++i) bkg_cell->width[i] = 100. / 8.;

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
 * @brief Run #zoom_mesh_init in a child process and report whether it
 * survived.
 *
 * #zoom_mesh_init calls #error() (i.e. exit(1)) when it rejects the
 * configuration, so we fork to observe that without killing the test runner.
 *
 * @param zoom_mesh_side_length The value of Gravity:zoom_mesh_side_length.
 * @return 1 if #zoom_mesh_init completed normally, 0 if it errored.
 */
static int zoom_mesh_init_survives(const int zoom_mesh_side_length) {

  struct swift_params params;
  struct gravity_props props;
  struct space s;
  struct zoom_region_properties zoom_props;
  struct cell bkg_cell, zoom_cell;
  setup_mock_engine(&params, &props, &s, &zoom_props, &bkg_cell, &zoom_cell,
                    zoom_mesh_side_length);

  const pid_t pid = fork();
  if (pid < 0) error("fork() failed");

  if (pid == 0) {
    /* Child: this either returns (success) or calls error()/exit(1)
     * (rejected). The diagnostic message printed to stderr on rejection is
     * expected and not a test failure. */
    struct zoom_pm_mesh zoom_mesh;
    bzero(&zoom_mesh, sizeof(struct zoom_pm_mesh));
    zoom_mesh_init(&zoom_mesh, &params, &props, &s, /*verbose=*/0);
    zoom_mesh_clean(&zoom_mesh);
    _exit(0);
  }

  int status;
  if (waitpid(pid, &status, 0) < 0) error("waitpid() failed");

  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

int main(int argc, char *argv[]) {

  (void)argc;
  (void)argv;

  /* A zoom mesh with side_length=128 over the ~35-wide padded zoom domain
   * (10 + 2 * 1 buffer cell of width 12.5) is far higher resolution than the
   * r_s=0.9765625 background mesh: r_s_zoom = 1.25 * 35/128 ~ 0.342. This
   * must be accepted. */
  if (!zoom_mesh_init_survives(128))
    error("A genuinely higher-resolution zoom mesh should not be rejected.");

  /* A zoom mesh with the minimum allowed side_length=8 over the same ~35-wide
   * domain is far *lower* resolution than the background mesh: r_s_zoom =
   * 1.25 * 35/8 ~ 5.47 >> r_s_background = 0.9765625. This must be rejected;
   * accepting it silently makes the mesh-tree splitting criterion *more*
   * conservative inside the zoom region, exploding the gravity task count
   * instead of reducing it. */
  if (zoom_mesh_init_survives(8))
    error(
        "A zoom mesh lower resolution than the background mesh should have "
        "been rejected by zoom_mesh_init.");

  message("Zoom mesh resolution guard behaves as expected.");
  return 0;
}
