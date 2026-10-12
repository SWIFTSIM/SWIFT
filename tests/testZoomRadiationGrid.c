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

/* Checks that the radiation grid helpers, the radiation cap and the fixed
 * radius floor use the zoom grid in a zoom run. The zoom geometry is built
 * by the real space_regrid() on a mock particle set, as in testZoomRegrid. */

/* Config parameters. */
#include <config.h>

/* Standard headers. */
#include <math.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"
#include "zoom_region/zoom.h"

/**
 * @brief Build a periodic zoom space and its zoom geometry.
 *
 * Background particles sit on a lattice through the box, high resolution
 * particles on a lattice in a cube of side @p zoom_extent at the centre.
 *
 * @param params The zoom parameters.
 * @param zoom_extent The side of the high resolution particle cube.
 * @param floor_width The fixed radius floor on the zoom cell width.
 */
static struct space *make_zoom_space(struct swift_params *params,
                                     const double zoom_extent,
                                     const double floor_width) {

  struct space *s = malloc(sizeof(struct space));
  bzero(s, sizeof(struct space));
  s->dim[0] = s->dim[1] = s->dim[2] = 1000.;
  s->periodic = 1;
  s->with_zoom_region = 1;

  const int n_bkg = 4, n_zoom = 4;
  s->nr_gparts = n_bkg * n_bkg * n_bkg + n_zoom * n_zoom * n_zoom;
  if (posix_memalign((void **)&s->gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0)
    error("Failed to allocate the gparts.");
  bzero(s->gparts, s->nr_gparts * sizeof(struct gpart));

  size_t k = 0;
  for (int i = 0; i < n_bkg; i++)
    for (int j = 0; j < n_bkg; j++)
      for (int l = 0; l < n_bkg; l++, k++) {
        s->gparts[k].x[0] = (i + 0.1) * s->dim[0] / n_bkg;
        s->gparts[k].x[1] = (j + 0.1) * s->dim[1] / n_bkg;
        s->gparts[k].x[2] = (l + 0.1) * s->dim[2] / n_bkg;
        s->gparts[k].type = swift_type_dark_matter_background;
        s->gparts[k].mass = 1.;
      }
  for (int i = 0; i < n_zoom; i++)
    for (int j = 0; j < n_zoom; j++)
      for (int l = 0; l < n_zoom; l++, k++) {
        const double step = zoom_extent / (n_zoom - 1);
        s->gparts[k].x[0] = 0.5 * (s->dim[0] - zoom_extent) + i * step;
        s->gparts[k].x[1] = 0.5 * (s->dim[1] - zoom_extent) + j * step;
        s->gparts[k].x[2] = 0.5 * (s->dim[2] - zoom_extent) + l * step;
        s->gparts[k].type = swift_type_dark_matter;
        s->gparts[k].mass = 1.;
      }

  zoom_props_init(params, s, /*verbose=*/0);
  s->zoom_props->cell_min_fixed_radius = floor_width;
  space_regrid(s, /*verbose=*/0);

  return s;
}

/**
 * @brief Run make_zoom_space() in a child process.
 *
 * @return 1 if the child stopped with an error, 0 if it returned normally.
 */
static int zoom_space_fails(struct swift_params *params,
                            const double zoom_extent,
                            const double floor_width) {

  fflush(stdout);
  fflush(stderr);
  const pid_t pid = fork();
  if (pid < 0) error("fork() failed.");
  if (pid == 0) {
    make_zoom_space(params, zoom_extent, floor_width);
    _exit(0);
  }
  int status = 0;
  if (waitpid(pid, &status, 0) != pid) error("waitpid() failed.");
  return !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
}

static void make_params(struct swift_params *params, const char *depth,
                        const char *bkg_cdim) {
  parser_init("testZoomRadiationGrid-in-memory", params);
  parser_set_param(params, "Gravity:mesh_side_length:64");
  parser_set_param(params, depth);
  if (bkg_cdim != NULL) parser_set_param(params, bkg_cdim);
}

int main(int argc, char *argv[]) {

  clocks_set_cpufreq(0);

  /* Zoom region inside one background cell (automatic background grid). */
  struct swift_params params;
  make_params(&params, "ZoomRegion:zoom_top_level_depth:2", NULL);
  struct space *s = make_zoom_space(&params, 40., 0.);
  const struct zoom_region_properties *zp = s->zoom_props;
  message("zoom cdim=%d depth=%d zoom dim=%g bkg cdim=%d", zp->cdim[0],
          zp->zoom_cell_depth, zp->dim[0], s->cdim[0]);

  if (zp->cdim[0] != 4) error("Expected 4 zoom cells per axis.");

  /* The zoom grid decides coverage, whatever the background grid is. */
  if (space_radiation_top_stencil_covers_box(s))
    error("A 4x4x4 zoom grid without wrapping cannot cover the zoom region.");
  const int bkg_cdim = s->cdim[0];
  s->cdim[0] = s->cdim[1] = s->cdim[2] = 3;
  if (space_radiation_top_stencil_covers_box(s))
    error("Coverage was taken from the background grid.");
  s->cdim[0] = s->cdim[1] = s->cdim[2] = bkg_cdim;
  if (!space_radiation_top_stencil_can_cover_box(s))
    error("A zoom region in one background cell can reach 2 cells/axis.");

  /* The cap on the star radiation term is set by the zoom region, and an
   * unbounded HII reach asks for exactly 2 zoom cells per axis. */
  const double cap_width = space_radiation_cell_max_width(s);
  if (cap_width != 0.5 * zp->dim[0])
    error("Radiation cap width %g is not half the zoom region %g.", cap_width,
          zp->dim[0]);
  const float term =
      space_regrid_star_radiation_term_for(1e30f, 0.f, cap_width);
  const double radius = space_regrid_search_radius_for(0.f, term);
  const int capped_cdim = (int)floor(zp->dim[0] / radius);
  if (capped_cdim != 2)
    error("An unbounded HII reach asks for %d zoom cells per axis, not 2.",
          capped_cdim);

  /* A small HII reach is not capped. */
  const float small_h = 0.01f * (float)zp->width[0];
  if (space_regrid_star_radiation_term_for(small_h, 0.f, cap_width) != small_h)
    error("A small HII reach was capped.");
  message("Zoom radiation grid helpers passed (cap width %g).", cap_width);

  /* Zoom region across several background cells: no 2 cells/axis grid. */
  struct swift_params params_wide;
  make_params(&params_wide, "ZoomRegion:zoom_top_level_depth:1",
              "ZoomRegion:bkg_top_level_cells:64");
  struct space *s_wide = make_zoom_space(&params_wide, 40., 0.);
  message("wide zoom cdim=%d depth=%d", s_wide->zoom_props->cdim[0],
          s_wide->zoom_props->zoom_cell_depth);
  if (space_radiation_top_stencil_can_cover_box(s_wide))
    error("A zoom region over several background cells cannot reach 2 "
          "cells/axis.");
  message("Wide zoom region correctly flagged.");

  /* The fixed radius floor: a floor at the zoom cell width is accepted, a
   * floor above it stops the run. */
  const double zoom_width = zp->width[0];
  if (zoom_space_fails(&params, 40., zoom_width))
    error("A floor equal to the zoom cell width was rejected.");
  if (!zoom_space_fails(&params, 40., 1.01 * zoom_width))
    error("A floor above the zoom cell width was accepted.");
  message("Zoom cell width floor check passed (width %g).", zoom_width);

#ifdef SINKS_WITH_FIXED_CUTOFF_RADIUS
  /* space_init() hands the sink aperture floor to the zoom grid. */
  {
    struct hydro_props hp;
    hydro_props_init_no_hydro(&hp);
    struct cosmology cosmo;
    cosmology_init_no_cosmo(&cosmo);
    double dim[3] = {1000., 1000., 1000.};

    struct sink_props sink_properties;
    memset(&sink_properties, 0, sizeof(struct sink_props));
    sink_properties.use_fixed_r_cut = 1;
    sink_properties.cut_off_radius = 3.f;

    struct space s_init;
    space_init(&s_init, &params, &cosmo, dim, &hp, &sink_properties,
               /*parts=*/NULL, /*gparts=*/NULL, /*sinks=*/NULL,
               /*sparts=*/NULL, /*bparts=*/NULL, /*Npart=*/0, /*Ngpart=*/0,
               /*Nsink=*/0, /*Nspart=*/0, /*Nbpart=*/0, /*Nnupart=*/0,
               /*periodic=*/1, /*replicate=*/1, /*remap_ids=*/0,
               /*generate_gas_in_ics=*/0, /*hydro=*/1, /*self_gravity=*/1,
               /*star_formation=*/0, /*with_sink=*/0, /*with_DM=*/1,
               /*with_DM_background=*/1, /*neutrinos=*/0,
               /*with_zoom_region=*/1, /*verbose=*/0, /*dry_run=*/1,
               /*nr_nodes=*/1);
    const double expected = space_stretch * sink_properties.cut_off_radius;
    if (fabs(s_init.zoom_props->cell_min_fixed_radius - expected) > 1e-6)
      error("Zoom floor %g does not match the sink aperture floor %g.",
            s_init.zoom_props->cell_min_fixed_radius, expected);
    message("space_init() zoom aperture floor passed (%g).", expected);
  }
#else
  message("space_init() zoom aperture floor SKIPPED (no fixed sink radius).");
#endif

  return 0;
}
