/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2024 Will J. Roper (w.roper@sussex.ac.uk).
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

/* Standard headers. */
#include <fenv.h>

/* Local headers. */
#include "swift.h"
#include "zoom_region/zoom.h"

double generate_gaussian_coordinate(const double mean, const double std,
                                    const double max_width) {
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;
  double z0 = (sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2)) * std + mean;

  /* Try again if we're out of bounds. */
  if (z0 < mean - max_width / 2 || z0 > mean + max_width / 2) {
    return generate_gaussian_coordinate(mean, std, max_width);
  }

  return z0;
}
void make_mock_cells(struct space *s) {

  /* Allocate memory for the cells. */
  s->cells_top = (struct cell *)malloc(s->nr_cells * sizeof(struct cell));

  /* Set up the top level cells. */
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        struct cell *c = &s->cells_top[i * 16 + j * 4 + k];
        c->width[0] = s->zoom_props->dim[0] / 4;
        c->width[1] = s->zoom_props->dim[1] / 4;
        c->width[2] = s->zoom_props->dim[2] / 4;
        c->loc[0] = s->zoom_props->region_lower_bounds[0] + i * c->width[0];
        c->loc[1] = s->zoom_props->region_lower_bounds[1] + j * c->width[1];
        c->loc[2] = s->zoom_props->region_lower_bounds[2] + k * c->width[2];
        c->type = cell_type_zoom;
        c->subtype = cell_subtype_regular;
        c->split = 0;
        c->void_parent = NULL;
        c->grav.count = 0;
        c->grav.parts = NULL;
        c->top = c;
      }
    }
  }

  /* Set up the background cell. */
  struct cell *c = &s->cells_top[64];
  c->width[0] = 100;
  c->width[1] = 100;
  c->width[2] = 100;
  c->loc[0] = 450;
  c->loc[1] = 450;
  c->loc[2] = 450;
  c->type = cell_type_bkg;
  c->subtype = cell_subtype_void;
  c->split = 0;
  c->void_parent = NULL;
  c->grav.count = 0;
  c->grav.parts = NULL;
  c->top = c;

  /* Flag the void cell indices. */
  s->zoom_props->void_cell_indices = (int *)malloc(64 * sizeof(int));
  s->zoom_props->void_cell_indices[0] = 64;
  s->zoom_props->nr_void_cells = 1;
}

void make_mock_space(struct space *s) {

  /* Define the boxsize. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;

  /* The simulation is periodic */
  s->periodic = 1;

  /* Define the gpart count (we only care about the zoom region) */
  s->nr_gparts = 10000;

  /* Allocate memory for the gparts. */
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0) {
    error("Failed to allocate memory for gparts");
  }
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Define the width of the zoom region (randomly). */
  double zoom_width = 50 * ((double)rand() / RAND_MAX) + 1;
  message("Zoom width = %f", zoom_width);

  /* Define the zoom particles by sampling from a normal distribution. */
  for (int i = 10000; i < 20000; i++) {
    gparts[i].x[0] =
        generate_gaussian_coordinate(s->dim[0] / 2, zoom_width, 100);
    gparts[i].x[1] =
        generate_gaussian_coordinate(s->dim[1] / 2, zoom_width, 100);
    gparts[i].x[2] =
        generate_gaussian_coordinate(s->dim[2] / 2, zoom_width, 100);
    gparts[i].type = swift_type_dark_matter;
    gparts[i].mass = 1.0;
  }

  s->gparts = gparts;

  /* Set up the zoom region properties. */
  s->zoom_props = (struct zoom_region_properties *)malloc(
      sizeof(struct zoom_region_properties));
  struct zoom_region_properties *zoom_props = s->zoom_props;
  zoom_props->region_lower_bounds[0] = s->dim[0] / 2 - 50;
  zoom_props->region_lower_bounds[1] = s->dim[1] / 2 - 50;
  zoom_props->region_lower_bounds[2] = s->dim[2] / 2 - 50;
  zoom_props->region_upper_bounds[0] = s->dim[0] / 2 + 50;
  zoom_props->region_upper_bounds[1] = s->dim[1] / 2 + 50;
  zoom_props->region_upper_bounds[2] = s->dim[2] / 2 + 50;
  zoom_props->dim[0] = 100;
  zoom_props->dim[1] = 100;
  zoom_props->dim[2] = 100;
  zoom_props->cdim[0] = 4;
  zoom_props->cdim[1] = 4;
  zoom_props->cdim[2] = 4;
  zoom_props->nr_zoom_cells =
      zoom_props->cdim[0] * zoom_props->cdim[1] * zoom_props->cdim[2];
  zoom_props->zoom_cell_depth = 2;
  zoom_props->nr_bkg_cells = 1;
  zoom_props->nr_buffer_cells = 0;
  zoom_props->with_buffer_cells = 0;
  zoom_props->bkg_cell_offset = zoom_props->nr_zoom_cells;

  /* We have nr_zoom_cells zoom cells and 1 background cell. */
  s->nr_cells = zoom_props->nr_zoom_cells + 1;

  /* Allocate the cells. */
  make_mock_cells(s);

  /* Allocate sub cells and multipoles. */
  s->cells_sub = (struct cell **)calloc(2, sizeof(struct cell *));
  s->multipoles_sub =
      (struct gravity_tensors **)calloc(2, sizeof(struct gravity_tensors *));
}

void test_cell_tree(struct cell *c, struct space *s) {

  /* Recurse if the the cell is split. */
  if (c->split) {
    for (int i = 0; i < 8; i++) {
      test_cell_tree(c->progeny[i], s);
    }
  }

  /* Otherwise we are in a void leaf which should have zoom cell progeny. */
  else {

    /* Check this void leaf is attached to zoom cells and the zoom cells are
     * correctly attached. */
    for (int i = 0; i < 8; i++) {
      assert(c->progeny[i]->void_parent != NULL);
      assert(c->progeny[i]->void_parent->subtype == cell_subtype_void);
      assert(c->progeny[i]->type == cell_type_zoom ||
             c->progeny[i]->type == cell_type_buffer);
    }

    /* NOTE: zoom_void_space_split contains a lot of its own checks so this
     * is sufficient here. */
  }
}

int main(int argc, char *argv[]) {

  /* NOTE: This is a modified version of testZoomRegrid. It uses the exact same
   * construction mechanism and then constructs the void tree and tests it.
   * This does mean there are redundant background particles  */

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(argv[1], &param_file);

  /* Create a space structure. (this creates fake particles to get the layout
   * of cells right but we will modify this shortly to get a particle per
   * cell) */
  struct space *s = malloc(sizeof(struct space));
  bzero(s, sizeof(struct space));
  make_mock_space(s);

  /* Construct the void cell tree. */
  zoom_void_space_split(s, /*verbose*/ 1);

  /* Test the cell tree. */
  for (int cid = 0; cid < s->nr_cells; cid++) {
    /* Only test void cells. */
    if (s->cells_top[cid].subtype == cell_subtype_void)
      test_cell_tree(&s->cells_top[cid], s);
  }

  /* Free the space. */
  free(s->local_cells_top);
  free(s->multipoles_top);
  free(s->local_cells_with_tasks_top);
  free(s->cells_with_particles_top);
  free(s->local_cells_with_particles_top);
  free(s->zoom_props->local_zoom_cells_top);
  free(s->zoom_props->void_cell_indices);
  free(s->zoom_props->neighbour_cells_top);
  free(s->zoom_props->local_bkg_cells_top);
  free(s->zoom_props->local_buffer_cells_top);
  free(s->zoom_props->local_zoom_cells_with_particles_top);
  free(s->zoom_props->local_bkg_cells_with_particles_top);
  free(s->zoom_props->local_buffer_cells_with_particles_top);
  free(s->cells_top);
  free(s->gparts);
  free(s->zoom_props);
  free(s->cells_sub);
  free(s->multipoles_sub);
  free(s);

  return 0;
}
