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

/**
 * @brief Generate a random coordinate from a gaussian distribution.
 *
 * @param mean The mean of the gaussian distribution.
 * @param std The standard deviation of the gaussian distribution.
 * @param max_width The maximum width of the cell.
 * @param id The ID of the particle for which to generate the number.
 * @param ti_current The time (on the time-line) for which to generate the
 * number.
 * @return A random number drawn from the gaussian distribution.
 */
double random_gaussian_coordinate(const double mean, const double std,
                                  const double max_width, const int id,
                                  const integertime_t ti_current,
                                  const enum random_number_type type) {

  /* Generate a random number from a normal distribution. */
  double z0 = random_gaussian(mean, std, id, ti_current, type);

  /* We only want to go out at most by the size of a cell. If we've got a
   * coordinate out too far we should try again changing the random number
   * seed. */
  if (z0 < mean - max_width / 2 || z0 > mean + max_width / 2) {
    return random_gaussian_coordinate(
        mean, std, max_width, id, ti_current,
        type + 3 % random_number_powerspectrum_split);
  }

  return z0;
}

void make_mock_space(struct space *s) {

  /* Define the boxsize. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;

  /* The simulation is periodic */
  s->periodic = 1;

  /* Define the gpart count (100 high and 100 low resolution) */
  s->nr_gparts = 100 + 100;

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Allocate memory for the gparts. */
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, gpart_align,
                     s->nr_gparts * sizeof(struct gpart)) != 0) {
    error("Failed to allocate memory for gparts");
  }
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Randomly place the background particles. */
  for (int i = 0; i < 100; i++) {
    gparts[i].x[0] = s->dim[0] * 0.99 * ((double)rand() / RAND_MAX) + 1;
    gparts[i].x[1] = s->dim[1] * 0.99 * ((double)rand() / RAND_MAX) + 1;
    gparts[i].x[2] = s->dim[2] * 0.99 * ((double)rand() / RAND_MAX) + 1;
    gparts[i].type = swift_type_dark_matter_background;
    gparts[i].mass = 1.0;
  }

  /* Define the width of the zoom region (randomly). */
  double zoom_width = 50;

  /* Get the "current time" */
  time_t ti_current = time(NULL);

  /* Define the zoom particles by sampling from a normal distribution. */
  for (int i = 100; i < 200; i++) {
    gparts[i].x[0] = random_gaussian_coordinate(s->dim[0] / 2, zoom_width, 100,
                                                i, ti_current, 0);
    gparts[i].x[1] = random_gaussian_coordinate(s->dim[1] / 2, zoom_width, 100,
                                                i, ti_current, 1);
    gparts[i].x[2] = random_gaussian_coordinate(s->dim[2] / 2, zoom_width, 100,
                                                i, ti_current, 2);
    gparts[i].type = swift_type_dark_matter;
    gparts[i].mass = 1.0;
  }

  s->gparts = gparts;

  /* Allocate sub cells and multipoles. */
  s->cells_sub = (struct cell **)calloc(2, sizeof(struct cell *));
  s->multipoles_sub =
      (struct gravity_tensors **)calloc(2, sizeof(struct gravity_tensors *));
}

void associate_gparts_to_cells(struct space *s) {
  for (size_t i = 0; i < s->nr_gparts; i++) {
    struct gpart *gpart = &s->gparts[i];

    int cid = cell_getid_from_pos(s, gpart->x[0], gpart->x[1], gpart->x[2]);

    struct cell *c = &s->cells_top[cid];
    if (c == NULL) {
      error("Failed to find cell for gpart.");
    }
    c->grav.count++;
  }
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
      assert(c->progeny[i]->type == cell_type_zoom);
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

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, s, /*verbose*/ 0);

  /* Run the regridding. */
  space_regrid(s, /*verbose*/ 0);

  /* Associate gparts. */
  associate_gparts_to_cells(s);

  /* Construct the void cell tree. */
  zoom_void_space_split(s, /*verbose*/ 0);

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
