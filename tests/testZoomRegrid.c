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

/* This script tests the process of regridding a zoom simulation. It'll first
 * calculate the zoom geometry, construct the cell structures as would be done
 * in a run during a regrid, populates the cells and then tests all particles
 * have ended up where they are supposed to be. */

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
  double zoom_width = 50 * ((double)rand() / RAND_MAX) + 1;
  message("Zoom width = %f", zoom_width);

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

int main(int argc, char *argv[]) {

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

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, s, /*verbose*/ 0);

  /* Run the regridding. */
  space_regrid(s, /*verbose*/ 0);

  /* Associate gparts. */
  associate_gparts_to_cells(s);

  /* Ensure void cells have 0 counts. */
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->subtype == cell_subtype_void) {
      if (c->grav.count != 0) {
        error(
            "Cell %d has %d particles (c->type = %s, c->subtype = %s, c->loc = "
            "[%f, %f, %f])",
            i, c->grav.count, cellID_names[c->type],
            subcellID_names[c->subtype], c->loc[0], c->loc[1], c->loc[2]);
      }
    }
  }

  /* Test each cell type contains the expected counts. (This ensures the zoom
   * and void cells make sense too.)*/
  int zoom_count = 0;
  int bkg_count = 0;
  int buffer_count = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    struct cell *c = &s->cells_top[i];
    if (c->type == cell_type_zoom) {
      zoom_count += c->grav.count;
    } else if (c->type == cell_type_bkg) {
      bkg_count += c->grav.count;
    } else if (c->type == cell_type_buffer) {
      buffer_count += c->grav.count;
    } else {
      error("Cell %d has unexpected type %s", i, cellID_names[c->type]);
    }
  }

  /* Because we are doing things randomly we need to check for both the buffer
   * cell case and no buffer cell case. */
  if (s->zoom_props->with_buffer_cells) {
    /* Zoom region has 100 high res particles and any low res interlopers. */
    assert(zoom_count == 100 + 100 - bkg_count - buffer_count);

    /* Buffer region + background region + zoom interlopers == 100 background
     * particles. */
    assert(bkg_count + buffer_count + zoom_count - 100 == 100);

  } else {
    /* Zoom region has 100 high res particles and any low res interlopers. */
    assert(zoom_count == 100 + 100 - bkg_count);

    /* Background region + zoom interlopers == 100 background particles. */
    assert(bkg_count + zoom_count - 100 == 100);
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
  free(s);

  return 0;
}
