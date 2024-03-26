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
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "parser.h"
#include "space.h"
#include "swift.h"
#include "zoom_region/zoom_init.h"
#include "zoom_region/zoom_regrid.h"

void make_mock_space(struct space *s) {

  /* Define the boxsize. */
  s->dim[0] = 1000;
  s->dim[1] = 1000;
  s->dim[2] = 1000;

  /* Define the gpart count (this is a particle per cell) */
  s->nr_gparts = 18;

  /* We need the engine to be NULL for the logic. */
  s->e = NULL;

  /* Allocate memory for the gparts. */
  struct gpart *gparts =
      (struct gpart *)malloc(s->nr_gparts * sizeof(struct gpart));
  bzero(gparts, s->nr_gparts * sizeof(struct gpart));

  /* Define the corners of the region */
  const double mid = 600;
  const double offset = 9;
  double cube_corners[8][3] = {{mid - offset, mid - offset, mid - offset},
                               {mid - offset, mid - offset, mid + offset},
                               {mid - offset, mid + offset, mid - offset},
                               {mid - offset, mid + offset, mid + offset},
                               {mid + offset, mid - offset, mid - offset},
                               {mid + offset, mid - offset, mid + offset},
                               {mid + offset, mid + offset, mid - offset},
                               {mid + offset, mid + offset, mid + offset}};

  /* Loop over the gparts and set up baxckground and zoom particles. */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    gparts[i].mass = 1.0;

    /* Handle background and zoom region particles differently. */
    if (i < 10) {
      /* Set background particles to be evenly spaced. */
      gparts[i].x[0] = s->dim[0] / s->nr_gparts * i;
      gparts[i].x[1] = s->dim[1] / s->nr_gparts * i;
      gparts[i].x[2] = s->dim[2] / s->nr_gparts * i;
      gparts[i].type = swift_type_dark_matter_background;

    } else {
      /* Set zoom region particles to be at the corners of the region. */
      gparts[i].x[0] = cube_corners[i - 10][0];
      gparts[i].x[1] = cube_corners[i - 10][1];
      gparts[i].x[2] = cube_corners[i - 10][2];
      gparts[i].type = swift_type_dark_matter;
    }
  }

  s->gparts = gparts;
}

void make_gparts_grid(struct space *s) {}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  const char *input_file = argv[1];

  /* Create a structure to read file into. */
  struct swift_params param_file;

  /* Read the parameter file. */
  parser_read_file(input_file, &param_file);

  /* Create a space structure. (this creates fake particles to get the layout
   * of cells right but we will modify this shortly to get a particle per
   * cell) */
  struct space *s = malloc(sizeof(struct space));
  bzero(s, sizeof(struct space));
  make_mock_space(s);

  /* Flag that we are running a zoom. */
  s->with_zoom_region = 1;

  /* Run the zoom_init function. */
  zoom_props_init(&param_file, s, 0);

  /* Run the regridding. */
  space_regrid(s, 1);

  free(s->cells_top);
  free(s->gparts);
  free(s->zoom_props);
  free(s);

  return 0;
}
