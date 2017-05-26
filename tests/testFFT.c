/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Some standard headers. */

#include "../config.h"

#ifndef HAVE_FFTW

int main() { return 0; }

#else

/* Some standard headers. */
#include <stdlib.h>
#include <string.h>

/* Includes. */
#include "swift.h"

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Make one particle */
  int nr_gparts = 1;
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&gparts, 64, nr_gparts * sizeof(struct gpart)) !=
      0)
    error("Impossible to allocate memory for gparts.");
  bzero(gparts, nr_gparts * sizeof(struct gpart));

  gparts[0].x[0] = 0.3;
  gparts[0].x[1] = 0.8;
  gparts[0].x[2] = 0.2;
  gparts[0].mass = 1.f;

  /* Read the parameter file */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  parser_read_file("fft_params.yml", params);

  /* Initialise the gravity properties */
  struct gravity_props gravity_properties;
  gravity_props_init(&gravity_properties, params);

  /* Build the infrastructure */
  struct space space;
  double dim[3] = {1., 1., 1.};
  space_init(&space, params, dim, NULL, gparts, NULL, 0, nr_gparts, 0, 1, 1, 1,
             0, 0);

  struct engine engine;
  engine.s = &space;
  space.e = &engine;
  engine.time = 0.1f;
  engine.ti_current = 0;
  engine.ti_old = 0;
  engine.max_active_bin = num_time_bins;
  engine.gravity_properties = &gravity_properties;
  engine.nr_threads = 1;

  struct runner runner;
  runner.e = &engine;

  /* Initialize the threadpool. */
  threadpool_init(&engine.threadpool, engine.nr_threads);

  space_rebuild(&space, 0);

  /* Run the FFT task */
  runner_do_grav_fft(&runner, 1);

  /* Now check that we got the right answer */
  int nr_cells = space.nr_cells;
  double *r = malloc(nr_cells * sizeof(double));
  double *pot = malloc(nr_cells * sizeof(double));
  double *pot_exact = malloc(nr_cells * sizeof(double));

  // FILE *file = fopen("potential.dat", "w");
  for (int i = 0; i < nr_cells; ++i) {
    pot[i] = space.multipoles_top[i].pot.F_000;
    double dx =
        nearest(space.multipoles_top[i].CoM[0] - gparts[0].x[0], dim[0]);
    double dy =
        nearest(space.multipoles_top[i].CoM[1] - gparts[0].x[1], dim[1]);
    double dz =
        nearest(space.multipoles_top[i].CoM[2] - gparts[0].x[2], dim[2]);
    r[i] = sqrt(dx * dx + dy * dy + dz * dz);
    if (r[i] > 0) pot_exact[i] = -1. / r[i];
    // fprintf(file, "%e %e %e\n", r[i], pot[i], pot_exact[i]);
  }
  // fclose(file);

  /* Clean up */
  free(r);
  free(pot);
  free(pot_exact);
  free(params);
  free(gparts);
  return 0;
}

#endif
