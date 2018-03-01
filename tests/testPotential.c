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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "swift.h"
#include "runner_doiact_grav.h"

const int num_tests = 100;

int main() {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);
  
  struct engine e;
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.timeBase = 1e-10;

  struct runner r;
  bzero(&r, sizeof(struct runner));
  r.e = &e;

  gravity_cache_init(&r.ci_gravity_cache, num_tests * 2);
  
  
  /* Let's create one cell with a massive particle and a bunch of test particles */
  struct cell c;
  bzero(&c, sizeof(struct cell));
  c.width[0] = 1.;
  c.width[1] = 1.;
  c.width[2] = 1.;
  c.gcount = 1 + num_tests;
  c.ti_old_gpart = 8;
  c.ti_gravity_end_min = 8;
  c.ti_gravity_end_max = 8;
  
  posix_memalign((void**) &c.gparts, gpart_align, c.gcount * sizeof(struct gpart));
  bzero(c.gparts, c.gcount * sizeof(struct gpart));

  /* Create the massive particle */
  c.gparts[0].x[0] = 0.;
  c.gparts[0].x[1] = 0.5;
  c.gparts[0].x[2] = 0.5;
  c.gparts[0].mass = 1.;
  c.gparts[0].epsilon = 0.;
  c.gparts[0].time_bin = 1;
  c.gparts[0].type = swift_type_dark_matter;
  c.gparts[0].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
  c.gparts[0].ti_drift = 8;
#endif
  
  /* Create the mass-less particles */
  for (int n = 1; n < num_tests + 1; ++n) {

    struct gpart *gp = &c.gparts[n];

    gp->x[0] = n / ((double) num_tests);
    gp->x[1] = 0.5;
    gp->x[2] = 0.5;
    gp->mass = 0.;
    gp->epsilon = 0.;
    gp->time_bin = 1;
    gp->type = swift_type_dark_matter;
    gp->id_or_neg_offset = n + 1;
#ifdef SWIFT_DEBUG_CHECKS
    gp->ti_drift = 8;
#endif

  }
  
  /* Now compute the forces */
  runner_doself_grav_pp_full(&r, &c);
  
  free(c.gparts);
  return 0;
}
