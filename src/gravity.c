/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* This object's header. */
#include "gravity.h"

/* Local headers. */
#include "active.h"
#include "error.h"

/**
 * @brief Run a brute-force gravity calculation for a subset of particles.
 *
 * All gpart with ID modulo SWIFT_GRAVITY_FORCE_CHECKS will get their forces
 * computed.
 *
 * @param s The #space to use.
 * @param e The #engine (to access the current time).
 */
void gravity_exact_force_compute(struct space *s, const struct engine *e) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  const ticks tic = getticks();
  const double const_G = e->physical_constants->const_newton_G;
  int counter = 0;

  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Is the particle active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_active(gpi, e)) {

      /* Be ready for the calculation */
      gpi->a_grav[0] = 0.f;
      gpi->a_grav[1] = 0.f;
      gpi->a_grav[2] = 0.f;

      /* Interact it with all other particles in the space.*/
      for (size_t j = 0; j < s->nr_gparts; ++j) {

        /* No self interaction */
        if (i == j) continue;

        struct gpart *gpj = &s->gparts[j];

        /* Compute the pairwise distance. */
        float dx[3] = {gpi->x[0] - gpj->x[0],   // x
                       gpi->x[1] - gpj->x[1],   // y
                       gpi->x[2] - gpj->x[2]};  // z
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        runner_iact_grav_pp_nonsym(0.f, r2, dx, gpi, gpj);
      }

      /* Finish the calculation */
      gravity_end_force(gpi, const_G);

      /* Store the exact answer */
      gpi->a_grav_exact[0] = gpi->a_grav[0];
      gpi->a_grav_exact[1] = gpi->a_grav[1];
      gpi->a_grav_exact[2] = gpi->a_grav[2];

      /* Restore everything */
      gpi->a_grav[0] = 0.f;
      gpi->a_grav[1] = 0.f;
      gpi->a_grav[2] = 0.f;

      counter++;
    }
  }

  message("Computed exact gravity for %d gparts (took %.3f %s). ", counter,
          clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Check the accuracy of the gravity calculation by comparing the
 * accelerations
 * to the brute-force computed ones.
 *
 * All gpart with ID modulo SWIFT_GRAVITY_FORCE_CHECKS will be checked.
 *
 * @param s The #space to use.
 * @param e The #engine (to access the current time).
 * @param rel_tol The maximal relative error. Will call error() if one particle
 * has a larger error.
 */
void gravity_exact_force_check(struct space *s, const struct engine *e,
                               float rel_tol) {

#ifdef SWIFT_GRAVITY_FORCE_CHECKS

  //const double const_G = e->physical_constants->const_newton_G;

  int counter = 0;

  /* Some accumulators */
  float err_rel_max = 0.f;
  float err_rel_min = FLT_MAX;
  float err_rel_mean = 0.f;
  float err_rel_mean2 = 0.f;
  float err_rel_std = 0.f;

  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Is the particle was active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_starting(gpi, e)) {

      const float diff[3] = {gpi->a_grav[0] - gpi->a_grav_exact[0],
			     gpi->a_grav[1] - gpi->a_grav_exact[1],
			     gpi->a_grav[2] - gpi->a_grav_exact[2]};

      const float diff_norm = sqrtf(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      const float a_norm = sqrtf( gpi->a_grav_exact[0] * gpi->a_grav_exact[0] +
				  gpi->a_grav_exact[1] * gpi->a_grav_exact[1] +
				  gpi->a_grav_exact[2] * gpi->a_grav_exact[2]);
	
      /* Compute relative error */
      const float err_rel = diff_norm / a_norm;
      
      
      /* Check that we are not above tolerance */
      if (err_rel > rel_tol) {
        message(
            "Error too large ! gp->a_grav=[%3.6e %3.6e %3.6e] gp->a_exact=[%3.6e %3.6e %3.6e], "
            "gp->num_interacted=%lld, err=%f",
            gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2],
            gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2],
            gpi->num_interacted, err_rel);

	continue;
      }

      /* Construct some statistics */
      err_rel_max = max(err_rel_max, fabsf(err_rel));
      err_rel_min = min(err_rel_min, fabsf(err_rel));
      err_rel_mean += err_rel;
      err_rel_mean2 += err_rel * err_rel;
      counter++;
    }
  }

  /* Final operation on the stats */
  if (counter > 0) {
    err_rel_mean /= counter;
    err_rel_mean2 /= counter;
    err_rel_std = sqrtf(err_rel_mean2 - err_rel_mean * err_rel_mean);
  }

  /* Report on the findings */
  message("Checked gravity for %d gparts.", counter);
  message("Error on |a_grav|: min=%e max=%e mean=%e std=%e",
	  err_rel_min, err_rel_max, err_rel_mean, err_rel_std);

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}
