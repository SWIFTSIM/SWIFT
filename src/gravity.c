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

  const double const_G = e->physical_constants->const_newton_G;

  int counter = 0;

  /* Some accumulators */
  float err_rel[3];
  float err_rel_max[3] = {0.f, 0.f, 0.f};
  float err_rel_min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
  float err_rel_mean[3] = {0.f, 0.f, 0.f};
  float err_rel_mean2[3] = {0.f, 0.f, 0.f};
  float err_rel_std[3] = {0.f, 0.f, 0.f};

  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Is the particle was active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_starting(gpi, e)) {

      /* Compute relative error */
      for (int k = 0; k < 3; ++k)
        if (fabsf(gpi->a_grav_exact[k]) > FLT_EPSILON * const_G)
          err_rel[k] = (gpi->a_grav[k] - gpi->a_grav_exact[k]) /
                       fabsf(gpi->a_grav_exact[k]);
        else
          err_rel[k] = 0.f;

      /* Check that we are not above tolerance */
      if (fabsf(err_rel[0]) > rel_tol || fabsf(err_rel[1]) > rel_tol ||
          fabsf(err_rel[2]) > rel_tol)
        error(
            "Error too large ! gp->a_grav=[%e %e %e] gp->a_exact=[%e %e %e], "
            "gp->mass_interacted=%e",
            gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2],
            gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2],
            gpi->mass_interacted);

      /* Construct some statistics */
      for (int k = 0; k < 3; ++k) {
        err_rel_max[k] = max(err_rel_max[k], fabsf(err_rel[k]));
        err_rel_min[k] = min(err_rel_min[k], fabsf(err_rel[k]));
        err_rel_mean[k] += err_rel[k];
        err_rel_mean2[k] += err_rel[k] * err_rel[k];
      }

      counter++;
    }
  }

  /* Final operation on the stats */
  if (counter > 0) {
    for (int k = 0; k < 3; ++k) {
      err_rel_mean[k] /= counter;
      err_rel_mean2[k] /= counter;
      err_rel_std[k] =
          sqrtf(err_rel_mean2[k] - err_rel_mean[k] * err_rel_mean[k]);
    }
  }

  /* Report on the findings */
  message("Checked gravity for %d gparts.", counter);
  for (int k = 0; k < 3; ++k)
    message("Error on a_grav[%d]: min=%e max=%e mean=%e std=%e", k,
            err_rel_min[k], err_rel_max[k], err_rel_mean[k], err_rel_std[k]);

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}
