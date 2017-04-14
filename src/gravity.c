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

/* Some standard headers. */
#include <stdio.h>

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
      double a_grav[3] = {0., 0., 0.};

      /* Interact it with all other particles in the space.*/
      for (size_t j = 0; j < s->nr_gparts; ++j) {

        /* No self interaction */
        if (i == j) continue;

        struct gpart *gpj = &s->gparts[j];

        /* Compute the pairwise distance. */
        const double dx[3] = {gpi->x[0] - gpj->x[0],   // x
                              gpi->x[1] - gpj->x[1],   // y
                              gpi->x[2] - gpj->x[2]};  // z
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        const double r = sqrt(r2);
        const double ir = 1. / r;
        const double mj = gpj->mass;
        const double hi = gpi->epsilon;
        double f;
        const double f_lr = 1.;

        if (r >= hi) {

          /* Get Newtonian gravity */
          f = mj * ir * ir * ir * f_lr;

        } else {

          const double hi_inv = 1. / hi;
          const double hi_inv3 = hi_inv * hi_inv * hi_inv;
          const double ui = r * hi_inv;
          double W;

          kernel_grav_eval_double(ui, &W);

          /* Get softened gravity */
          f = mj * hi_inv3 * W * f_lr;

          // printf("r=%e hi=%e W=%e fac=%e\n", r, hi, W, f);
        }

        const double fdx[3] = {f * dx[0], f * dx[1], f * dx[2]};

        a_grav[0] -= fdx[0];
        a_grav[1] -= fdx[1];
        a_grav[2] -= fdx[2];
      }

      /* Store the exact answer */
      gpi->a_grav_exact[0] = a_grav[0] * const_G;
      gpi->a_grav_exact[1] = a_grav[1] * const_G;
      gpi->a_grav_exact[2] = a_grav[2] * const_G;

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

  int counter = 0;

  /* Some accumulators */
  float err_rel_max = 0.f;
  float err_rel_min = FLT_MAX;
  float err_rel_mean = 0.f;
  float err_rel_mean2 = 0.f;
  float err_rel_std = 0.f;

  char file_name[100];
  sprintf(file_name, "gravity_checks_step%d_order%d.dat", e->step,
          SELF_GRAVITY_MULTIPOLE_ORDER);
  FILE *file = fopen(file_name, "w");
  fprintf(file, "# Gravity accuracy test G = %16.8e\n",
          e->physical_constants->const_newton_G);
  fprintf(file, "# %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\n", "id",
          "pos[0]", "pos[1]", "pos[2]", "a_exact[0]", "a_exact[1]",
          "a_exact[2]", "a_grav[0]", "a_grav[1]", "a_grav[2]");

  for (size_t i = 0; i < s->nr_gparts; ++i) {

    struct gpart *gpi = &s->gparts[i];

    /* Is the particle was active and part of the subset to be tested ? */
    if (gpi->id_or_neg_offset % SWIFT_GRAVITY_FORCE_CHECKS == 0 &&
        gpart_is_starting(gpi, e)) {

      fprintf(file,
              "%18lld %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e "
              "%16.8e \n",
              gpi->id_or_neg_offset, gpi->x[0], gpi->x[1], gpi->x[2],
              gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2],
              gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2]);

      const float diff[3] = {gpi->a_grav[0] - gpi->a_grav_exact[0],
                             gpi->a_grav[1] - gpi->a_grav_exact[1],
                             gpi->a_grav[2] - gpi->a_grav_exact[2]};

      const float diff_norm =
          sqrtf(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      const float a_norm = sqrtf(gpi->a_grav_exact[0] * gpi->a_grav_exact[0] +
                                 gpi->a_grav_exact[1] * gpi->a_grav_exact[1] +
                                 gpi->a_grav_exact[2] * gpi->a_grav_exact[2]);

      /* Compute relative error */
      const float err_rel = diff_norm / a_norm;

      /* Check that we are not above tolerance */
      if (err_rel > rel_tol) {
        message(
            "Error too large !"
            "gp->a_grav=[%3.6e %3.6e %3.6e] "
            "gp->a_exact=[%3.6e %3.6e %3.6e], "
            "gp->num_interacted=%lld, err=%f",
            gpi->a_grav_exact[0], gpi->a_grav_exact[1], gpi->a_grav_exact[2],
            gpi->a_grav[0], gpi->a_grav[1], gpi->a_grav[2], gpi->num_interacted,
            err_rel);

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

  /* Be nice */
  fclose(file);

  /* Final operation on the stats */
  if (counter > 0) {
    err_rel_mean /= counter;
    err_rel_mean2 /= counter;
    err_rel_std = sqrtf(err_rel_mean2 - err_rel_mean * err_rel_mean);
  }

  /* Report on the findings */
  message("Checked gravity for %d gparts.", counter);
  message("Error on |a_grav|: min=%e max=%e mean=%e std=%e", err_rel_min,
          err_rel_max, err_rel_mean, err_rel_std);

#else
  error("Gravity checking function called without the corresponding flag.");
#endif
}
