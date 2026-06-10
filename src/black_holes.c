/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Robert McGibbon (robjmcgibbon@gmail.com)
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

/* Some standard headers. */
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Local headers. */
#include "active.h"
#include "black_holes.h"
#include "error.h"
#include "version.h"

#if defined(SWIFT_BH_STARS_DENSITY_CHECKS) && \
    !defined(BLACK_HOLES_HAVE_STAR_DENSITY)
#error \
    "Cannot use --enable-bh-stars-density-checks with a black hole model " \
    "that does not support the star-density loops."
#endif

struct exact_bh_stars_density_data {
  const struct engine *e;
  const struct space *s;
  int counter_global;
};

/**
 * @brief Mapper function for the exact BH star-density checks.
 *
 * @brief map_data The #bpart's.
 * @brief nr_bparts The number of black holes.
 * @brief extra_data Pointers to the structure containing global interaction
 * counters.
 */
void black_holes_exact_stars_density_compute_mapper(void *map_data,
                                                    int nr_bparts,
                                                    void *extra_data) {
#ifdef SWIFT_BH_STARS_DENSITY_CHECKS

  /* Unpack the data */
  struct bpart *restrict bparts = (struct bpart *)map_data;
  struct exact_bh_stars_density_data *data =
      (struct exact_bh_stars_density_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int counter = 0;

  for (int i = 0; i < nr_bparts; ++i) {

    struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;

    /* Is the particle active and part of the subset to be tested ? */
    if (id % SWIFT_BH_STARS_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      /* Get some information about the particle */
      const double pix[3] = {bpi->x[0], bpi->x[1], bpi->x[2]};
      const double hi = bpi->h_star;
      const float hi_inv = 1.f / hi;
      const float hig2 = hi * hi * kernel_gamma2;

      /* Be ready for the calculation */
      int N_density_exact = 0;
      double rho_exact = 0.;
      double n_exact = 0.;
      bpi->stars_inhibited_exact = 0;

      /* Interact it with all star particles in the space.*/
      for (int j = 0; j < (int)s->nr_sparts; ++j) {

        const struct spart *sj = &s->sparts[j];

        /* Skip the placeholder particles used by on-the-fly star formation:
         * they sit at the cell locations with zero mass and are not in the
         * cell counts the tree loops see. */
        if (sj->time_bin == time_bin_not_created) continue;

        /* Compute the pairwise distance. */
        double dx = sj->x[0] - pix[0];
        double dy = sj->x[1] - pix[1];
        double dz = sj->x[2] - pix[2];

        /* Now apply periodic BC */
        if (periodic) {
          dx = nearest(dx, dim[0]);
          dy = nearest(dy, dim[1]);
          dz = nearest(dz, dim[2]);
        }

        const double r2 = dx * dx + dy * dy + dz * dz;

        /* Within the kernel? */
        if (r2 < hig2) {

          const float mj = sj->mass;

          float wi, wi_dx;

          /* Kernel function */
          const float r = sqrtf(r2);
          const float ui = r * hi_inv;
          kernel_deval(ui, &wi, &wi_dx);

          /* Flag that we found an inhibited neighbour */
          if (spart_is_inhibited(sj, e)) {
            bpi->stars_inhibited_exact = 1;
          } else {

            /* Density */
            rho_exact += mj * wi;

            /* Number density */
            n_exact += wi;

            /* Number of neighbours */
            N_density_exact++;
          }
        }
      }

      /* Store the exact answer */
      bpi->N_stars_density_exact = N_density_exact;
      bpi->stars_rho_exact = rho_exact * pow_dimension(hi_inv);
      bpi->stars_n_exact = n_exact * pow_dimension(hi_inv);

      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error("BH star-density checking function called without the right flag.");
#endif
}

/**
 * @brief Compute the exact star-density interactions for a selection of black
 * holes by running a brute force loop over all the star particles in the
 * simulation.
 *
 * Will be incorrect over MPI.
 *
 * @param s The #space.
 * @param e The #engine.
 */
void black_holes_exact_stars_density_compute(struct space *s,
                                             const struct engine *e) {

#ifdef SWIFT_BH_STARS_DENSITY_CHECKS

  const ticks tic = getticks();

  struct exact_bh_stars_density_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;

  threadpool_map(&s->e->threadpool,
                 black_holes_exact_stars_density_compute_mapper, s->bparts,
                 s->nr_bparts, sizeof(struct bpart), 0, &data);

  if (e->verbose)
    message("Computed exact star densities for %d bparts (took %.3f %s). ",
            data.counter_global, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("BH star-density checking function called without the right flag.");
#endif
}

/**
 * @brief Check the black holes' star-density calculation against the
 * values obtained via the brute-force summation.
 *
 * @param s The #space.
 * @param e The #engine.
 * @param rel_tol Relative tolerance for the checks
 */
void black_holes_exact_stars_density_check(struct space *s,
                                           const struct engine *e,
                                           const double rel_tol) {

#ifdef SWIFT_BH_STARS_DENSITY_CHECKS

  const ticks tic = getticks();

  const struct bpart *bparts = s->bparts;
  const size_t nr_bparts = s->nr_bparts;

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "bh_stars_checks_swift_step%.4d.dat", e->step);

  /* Create files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  if (file_swift == NULL) error("Could not create file '%s'.", file_name_swift);
  fprintf(file_swift, "# BH star-density accuracy test - SWIFT DENSITIES\n");
  fprintf(file_swift, "# N= %d\n", SWIFT_BH_STARS_DENSITY_CHECKS);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift, "# %16s %16s %16s %16s %16s %7s %16s %16s\n", "id",
          "pos[0]", "pos[1]", "pos[2]", "h_star", "Nd", "rho", "n_rho");

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_bparts; ++i) {

    const struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;

    if (id % SWIFT_BH_STARS_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      fprintf(file_swift, "%18lld %16.8e %16.8e %16.8e %16.8e %7d %16.8e %16.8e\n",
              id, bpi->x[0], bpi->x[1], bpi->x[2], bpi->h_star,
              bpi->N_stars_density, bpi->stars_density.rho,
              bpi->stars_density.wcount);
    }
  }

  if (e->verbose)
    message("Written SWIFT densities in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  /* File name */
  char file_name_exact[100];
  sprintf(file_name_exact, "bh_stars_checks_exact_step%.4d.dat", e->step);

  /* Create files and write header */
  FILE *file_exact = fopen(file_name_exact, "w");
  if (file_exact == NULL) error("Could not create file '%s'.", file_name_exact);
  fprintf(file_exact, "# BH star-density accuracy test - EXACT DENSITIES\n");
  fprintf(file_exact, "# N= %d\n", SWIFT_BH_STARS_DENSITY_CHECKS);
  fprintf(file_exact, "# periodic= %d\n", s->periodic);
  fprintf(file_exact, "# Git Branch: %s\n", git_branch());
  fprintf(file_exact, "# Git Revision: %s\n", git_revision());
  fprintf(file_exact, "# %16s %16s %16s %16s %16s %7s %16s %16s\n", "id",
          "pos[0]", "pos[1]", "pos[2]", "h_star", "Nd", "rho_exact",
          "n_rho_exact");

  int wrong_rho = 0;
  int counter = 0;

  /* Output the exact densities and check against the SWIFT ones */
  for (size_t i = 0; i < nr_bparts; ++i) {

    const struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;
    const int found_inhibited = bpi->stars_inhibited_exact;

    if (id % SWIFT_BH_STARS_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      counter++;

      fprintf(file_exact,
              "%18lld %16.8e %16.8e %16.8e %16.8e %7d %16.8e %16.8e\n", id,
              bpi->x[0], bpi->x[1], bpi->x[2], bpi->h_star,
              bpi->N_stars_density_exact, bpi->stars_rho_exact,
              bpi->stars_n_exact);

      /* Check that we did not go above the threshold.
       * Note that we ignore particles that saw an inhibited particle as a
       * neighbour as we don't know whether that neighbour became inhibited in
       * that step or not.
       * Black holes with zero star neighbours are also ignored (their h_star
       * has hit the maximal allowed search radius). */
      if (!found_inhibited && bpi->N_stars_density_exact > 0 &&
          (fabsf(bpi->stars_n_exact) > 0.f) &&
          (fabsf(bpi->stars_density.wcount / bpi->stars_n_exact - 1.f) >
               rel_tol ||
           fabsf(bpi->stars_n_exact / bpi->stars_density.wcount - 1.f) >
               rel_tol)) {
        message("N_DENSITY: id=%lld swift=%e exact=%e N_true=%d N_swift=%d", id,
                bpi->stars_density.wcount, bpi->stars_n_exact,
                bpi->N_stars_density_exact, bpi->N_stars_density);
        wrong_rho++;
      }
    }
  }

  if (e->verbose)
    message("Written exact densities in file '%s'.", file_name_exact);

  /* Be nice */
  fclose(file_exact);

  if (wrong_rho)
    error(
        "Star density difference larger than the allowed tolerance for %d "
        "black holes! (out of %d black holes)",
        wrong_rho, counter);
  else
    message("Verified the star densities of %d black holes", counter);

  if (e->verbose)
    message("Writing brute-force density files took %.3f %s. ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("BH star-density checking function called without the right flag.");
#endif
}
