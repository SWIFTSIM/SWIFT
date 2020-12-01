/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Local headers. */
#include "active.h"
#include "black_holes.h"
#include "error.h"
#include "version.h"

struct exact_density_data {
  const struct engine *e;
  const struct space *s;
  int counter_global;
};

/**
 * @brief Mapper function for the exact BH checks.
 *
 * @brief map_data The #bpart's.
 * @brief nr_bparts The number of BH particles.
 * @brief extra_data Pointers to the structure containing global interaction
 * counters.
 */
void black_holes_exact_density_compute_mapper(void *map_data, int nr_bparts,
                                              void *extra_data) {
#ifdef SWIFT_BH_DENSITY_CHECKS

  /* Unpack the data */
  struct bpart *restrict bparts = (struct bpart *)map_data;
  struct exact_density_data *data = (struct exact_density_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int counter = 0;

  for (int i = 0; i < nr_bparts; ++i) {

    struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;

    /* Is the particle active and part of the subset to be tested ? */
    if (id % SWIFT_BH_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      /* Get some information about the particle */
      const double pix[3] = {bpi->x[0], bpi->x[1], bpi->x[2]};
      const double hi = bpi->h;
      const float hi_inv = 1.f / hi;
      const float hig2 = hi * hi * kernel_gamma2;

      /* Be ready for the calculation */
      int N_density_exact = 0;
      double rho_exact = 0.;
      double n_exact = 0.;

      /* Interact it with all other particles in the space.*/
      for (int j = 0; j < (int)s->nr_parts; ++j) {

        const struct part *pj = &s->parts[j];

        /* Compute the pairwise distance. */
        double dx = pj->x[0] - pix[0];
        double dy = pj->x[1] - pix[1];
        double dz = pj->x[2] - pix[2];

        /* Now apply periodic BC */
        if (periodic) {
          dx = nearest(dx, dim[0]);
          dy = nearest(dy, dim[1]);
          dz = nearest(dz, dim[2]);
        }

        const double r2 = dx * dx + dy * dy + dz * dz;

        /* Interact loop of type 1? */
        if (r2 < hig2) {

          const float mj = pj->mass;

          float wi, wi_dx;

          /* Kernel function */
          const float r = sqrtf(r2);
          const float ui = r * hi_inv;
          kernel_deval(ui, &wi, &wi_dx);

          /* Flag that we found an inhibited neighbour */
          if (part_is_inhibited(pj, e)) {
            bpi->inhibited_exact = 1;
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
      bpi->N_density_exact = N_density_exact;
      bpi->rho_exact = rho_exact * pow_dimension(hi_inv);
      bpi->n_exact = n_exact * pow_dimension(hi_inv);

      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error("Black_Holes checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Compute the exact interactions for a selection of BH particles
 * by running a brute force loop over all the particles in the simulation.
 *
 * Will be incorrect over MPI.
 *
 * @param s The #space.
 * @param e The #engine.
 */
void black_holes_exact_density_compute(struct space *s, const struct engine *e) {

#ifdef SWIFT_BH_DENSITY_CHECKS

  const ticks tic = getticks();

  struct exact_density_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;

  threadpool_map(&s->e->threadpool, black_holes_exact_density_compute_mapper,
                 s->bparts, s->nr_bparts, sizeof(struct bpart), 0, &data);

  if (e->verbose)
    message("Computed exact densities for %d parts (took %.3f %s). ",
            data.counter_global, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("Black_Holes checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Check the BH particles' density and force calculations against the
 * values obtained via the brute-force summation.
 *
 * @param s The #space.
 * @param e The #engine.
 * @param rel_tol Relative tolerance for the checks
 */
void black_holes_exact_density_check(struct space *s, const struct engine *e,
                                     const double rel_tol) {

#ifdef SWIFT_BH_DENSITY_CHECKS

  const ticks tic = getticks();

  const struct bpart *bparts = s->bparts;
  const size_t nr_bparts = s->nr_bparts;

  const double eta = e->black_holes_properties->eta_neighbours;
  const double N_ngb_target =
      (4. / 3.) * M_PI * pow_dimension(kernel_gamma * eta);
  const double N_ngb_max =
      N_ngb_target + 5. * e->black_holes_properties->delta_neighbours;
  const double N_ngb_min =
      N_ngb_target - 5. * e->black_holes_properties->delta_neighbours;

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "black_holes_checks_swift_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  fprintf(file_swift, "# Black_Holes accuracy test - SWIFT DENSITIES\n");
  fprintf(file_swift, "# N= %d\n", SWIFT_BH_DENSITY_CHECKS);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->black_holes_properties->delta_neighbours);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift, "# %16s %16s %16s %16s %16s %7s %7s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Nf", "rho", "n_rho",
          "N_ngb");

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_bparts; ++i) {

    const struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;

    const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                         kernel_gamma * bpi->h * bpi->h * bpi->h * bpi->n;

    if (id % SWIFT_BH_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      fprintf(
          file_swift,
          "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %16.8e %16.8e %16.8e\n",
          id, bpi->x[0], bpi->x[1], bpi->x[2], bpi->h, bpi->N_density, 0,
          bpi->rho, bpi->n, N_ngb);
    }
  }

  if (e->verbose)
    message("Written SWIFT densities in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  /* File name */
  char file_name_exact[100];
  sprintf(file_name_exact, "black_holes_checks_exact_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_exact = fopen(file_name_exact, "w");
  fprintf(file_exact, "# Black_Holes accuracy test - EXACT DENSITIES\n");
  fprintf(file_exact, "# N= %d\n", SWIFT_BH_DENSITY_CHECKS);
  fprintf(file_exact, "# periodic= %d\n", s->periodic);
  fprintf(file_exact, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->black_holes_properties->delta_neighbours);
  fprintf(file_exact, "# Git Branch: %s\n", git_branch());
  fprintf(file_exact, "# Git Revision: %s\n", git_revision());
  fprintf(file_exact, "# %16s %16s %16s %16s %16s %7s %7s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Nf", "rho_exact",
          "n_rho_exact", "N_ngb");

  int wrong_rho = 0;
  int wrong_n_ngb = 0;
  int counter = 0;

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_bparts; ++i) {

    const struct bpart *bpi = &bparts[i];
    const long long id = bpi->id;
    const int found_inhibited = bpi->inhibited_exact;

    const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                         kernel_gamma * bpi->h * bpi->h * bpi->h * bpi->n_exact;

    if (id % SWIFT_BH_DENSITY_CHECKS == 0 && bpart_is_starting(bpi, e)) {

      counter++;

      fprintf(
          file_exact,
          "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %16.8e %16.8e %16.8e\n",
          id, bpi->x[0], bpi->x[1], bpi->x[2], bpi->h, bpi->N_density_exact, 0,
          bpi->rho_exact, bpi->n_exact, N_ngb);

      /* Check that we did not go above the threshold.
       * Note that we ignore particles that saw an inhibted particle as a
       * neighbour as we don't know whether that neighbour became inhibited in
       * that step or not. */
      if (!found_inhibited && bpi->N_density_exact != bpi->N_density &&
          (fabsf(bpi->rho / bpi->rho_exact - 1.f) > rel_tol ||
           fabsf(bpi->rho_exact / bpi->rho - 1.f) > rel_tol)) {
        message("RHO: id=%lld swift=%e exact=%e N_true=%d N_swift=%d", id,
                bpi->rho, bpi->rho_exact, bpi->N_density_exact, bpi->N_density);
        wrong_rho++;
      }

      if (!found_inhibited && (N_ngb > N_ngb_max || N_ngb < N_ngb_min)) {

        message("N_NGB: id=%lld exact=%f N_true=%d N_swift=%d", id, N_ngb,
                bpi->N_density_exact, bpi->N_density);

        wrong_n_ngb++;
      }
    }
  }

  if (e->verbose)
    message("Written exact densities in file '%s'.", file_name_exact);

  /* Be nice */
  fclose(file_exact);

  if (wrong_rho)
    error(
        "Density difference larger than the allowed tolerance for %d "
        "BH particles! (out of %d particles)",
        wrong_rho, counter);
  else
    message("Verified %d BH particles", counter);

  if (wrong_n_ngb)
    error(
        "N_ngb difference larger than the allowed tolerance for %d "
        "BH particles! (out of %d particles)",
        wrong_n_ngb, counter);
  else
    message("Verified %d BH particles", counter);

  if (e->verbose)
    message("Writting brute-force density files took %.3f %s. ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("Black_Holes checking function called without the corresponding flag.");
#endif
}
