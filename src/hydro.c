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
#include "error.h"
#include "hydro_properties.h"
#include "version.h"

struct exact_density_data {
  const struct engine *e;
  const struct space *s;
  int counter_global;
  int check_force;
};

/**
 * @brief Mapper function for the exact hydro checks.
 *
 * @brief map_data The #part's.
 * @brief nr_parts The number of particles.
 * @brief extra_data Pointers to the structure containing global interaction
 * counters.
 */
void hydro_exact_density_compute_mapper(void *map_data, int nr_parts,
                                        void *extra_data) {
#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  /* Unpack the data */
  struct part *restrict parts = (struct part *)map_data;
  struct exact_density_data *data = (struct exact_density_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int check_force = data->check_force;
  int counter = 0;

  for (int i = 0; i < nr_parts; ++i) {

    struct part *pi = &parts[i];
    const long long id = pi->id;

    /* Is the particle active and part of the subset to be tested ? */
    if (id % SWIFT_HYDRO_DENSITY_CHECKS == 0 && part_is_starting(pi, e)) {

      /* Get some information about the particle */
      const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
      const double hi = pi->h;
      const float hi_inv = 1.f / hi;
      const float hig2 = hi * hi * kernel_gamma2;

      /* Be ready for the calculation */
      int N_density_exact = 0;
      int N_gradient_exact = 0;
      int N_force_exact = 0;
      int N_limiter_exact = 0;
      double rho_exact = 0.;
      double n_density_exact = 0.;
      double n_gradient_exact = 0.;
      double n_force_exact = 0.;
      double n_limiter_exact = 0.;

      /* Interact it with all other particles in the space.*/
      for (int j = 0; j < (int)s->nr_parts; ++j) {

        const struct part *pj = &s->parts[j];
        const double hj = pj->h;
        const float hj_inv = 1.f / hj;
        const float hjg2 = hj * hj * kernel_gamma2;

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

          // if (id != pj->id && r2 < 1e-10 && e->step == SCHECK) continue;

          const float mj = pj->mass;

          float wi, wi_dx;

          /* Kernel function */
          const float r = sqrtf(r2);
          const float ui = r * hi_inv;
          kernel_deval(ui, &wi, &wi_dx);

          if (id == ICHECK && e->step == SCHECK)
            message("pi=%lld interacting with pj=%lld inhibited=%d r2=%e", id,
                    pj->id, pj->time_bin == time_bin_inhibited, r2);

          /* Flag that we found an inhibited neighbour */
          if (part_is_inhibited(pj, e)) {
            pi->inhibited_exact = 1;
          } else {

            /* Density */
            rho_exact += mj * wi;

            n_density_exact += wi;

            /* Number of neighbours */
            N_density_exact++;

            n_limiter_exact += wi;

            /* Number of neighbours */
            N_limiter_exact++;

            /* Gradient */
            n_gradient_exact += wi;

            /* Number of neighbours */
            N_gradient_exact++;
          }
        }

        /* Interact loop of type 2? */
        if (check_force && (pi != pj) && (r2 < hig2 || r2 < hjg2)) {

          float wi, wi_dx;
          float wj, wj_dx;

          /* Kernel function */
          const float r = sqrtf(r2);
          const float ui = r * hi_inv;
          kernel_deval(ui, &wi, &wi_dx);
          const float uj = r * hj_inv;
          kernel_deval(uj, &wj, &wj_dx);

          /* Flag that we found an inhibited neighbour */
          if (part_is_inhibited(pj, e)) {
            pi->inhibited_exact = 1;
          } else {
            /* Force count */
            n_force_exact += wi + wj;
          }

          /* Number of neighbours */
          N_force_exact++;
        }
      }

      /* Store the exact answer */
      pi->N_density_exact = N_density_exact;
      pi->N_gradient_exact = N_gradient_exact;
      pi->N_force_exact = N_force_exact;
      pi->limiter_data.N_limiter_exact = N_limiter_exact;

      pi->rho_exact = rho_exact * pow_dimension(hi_inv);
      pi->n_density_exact = n_density_exact * pow_dimension(hi_inv);
      pi->n_gradient_exact = n_gradient_exact;
      pi->n_force_exact = n_force_exact;
      pi->limiter_data.n_limiter_exact =
          n_limiter_exact * pow_dimension(hi_inv);

      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error("Hydro checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Compute the exact interactions for a selection of gas particles
 * by running a brute force loop over all the particles in the simulation.
 *
 * Will be incorrect over MPI.
 *
 * @param s The #space.
 * @param e The #engine.
 * @param check_force Whether or not to run checks also for the force loop.
 */
void hydro_exact_density_compute(struct space *s, const struct engine *e,
                                 const int check_force) {

#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const ticks tic = getticks();

  struct exact_density_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;
  data.check_force = check_force;

  threadpool_map(&s->e->threadpool, hydro_exact_density_compute_mapper,
                 s->parts, s->nr_parts, sizeof(struct part), 0, &data);

  if (e->verbose)
    message("Computed exact densities for %d parts (took %.3f %s). ",
            data.counter_global, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("Hydro checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Check the gas particles' density and force calculations against the
 * values obtained via the brute-force summation.
 *
 * @param s The #space.
 * @param e The #engine.
 * @param rel_tol Relative tolerance for the checks
 * @param check_force Whether or not to run checks also for the force loop.
 */
void hydro_exact_density_check(struct space *s, const struct engine *e,
                               const float rel_tol, const int check_force) {

#ifdef SWIFT_HYDRO_DENSITY_CHECKS

  const ticks tic = getticks();

  const struct part *parts = s->parts;
  const size_t nr_parts = s->nr_parts;

  const int with_limiter = e->policy & engine_policy_timestep_limiter;
  const double eta = e->hydro_properties->eta_neighbours;
  const float h_max = e->hydro_properties->h_max;
  const float h_min = e->hydro_properties->h_min;
  const double N_ngb_target =
      (4. / 3.) * M_PI * pow_dimension(kernel_gamma * eta);
  const double N_ngb_max =
      N_ngb_target + 5. * e->hydro_properties->delta_neighbours;
  const double N_ngb_min =
      N_ngb_target - 5. * e->hydro_properties->delta_neighbours;

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "hydro_checks_swift_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  fprintf(file_swift, "# Hydro accuracy test - SWIFT DENSITIES\n");
  fprintf(file_swift, "# N= %d\n", SWIFT_HYDRO_DENSITY_CHECKS);
  fprintf(file_swift, "# h_max= %e\n", h_max);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->hydro_properties->delta_neighbours);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift,
          "# %16s %16s %16s %16s %16s %7s %7s %7s %16s %16s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Ng", "Nf", "rho",
          "n_gradient", "n_force", "n_limiter", "N_ngb");

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_parts; ++i) {

    const struct part *pi = &parts[i];
    const long long id = pi->id;
    if (pi->limited_part) continue;

    const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                         kernel_gamma * pi->h * pi->h * pi->h * pi->n_density;

    if (id % SWIFT_HYDRO_DENSITY_CHECKS == 0 && part_is_starting(pi, e)) {

      fprintf(file_swift,
              "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %7d %16.8e %16.8e "
              "%16.8e %16.8e %16.8e\n",
              id, pi->x[0], pi->x[1], pi->x[2], pi->h, pi->N_density,
              pi->N_gradient, pi->N_force, pi->rho, pi->n_gradient, pi->n_force,
              /* pi->limiter_data.n_limiter, */ 0., N_ngb);
    }
  }

  if (e->verbose)
    message("Written SWIFT densities in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  /* File name */
  char file_name_exact[100];
  sprintf(file_name_exact, "hydro_checks_exact_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_exact = fopen(file_name_exact, "w");
  fprintf(file_exact, "# Hydro accuracy test - EXACT DENSITIES\n");
  fprintf(file_exact, "# N= %d\n", SWIFT_HYDRO_DENSITY_CHECKS);
  fprintf(file_exact, "# h_max= %e\n", h_max);
  fprintf(file_exact, "# periodic= %d\n", s->periodic);
  fprintf(file_exact, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->hydro_properties->delta_neighbours);
  fprintf(file_exact, "# Git Branch: %s\n", git_branch());
  fprintf(file_exact, "# Git Revision: %s\n", git_revision());
  fprintf(file_exact,
          "# %16s %16s %16s %16s %16s %7s %7s %7s %16s %16s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Ng", "Nf",
          "rho_exact", "n_gradient", "n_force_exact", "n_limiter_exact",
          "N_ngb");

  int wrong_rho = 0;
  int wrong_gradient = 0;
  int wrong_n_ngb = 0;
  int wrong_n_force = 0;
  int wrong_limiter = 0;
  int counter = 0;

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_parts; ++i) {

    const struct part *pi = &parts[i];
    const long long id = pi->id;
    const int found_inhibited = pi->inhibited_exact;
    const int h_max_limited = pi->h >= 0.99 * h_max;
    const int h_min_limited = pi->h <= 1.01 * h_min;
    if (pi->limited_part) continue;

    if (id % SWIFT_HYDRO_DENSITY_CHECKS == 0 && part_is_starting(pi, e)) {

      counter++;

      const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                           kernel_gamma * pi->h * pi->h * pi->h *
                           pi->n_density_exact;

      fprintf(file_exact,
              "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %7d %16.8e %16.8e "
              "%16.8e %16.8e %16.8e\n",
              id, pi->x[0], pi->x[1], pi->x[2], pi->h, pi->N_density_exact,
              pi->N_gradient_exact, pi->N_force_exact, pi->rho_exact,
              pi->n_gradient_exact, pi->n_force_exact,
              /* pi->limiter_data.n_limiter_exact, */ 0., N_ngb);

      /* Check that we did not go above the threshold.
       * Note that we ignore particles that saw an inhibited particle as a
       * neighbour as we don't know whether that neighbour became inhibited in
       * that step or not. */
      if (!found_inhibited && pi->N_density != pi->N_density_exact &&
          (fabsf(pi->rho / pi->rho_exact - 1.f) > rel_tol ||
           fabsf(pi->rho_exact / pi->rho - 1.f) > rel_tol)) {
        message("RHO: id=%lld swift=%e exact=%e N=%d N_exact=%d h_max=%d", id,
                pi->rho, pi->rho_exact, pi->N_density, pi->N_density_exact,
                h_max_limited);
        wrong_rho++;
      }
#if defined(SPHENIX_SPH)
      if (check_force && !found_inhibited &&
          (fabsf(pi->n_gradient / pi->n_gradient_exact - 1.f) > rel_tol ||
           fabsf(pi->n_gradient_exact / pi->n_gradient - 1.f) > rel_tol)) {
        message("GRADIENT: id=%lld swift=%e exact=%e", id, pi->n_gradient,
                pi->n_gradient_exact);
        wrong_gradient++;
      }
#endif
      if (check_force && !found_inhibited &&
          (fabsf(pi->n_force / pi->n_force_exact - 1.f) > 10. * rel_tol ||
           fabsf(pi->n_force_exact / pi->n_force - 1.f) > 10. * rel_tol)) {
        message("N_FORCE: id=%lld swift=%e exact=%e", id, pi->n_force,
                pi->n_force_exact);
        wrong_n_force++;
      }
      if (with_limiter && check_force && !found_inhibited && !h_max_limited &&
          !h_min_limited &&
          (fabsf(pi->limiter_data.n_limiter / pi->limiter_data.n_limiter_exact -
                 1.f) > rel_tol ||
           fabsf(pi->limiter_data.n_limiter_exact / pi->limiter_data.n_limiter -
                 1.f) > rel_tol)) {
        message(
            "LIMITER: id=%lld swift=%e exact=%e N_limiter=%d "
            "N_limiter_exact=%d",
            id, pi->limiter_data.n_limiter, pi->limiter_data.n_limiter_exact,
            pi->limiter_data.N_limiter, pi->limiter_data.N_limiter_exact);
        wrong_limiter++;
      }

      if (!found_inhibited && !h_max_limited && !h_min_limited &&
          (N_ngb > N_ngb_max || N_ngb < N_ngb_min)) {

        message(
            "N_NGB: id=%lld exact=%f expected=%f/%f N_true=%d N_swift=%d h=%e "
            "time_bin=%d",
            id, N_ngb, N_ngb_target, N_ngb_max - N_ngb_target,
            pi->N_density_exact, pi->N_density, pi->h, pi->time_bin);

        wrong_n_ngb++;
      }
    }
  }

  if (e->verbose)
    message("Written exact densities in file '%s'.", file_name_exact);

  /* Be nice */
  fclose(file_exact);

  if (wrong_n_ngb)
    error(
        "N_ngb difference larger than the allowed tolerance for %d "
        "gas particles! (out of %d particles)",
        wrong_n_ngb, counter);
  else
    message("Verified %d gas particles", counter);

  if (wrong_rho)
    error(
        "Density difference larger than the allowed tolerance for %d "
        "particles!",
        wrong_rho);

  if (wrong_gradient)
    error(
        "Gradient difference larger than the allowed tolerance for %d "
        "particles!",
        wrong_gradient);

  if (wrong_n_force)
    message(
        "N_force difference larger than the allowed tolerance for %d "
        "particles!",
        wrong_n_force);

  if (wrong_limiter)
    error(
        "Limiter difference larger than the allowed tolerance for %d "
        "particles! (out of %d particles)",
        wrong_limiter, counter);

  if (e->verbose)
    message("Writing brute-force density files took %.3f %s. ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("Hydro checking function called without the corresponding flag.");
#endif
}
