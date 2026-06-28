/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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
#include "error.h"
#include "sink_properties.h"
#include "version.h"

struct exact_density_data {
  const struct engine *e;
  const struct space *s;
  int counter_global;
};

/**
 * @brief Mapper function for the exact sink checks.
 *
 * @brief map_data The #sink.
 * @brief nr_sinks The number of star particles.
 * @brief extra_data Pointers to the structure containing global interaction
 * counters.
 */
void sink_exact_density_compute_mapper(void *map_data, int nr_sinks,
                                       void *extra_data) {
#ifdef SWIFT_SINK_DENSITY_CHECKS

  /* Unpack the data */
  struct sink *restrict sinks = (struct sink *)map_data;
  struct exact_density_data *data = (struct exact_density_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  int counter = 0;

  for (int i = 0; i < nr_sinks; ++i) {

    struct sink *si = &sinks[i];
    const long long id = si->id;

    /* Is the particle active and part of the subset to be tested ? */
    if (id % SWIFT_SINK_DENSITY_CHECKS == 0 && sink_is_starting(si, e)) {

      /* Get some information about the particle */
      const double pix[3] = {si->x[0], si->x[1], si->x[2]};
      const double hi = si->h;
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
            si->inhibited_check_exact = 1;
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
      si->N_check_density_exact = N_density_exact;
      si->rho_check_exact = rho_exact * pow_dimension(hi_inv);
      si->n_check_exact = n_exact * pow_dimension(hi_inv);

      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error("Sink checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Compute the exact interactions for a selection of star particles
 * by running a brute force loop over all the particles in the simulation.
 *
 * Will be incorrect over MPI.
 *
 * @param s The #space.
 * @param e The #engine.
 */
void sink_exact_density_compute(struct space *s, const struct engine *e) {

#ifdef SWIFT_SINK_DENSITY_CHECKS

  const ticks tic = getticks();

  struct exact_density_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;

  threadpool_map(&s->e->threadpool, sink_exact_density_compute_mapper, s->sinks,
                 s->nr_sinks, sizeof(struct sink), 0, &data);

  if (e->verbose)
    message("Computed exact densities for %d sinks (took %.3f %s). ",
            data.counter_global, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("Sink checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Check the star particles' density and force calculations against the
 * values obtained via the brute-force summation.
 *
 * @param s The #space.
 * @param e The #engine.
 * @param rel_tol Relative tolerance for the checks
 */
void sink_exact_density_check(struct space *s, const struct engine *e,
                              const double rel_tol) {

#ifdef SWIFT_SINK_DENSITY_CHECKS

  const ticks tic = getticks();

  const struct sink *sinks = s->sinks;
  const size_t nr_sinks = s->nr_sinks;

  const double eta = e->sink_properties->eta_neighbours;
  const double N_ngb_target =
      (4. / 3.) * M_PI * pow_dimension(kernel_gamma * eta);
  const double N_ngb_max =
      N_ngb_target + 2. * e->sink_properties->delta_neighbours;
  const double N_ngb_min =
      N_ngb_target - 2. * e->sink_properties->delta_neighbours;

  /* File name */
  char file_name_swift[100];
  sprintf(file_name_swift, "sink_checks_swift_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_swift = fopen(file_name_swift, "w");
  if (file_swift == NULL) error("Could not create file '%s'.", file_name_swift);
  fprintf(file_swift, "# Sink accuracy test - SWIFT DENSITIES\n");
  fprintf(file_swift, "# N= %d\n", SWIFT_SINK_DENSITY_CHECKS);
  fprintf(file_swift, "# periodic= %d\n", s->periodic);
  fprintf(file_swift, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->sink_properties->delta_neighbours);
  fprintf(file_swift, "# Git Branch: %s\n", git_branch());
  fprintf(file_swift, "# Git Revision: %s\n", git_revision());
  fprintf(file_swift, "# %16s %16s %16s %16s %16s %7s %7s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Nf", "rho", "n_rho",
          "N_ngb");

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_sinks; ++i) {

    const struct sink *si = &sinks[i];
    const long long id = si->id;

    const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                         kernel_gamma * si->h * si->h * si->h * si->n_check;

    if (id % SWIFT_SINK_DENSITY_CHECKS == 0 && sink_is_starting(si, e)) {

      fprintf(
          file_swift,
          "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %16.8e %16.8e %16.8e\n",
          id, si->x[0], si->x[1], si->x[2], si->h, si->N_check_density, 0,
          si->rho_check, si->n_check, N_ngb);
    }
  }

  if (e->verbose)
    message("Written SWIFT densities in file '%s'.", file_name_swift);

  /* Be nice */
  fclose(file_swift);

  /* File name */
  char file_name_exact[100];
  sprintf(file_name_exact, "sink_checks_exact_step%.4d.dat", e->step);

  /* Creare files and write header */
  FILE *file_exact = fopen(file_name_exact, "w");
  if (file_exact == NULL) error("Could not create file '%s'.", file_name_exact);
  fprintf(file_exact, "# Sink accuracy test - EXACT DENSITIES\n");
  fprintf(file_exact, "# N= %d\n", SWIFT_SINK_DENSITY_CHECKS);
  fprintf(file_exact, "# periodic= %d\n", s->periodic);
  fprintf(file_exact, "# N_ngb_target= %f +/- %f\n", N_ngb_target,
          e->sink_properties->delta_neighbours);
  fprintf(file_exact, "# Git Branch: %s\n", git_branch());
  fprintf(file_exact, "# Git Revision: %s\n", git_revision());
  fprintf(file_exact, "# %16s %16s %16s %16s %16s %7s %7s %16s %16s %16s\n",
          "id", "pos[0]", "pos[1]", "pos[2]", "h", "Nd", "Nf", "rho_exact",
          "n_rho_exact", "N_ngb");

  int wrong_rho = 0;
  int wrong_n_ngb = 0;
  int counter = 0;

  /* Output particle SWIFT densities */
  for (size_t i = 0; i < nr_sinks; ++i) {

    const struct sink *si = &sinks[i];
    const long long id = si->id;
    const int found_inhibited = si->inhibited_check_exact;

    const double N_ngb = (4. / 3.) * M_PI * kernel_gamma * kernel_gamma *
                         kernel_gamma * si->h * si->h * si->h *
                         si->n_check_exact;

    if (id % SWIFT_SINK_DENSITY_CHECKS == 0 && sink_is_starting(si, e)) {

      counter++;

      fprintf(
          file_exact,
          "%18lld %16.8e %16.8e %16.8e %16.8e %7d %7d %16.8e %16.8e %16.8e\n",
          id, si->x[0], si->x[1], si->x[2], si->h, si->N_check_density_exact, 0,
          si->rho_check_exact, si->n_check_exact, N_ngb);

      /* Check that we did not go above the threshold.
       * Note that we ignore particles that saw an inhibted particle as a
       * neighbour as we don't know whether that neighbour became inhibited in
       * that step or not. */
      if (!found_inhibited &&
          si->N_check_density_exact != si->N_check_density &&
          (fabsf(si->rho_check / si->rho_check_exact - 1.f) > rel_tol ||
           fabsf(si->rho_check_exact / si->rho_check - 1.f) > rel_tol)) {
        message("RHO: id=%lld swift=%e exact=%e N_swift=%d N_true=%d", id,
                si->rho_check, si->rho_check_exact, si->N_check_density,
                si->N_check_density_exact);
        wrong_rho++;
      }

      if (!found_inhibited && (N_ngb > N_ngb_max || N_ngb < N_ngb_min)) {

        message("N_NGB: id=%lld exact=%f N_true=%d N_swift=%d", id, N_ngb,
                si->N_check_density_exact, si->N_check_density);

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
        "sink particles! (out of %d particles)",
        wrong_rho, counter);
  else
    message("Verified %d sink particles", counter);

  /* if (wrong_n_ngb) */
  /*   error( */
  /*       "N_ngb difference larger than the allowed tolerance for %d " */
  /*       "star particles! (out of %d particles)", */
  /*       wrong_n_ngb, counter); */
  /* else */
  /*   message("Verified %d star particles", counter); */

  if (e->verbose)
    message("Writting brute-force density files took %.3f %s. ",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("Sink checking function called without the corresponding flag.");
#endif
}

/**
 * @brief Mapper function for the exact sink formation count checks.
 *
 * @param map_data Pointer to gas particles array.
 * @param nr_parts Number of gas particles.
 * @param extra_data Pointers to space and engine.
 */
void sink_exact_formation_count_compute_mapper(void *map_data, int nr_parts,
                                               void *extra_data) {
#ifdef SWIFT_DEBUG_CHECKS_HYDRO_SINKS_FORMATION_COUNT_CHECKS

  struct part *restrict parts = (struct part *)map_data;
  struct exact_density_data *data = (struct exact_density_data *)extra_data;
  const struct space *s = data->s;
  const struct engine *e = data->e;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  /* Use float precision to match the formation loop, which computes
   * r_cut2 = r_cut * r_cut and dx/r2 all in float.  Using double here
   * would give a more precise r_cut2 and could exclude particles that
   * the formation loop (at float precision) correctly included. */
  const float r_cut = (float)e->sink_properties->cut_off_radius;
  const float r_cut2 = r_cut * r_cut;
  int counter = 0;

  for (int i = 0; i < nr_parts; ++i) {

    struct part *pi = &parts[i];
    const long long id = pi->id;

    /* Is the particle part of the subset to be tested? */
    if (id % SWIFT_DEBUG_CHECKS_HYDRO_SINKS_FORMATION_COUNT_CHECKS == 0 &&
        part_is_starting(pi, e)) {

      /* Only process particles for which sink_init_part ran this step.
       * sink_init_part sets N_check_formation_exact = -2 as a sentinel.
       * If that sentinel is absent the particle's cell was not in the drift
       * task graph (e.g. it was woken by the limiter after the graph was
       * built) so N_check_formation is stale — skip it. */
      if (pi->sink_data.N_check_formation_exact != -2) {
        counter++;
        continue;
      }

      /* Get position of gas particle i (in double for nearest() accuracy) */
      const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};

      /* Brute-force count: loop over all gas particles */
      int N_formation_exact = 0;

      for (int j = 0; j < (int)s->nr_parts; ++j) {

        const struct part *pj = &s->parts[j];

        /* Skip self-interaction: compare pointers, not chunk-local index,
         * because map_data may start mid-array when threadpool chunks work. */
        if (pj == pi) continue;

        /* Compute pairwise distance.  Subtract in double (like the self-cell
         * formation loop: `(float)(pi->x[k] - pj->x[k])`), apply periodic BC
         * in double for accuracy, then cast to float before squaring so that
         * the r2 < r_cut2 comparison uses the same float precision as the
         * formation loop.  This prevents spurious mismatches from particles
         * sitting exactly on the r_cut boundary. */
        double ddx = pj->x[0] - pix[0];
        double ddy = pj->x[1] - pix[1];
        double ddz = pj->x[2] - pix[2];

        /* Apply periodic BC in double */
        if (periodic) {
          ddx = nearest(ddx, dim[0]);
          ddy = nearest(ddy, dim[1]);
          ddz = nearest(ddz, dim[2]);
        }

        const float dx = (float)ddx;
        const float dy = (float)ddy;
        const float dz = (float)ddz;
        const float r2 = dx * dx + dy * dy + dz * dz;

        /* Count if within fixed aperture.
         * Particles swallowed *during this step* (by the swallow task,
         * which runs after the formation loop) must be counted: the
         * formation loop saw them as live.  Such particles were drifted
         * to ti_current at the start of the step (before being inhibited),
         * so pj->ti_drift == e->ti_current distinguishes them from
         * particles that were already inhibited before this step
         * (pj->ti_drift < e->ti_current, drift was skipped for them).
         * In non-debug builds ti_drift is unavailable; fall back to
         * excluding all inhibited particles (may produce false alarms
         * when swallowing occurs between rebuilds). */
        if (r2 < r_cut2) {
#ifdef SWIFT_DEBUG_CHECKS
          const int skip =
              part_is_inhibited(pj, e) && (pj->ti_drift != e->ti_current);
#else
          const int skip = part_is_inhibited(pj, e);
#endif
          if (!skip) N_formation_exact++;
        }
      }

      /* Store the exact count */
      pi->sink_data.N_check_formation_exact = N_formation_exact;
      counter++;
    }
  }
  atomic_add(&data->counter_global, counter);

#else
  error(
      "Formation count checking function called without the corresponding "
      "flag.");
#endif
}

/**
 * @brief Compute exact gas-gas neighbor counts for a selection of gas particles
 * by running a brute-force loop over all particles in the simulation.
 *
 * @param s The space.
 * @param e The engine.
 */
void sink_exact_formation_count_compute(struct space *s,
                                        const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS_HYDRO_SINKS_FORMATION_COUNT_CHECKS

  const ticks tic = getticks();

  struct exact_density_data data;
  data.e = e;
  data.s = s;
  data.counter_global = 0;

  threadpool_map(&s->e->threadpool, sink_exact_formation_count_compute_mapper,
                 s->parts, s->nr_parts, sizeof(struct part), 0, &data);

  if (e->verbose)
    message(
        "Computed exact formation neighbor counts for %d gas particles "
        "(took %.3f %s).",
        data.counter_global, clocks_from_ticks(getticks() - tic),
        clocks_getunit());

#else
  error(
      "Formation count checking function called without the corresponding "
      "flag.");
#endif
}

/**
 * @brief Check gas particles' gas-gas neighbor counts (formation loop) against
 * values obtained via brute-force summation.
 *
 * @param s The space.
 * @param e The engine.
 */
void sink_exact_formation_count_check(struct space *s, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS_HYDRO_SINKS_FORMATION_COUNT_CHECKS

  const ticks tic = getticks();

  const struct part *parts = s->parts;
  const size_t nr_parts = s->nr_parts;

  int wrong_count = 0;
  int counter = 0;

  for (size_t i = 0; i < nr_parts; ++i) {

    const struct part *pi = &parts[i];
    const long long id = pi->id;

    if (id % SWIFT_DEBUG_CHECKS_HYDRO_SINKS_FORMATION_COUNT_CHECKS == 0 &&
        part_is_starting(pi, e)) {

      counter++;

      const int N_formation = pi->sink_data.N_check_formation;
      const int N_formation_exact = pi->sink_data.N_check_formation_exact;

      /* Skip particles with a negative N_check_formation_exact:
       *   -2  brute-force skipped this particle (sink_init_part didn't run)
       *   -1  limiter-woken sentinel (timestep_limit_part)
       * In both cases N_check_formation is stale — skip the comparison. */
      if (N_formation_exact < 0) continue;

      if (N_formation != N_formation_exact) {
        message("FORMATION_COUNT: id=%lld optimised=%d exact=%d", id,
                N_formation, N_formation_exact);
        wrong_count++;
      }
    }
  }

  if (wrong_count)
    error(
        "Gas-gas formation neighbor count mismatch for %d particles "
        "(out of %d checked).",
        wrong_count, counter);
  else if (counter > 0)
    message("Verified formation neighbor counts for %d gas particles.",
            counter);

  if (e->verbose)
    message("Formation count checks took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error(
      "Formation count checking function called without the corresponding "
      "flag.");
#endif
}
