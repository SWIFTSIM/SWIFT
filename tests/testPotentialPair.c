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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Local headers. */
#include "runner_doiact_grav.h"
#include "swift.h"

const int num_tests = 100;
const double eps = 0.1;

/**
 * @brief Check that a and b are consistent (up to some relative error)
 *
 * @param a First value
 * @param b Second value
 * @param s String used to identify this check in messages
 * @param rel_tol Maximal relative error
 * @param limit Minimal value to consider in the tests
 */
void check_value_backend(double a, double b, const char *s, double rel_tol,
                         double limit) {
  if (fabs(a - b) / fabs(a + b) > rel_tol && fabs(a - b) > limit)
    error("Values are inconsistent: SWIFT:%12.15e true:%12.15e (%s)!", a, b, s);
}

void check_value(double a, double b, const char *s) {
  check_value_backend(a, b, s, 2e-6, 1e-6);
}

/* Definitions of the potential and force that match
   exactly the theory document */
double S(double x) { return exp(x) / (1. + exp(x)); }

double S_prime(double x) { return exp(x) / ((1. + exp(x)) * (1. + exp(x))); }

double potential(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double pot;
  if (u > 1.)
    pot = -mass / r;
  else
    pot = -mass *
          (-3. * u * u * u * u * u * u * u + 15. * u * u * u * u * u * u -
           28. * u * u * u * u * u + 21. * u * u * u * u - 7. * u * u + 3.) /
          H;

  return pot * (2. - 2. * S(2. * x));
}

double acceleration(double mass, double r, double H, double rlr) {

  const double u = r / H;
  const double x = r / rlr;
  double acc;
  if (u > 1.)
    acc = -mass / (r * r * r);
  else
    acc = -mass *
          (21. * u * u * u * u * u - 90. * u * u * u * u + 140. * u * u * u -
           84. * u * u + 14.) /
          (H * H * H);

  return r * acc * (4. * x * S_prime(2 * x) - 2. * S(2. * x) + 2.);
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Initialise a few things to get us going */

  /* Non-truncated forces first */
  double rlr = FLT_MAX;

  struct engine e;
  e.max_active_bin = num_time_bins;
  e.time = 0.1f;
  e.ti_current = 8;
  e.time_base = 1e-10;
  e.nodeID = 0;

  struct space s;
  s.periodic = 0;
  e.s = &s;

  struct pm_mesh mesh;
  mesh.periodic = 0;
  mesh.dim[0] = 10.;
  mesh.dim[1] = 10.;
  mesh.dim[2] = 10.;
  mesh.r_s = rlr;
  mesh.r_s_inv = 1. / rlr;
  mesh.r_cut_min = 0.;
  mesh.r_cut_max = FLT_MAX;
  e.mesh = &mesh;

  struct gravity_props props;
  props.theta_crit = 0.;
  props.epsilon_DM_cur = eps;
  props.epsilon_baryon_cur = eps;
  e.gravity_properties = &props;

  struct runner r;
  bzero(&r, sizeof(struct runner));
  r.e = &e;

  /* Init the cache for gravity interaction */
  gravity_cache_init(&r.ci_gravity_cache, num_tests);
  gravity_cache_init(&r.cj_gravity_cache, num_tests);

  /* Let's create one cell with a massive particle and a bunch of test particles
   */
  struct cell ci, cj;
  bzero(&ci, sizeof(struct cell));
  bzero(&cj, sizeof(struct cell));

  ci.nodeID = 0;
  ci.width[0] = 1.;
  ci.width[1] = 1.;
  ci.width[2] = 1.;
  ci.loc[0] = 0.;
  ci.loc[1] = 0.;
  ci.loc[2] = 0.;
  ci.grav.count = 1;
  ci.grav.ti_old_part = 8;
  ci.grav.ti_old_multipole = 8;
  ci.grav.ti_end_min = 8;
  ci.grav.ti_end_max = 8;

  cj.nodeID = 0;
  cj.width[0] = 1.;
  cj.width[1] = 1.;
  cj.width[2] = 1.;
  cj.loc[0] = 1.;
  cj.loc[1] = 0.;
  cj.loc[2] = 0.;
  cj.grav.count = num_tests;
  cj.grav.ti_old_part = 8;
  cj.grav.ti_old_multipole = 8;
  cj.grav.ti_end_min = 8;
  cj.grav.ti_end_max = 8;

  /* Allocate multipoles */
  ci.grav.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  cj.grav.multipole =
      (struct gravity_tensors *)malloc(sizeof(struct gravity_tensors));
  bzero(ci.grav.multipole, sizeof(struct gravity_tensors));
  bzero(cj.grav.multipole, sizeof(struct gravity_tensors));

  /* Set the multipoles */
  ci.grav.multipole->r_max = 0.1;
  cj.grav.multipole->r_max = 0.1;

  /* Allocate the particles */
  if (posix_memalign((void **)&ci.grav.parts, gpart_align,
                     ci.grav.count * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(ci.grav.parts, ci.grav.count * sizeof(struct gpart));

  if (posix_memalign((void **)&cj.grav.parts, gpart_align,
                     cj.grav.count * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(cj.grav.parts, cj.grav.count * sizeof(struct gpart));

  /* Create the mass-less test particles */
  for (int n = 0; n < num_tests; ++n) {

    struct gpart *gp = &cj.grav.parts[n];

    gp->x[0] = 1. + (n + 1) / ((double)num_tests);
    gp->x[1] = 0.5;
    gp->x[2] = 0.5;
    gp->mass = 0.;
    gp->time_bin = 1;
    gp->type = swift_type_dark_matter;
    gp->id_or_neg_offset = n + 1;
#ifdef SWIFT_DEBUG_CHECKS
    gp->ti_drift = 8;
    gp->initialised = 1;
#endif
  }

  /***********************************************/
  /* Let's start by testing the P-P interactions */
  /***********************************************/

  /* Create the massive particle */
  ci.grav.parts[0].x[0] = 0.;
  ci.grav.parts[0].x[1] = 0.5;
  ci.grav.parts[0].x[2] = 0.5;
  ci.grav.parts[0].mass = 1.;
  ci.grav.parts[0].time_bin = 1;
  ci.grav.parts[0].type = swift_type_dark_matter;
  ci.grav.parts[0].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
  ci.grav.parts[0].ti_drift = 8;
  ci.grav.parts[0].initialised = 1;
#endif

  /* Now compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj, 1, 1);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.grav.parts[n];
    const struct gpart *gp2 = &ci.grav.parts[0];
    const double epsilon = gravity_get_softening(gp, &props);

#if defined(POTENTIAL_GRAVITY)
    double pot_true =
        potential(ci.grav.parts[0].mass, gp->x[0] - gp2->x[0], epsilon, rlr);
    check_value(gp->potential, pot_true, "potential");
#endif

    double acc_true =
        acceleration(ci.grav.parts[0].mass, gp->x[0] - gp2->x[0], epsilon, rlr);

    /* message("x=%e f=%e f_true=%e pot=%e pot_true=%e", gp->x[0] - gp2->x[0],
       gp->a_grav[0], acc_true, gp->potential, pot_true); */

    check_value(gp->a_grav[0], acc_true, "acceleration");
  }

  message("\n\t\t P-P interactions all good\n");

  /* Reset the accelerations */
  for (int n = 0; n < num_tests; ++n) gravity_init_gpart(&cj.grav.parts[n]);

  /**********************************/
  /* Test the basic PM interactions */
  /**********************************/

  /* Set an opening angle that allows P-M interactions */
  props.theta_crit = 1.;

  ci.grav.parts[0].mass = 0.;
  ci.grav.multipole->CoM[0] = 0.;
  ci.grav.multipole->CoM[1] = 0.5;
  ci.grav.multipole->CoM[2] = 0.5;

  bzero(&ci.grav.multipole->m_pole, sizeof(struct multipole));
  bzero(&cj.grav.multipole->m_pole, sizeof(struct multipole));
  ci.grav.multipole->m_pole.M_000 = 1.;

  /* Now compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj, 1, 1);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.grav.parts[n];
    const struct gravity_tensors *mpole = ci.grav.multipole;
    const double epsilon = gravity_get_softening(gp, &props);

#if defined(POTENTIAL_GRAVITY)
    double pot_true =
        potential(mpole->m_pole.M_000, gp->x[0] - mpole->CoM[0], epsilon, rlr);
    check_value(gp->potential, pot_true, "potential");
#endif

    double acc_true = acceleration(mpole->m_pole.M_000,
                                   gp->x[0] - mpole->CoM[0], epsilon, rlr);
    check_value(gp->a_grav[0], acc_true, "acceleration");

    /* message("x=%e f=%e f_true=%e pot=%e pot_true=%e", gp->x[0] -
     * mpole->CoM[0], gp->a_grav[0], acc_true, gp->potential, pot_true); */
  }

  message("\n\t\t basic P-M interactions all good\n");

#ifndef GADGET2_LONG_RANGE_CORRECTION

  /* Reset the accelerations */
  for (int n = 0; n < num_tests; ++n) gravity_init_gpart(&cj.grav.parts[n]);

  /***************************************/
  /* Test the truncated PM interactions  */
  /***************************************/
  rlr = 2.;
  mesh.r_s = rlr;
  mesh.r_s_inv = 1. / rlr;
  mesh.periodic = 1;
  s.periodic = 1;
  props.epsilon_cur = FLT_MIN; /* No softening */

  /* Now compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj, 1, 1);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.grav.parts[n];
    const struct gravity_tensors *mpole = ci.grav.multipole;
    const double epsilon = gravity_get_softening(gp, &props);

#if defined(POTENTIAL_GRAVITY)
    double pot_true =
        potential(mpole->m_pole.M_000, gp->x[0] - mpole->CoM[0], epsilon, rlr);
    check_value(gp->potential, pot_true, "potential");
#endif

    double acc_true = acceleration(mpole->m_pole.M_000,
                                   gp->x[0] - mpole->CoM[0], epsilon, rlr);
    check_value(gp->a_grav[0], acc_true, "acceleration");

    /* message("x=%e f=%e f_true=%e pot=%e pot_true=%e", gp->x[0] -
     * mpole->CoM[0], */
    /*         gp->a_grav[0], acc_true, gp->potential, pot_true); */
  }

  message("\n\t\t truncated P-M interactions all good\n");

#endif

  /************************************************/
  /* Test the high-order periodic PM interactions */
  /************************************************/

  /* Reset the accelerations */
  for (int n = 0; n < num_tests; ++n) gravity_init_gpart(&cj.grav.parts[n]);

#if SELF_GRAVITY_MULTIPOLE_ORDER >= 3

  /* Let's make ci more interesting */
  free(ci.grav.parts);
  ci.grav.count = 8;
  if (posix_memalign((void **)&ci.grav.parts, gpart_align,
                     ci.grav.count * sizeof(struct gpart)) != 0)
    error("Error allocating gparts for cell ci");
  bzero(ci.grav.parts, ci.grav.count * sizeof(struct gpart));

  /* Place particles on a simple cube of side-length 0.2 */
  for (int n = 0; n < 8; ++n) {
    if (n & 1)
      ci.grav.parts[n].x[0] = 0.0 - 0.1;
    else
      ci.grav.parts[n].x[0] = 0.0 + 0.1;

    if (n & 2)
      ci.grav.parts[n].x[1] = 0.5 - 0.1;
    else
      ci.grav.parts[n].x[1] = 0.5 + 0.1;

    if (n & 2)
      ci.grav.parts[n].x[2] = 0.5 - 0.1;
    else
      ci.grav.parts[n].x[2] = 0.5 + 0.1;

    ci.grav.parts[n].mass = 1. / 8.;

    ci.grav.parts[n].time_bin = 1;
    ci.grav.parts[n].type = swift_type_dark_matter;
    ci.grav.parts[n].id_or_neg_offset = 1;
#ifdef SWIFT_DEBUG_CHECKS
    ci.grav.parts[n].ti_drift = 8;
    ci.grav.parts[n].initialised = 1;
#endif
  }

  /* Now let's make a multipole out of it. */
  gravity_reset(ci.grav.multipole);
  gravity_P2M(ci.grav.multipole, ci.grav.parts, ci.grav.count, &props);

  gravity_multipole_print(&ci.grav.multipole->m_pole);

  /* Compute the forces */
  runner_dopair_grav_pp(&r, &ci, &cj, 1, 1);

  /* Verify everything */
  for (int n = 0; n < num_tests; ++n) {
    const struct gpart *gp = &cj.grav.parts[n];

#if defined(POTENTIAL_GRAVITY)
    double pot_true = 0;
#endif
    double acc_true[3] = {0., 0., 0.};

    for (int i = 0; i < 8; ++i) {
      const struct gpart *gp2 = &ci.grav.parts[i];
      const double epsilon = gravity_get_softening(gp, &props);

      const double dx[3] = {gp2->x[0] - gp->x[0], gp2->x[1] - gp->x[1],
                            gp2->x[2] - gp->x[2]};
      const double d = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

#if defined(POTENTIAL_GRAVITY)
      pot_true += potential(gp2->mass, d, epsilon, rlr);
#endif

      acc_true[0] -= acceleration(gp2->mass, d, epsilon, rlr) * dx[0] / d;
      acc_true[1] -= acceleration(gp2->mass, d, epsilon, rlr) * dx[1] / d;
      acc_true[2] -= acceleration(gp2->mass, d, epsilon, rlr) * dx[2] / d;
    }

#if defined(POTENTIAL_GRAVITY)
    check_value_backend(gp->potential, pot_true, "potential", 1e-2, 1e-6);
#endif
    check_value_backend(gp->a_grav[0], acc_true[0], "acceleration", 1e-2, 1e-6);

    /* const struct gravity_tensors *mpole = ci.grav.multipole; */
    /* message("x=%e f=%e f_true=%e pot=%e pot_true=%e %e %e", */
    /*         gp->x[0] - mpole->CoM[0], gp->a_grav[0], acc_true[0],
     * gp->potential, */
    /*         pot_true, acc_true[1], acc_true[2]); */
  }

  message("\n\t\t high-order P-M interactions all good\n");

#endif

  free(ci.grav.multipole);
  free(cj.grav.multipole);
  free(ci.grav.parts);
  free(cj.grav.parts);

  /* Clean up the caches */
  gravity_cache_clean(&r.ci_gravity_cache);
  gravity_cache_clean(&r.cj_gravity_cache);

  return 0;
}
