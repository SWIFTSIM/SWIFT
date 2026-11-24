/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#include <config.h>

/* Local headers. */
#include "swift.h"

/* System headers. */
#include <float.h>
#include <stdio.h>
#include <string.h>

/**
 * @file testChemistryPositivityLimiter.c
 * @brief Check that the positivity limiter of the FVPM chemistry keeps a
 * nearly empty donor non-negative and conserves the metal mass.
 *
 * One donor with little metal mass sits among 12 metal-free neighbours. The
 * test runs the prep pass, the theta computation and the apply pass through
 * the real interaction hooks.
 */

#if defined(CHEMISTRY_GEAR_FVPM_DIFFUSION) || \
    defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)

#define N_NEIGHBOURS 12

/* A long flux interval (1e3 at a time base of 1e-5) makes the donor's
   outflow exceed its mass. */
#define TI_CURRENT 100000000LL
#define MAX_ACTIVE_BIN 7

/**
 * @brief Is x finite?
 *
 * -ffast-math lets the compiler assume that x is finite and fold the test,
 * even a test of its bits. A load from a volatile copy carries no such
 * assumption.
 *
 * @param x The value.
 */
static int is_finite_double(const double x) {
  const volatile double stored = x;
  const double loaded = stored;
  uint64_t bits;
  memcpy(&bits, &loaded, sizeof(bits));
  return ((bits >> 52) & 0x7ff) != 0x7ff;
}

static void init_part(struct part *p, const long long id, const double x[3],
                      const double metal_mass_fraction) {
  bzero(p, sizeof(struct part));
  p->id = id;
  for (int k = 0; k < 3; k++) p->x[k] = x[k];
  p->h = 1.f;
  hydro_set_mass(p, 1.f);
  p->rho = 1.f;
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH)
  p->P = 1.f;
#else
  p->u = 1.f;
#endif
  p->time_bin = 1;

  p->geometry.volume = 1.f;
  p->geometry.matrix_E[0][0] = 1.f;
  p->geometry.matrix_E[1][1] = 1.f;
  p->geometry.matrix_E[2][2] = 1.f;
  p->geometry.wcorr = 1.f;
  p->geometry.condition_number = 1.f;

  struct chemistry_part_data *ch = &p->chemistry_data;
  ch->kappa = 1.0;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  ch->tau = 1.0;
#endif
  ch->limiter.maxr = 1.f;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    ch->limiter.rhoZ[m][0] = -1e10;
    ch->limiter.rhoZ[m][1] = 1e10;
    ch->limiter.Z[m][0] = -1e10;
    ch->limiter.Z[m][1] = 1e10;
    ch->fct_theta[m] = 1.0;
    /* A different mass per element, so that an element index mixed up in the
       limiter changes the result. */
    ch->metal_mass[m] = hydro_get_mass(p) * metal_mass_fraction * (1. + m);
  }

  /* The particles last applied their fluxes at ti = 0. */
  ch->flux.ti_flux_last = 0;
}

/**
 * @brief Call the flux exchange of the donor with each neighbour.
 */
static void exchange_all(struct part *donor, struct part *ngb, const int prep,
                         const struct cosmology *cosmo,
                         const struct chemistry_global_data *cd) {
  for (int n = 0; n < N_NEIGHBOURS; n++) {
    float dx[3];
    for (int k = 0; k < 3; k++) dx[k] = donor->x[k] - ngb[n].x[k];
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    runner_iact_chemistry_flux_exchange(
        r2, dx, donor->h, ngb[n].h, donor, &ngb[n],
        /*both_updatable_here=*/1, /*local_first=*/1, /*local_second=*/1,
        /*a=*/1.f, /*H=*/0.f, /*time_base=*/1e-5, TI_CURRENT, cosmo,
        /*with_cosmology=*/0, cd,
        prep ? chemistry_flux_pass_prep : chemistry_flux_pass_apply,
        /*max_active_bin=*/MAX_ACTIVE_BIN);
  }
}

/**
 * @brief Run prep, theta and apply on the fixture.
 *
 * @param positive (return) Is the donor mass non-negative for all elements?
 * @param theta_ok (return) Does the theta of the donor match the formula of
 * the limiter for all elements?
 * @return Is the metal mass conserved for all elements?
 */
static int run_fixture(const struct chemistry_global_data *cd,
                       const struct cosmology *cosmo, int *positive,
                       int *theta_ok) {

  /* Icosahedron vertices at distance 0.5 from the donor. */
  const double phi = 0.5 * (1. + sqrt(5.));
  const double v[N_NEIGHBOURS][3] = {
      {0, 1, phi}, {0, -1, phi}, {0, 1, -phi}, {0, -1, -phi},
      {1, phi, 0}, {-1, phi, 0}, {1, -phi, 0}, {-1, -phi, 0},
      {phi, 0, 1}, {-phi, 0, 1}, {phi, 0, -1}, {-phi, 0, -1}};
  const double scale = 0.5 / sqrt(1. + phi * phi);

  struct part donor;
  struct part ngb[N_NEIGHBOURS];
  const double origin[3] = {0., 0., 0.};
  init_part(&donor, 1, origin, 1e-6);
  for (int n = 0; n < N_NEIGHBOURS; n++) {
    const double x[3] = {scale * v[n][0], scale * v[n][1], scale * v[n][2]};
    init_part(&ngb[n], n + 2, x, 0.);

    /* Metal density falling away from the donor drives an outward flux. The
       neighbours' metal mass and gradient differ per element, so the outflow
       per unit of donor mass, and thus theta, differs per element. */
    for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
      ngb[n].chemistry_data.metal_mass[m] =
          donor.chemistry_data.metal_mass[m] * (0.4 + 0.05 * m);
      for (int k = 0; k < 3; k++)
        ngb[n].chemistry_data.gradients.rhoZ[m][k] =
            -4e-6 * x[k] * (1. + 2. * m);
    }
  }

  exchange_all(&donor, ngb, /*prep=*/1, cosmo, cd);

  double expected_theta[GEAR_CHEMISTRY_ELEMENT_COUNT];
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    const double available = donor.chemistry_data.metal_mass[m];
    const double sum_out = donor.chemistry_data.fct_sum_out[m];
    message("Element %d: donor metal mass %e, prospective outflow %e", m,
            available, sum_out);
    if (!(sum_out > available))
      error("The fixture does not over-drain the donor: sum_out=%e <= %e",
            sum_out, available);

    /* The donor keeps 1e-6 of its mass, as in chemistry_end_flux_prep(). */
    expected_theta[m] = available * (1.0 - 1e-6) / sum_out;
  }

  chemistry_end_flux_prep(&donor, cd);
  for (int n = 0; n < N_NEIGHBOURS; n++) chemistry_end_flux_prep(&ngb[n], cd);

  /* Three roundings (the factor, the product, the quotient) of at most
     DBL_EPSILON / 2 each, so 4 * DBL_EPSILON bounds the difference. */
  *theta_ok = 1;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    const double theta = donor.chemistry_data.fct_theta[m];
    const int ok =
        is_finite_double(theta) &&
        fabs(theta - expected_theta[m]) <= 4. * DBL_EPSILON * expected_theta[m];
    if (!ok)
      message("element %d: theta %.17e, expected %.17e", m, theta,
              expected_theta[m]);
    *theta_ok &= ok;
  }

  exchange_all(&donor, ngb, /*prep=*/0, cosmo, cd);

  int conserved_all = 1;
  *positive = 1;
  for (int m = 0; m < GEAR_CHEMISTRY_ELEMENT_COUNT; m++) {
    const double donor_final = donor.chemistry_data.metal_mass[m] +
                               donor.chemistry_data.flux.metal_mass[m];
    double net = donor.chemistry_data.flux.metal_mass[m];
    double traffic = fabs(net);

    /* Test the stored values, not the sums: under -ffast-math a sum can
       lose a NaN. */
    int finite = is_finite_double(donor.chemistry_data.metal_mass[m]) &&
                 is_finite_double(donor.chemistry_data.flux.metal_mass[m]);
    for (int n = 0; n < N_NEIGHBOURS; n++) {
      net += ngb[n].chemistry_data.flux.metal_mass[m];
      traffic += fabs(ngb[n].chemistry_data.flux.metal_mass[m]);
      finite &= is_finite_double(ngb[n].chemistry_data.flux.metal_mass[m]);
    }

    /* Each pair adds +F and -F, so only the summation rounding remains. */
    const double conservation_bar = 2. * (N_NEIGHBOURS + 1) * DBL_EPSILON;
    const int positive_m = finite && donor_final >= 0.;
    const int conserved =
        finite && fabs(net) <= conservation_bar * fmax(traffic, DBL_MIN);
    if (m == 0 || !conserved)
      message(
          "element %d: donor final mass %e (theta %e), net flux %e, "
          "traffic %e",
          m, donor_final, donor.chemistry_data.fct_theta[m], net, traffic);
    *positive &= positive_m;
    conserved_all &= conserved;
  }
  return conserved_all;
}

int main(int argc, char *argv[]) {

  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  struct chemistry_global_data cd;
  bzero(&cd, sizeof(struct chemistry_global_data));
  cd.diffusion_coefficient = 1.f;
  cd.diffusion_mode = isotropic_constant;
#if defined(CHEMISTRY_GEAR_FVPM_HYPERBOLIC_DIFFUSION)
  cd.tau = 1.;
  cd.relaxation_time_mode = constant_mode;
  cd.hyperbolic_limiter_scope = limiter_all_components;
#endif
  cd.hll_riemann_solver_psi = 0.1f;
  cd.C_CFL_chemistry = 0.5f;
  cd.flux_limiter_noise_gate = 1e-15;
  cd.flux_limiter_safety = 0.5;
  cd.flux_limiter_sink_stability = 0.25;
  cd.flux_limiter_startup = 0.1;

  int positive_on, positive_off, theta_on, theta_off;
  cd.use_positivity_limiter = 1;
  const int conserved_on = run_fixture(&cd, &cosmo, &positive_on, &theta_on);
  message("Limiter on:  positive=%d conserved=%d theta=%d", positive_on,
          conserved_on, theta_on);

  /* With the limiter off, the donor must go negative and the mass stay
     conserved. */
  cd.use_positivity_limiter = 0;
  const int conserved_off = run_fixture(&cd, &cosmo, &positive_off, &theta_off);
  message("Limiter off: positive=%d conserved=%d", positive_off, conserved_off);

  /* The donor stays positive with theta = 0 too, so the theta is required. */
  if (!positive_on || !conserved_on || !theta_on || positive_off ||
      !conserved_off)
    error("testChemistryPositivityLimiter: FAILED");
  message("testChemistryPositivityLimiter: PASSED");
  return 0;
}

#else

int main(int argc, char *argv[]) { return 0; }

#endif
