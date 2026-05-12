/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2025 Katy Proctor (katy.proctor@fysik.su.se).
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

/* Local includes. */
#include "swift.h"
#include "timestep_limiter_iact.h"

/* System includes. */
#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void print_bytes(void *p, size_t len) {
  printf("(");
  for (size_t i = 0; i < len; ++i) {
    printf("%02x", ((unsigned char *)p)[i]);
    if (i % 4 == 3) printf("|");
  }
  printf(")\n");
}

void test(const char *param_filename) {

  /* Start with some values for the cosmological parameters and unused variables
   */
  const float a = (float)random_uniform(0.8, 1.);
  const float H = 1.f;
  const int with_cosmology = 0;
  const float ti_current = 0;
  const float time_base = 1;

  /* parse parameters */
  struct swift_params param_file;
  parser_read_file(param_filename, &param_file);

  /* Default unit system */
  struct unit_system us;
  units_init_cgs(&us);

  /* Default physical constants */
  struct phys_const prog_const;
  phys_const_init(&us, &param_file, &prog_const);

  struct cosmology cosmo;
  cosmology_init_no_cosmo(&cosmo);

  struct hydro_props hydro_properties = NULL;

  struct sidm_props sidm_props;
  sidm_props_init(&sidm_props, &prog_const, &us, &param_file, &hydro_properties,
                  &cosmo);
  sidm_props.eta_neighbours = 1.2348f;
  sidm_props.sigma_over_m = 0.f;
  sidm_props.h_tolerance = 1e0;
  sidm_props.h_max = FLT_MAX;
  sidm_props.h_min = 0.f;
  sidm_props.h_min_ratio = 0.f;
  sidm_props.max_smoothing_iterations = 10;

  /* Create two random particles (don't do this at home !) */
  struct sipart pi, pj;
  for (size_t i = 0; i < sizeof(struct sipart) / sizeof(float); ++i) {
    *(((float *)&pi) + i) = (float)random_uniform(0., 2.);
    *(((float *)&pj) + i) = (float)random_uniform(0., 2.);
  }

  /* Make the particle smoothing length and position reasonable */
  for (size_t i = 0; i < 3; ++i) pi.x[i] = random_uniform(-1., 1.);
  for (size_t i = 0; i < 3; ++i) pj.x[i] = random_uniform(-1., 1.);
  pi.h = 2.f;
  pj.h = 2.f;
  pi.id = 1ll;
  pj.id = 2ll;
  pi.time_bin = 1;
  pj.time_bin = 1;
  pi.v[0] = random_uniform(-10.0f, 10.0f);
  pi.v[1] = random_uniform(-10.0f, 10.0f);
  pi.v[2] = random_uniform(-10.0f, 10.0f);
  pj.v[0] = random_uniform(-10.0f, 10.0f);
  pj.v[1] = random_uniform(-10.0f, 10.0f);
  pj.v[2] = random_uniform(-10.0f, 10.0f);

  /* Make some copies */
  struct sipart pi2, pj2;
  memcpy(&pi2, &pi, sizeof(struct sipart));
  memcpy(&pj2, &pj, sizeof(struct sipart));

  int i_not_ok = memcmp(&pi, &pi2, sizeof(struct sipart));
  int j_not_ok = memcmp(&pj, &pj2, sizeof(struct sipart));

  if (i_not_ok) error("Particles 'sipi' do not match after copy");
  if (j_not_ok) error("Particles 'sipj' do not match after copy");

  /* Compute distance vector */
  float dx[3];
  dx[0] = pi.x[0] - pj.x[0];
  dx[1] = pi.x[1] - pj.x[1];
  dx[2] = pi.x[2] - pj.x[2];
  float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* --- Test the density loop --- */

  /* Call the symmetric version */
  runner_iact_sidm_density(r2, dx, pi.h, pj.h, &pi, &pj, a, H, with_cosmology,
                           &cosmo, &sidm_props, ti_current, time_base);

  /* Call the non-symmetric version */
  runner_iact_nonsym_sidm_density(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H,
                                  with_cosmology, &cosmo, &sidm_props,
                                  ti_current, time_base);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_sidm_density(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H,
                                  with_cosmology, &cosmo, &sidm_props,
                                  ti_current, time_base);

  /* Check that the particles are the same */
  i_not_ok = memcmp(&pi, &pi2, sizeof(struct sipart));
  j_not_ok = memcmp(&pj, &pj2, sizeof(struct sipart));

  if (i_not_ok) {
    printSIDMParticle_single(&pi);
    printSIDMParticle_single(&pi2);
    print_bytes(&pi, sizeof(struct sipart));
    print_bytes(&pi2, sizeof(struct sipart));
    error("Particles 'sipi' do not match after density (byte = %d)", i_not_ok);
  }
  if (j_not_ok) {
    printSIDMParticle_single(&pj);
    printSIDMParticle_single(&pj2);
    print_bytes(&pj, sizeof(struct sipart));
    print_bytes(&pj2, sizeof(struct sipart));
    error("Particles 'sipj' do not match after density (byte = %d)", j_not_ok);
  }

  /* --- Test the force loop --- */

  /* Call the symmetric version */
  runner_iact_sidm_force(r2, dx, pi.h, pj.h, &pi, &pj, a, H, with_cosmology,
                         &cosmo, &sidm_props, ti_current, time_base);

  /* Call the non-symmetric version */
  runner_iact_nonsym_sidm_force(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H,
                                with_cosmology, &cosmo, &sidm_props, ti_current,
                                time_base);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_sidm_force(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H,
                                with_cosmology, &cosmo, &sidm_props, ti_current,
                                time_base);

  /* Check that the particles have the same scattering rates
  Only meaningful for sigma_over_m=0 since particle velocities change if a kick
  is triggered */
  const float tol = 1e-6f;
  i_not_ok = fabs(pi.SIDM_rate - pi2.SIDM_rate) > tol;
  j_not_ok = fabs(pj.SIDM_rate - pj2.SIDM_rate) > tol;

  if (i_not_ok) {
    error("Particle i SIDM_rate mismatch: %.8e != %.8e", pi.SIDM_rate,
          pi2.SIDM_rate);
  }

  if (j_not_ok) {
    error("Particle j SIDM_rate mismatch: %.8e != %.8e", pj.SIDM_rate,
          pj2.SIDM_rate);
  }
}

int main(void) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  const char *param_filename = "testSIDM.yml";

#ifdef SIDM_NONE
  message("SIDM disabled: skipping unit test.");
  message("Test automatically passes.");
  return 0;
#endif

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  for (int i = 0; i < 100; ++i) {
    message("Random test %d/100", i + 1);
    test(param_filename);
    message("Passed test %d/100", i + 1);
  }
  message("All good");

  return 0;
}
