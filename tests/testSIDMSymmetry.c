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

void test(void) {

  /* Start with some values for the cosmological parameters */
  const float a = (float)random_uniform(0.8, 1.);
  const float H = 1.f;

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

  /* Make an xpart companion */
  struct xpart xpi, xpj;
  bzero(&xpi, sizeof(struct xpart));
  bzero(&xpj, sizeof(struct xpart));

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
  runner_iact_sidm_density(r2, dx, pi.h, pj.h, &pi, &pj, a, H);

  /* Call the non-symmetric version */
  runner_iact_nonsym_sidm_density(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_sidm_density(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H);

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
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /* Get some randomness going */
  const int seed = time(NULL);
  message("Seed = %d", seed);
  srand(seed);

  for (int i = 0; i < 100; ++i) {
    message("Random test %d/100", i);
    test();
    message("Passed test %d/100", i);
  }
  message("All good");

  return 0;
}
