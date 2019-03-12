/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "swift.h"

void print_bytes(void *p, size_t len) {
  printf("(");
  for (size_t i = 0; i < len; ++i) {
    printf("%02x", ((unsigned char *)p)[i]);
    if (i % 4 == 3) printf("|");
  }
  printf(")\n");
}

void test(void) {

#if defined(SHADOWFAX_SPH)
  /* Initialize the Voronoi simulation box */
  double box_anchor[3] = {-2.0f, -2.0f, -2.0f};
  double box_side[3] = {6.0f, 6.0f, 6.0f};
/*  voronoi_set_box(box_anchor, box_side);*/
#endif

  /* Start with some values for the cosmological parameters */
  const float a = (float)random_uniform(0.8, 1.);
  const float H = 1.f;

  /* Create two random particles (don't do this at home !) */
  struct part pi, pj;
  for (size_t i = 0; i < sizeof(struct part) / sizeof(float); ++i) {
    *(((float *)&pi) + i) = (float)random_uniform(0., 2.);
    *(((float *)&pj) + i) = (float)random_uniform(0., 2.);
  }

  /* Make the particle smoothing length and position reasonable */
  for (size_t i = 0; i < 3; ++i) pi.x[0] = random_uniform(-1., 1.);
  for (size_t i = 0; i < 3; ++i) pj.x[0] = random_uniform(-1., 1.);
  pi.h = 2.f;
  pj.h = 2.f;
  pi.id = 1ll;
  pj.id = 2ll;
  pi.time_bin = 1;
  pj.time_bin = 1;

#if defined(GIZMO_MFV_SPH) || defined(SHADOWFAX_SPH)
  /* Give the primitive variables sensible values, since the Riemann solver does
     not like negative densities and pressures */
  pi.primitives.rho = random_uniform(0.1f, 1.0f);
  pi.primitives.v[0] = random_uniform(-10.0f, 10.0f);
  pi.primitives.v[1] = random_uniform(-10.0f, 10.0f);
  pi.primitives.v[2] = random_uniform(-10.0f, 10.0f);
  pi.primitives.P = random_uniform(0.1f, 1.0f);
  pj.primitives.rho = random_uniform(0.1f, 1.0f);
  pj.primitives.v[0] = random_uniform(-10.0f, 10.0f);
  pj.primitives.v[1] = random_uniform(-10.0f, 10.0f);
  pj.primitives.v[2] = random_uniform(-10.0f, 10.0f);
  pj.primitives.P = random_uniform(0.1f, 1.0f);
  /* make gradients zero */
  pi.primitives.gradients.rho[0] = 0.0f;
  pi.primitives.gradients.rho[1] = 0.0f;
  pi.primitives.gradients.rho[2] = 0.0f;
  pi.primitives.gradients.v[0][0] = 0.0f;
  pi.primitives.gradients.v[0][1] = 0.0f;
  pi.primitives.gradients.v[0][2] = 0.0f;
  pi.primitives.gradients.v[1][0] = 0.0f;
  pi.primitives.gradients.v[1][1] = 0.0f;
  pi.primitives.gradients.v[1][2] = 0.0f;
  pi.primitives.gradients.v[2][0] = 0.0f;
  pi.primitives.gradients.v[2][1] = 0.0f;
  pi.primitives.gradients.v[2][2] = 0.0f;
  pi.primitives.gradients.P[0] = 0.0f;
  pi.primitives.gradients.P[1] = 0.0f;
  pi.primitives.gradients.P[2] = 0.0f;
  pj.primitives.gradients.rho[0] = 0.0f;
  pj.primitives.gradients.rho[1] = 0.0f;
  pj.primitives.gradients.rho[2] = 0.0f;
  pj.primitives.gradients.v[0][0] = 0.0f;
  pj.primitives.gradients.v[0][1] = 0.0f;
  pj.primitives.gradients.v[0][2] = 0.0f;
  pj.primitives.gradients.v[1][0] = 0.0f;
  pj.primitives.gradients.v[1][1] = 0.0f;
  pj.primitives.gradients.v[1][2] = 0.0f;
  pj.primitives.gradients.v[2][0] = 0.0f;
  pj.primitives.gradients.v[2][1] = 0.0f;
  pj.primitives.gradients.v[2][2] = 0.0f;
  pj.primitives.gradients.P[0] = 0.0f;
  pj.primitives.gradients.P[1] = 0.0f;
  pj.primitives.gradients.P[2] = 0.0f;

#ifdef SHADOWFAX_SPH
  /* set time step to reasonable value */
  pi.force.dt = 0.001;
  pj.force.dt = 0.001;

  voronoi_cell_init(&pi.cell, pi.x, box_anchor, box_side);
  voronoi_cell_init(&pj.cell, pj.x, box_anchor, box_side);
#endif

#endif

  /* Make an xpart companion */
  struct xpart xpi, xpj;
  bzero(&xpi, sizeof(struct xpart));
  bzero(&xpj, sizeof(struct xpart));

  /* Make some copies */
  struct part pi2, pj2;
  memcpy(&pi2, &pi, sizeof(struct part));
  memcpy(&pj2, &pj, sizeof(struct part));

  int i_not_ok = memcmp(&pi, &pi2, sizeof(struct part));
  int j_not_ok = memcmp(&pj, &pj2, sizeof(struct part));

  if (i_not_ok) error("Particles 'pi' do not match after copy");
  if (j_not_ok) error("Particles 'pj' do not match after copy");

  /* Compute distance vector */
  float dx[3];
  dx[0] = pi.x[0] - pj.x[0];
  dx[1] = pi.x[1] - pj.x[1];
  dx[2] = pi.x[2] - pj.x[2];
  float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* --- Test the density loop --- */

  /* Call the symmetric version */
  runner_iact_density(r2, dx, pi.h, pj.h, &pi, &pj, a, H);
  runner_iact_chemistry(r2, dx, pi.h, pj.h, &pi, &pj, a, H);

  /* Call the non-symmetric version */
  runner_iact_nonsym_density(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H);
  runner_iact_nonsym_chemistry(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_density(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H);
  runner_iact_nonsym_chemistry(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H);

  /* Check that the particles are the same */
  i_not_ok = memcmp(&pi, &pi2, sizeof(struct part));
  j_not_ok = memcmp(&pj, &pj2, sizeof(struct part));

  if (i_not_ok) {
    printParticle_single(&pi, &xpi);
    printParticle_single(&pi2, &xpi);
    print_bytes(&pi, sizeof(struct part));
    print_bytes(&pi2, sizeof(struct part));
    error("Particles 'pi' do not match after density (byte = %d)", i_not_ok);
  }
  if (j_not_ok) {
    printParticle_single(&pj, &xpj);
    printParticle_single(&pj2, &xpj);
    print_bytes(&pj, sizeof(struct part));
    print_bytes(&pj2, sizeof(struct part));
    error("Particles 'pj' do not match after density (byte = %d)", j_not_ok);
  }

    /* --- Test the gradient loop --- */
#ifdef EXTRA_HYDRO_LOOP

  /* Call the symmetric version */
  runner_iact_gradient(r2, dx, pi.h, pj.h, &pi, &pj, a, H);

  /* Call the non-symmetric version */
  runner_iact_nonsym_gradient(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_gradient(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H);

  i_not_ok = memcmp((char *)&pi, (char *)&pi2, sizeof(struct part));
  j_not_ok = memcmp((char *)&pj, (char *)&pj2, sizeof(struct part));

  if (i_not_ok) {
    printParticle_single(&pi, &xpi);
    printParticle_single(&pi2, &xpi);
    print_bytes(&pi, sizeof(struct part));
    print_bytes(&pi2, sizeof(struct part));
    error("Particles 'pi' do not match after gradient (byte = %d)", i_not_ok);
  }
  if (j_not_ok) {
    printParticle_single(&pj, &xpj);
    printParticle_single(&pj2, &xpj);
    print_bytes(&pj, sizeof(struct part));
    print_bytes(&pj2, sizeof(struct part));
    error("Particles 'pj' do not match after gradient (byte = %d)", j_not_ok);
  }
#endif

  /* --- Test the force loop --- */

  /* Call the symmetric version */
  runner_iact_force(r2, dx, pi.h, pj.h, &pi, &pj, a, H);

  /* Call the non-symmetric version */
  runner_iact_nonsym_force(r2, dx, pi2.h, pj2.h, &pi2, &pj2, a, H);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_force(r2, dx, pj2.h, pi2.h, &pj2, &pi2, a, H);

/* Check that the particles are the same */
#if defined(GIZMO_MFV_SPH)
  i_not_ok = 0;
  j_not_ok = 0;
  for (size_t i = 0; i < sizeof(struct part) / sizeof(float); ++i) {
    float aa = *(((float *)&pi) + i);
    float bb = *(((float *)&pi2) + i);
    float cc = *(((float *)&pj) + i);
    float dd = *(((float *)&pj2) + i);

    int a_is_b;
    if ((aa + bb)) {
      a_is_b = (fabs((aa - bb) / (aa + bb)) > 1.e-4);
    } else {
      a_is_b = !(aa == 0.0f);
    }
    int c_is_d;
    if ((cc + dd)) {
      c_is_d = (fabs((cc - dd) / (cc + dd)) > 1.e-4);
    } else {
      c_is_d = !(cc == 0.0f);
    }

    if (a_is_b) {
      message("%.8e, %.8e, %lu", aa, bb, i);
    }
    if (c_is_d) {
      message("%.8e, %.8e, %lu", cc, dd, i);
    }

    i_not_ok |= a_is_b;
    j_not_ok |= c_is_d;
  }
#else
  i_not_ok = memcmp((char *)&pi, (char *)&pi2, sizeof(struct part));
  j_not_ok = memcmp((char *)&pj, (char *)&pj2, sizeof(struct part));
#endif

  if (i_not_ok) {
    printParticle_single(&pi, &xpi);
    printParticle_single(&pi2, &xpi);
    print_bytes(&pi, sizeof(struct part));
    print_bytes(&pi2, sizeof(struct part));
    error("Particles 'pi' do not match after force (byte = %d)", i_not_ok);
  }
  if (j_not_ok) {
    printParticle_single(&pj, &xpj);
    printParticle_single(&pj2, &xpj);
    print_bytes(&pj, sizeof(struct part));
    print_bytes(&pj2, sizeof(struct part));
    error("Particles 'pj' do not match after force (byte = %d)", j_not_ok);
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
  }
  message("All good");

  return 0;
}
