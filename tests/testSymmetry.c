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

int main(int argc, char *argv[]) {

/* Choke on FPEs */
#ifdef HAVE_FE_ENABLE_EXCEPT
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

#if defined(SHADOWFAX_SPH)
  /* Initialize the Voronoi simulation box */
  double box_anchor[3] = {-2.0f, -2.0f, -2.0f};
  double box_side[3] = {6.0f, 6.0f, 6.0f};
/*  voronoi_set_box(box_anchor, box_side);*/
#endif

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
  pi.id = 1;
  pj.id = 2;

#if defined(GIZMO_SPH) || defined(SHADOWFAX_SPH)
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
  /* set time step to reasonable value */
  pi.force.dt = 0.001;
  pj.force.dt = 0.001;

#ifdef SHADOWFAX_SPH
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

  int i_ok = memcmp(&pi, &pi2, sizeof(struct part));
  int j_ok = memcmp(&pj, &pj2, sizeof(struct part));

  if (i_ok != 0) error("Particles 'pi' do not match after copy");
  if (j_ok != 0) error("Particles 'pj' do not match after copy");

  /* Compute distance vector */
  float dx[3];
  dx[0] = pi.x[0] - pj.x[0];
  dx[1] = pi.x[1] - pj.x[1];
  dx[2] = pi.x[2] - pj.x[2];
  float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  /* --- Test the density loop --- */

  /* Call the symmetric version */
  runner_iact_density(r2, dx, pi.h, pj.h, &pi, &pj);

  /* Call the non-symmetric version */
  runner_iact_nonsym_density(r2, dx, pi2.h, pj2.h, &pi2, &pj2);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_density(r2, dx, pj2.h, pi2.h, &pj2, &pi2);

  /* Check that the particles are the same */
  i_ok = memcmp(&pi, &pi2, sizeof(struct part));
  j_ok = memcmp(&pj, &pj2, sizeof(struct part));

  if (i_ok) error("Particles 'pi' do not match after density");
  if (j_ok) error("Particles 'pj' do not match after density");

  /* --- Test the force loop --- */

  /* Call the symmetric version */
  runner_iact_force(r2, dx, pi.h, pj.h, &pi, &pj);

  /* Call the non-symmetric version */
  runner_iact_nonsym_force(r2, dx, pi2.h, pj2.h, &pi2, &pj2);
  dx[0] = -dx[0];
  dx[1] = -dx[1];
  dx[2] = -dx[2];
  runner_iact_nonsym_force(r2, dx, pj2.h, pi2.h, &pj2, &pi2);

/* Check that the particles are the same */
#if defined(GIZMO_SPH)
  i_ok = 0;
  j_ok = 0;
  for (size_t i = 0; i < sizeof(struct part) / sizeof(float); ++i) {
    float a = *(((float *)&pi) + i);
    float b = *(((float *)&pi2) + i);
    float c = *(((float *)&pj) + i);
    float d = *(((float *)&pj2) + i);

    int a_is_b;
    if ((a + b)) {
      a_is_b = (fabs((a - b) / (a + b)) > 1.e-4);
    } else {
      a_is_b = !(a == 0.0f);
    }
    int c_is_d;
    if ((c + d)) {
      c_is_d = (fabs((c - d) / (c + d)) > 1.e-4);
    } else {
      c_is_d = !(c == 0.0f);
    }

    if (a_is_b) {
      message("%.8e, %.8e, %lu", a, b, i);
    }
    if (c_is_d) {
      message("%.8e, %.8e, %lu", c, d, i);
    }

    i_ok |= a_is_b;
    j_ok |= c_is_d;
  }
#else
  i_ok = memcmp(&pi, &pi2, sizeof(struct part));
  j_ok = memcmp(&pj, &pj2, sizeof(struct part));
#endif

  if (i_ok) {
    printParticle_single(&pi, &xpi);
    printParticle_single(&pi2, &xpi);
    error("Particles 'pi' do not match after force");
  }
  if (j_ok) {
    printParticle_single(&pj, &xpj);
    printParticle_single(&pj2, &xpj);
    error("Particles 'pj' do not match after force");
  }

  return 0;
}
