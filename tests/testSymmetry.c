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

  /* Choke if need be */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

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
  i_ok = memcmp(&pi, &pi2, sizeof(struct part));
  j_ok = memcmp(&pj, &pj2, sizeof(struct part));

  if (i_ok) error("Particles 'pi' do not match after force");
  if (j_ok) error("Particles 'pj' do not match after force");

  return 0;
}
