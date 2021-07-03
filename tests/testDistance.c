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
#include "swift.h"

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void compute_interaction(struct part *pi, struct part *pj, float a, float H) {

  /* Compute the distance between the two particles */
  const float dx[3] = {pi->x[0] - pj->x[0], pi->x[1] - pj->x[1],
                       pi->x[2] - pj->x[2]};
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  if (r2 < pi->h * pi->h * kernel_gamma2) {

    /* And interact them (density) */
    runner_iact_density(r2, dx, pi->h, pj->h, pi, pj, a, H);
    runner_iact_chemistry(r2, dx, pi->h, pj->h, pi, pj, a, H);
    runner_iact_pressure_floor(r2, dx, pi->h, pj->h, pi, pj, a, H);
    runner_iact_star_formation(r2, dx, pi->h, pj->h, pi, pj, a, H);

#ifdef EXTRA_HYDRO_LOOP

    /* And interact them (gradient) */
    runner_iact_gradient(r2, dx, pi->h, pj->h, pi, pj, a, H);
#endif

    /* And interact them (force) */
    runner_iact_force(r2, dx, pi->h, pj->h, pi, pj, a, H);
  }
}

void test(void) {

  /* Start with some values for the cosmological parameters */
  const float a = (float)random_uniform(0.8, 1.);
  const float H = 1.f;

  /* Create two random particles (don't do this at home !) */
  struct part pi, pj;
  for (size_t i = 0; i < sizeof(struct part) / sizeof(float); ++i) {
    *(((float *)&pi) + i) = (float)random_uniform(0., 2.);
    *(((float *)&pj) + i) = (float)random_uniform(0., 2.);
  }

  /* Make the particle smoothing length, id and time-bin reasonable */
  pi.h = 1.f;
  pj.h = 1.f;
  pi.id = 1ll;
  pj.id = 2ll;
  pi.time_bin = 1;
  pj.time_bin = 1;

  /* Place the first particle at (1, 1, 1) */
  pi.x[0] = 1.;
  pi.x[1] = 1.;
  pi.x[2] = 1.;

  /* Move the second particle at various distances from the first */
  for (double dist = 1.0f; 1.0 + dist > 1.0; dist /= 2.) {

    pj.x[0] = pi.x[0] + random_uniform(0., dist * pi.h);
    pj.x[1] = pi.x[1] + random_uniform(0., dist * pi.h);
    pj.x[2] = pi.x[2] + random_uniform(0., dist * pi.h);

    compute_interaction(&pi, &pj, a, H);
  }

  /* Also test 0 distance */
  pj.x[0] = pi.x[0];
  pj.x[1] = pi.x[1];
  pj.x[2] = pi.x[2];

  compute_interaction(&pi, &pj, a, H);
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
