/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#include <fenv.h>

/* Local headers. */
#include "swift.h"

const int N = 1000000;

int main(int argc, char* argv[]) {

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

  /* Log the swift random seed */
  message("SWIFT random seed = %d", SWIFT_RANDOM_SEED_XOR);

  /* Get a random poisson parameter and ID */
  const double lambda = exp10(random_uniform(-1, 3));
  const long long id = random_uniform(0, 1e10);
  const int offset = random_uniform(0, 1e6);
  const long long type = random_uniform(0, 1e10);

  message("Testing the generator for Lambda=%e", lambda);

  /* Verify that the mean and std. dev. are as expected */
  double mean = 0., mean2 = 0.;
  for (int i = 0; i < N; ++i) {

    const double r =
        random_poisson(id, lambda, offset + i, (enum random_number_type)type);

    mean += r;
    mean2 += r * r;
  }
  mean /= (double)N;
  mean2 /= (double)N;
  const double var = mean2 - mean * mean;

  /* Verify that the mean and variance are correct */
  message("Mean = %e expected = %e", mean, lambda);
  message("Vari = %e expected = %e", var, lambda);

  if (fabs(mean / lambda - 1.) > 0.01) {
    error("Incorrect mean!");
    return 1;
  }
  if (fabs(var / lambda - 1.) > 0.01) {
    error("Incorrect mean!");
    return 1;
  }

  return 0;
}
