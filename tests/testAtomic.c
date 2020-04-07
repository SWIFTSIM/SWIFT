/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2020 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Standard includes. */
#include <fenv.h>

/* Local includes */
#include "swift.h"

const int array_size = 2048 * 2048;
const int num_threads = 64;
const int chunk_size = 64;

void map_function_sum_f(void *data, int num_elements, void *extra_data) {

  float *array = (float *)data;
  float *sum = (float *)extra_data;

  for (int i = 0; i < num_elements; ++i) atomic_add_f(sum, array[i]);
}

void map_function_sum_ll(void *data, int num_elements, void *extra_data) {

  long long *array = (long long *)data;
  long long *sum = (long long *)extra_data;

  for (int i = 0; i < num_elements; ++i) atomic_add(sum, array[i]);
}

void map_function_inc_ll(void *data, int num_elements, void *extra_data) {

  long long *sum = (long long *)extra_data;

  for (int i = 0; i < num_elements; ++i) atomic_inc(sum);
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

  /* Start a bunch of threads */
  printf("# Creating threadpool with %d threads\n", num_threads);
  struct threadpool tp;
  threadpool_init(&tp, num_threads);

  /* Create some random data */
  float *array_f = malloc(array_size * sizeof(float));
  long long *array_ll = malloc(array_size * sizeof(long long));

  for (int i = 0; i < array_size; ++i) {
    array_f[i] = rand() / ((float)RAND_MAX);
    array_ll[i] = rand();
  }

  /*** Test the addition atomic ops *******************************/

  /* float case */

  /* Compute the real answer */
  float real_sum_f = 0.f;
  for (int i = 0; i < array_size; ++i) {
    real_sum_f += array_f[i];
  }

  /* Compute the answer via threads and atomic */
  float atomic_sum_f = 0.f;
  threadpool_map(&tp, map_function_sum_f, array_f, array_size, sizeof(float),
                 chunk_size, &atomic_sum_f);

  const double diff_sum_f = (double)real_sum_f - (double)atomic_sum_f;
  const double sum_sum_f = (double)real_sum_f + (double)atomic_sum_f;
  const double rel_sum_f = 0.5 * fabs(diff_sum_f) / sum_sum_f;
  message("Real sum = %.7e -- atomic sum = %.7e rel=%e", real_sum_f,
          atomic_sum_f, rel_sum_f);

  /* long long case */

  /* Compute the real answer */
  long long real_sum_ll = 0.f;
  for (int i = 0; i < array_size; ++i) {
    real_sum_ll += array_ll[i];
  }

  /* Compute the answer via threads and atomic */
  long long atomic_sum_ll = 0LL;
  threadpool_map(&tp, map_function_sum_ll, array_ll, array_size,
                 sizeof(long long), chunk_size, &atomic_sum_ll);

  const double diff_sum_ll = (double)real_sum_ll - (double)atomic_sum_ll;
  const double sum_sum_ll = (double)real_sum_ll + (double)atomic_sum_ll;
  const double rel_sum_ll = 0.5 * fabs(diff_sum_ll) / sum_sum_ll;
  message("Real sum = %lld -- atomic sum = %lld rel=%e", real_sum_ll,
          atomic_sum_ll, rel_sum_ll);

  /*** Test the inc atomic ops *******************************/

  long long real_inc_ll = array_size;

  /* Compute the answer via threads and atomic */
  long long atomic_inc_ll = 0LL;
  threadpool_map(&tp, map_function_inc_ll, array_ll, array_size,
                 sizeof(long long), chunk_size, &atomic_inc_ll);

  const double diff_inc_ll = (double)real_inc_ll - (double)atomic_inc_ll;
  const double sum_inc_ll = (double)real_inc_ll + (double)atomic_inc_ll;
  const double rel_inc_ll = 0.5 * fabs(diff_inc_ll) / sum_inc_ll;
  message("Real inc = %lld -- atomic inc = %lld rel=%e", real_inc_ll,
          atomic_inc_ll, rel_inc_ll);

  /* Be clean */
  threadpool_clean(&tp);
  free(array_f);
  free(array_ll);
  return 0;
}
