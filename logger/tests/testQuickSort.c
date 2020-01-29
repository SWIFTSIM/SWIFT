
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

#include "memswap.h"
#include "quick_sort.h"

#define N 10000

int64_t getRandom(void) {
  return llabs(((int64_t)lrand48() << 32) + (int64_t)lrand48());
}
/**
 * @brief Initialize the array.
 */
void init_array(struct index_data *data) {
  /* Assign the ids */
  for (int i = 0; i < N; i++) {
    data[i].id = getRandom();
  }

  /* randomize the array */
  for (int i = 0; i < N; i++) {
    /* Select an index randomly */
    int j = rand() % N;

    /* Swap the two elements */
    memswap(&data[i], &data[j], sizeof(struct index_data));
  }
}

/**
 * @brief Ensure that the array is sorted
 */
void check_sort(struct index_data *data) {
  for (size_t i = 1; i < N; i++) {
    if (data[i].id < data[i - 1].id) {
      error("The array is not sorted index=%zi, prev=%li, cur=%li", i,
            data[i - 1].id, data[i].id);
    }
  }
}

/**
 * @brief print the array
 */
void print_array(struct index_data *data) {
  for (int i = 0; i < N; i++) {
    printf("%li \t %zi\n", data[i].id, data[i].offset);
  }
}

int main(int argc, char *argv[]) {

  /* Create the array */
  struct index_data *data =
      (struct index_data *)malloc(N * sizeof(struct index_data));

  if (data == NULL) {
    error("Failed to allocate the memory");
  }

  /* Initialize the array */
  init_array(data);

  /* Print the array */
  message("\nArray before sort\n");
  print_array(data);

  /* Sort the array */
  quick_sort(data, N);

  /* Print the array */
  message("\nArray after sort\n");
  print_array(data);

  /* Check the results */
  check_sort(data);

  return 0;
}
