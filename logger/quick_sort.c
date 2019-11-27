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

/* Include corresponding header */
#include "quick_sort.h"

/* Include local headers */
#include "logger_index.h"

struct qstack {
  int64_t lo;
  int64_t hi;
};

/**
 * @brief Sort the data with the quicksort according to the ids.
 */
void quick_sort(struct index_data *data, size_t N) {
  int64_t qpos, i, j, lo, hi, imin, pivot;
  struct index_data temp;

  /* Allocate a stack of operations */
  int stack_size = log(N) + 1;
  struct qstack *qstack =
      (struct qstack *)malloc(sizeof(struct qstack) * stack_size);

  /* Sort parts in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    /* Do we have a low number of element to sort? */
    if (hi - lo < 15) {
      /* Sort the last elements. */
      for (i = lo; i < hi; i++) {
        imin = i;
        /* Find the minimal value. */
        for (j = i + 1; j <= hi; j++)
          if (data[j].id < data[imin].id) imin = j;
        /* Did we find the minimal value? */
        if (imin != i) {
          /* Swap the elements. */
          temp = data[imin];
          data[imin] = data[i];
          data[i] = temp;
        }
      }
    } else {
      /* Select a pivot */
      pivot = data[(lo + hi) / 2].id;
      i = lo;
      j = hi;
      /* Ensure that the elements before/after the pivot
         are smaller/larger than the pivot. */
      while (i <= j) {
        /* Find the first elements that do not respect
           the order. */
        while (data[i].id < pivot) i++;
        while (data[j].id > pivot) j--;

        /* Did we get two elements */
        if (i <= j) {
          /* Are they different? */
          if (i < j) {
            /* Swap them */
            temp = data[i];
            data[i] = data[j];
            data[j] = temp;
          }
          /* The two elements are good now */
          i += 1;
          j -= 1;
        }
      }
      /* Add the next operations to the stack.
       * The order is important in order to decrease the stack size.
       */
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }

  /* Free the allocated memory */
  free(qstack);
}
