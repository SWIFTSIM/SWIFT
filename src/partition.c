/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/**
 *  @file partition.c
 *  @brief Pick sample cells to seed some partition of the space and
 *         apply the partition to the space.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

/* Local headers. */
#include "space.h"
#include "error.h"

/**
 *  @brief Pick a number of cell positions from a vectorized list.
 *
 *  Vectorize the space and pick positions in it for the number of expected
 *  partitions using a single step.
 *
 *  @param s the space.
 *  @param nparts the number of partitions.
 *  @param samplecells the list of sample cell positions, size of 3*nparts
 */
void part_vectorize(struct space *s, int nparts, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  if (nparts > length) {
    error("Too few cells (%d) for this number of partitions (%d)", length,
          nparts);
  }

  int step = length / nparts;
  int n = 0;
  int m = 0;
  int l = 0;

  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        if (n == 0 && l < nparts) {
          samplecells[m++] = i;
          samplecells[m++] = j;
          samplecells[m++] = k;
          l++;
        }
        n++;
        if (n == step) n = 0;
      }
    }
  }
}

/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometry closest
 * to each and return the counts per cells and optionally apply the partition
 * to the space.
 */
void part_apply(struct space *s, int nparts, int *samplecells,
                int apply, unsigned int *counts) {

  for (int i = 0; i < nparts; i++) {
    counts[i] = 0;
  }

  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nparts; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        if ( apply )
          s->cells[n++].nodeID = select;
        //message("@ %d %d %d %d", i, j, k, select);
        counts[select]++;
      }
    }
  }

  /* Test section */
  message("# counts:");
  unsigned int total = 0;
  for (int i = 0; i < nparts; i++) {
    message("#  %d %d", i, counts[i]);
    if ( counts[i] == 0 ) {
      message( "sampling failed" );
    }
    total += counts[i];
  }
  message("# total = %d", total);
}

int main(int argc, char *argv[]) {

  int N = 0;
  int D1 = 0;
  int D2 = 0;
  int D3 = 0;
  if ( argc > 4 ) {
    D1 = atoi( argv[1] );
    D2 = atoi( argv[2] );
    D3 = atoi( argv[3] );
    N = atoi( argv[4] );

    struct space s;
    int samplecells[3*N];
    unsigned int counts[N];

    s.cdim[0] = D1;
    s.cdim[1] = D2;
    s.cdim[2] = D3;

    part_vectorize(&s, N, samplecells);
    part_apply(&s, N, samplecells, 0, counts);

  } else {
    message( "no parts supplied" );
  }
}
