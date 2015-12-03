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
 *  @file metispartition.c
 *  @brief Use METIS to partition a space based on weighted cells.
 *         The weights should be the particle counts so that we get
 *         these shared equally around the partitions.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

/* METIS headers. */
#ifdef HAVE_METIS
#include <metis.h>
#endif

/* Local headers. */
#include "space.h"
#include "error.h"
#include "debug.h"

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/**
 * @brief Split the space into a number of partitions using METIS.
 *
 * Split the space using METIS to derive a partitions using the
 *       cell particle counts as weights.
 */

void part_metis(struct space *s, int nparts) {

#if defined(HAVE_METIS)


  /* Nothing to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if ( nparts == 1 ) return;

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Allocate weights and adjacency arrays . */
  idx_t *xadj;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
    error("Failed to allocate xadj buffer.");
  idx_t *adjncy;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
    error("Failed to allocate adjncy array.");
  idx_t *weights_v;
  if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
    error("Failed to allocate weights array");
  idx_t *nodeIDs;
  if ((nodeIDs = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
    error("Failed to allocate nodeIDs array");

  /* Fill the xadj and adjncy array to define the graph of cells. */
  /* Loop over all cells. */
  int cid = 0;
  for (int l = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {

        /* Visit all neighbours of this cell, wrapping space at edges. */
        int p = 0;
        for (int i = -1; i <= 1; i++) {
          int ii = l + i;
          if (ii < 0)
            ii += s->cdim[0];
          else if (ii >= s->cdim[0])
            ii -= s->cdim[0];
          for (int j = -1; j <= 1; j++) {
            int jj = m + j;
            if (jj < 0)
              jj += s->cdim[1];
            else if (jj >= s->cdim[1])
              jj -= s->cdim[1];
            for (int k = -1; k <= 1; k++) {
              int kk = n + k;
              if (kk < 0)
                kk += s->cdim[2];
              else if (kk >= s->cdim[2])
                kk -= s->cdim[2];

              /* If not self, record id of neighbour. */
              if (i || j || k) {
                adjncy[cid * 26 + p] = cell_getid(s->cdim, ii, jj, kk);
                p++;
              }
            }
          }
        }

        /* Next cell. */
        cid++;
      }
    }
  }
  xadj[0] = 0;
  for (int k = 0; k < ncells; k++) xadj[k + 1] = xadj[k] + 26;

  /* Init the weights array. XXX need particle counts. */
  for (int l = 0, k = 0; l < s->cdim[0]; l++) {
    int ll = l;
    if (ll > s->cdim[0]/2) ll -= s->cdim[0]/2;
    for (int m = 0; m < s->cdim[1]; m++) {
      int mm = l;
      if (mm > s->cdim[1]/2) mm -= s->cdim[1]/2;
      for (int n = 0; n < s->cdim[2]; n++) {
        int nn = l;
        if (nn > s->cdim[2]/2) nn -= s->cdim[2]/2;
        weights_v[k] = ll + mm + nn;
        //weights_v[k] = 1;
        message("weight: %d %d %d %d", l, m, n, weights_v[k]);
        k++;
      }
    }
  }

  //for (int k = 0; k < ncells; k++)
  //  weights_v[k] = (idx_t) 10000.0 * ((double)rand()/(double)RAND_MAX);

  /* Set the METIS options. */
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_NCUTS] = 10;
  options[METIS_OPTION_NITER] = 20;

  /* Call METIS. */
  idx_t one = 1;
  idx_t idx_ncells = ncells;
  idx_t idx_nparts = nparts;
  idx_t objval;

  /* Dump graph in METIS format */
  dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy,
                 weights_v, NULL, NULL);

  if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v,
                          NULL, NULL, &idx_nparts, NULL, NULL,
                          options, &objval, nodeIDs) != METIS_OK)
    error("Call to METIS_PartGraphKway failed.");

  /* Check that the nodeIDs are ok. */
  for (int k = 0; k < ncells; k++)
    if (nodeIDs[k] < 0 || nodeIDs[k] >= nparts)
      error("Got bad nodeID %"PRIDX" for cell %i.", nodeIDs[k], k);


  /* Check that the partition is complete and all nodes have some work. */
  int present[nparts];
  int failed = 0;
  for (int i = 0; i < nparts; i++) present[i] = 0;
  for (int i = 0; i < ncells; i++) present[nodeIDs[i]]++;
  for (int i = 0; i < nparts; i++) {
    message("count %d %d", i, present[i]);
    if (! present[i]) {
      failed = 1;
      message("Node %d is not present after repartition", i);
    }
  }
  if ( failed ) {
    message("Partitioning failed");
  }

  /* Set the cell nodeIDs and clear any non-local parts. */
  //for (int k = 0; k < ncells; k++) {
  //  s->cells[k].nodeID = nodeIDs[k];
  //}

  for (int l = 0, k = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {
        message("node %d %d %d -> %d", l, m, n, nodeIDs[k]);
        k++;
      }
    }
  }

  /* Clean up. */
  free(adjncy);
  free(weights_v);
  free(nodeIDs);

#else
  error("SWIFT was not compiled with METIS support.");
#endif
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

    s.cdim[0] = D1;
    s.cdim[1] = D2;
    s.cdim[2] = D3;

    message("# Partition space of %d,%d,%d cells into %d", s.cdim[0], s.cdim[1],
            s.cdim[2], N);

    part_metis(&s, N);

  } else {
    message( "no parts supplied" );
  }
}
