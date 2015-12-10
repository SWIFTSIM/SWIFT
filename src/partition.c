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
 *  @brief file of various techniques for partitioning a grid of cells
 *         into geometrically connected regions.
 *
 *  Currently supported types, grid, vectorise and metis.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS headers only used when MPI is also available. */
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "space.h"
#include "partition.h"
#include "const.h"
#include "error.h"
#include "debug.h"

/* Useful defines. */
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define CHUNK 512

/* Simple descriptions of initial partition types for reports. */
const char *initpart_name[] = {
    "gridded cells",
    "vectorized point associated cells",
    "METIS particle weighted cells", 
    "METIS unweighted cells"
};

/* Simple descriptions of repartition types for reports. */
const char *repart_name[] = {
    "no", 
    "METIS edge and vertex time weighted cells",
    "METIS particle count vertex weighted cells",
    "METIS time edge weighted cells",
    "METIS particle count vertex and time edge cells"
};

/*  Vectorisation support */
/*  ===================== */

/**
 *  @brief Pick a number of cell positions from a vectorised list.
 *
 *  Vectorise the cell space and pick positions in it for the number of
 *  expected regions using a single step. Vectorisation is guaranteed
 *  to work, providing there are more cells than regions.
 *
 *  @param s the space.
 *  @param nregions the number of regions
 *  @param samplecells the list of sample cell positions, size of 3*nregions
 */
void part_pick_vector(struct space *s, int nregions, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }

  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;

  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        if (n == 0 && l < nregions) {
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
 * to each and apply the partition to the space.
 */
void part_split_vector(struct space *s, int nregions, int *samplecells) {
  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells[n++].nodeID = select;
      }
    }
  }
}

/* METIS support
 * =============
 *
 * METIS partitions using a multi-level k-way scheme. We support using this in
 * a unweighted scheme, which works well and seems to be guaranteed, and a
 * weighted by the number of particles scheme. Note METIS is optional.
 */

/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the
 * cell particle counts as weights.
 *
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjcny array.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 */
void part_pick_metis(struct space *s, int nregions, int *vertexw, int *edgew,
                     int *celllist) {
#if defined(HAVE_METIS)

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* Allocate weights and adjacency arrays . */
  idx_t *xadj;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
    error("Failed to allocate xadj buffer.");
  idx_t *adjncy;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
    error("Failed to allocate adjncy array.");
  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate vertex weights array");
  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)malloc(26 * sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate edge weights array");
  idx_t *regionid;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
    error("Failed to allocate regionid array");

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

  /* Init the vertex weights array. */
  if (vertexw != NULL) {
    for (int k = 0; k < ncells; k++) {
      if (vertexw[k] > 0) {
        weights_v[k] = vertexw[k];
      } else {
        weights_v[k] = 1;
      }
    }
  }

  /* Init the edges weights array. */
  if (edgew != NULL) {
    for (int k = 0; k < 26 * ncells; k++) {
      if (edgew[k] > 0) {
        weights_e[k] = edgew[k];
      } else {
        weights_e[k] = 1;
      }
    }
  }

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
  idx_t idx_nregions = nregions;
  idx_t objval;

  /* Dump graph in METIS format */
  /* dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy,
   *                weights_v, weights_e, NULL);
   */

  if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v, weights_e,
                          NULL, &idx_nregions, NULL, NULL, options, &objval,
                          regionid) != METIS_OK)
    error("Call to METIS_PartGraphKway failed.");

  /* Check that the regionids are ok. */
  for (int k = 0; k < ncells; k++)
    if (regionid[k] < 0 || regionid[k] >= nregions)
      error("Got bad nodeID %" PRIDX " for cell %i.", regionid[k], k);

  /* Set the cell list to the region index. */
  for (int k = 0; k < ncells; k++) {
    celllist[k] = regionid[k];
  }

  /* Clean up. */
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(xadj);
  free(adjncy);
  free(regionid);

#else
  error("SWIFT was not compiled with METIS support.");
#endif
}

/**
 * @brief Apply METIS cell list partitioning to a cell structure.
 *
 * Uses the results of part_metis_pick to assign each cell's nodeID to the
 * picked region index, thus partitioning the space into regions.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
void part_split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) 
      s->cells[i].nodeID = celllist[i];
}

/*  General support */
/*  =============== */

/**
 * @brief Check if all regions have been assigned a node in the
 *        cells of a space.
 *
 * @param s the space containing the cells to check.
 * @param nregions number of regions expected.
 * @param verbose if true report the missing regions.
 * @return true if all regions have been found, false otherwise.
 */
int part_check_complete(struct space *s, int verbose, int nregions) {

  int *present = (int *)malloc(sizeof(int) * nregions);
  if (present == NULL) error("Failed to allocate present array");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    present[s->cells[i].nodeID]++;
  }
  for (int i = 0; i < nregions; i++) {
    if (!present[i]) {
      failed = 1;
      if (verbose) message("Region %d is not present in partition", i);
    }
  }
  free(present);
  return (!failed);
}
