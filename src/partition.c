/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
 *  @brief file of various techniques for partitioning and repartitioning
 *  a grid of cells into geometrically connected regions and distributing
 *  these around a number of MPI nodes.
 *
 *  Currently supported partitioning types: grid, vectorise and ParMETIS.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#endif

/* Local headers. */
#include "debug.h"
#include "error.h"
#include "partition.h"
#include "restart.h"
#include "space.h"
#include "tools.h"

/* Maximum weight used for ParMETIS. */
#define parmetis_maxweight 10000.0f

/* Simple descriptions of initial partition types for reports. */
const char *initial_partition_name[] = {
    "axis aligned grids of cells", "vectorized point associated cells",
    "memory balanced, using ParMETIS particle weighted cells",
    "similar sized regions, using ParMETIS unweighted cells"};

/* Simple descriptions of repartition types for reports. */
const char *repartition_name[] = {
    "no",
    "ParMETIS edge and vertex task cost weights",
    "ParMETIS task cost edge weights",
    "ParMETIS vertex task costs and edge delta timebin weights",
    "ParMETIS edge delta timebin weights",
};

/* Local functions, if needed. */
static int check_complete(struct space *s, int verbose, int nregions);

/*  Vectorisation support */
/*  ===================== */

#if defined(WITH_MPI)
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
static void pick_vector(struct space *s, int nregions, int *samplecells) {

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
#endif

#if defined(WITH_MPI)
/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometrically
 * closest and apply the partition to the space.
 */
static void split_vector(struct space *s, int nregions, int *samplecells) {
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
        s->cells_top[n++].nodeID = select;
      }
    }
  }
}
#endif

/* ParMETIS support
 * ================
 *
 * ParMETIS partitions using a multi-level k-way scheme. We support using this
 * in a unweighted scheme, which works well and seems to be guaranteed, and a
 * weighted by the number of particles scheme. Note ParMETIS is optional.
 *
 * Repartitioning is based on ParMETIS and uses weights determined from the
 * estimated costs that a cells tasks will take or the relative time bins of
 * the cells next updates.
 */

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
/**
 * @brief Fill the adjncy array defining the graph of cells in a space.
 *
 * See the ParMETIS and METIS manuals if you want to understand this
 * format. The cell graph consists of all nodes as vertices with edges as the
 * connections to all neighbours, so we have 26 per vertex. Note you will
 * also need an xadj array, for a single rank that would be:
 *
 *   xadj[0] = 0;
 *   for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
 *
 * but each rank needs a different xadj.
 *
 * @param s the space of cells.
 * @param adjncy the adjncy array to fill, must be of size 26 * the number of
 *               cells in the space.
 */
static void graph_init_adjncy(struct space *s, idx_t *adjncy) {

  /* Loop over all cells in the space. */
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
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)

struct counts_mapper {
    int cdim[3];
    double iwidth[3];
    double dim[3];
    int64_t *counts;
};

/**
 * @brief Threadpool helper for accumulating the counts of particles per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of particles per cell. Should be
 *               allocated as size s->nr_parts.
 */
static void accumulate_counts_mapper(void *map_data, int num_elements,
                                     void *extra_data) {

  struct part *parts = (struct part *) map_data;
  struct counts_mapper *mydata = (struct counts_mapper *)extra_data;

  for (int k = 0; k < num_elements; k++) {
      for (int j = 0; j < 3; j++) {
          if (parts[k].x[j] < 0.0)
              parts[k].x[j] += mydata->dim[j];
          else if (parts[k].x[j] >= mydata->dim[j])
              parts[k].x[j] -= mydata->dim[j];
      }
      const int cid =
          cell_getid(mydata->cdim, parts[k].x[0] * mydata->iwidth[0],
                     parts[k].x[1] * mydata->iwidth[1],
                     parts[k].x[2] * mydata->iwidth[2]);
      atomic_inc(&(mydata->counts[cid]));
  }
}

/**
 * @brief Accumulate the counts of particles per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of particles per cell. Should be
 *               allocated as size s->nr_cells.
 */
static void accumulate_counts(struct space *s, double *counts) {

  /* Need integer counts for atomic access in mapper. */
  int64_t *icounts = NULL;
  if ((icounts = (int64_t *)malloc(sizeof(int64_t) * s->nr_cells)) == NULL)
    error("Failed to allocate icounts buffer");
  bzero(icounts, sizeof(int64_t) * s->nr_cells);

  struct counts_mapper extra_data;
  extra_data.cdim[0] = s->cdim[0];
  extra_data.cdim[1] = s->cdim[1];
  extra_data.cdim[2] = s->cdim[2];
  extra_data.iwidth[0] = s->iwidth[0], 
  extra_data.iwidth[1] = s->iwidth[1], 
  extra_data.iwidth[2] = s->iwidth[2];
  extra_data.dim[0] = s->dim[0];
  extra_data.dim[1] = s->dim[1];
  extra_data.dim[2] = s->dim[2];
  extra_data.counts = icounts;

  threadpool_map(&s->e->threadpool, accumulate_counts_mapper, s->parts,
                 s->nr_parts, sizeof(struct part), 0, &extra_data);

  for (int k = 0; k < s->nr_cells; k++) counts[k] = icounts[k]; // XXX use
                                                                // int64_t for
                                                                // all weights
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
/**
 * @brief Apply ParMETIS cell-list partitioning to a cell structure.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
static void split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) s->cells_top[i].nodeID = celllist[i];

  /* To check or visualise the partition dump all the cells. */
  dumpCellRanks("metis_partition", s->cells_top, s->nr_cells);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)

/* qsort support. */
struct indexval {
  int index;
  int count;
  int old;
  int new;
};
static int indexvalcmp(const void *p1, const void *p2) {
  const struct indexval *iv1 = (const struct indexval *)p1;
  const struct indexval *iv2 = (const struct indexval *)p2;
  return iv2->count - iv1->count;
}

/**
 * @brief Check if there is a permutation of the region indices of our cells
 *        that will reduce the amount of particle movement and return it.
 *
 * @param newlist the new list of regions for our cells.
 * @param oldlist the old list of regions for our cells.
 * @param nregions the number of regions.
 * @param ncells the number of cells.
 * @param permlist the permutation of the newlist.
 */
void permute_regions(int *newlist, int *oldlist, int nregions, int ncells,
                     int *permlist) {

  /* We want a solution in which the current region assignments of the cells
   * are preserved when possible, to avoid unneccesary particle movement.  So
   * create a 2d-array of counts of cells that are common to all pairs of old
   * and new lists. Each element of the array has a count of cells and an
   * unique index so we can sort into decreasing counts.
   */
  int indmax = nregions * nregions;
  struct indexval *ivs = NULL;
  if ((ivs = (struct indexval *)malloc(sizeof(struct indexval) * indmax)) ==
      NULL)
    error("Failed to allocate ivs structs");
  bzero(ivs, sizeof(struct indexval) * indmax);

  for (int k = 0; k < ncells; k++) {
    int index = newlist[k] + nregions * oldlist[k];
    ivs[index].count++;
    ivs[index].index = index;
    ivs[index].old = oldlist[k];
    ivs[index].new = newlist[k];
  }
  qsort(ivs, indmax, sizeof(struct indexval), indexvalcmp);

  /* Go through the ivs using the largest counts first, these are the
   * regions with the most cells in common, old partition to new. If not
   * returning the permutation, avoid the associated work. */
  int *oldmap = NULL;
  int *newmap = NULL;
  oldmap = permlist; /* Reuse this */
  if ((newmap = malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate newmap array");

  for (int k = 0; k < nregions; k++) {
    oldmap[k] = -1;
    newmap[k] = -1;
  }

  for (int k = 0; k < indmax; k++) {

    /* Stop when all regions with common cells have been considered. */
    if (ivs[k].count == 0) break;

    /* Store old and new IDs, if not already used. */
    if (newmap[ivs[k].new] == -1 && oldmap[ivs[k].old] == -1) {
      newmap[ivs[k].new] = ivs[k].old;
      oldmap[ivs[k].old] = ivs[k].new;
    }
  }

  /* Handle any regions that did not get selected by picking an unused rank
   * from oldmap and assigning to newmap. */
  int spare = 0;
  for (int k = 0; k < nregions; k++) {
    if (newmap[k] == -1) {
      for (int j = spare; j < nregions; j++) {
        if (oldmap[j] == -1) {
          newmap[k] = j;
          oldmap[j] = j;
          spare = j;
          break;
        }
      }
    }
  }

  /* Permute the newlist into this order. */
  for (int k = 0; k < ncells; k++) {
    permlist[k] = newmap[newlist[k]];
  }
  free(newmap);
  free(ivs);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
/**
 * @brief Partition the given space into a number of connected regions using
 *        ParMETIS.
 *
 * Split the space using PARMETIS to derive a partitions using the
 * given edge and vertex weights. If no weights are given then an
 * unweighted partition is performed. If refine is set then an existing
 * partition is assumed to be present from the last call to this routine
 * in the celllist argument, that will get a refined partition, not a new
 * one.
 *
 * Assumes MPI is up and running and the number of ranks is the same as the
 * number of regions.
 *
 * @param nodeID our nodeID.
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights. Need to be in the range of idx_t.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array. Need to be in the range of
 *        idx_t.
 * @param refine whether to refine an existing partition, or create a new one.
 * @param seed seed for random numbers, usually when refining this is best
 *             kept to the same value. Use -1 for the ParMETIS default.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells. If refine is 1, then this should contain
 *        the old partition on entry.
 */
static void pick_parmetis(int nodeID, struct space *s, int nregions,
                          double *vertexw, double *edgew, int refine, int seed,
                          int *celllist) {
  int res;
  MPI_Comm comm;
  MPI_Status status;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Nothing much to do if only using a single MPI rank. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* We all get one of these with the same content. It defines the ranges of
   * vertices that are found on each rank. This contiguity constraint seems to
   * stop efficient local processing, since our cell distributions do not
   * meet this requirement. That means the graph and related information needs
   * to be all brought to one node and redistributed for processing in
   * approproiate batches. */
  idx_t *vtxdist;
  if ((vtxdist = (idx_t *)malloc(sizeof(idx_t) * (nregions + 1))) == NULL)
    error("Failed to allocate vtxdist buffer.");

  if (nodeID == 0) {

    /* Construct vtxdist and send it to all ranks. Each rank gets an equal
     * number of vertices. */
    vtxdist[0] = 0;
    int k = ncells;
    for (int i = 0; i < nregions; i++) {
      int l = k / (nregions - i);
      vtxdist[i + 1] = vtxdist[i] + l;
      k -= l;
    }
    res = MPI_Bcast((void *)vtxdist, nregions + 1, IDX_T, 0, comm);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast vtxdist.");

  } else {
    res = MPI_Bcast((void *)vtxdist, nregions + 1, IDX_T, 0, comm);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast vtxdist.");
  }

  /* Number of cells on this node and space for the expected arrays. */
  int nverts = vtxdist[nodeID + 1] - vtxdist[nodeID];

  idx_t *xadj = NULL;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
    error("Failed to allocate xadj buffer.");

  idx_t *adjncy = NULL;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * nverts)) == NULL)
    error("Failed to allocate adjncy array.");

  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate vertex weights array");

  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)malloc(26 * sizeof(idx_t) * nverts)) == NULL)
      error("Failed to allocate edge weights array");

  idx_t *regionid = NULL;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
    error("Failed to allocate regionid array");

  /* Only use one rank to organize everything. */
  if (nodeID == 0) {

    /* Space for largest lists. */
    idx_t *full_xadj = NULL;
    if ((full_xadj =
             (idx_t *)malloc(sizeof(idx_t) * (ncells + nregions + 1))) == NULL)
      error("Failed to allocate xadj buffer.");
    idx_t *full_adjncy = NULL;
    if ((full_adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
      error("Failed to allocate adjncy array.");
    idx_t *full_weights_v = NULL;
    if (weights_v != NULL)
      if ((full_weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate vertex weights array");
    idx_t *full_weights_e = NULL;
    if (weights_e != NULL)
      if ((full_weights_e = (idx_t *)malloc(26 * sizeof(idx_t) * ncells)) ==
          NULL)
        error("Failed to allocate edge weights array");

    idx_t *full_regionid = NULL;
    if (refine) {
      if ((full_regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate regionid array");
    }

    /* Define the cell graph. */
    graph_init_adjncy(s, full_adjncy);

    /* xadj is set for each rank, different to serial version in that each
     * rank starts with 0 */
    for (int rank = 0, j = 0; rank < nregions; rank++) {

      /* Number of vertices for this rank. */
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      /* Start from 0, and step forward 26 edges each value. */
      full_xadj[j] = 0;
      for (int k = 0; k <= nvt; k++) {
        full_xadj[j + 1] = full_xadj[j] + 26;
        j++;
      }
    }

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          full_weights_v[k] = vertexw[k];
        } else {
          full_weights_v[k] = 1;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < ncells; k++) {
        if ((idx_t)vertexw[k] < 0) {
          message("Input vertex weight out of range: %ld", (long)vertexw[k]);
          failed++;
        }
        if (full_weights_v[k] < 1) {
          message("Used vertex weight  out of range: %" PRIDX,
                  full_weights_v[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d vertex weights are out of range", failed);
#endif
    }

    /* Init the edges weights array. */
    if (edgew != NULL) {
      for (int k = 0; k < 26 * ncells; k++) {
        if (edgew[k] > 1) {
          full_weights_e[k] = edgew[k];
        } else {
          full_weights_e[k] = 1;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < 26 * ncells; k++) {

        if ((idx_t)edgew[k] < 0) {
          message("Input edge weight out of range: %ld", (long)edgew[k]);
          failed++;
        }
        if (full_weights_e[k] < 1) {
          message("Used edge weight out of range: %" PRIDX, full_weights_e[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d edge weights are out of range", failed);
#endif
    }

    /* Send ranges to the other ranks and keep our own. XXX async version XXX */
    for (int rank = 0, j1 = 0, j2 = 0, j3 = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      if (refine)
        for (int i = 0; i < nvt; i++) full_regionid[i] = celllist[j3 + i];

      if (rank == 0) {
        memcpy(xadj, &full_xadj[j1], sizeof(idx_t) * (nvt + 1));
        memcpy(adjncy, &full_adjncy[j2], sizeof(idx_t) * nvt * 26);
        if (weights_e != NULL)
          memcpy(weights_e, &full_weights_e[j2], sizeof(idx_t) * nvt * 26);
        if (weights_v != NULL)
          memcpy(weights_v, &full_weights_v[j3], sizeof(idx_t) * nvt);
        if (refine) memcpy(regionid, full_regionid, sizeof(idx_t) * nvt);

      } else {
        res = MPI_Send(&full_xadj[j1], nvt + 1, IDX_T, rank, 0, comm);
        if (res == MPI_SUCCESS)
          res = MPI_Send(&full_adjncy[j2], nvt * 26, IDX_T, rank, 1, comm);
        if (res == MPI_SUCCESS && weights_e != NULL)
          res = MPI_Send(&full_weights_e[j2], nvt * 26, IDX_T, rank, 2, comm);
        if (res == MPI_SUCCESS && weights_v != NULL)
          res = MPI_Send(&full_weights_v[j3], nvt, IDX_T, rank, 3, comm);
        if (refine && res == MPI_SUCCESS)
          res = MPI_Send(full_regionid, nvt, IDX_T, rank, 4, comm);
        if (res != MPI_SUCCESS) mpi_error(res, "Failed to send graph data");
      }
      j1 += nvt + 1;
      j2 += nvt * 26;
      j3 += nvt;
    }

    /* Clean up. */
    if (weights_v != NULL) free(full_weights_v);
    if (weights_e != NULL) free(full_weights_e);
    free(full_xadj);
    free(full_adjncy);
    if (refine) free(full_regionid);

  } else {

    /* Receive stuff from rank 0. */
    res = MPI_Recv(xadj, nverts + 1, IDX_T, 0, 0, comm, &status);
    if (res == MPI_SUCCESS)
      res = MPI_Recv(adjncy, nverts * 26, IDX_T, 0, 1, comm, &status);
    if (res == MPI_SUCCESS && weights_e != NULL)
      res = MPI_Recv(weights_e, nverts * 26, IDX_T, 0, 2, comm, &status);
    if (res == MPI_SUCCESS && weights_v != NULL)
      res = MPI_Recv(weights_v, nverts, IDX_T, 0, 3, comm, &status);
    if (refine && res == MPI_SUCCESS)
      res += MPI_Recv((void *)regionid, nverts, IDX_T, 0, 4, comm, &status);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive graph data");
  }

  /* Set up the tpwgts array. This is just 1/nregions. */
  real_t *tpwgts;
  if ((tpwgts = (real_t *)malloc(sizeof(real_t) * nregions)) == NULL)
    error("Failed to allocate tpwgts array");
  for (int i = 0; i < nregions; i++) tpwgts[i] = 1.0 / (real_t)nregions;

  /* Common parameters. */
  idx_t options[10];
  options[0] = 1;
  options[1] = 0;
  options[2] = seed;

  idx_t edgecut;
  idx_t ncon = 1;
  idx_t nparts = nregions;
  idx_t numflag = 0;
  idx_t wgtflag = 0;
  if (edgew != NULL) wgtflag += 1;
  if (vertexw != NULL) wgtflag += 2;
  int permute = 1;

  real_t ubvec[1];
  ubvec[0] = 1.05;

  if (refine) {
    /* Refine an existing partition, uncouple as we do not have the cells
     * present on their expected ranks. */
    options[3] = PARMETIS_PSR_UNCOUPLED;

    if (ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                               &wgtflag, &numflag, &ncon, &nparts, tpwgts,
                               ubvec, options, &edgecut, regionid,
                               &comm) != METIS_OK)
      error("Call to ParMETIS_V3_RefineKway failed.");

  } else {

    /* Create a new partition. */
    if (ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                             &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
                             options, &edgecut, regionid, &comm) != METIS_OK)
      error("Call to ParMETIS_V3_PartKway failed.");

    /* We need the existing partition so we can check for permutations. */
    int nsum = 0;
    for (int i = 0; i < s->nr_cells; i++) {
      celllist[i] = s->cells_top[i].nodeID;
      nsum += celllist[i];
    }

    /* If no previous partition all nodeIDs will be set to 0. */
    if (nsum == 0) permute = 0;
  }

  /* Gather the regionids from the other ranks. XXX async version XXX */
  int *newcelllist = NULL;
  if ((newcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
    error("Failed to allocate new celllist");
  if (nodeID == 0) {

    idx_t *remoteids = NULL;
    size_t sizeids = 0;
    for (int rank = 0, j = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      if (rank == 0) {

        /* Locals. */
        for (int k = 0; k < nvt; k++) newcelllist[k] = regionid[k];

      } else {

        /* Type mismatch, so we need to buffer. */
        if (sizeof(idx_t) * nvt > sizeids) {
          sizeids = sizeof(idx_t) * nvt;
          if (sizeids > 0) free(remoteids);
          if ((remoteids = (idx_t *)malloc(sizeids)) == NULL)
            error("Failed to (re)allocate remoteids array");
        }
        res = MPI_Recv((void *)remoteids, nvt, IDX_T, rank, 1, comm, &status);
        if (res != MPI_SUCCESS)
          mpi_error(res, "Failed to receive new regionids");
        for (int k = 0; k < nvt; k++) newcelllist[j + k] = remoteids[k];
      }

      j += nvt;
    }
    if (remoteids != NULL) free(remoteids);

    /* Check that the regionids are ok. */
    for (int k = 0; k < ncells; k++)
      if (newcelllist[k] < 0 || newcelllist[k] >= nregions)
        error("Got bad nodeID %" PRIDX " for cell %i.", newcelllist[k], k);

  } else {
    res = MPI_Send(regionid, vtxdist[nodeID + 1] - vtxdist[nodeID], IDX_T, 0, 1,
                   comm);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to send new regionids");
  }

  /* And everyone gets a copy. */
  res = MPI_Bcast(newcelllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");

  /* Now check the similarity to the old partition and permute if necessary.
   * Checks show that refinement can return a permutation of the partition, we
   * need to check that and correct as necessary. */
  if (permute) {
    int *permcelllist = NULL;
    if ((permcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
      error("Failed to allocate perm celllist array");
    permute_regions(newcelllist, celllist, nregions, ncells, permcelllist);

    /* And keep. */
    memcpy(celllist, permcelllist, sizeof(int) * ncells);
    free(permcelllist);
  } else {
    memcpy(celllist, newcelllist, sizeof(int) * ncells);
  }

  /* Check that the regionids are ok. */
  for (int k = 0; k < ncells; k++)
    if (celllist[k] < 0 || celllist[k] >= nregions)
      error("Got bad nodeID %" PRIDX " for cell %i.", celllist[k], k);

  /* Clean up. */
  free(newcelllist);
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(vtxdist);
  free(tpwgts);
  free(xadj);
  free(adjncy);
  free(regionid);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
/**
 * @brief Repartition the cells amongst the nodes using weights of
 *        various kinds.
 *
 * @param bothweights whether vertex and edge weights will be used, otherwise
 *                    only edge weights will be used.
 * @param timebins use timebins as edge weights.
 * @param repartition the partition struct of the local engine.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
static void repart_edge_parmetis(int bothweights, int timebins,
                                 struct repartition *repartition, int nodeID,
                                 int nr_nodes, struct space *s,
                                 struct task *tasks, int nr_tasks) {

  /* Create weight arrays using task ticks for vertices and edges (edges
   * assume the same graph structure as used in the part_ calls). */
  int nr_cells = s->nr_cells;
  struct cell *cells = s->cells_top;

  /* Allocate and fill the adjncy indexing array defining the graph of
   * cells. */
  idx_t *inds;
  if ((inds = (idx_t *)malloc(sizeof(idx_t) * 26 * nr_cells)) == NULL)
    error("Failed to allocate the inds array");
  graph_init_adjncy(s, inds);

  /* Allocate and init weights. */
  double *weights_v = NULL;
  double *weights_e = NULL;
  if (bothweights) {
    if ((weights_v = (double *)malloc(sizeof(double) * nr_cells)) == NULL)
      error("Failed to allocate vertex weights arrays.");
    bzero(weights_v, sizeof(double) * nr_cells);
  }
  if ((weights_e = (double *)malloc(sizeof(double) * 26 * nr_cells)) == NULL)
    error("Failed to allocate edge weights arrays.");
  bzero(weights_e, sizeof(double) * 26 * nr_cells);

  /* Loop over the tasks... */
  for (int j = 0; j < nr_tasks; j++) {
    /* Get a pointer to the kth task. */
    struct task *t = &tasks[j];

    /* Skip un-interesting tasks. */
    if (t->cost == 0) continue;

    /* Get the task weight based on costs. */
    double w = (double)t->cost;

    /* Get the top-level cells involved. */
    struct cell *ci, *cj;
    for (ci = t->ci; ci->parent != NULL; ci = ci->parent)
      ;
    if (t->cj != NULL)
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent)
        ;
    else
      cj = NULL;

    /* Get the cell IDs. */
    int cid = ci - cells;

    /* Different weights for different tasks. */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_ghost || t->type == task_type_extra_ghost ||
        t->type == task_type_kick1 || t->type == task_type_kick2 ||
        t->type == task_type_end_force || t->type == task_type_cooling ||
        t->type == task_type_timestep || t->type == task_type_init_grav ||
        t->type == task_type_grav_down ||
        t->type == task_type_grav_long_range) {

      /* Particle updates add only to vertex weight. */
      if (bothweights) weights_v[cid] += w;
    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub_self && cj == NULL &&
              ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (bothweights) weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair || (t->type == task_type_sub_pair)) {
      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (bothweights) weights_v[cid] += w;

      }

      /* Distinct cells. */
      else {
        /* Index of the jth cell. */
        int cjd = cj - cells;

        /* Local cells add weight to vertices. */
        if (bothweights && ci->nodeID == nodeID) {
          weights_v[cid] += 0.5 * w;
          if (cj->nodeID == nodeID) weights_v[cjd] += 0.5 * w;
        }

        /* Find indices of ci/cj neighbours. Note with gravity these cells may
         * not be neighbours, in that case we ignore any edge weight for that
         * pair. */
        int ik = -1;
        for (int k = 26 * cid; k < 26 * nr_cells; k++) {
          if (inds[k] == cjd) {
            ik = k;
            break;
          }
        }

        /* cj */
        int jk = -1;
        for (int k = 26 * cjd; k < 26 * nr_cells; k++) {
          if (inds[k] == cid) {
            jk = k;
            break;
          }
        }
        if (ik != -1 && jk != -1) {

          if (timebins) {
            /* Add weights to edge for all cells based on the expected
             * interaction time (calculated as the time to the last expected
             * time) as we want to avoid having active cells on the edges, so
             * we cut for that. Note that weight is added to the local and
             * remote cells, as we want to keep both away from any cuts, this
             * can overflow int, so take care. */
            int dti = num_time_bins - get_time_bin(ci->ti_hydro_end_min);
            int dtj = num_time_bins - get_time_bin(cj->ti_hydro_end_min);
            double dt = (double)(1 << dti) + (double)(1 << dtj);
            weights_e[ik] += dt;
            weights_e[jk] += dt;

          } else {

            /* Add weights from task costs to the edge. */
            weights_e[ik] += w;
            weights_e[jk] += w;
          }
        }
      }
    }
  }

  /* Merge the weights arrays across all nodes. */
  int res;
  if (bothweights) {
    res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_v, weights_v,
                     nr_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to allreduce vertex weights.");
  }

  res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_e, weights_e,
                   26 * nr_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to allreduce edge weights.");

  /* Allocate cell list for the partition. If not already done. */
  int refine = 1;
  if (repartition->ncelllist != s->nr_cells) {
    refine = 0;
    free(repartition->celllist);
    repartition->ncelllist = 0;
    if ((repartition->celllist = (int *)malloc(sizeof(int) * s->nr_cells)) ==
        NULL)
      error("Failed to allocate celllist");
    repartition->ncelllist = s->nr_cells;
  }

  /* We need to rescale the weights into the range of an integer for
   * ParMETIS (that is the range of idx_t). Also we would like the range of
   * vertex and edges weights to be similar so they balance. */
  double wminv = 0.0;
  double wmaxv = 0.0;
  if (bothweights) {
    wminv = weights_v[0];
    wmaxv = weights_v[0];
    for (int k = 0; k < nr_cells; k++) {
      wmaxv = weights_v[k] > wmaxv ? weights_v[k] : wmaxv;
      wminv = weights_v[k] < wminv ? weights_v[k] : wminv;
    }
  }

  double wmine = weights_e[0];
  double wmaxe = weights_e[0];
  for (int k = 0; k < 26 * nr_cells; k++) {
    wmaxe = weights_e[k] > wmaxe ? weights_e[k] : wmaxe;
    wmine = weights_e[k] < wmine ? weights_e[k] : wmine;
  }

  if (bothweights) {

    /* Make range the same in both weights systems. */
    if ((wmaxv - wminv) > (wmaxe - wmine)) {
      double wscale = 1.0;
      if ((wmaxe - wmine) > 0.0) {
        wscale = (wmaxv - wminv) / (wmaxe - wmine);
      }
      for (int k = 0; k < 26 * nr_cells; k++) {
        weights_e[k] = (weights_e[k] - wmine) * wscale + wminv;
      }
      wmine = wminv;
      wmaxe = wmaxv;

    } else {
      double wscale = 1.0;
      if ((wmaxv - wminv) > 0.0) {
        wscale = (wmaxe - wmine) / (wmaxv - wminv);
      }
      for (int k = 0; k < nr_cells; k++) {
        weights_v[k] = (weights_v[k] - wminv) * wscale + wmine;
      }
      wminv = wmine;
      wmaxv = wmaxe;
    }

    /* Scale to the ParMETIS range. */
    double wscale = 1.0;
    if ((wmaxv - wminv) > 0.0) {
      wscale = (parmetis_maxweight - 1.0) / (wmaxv - wminv);
    }
    for (int k = 0; k < nr_cells; k++) {
      weights_v[k] = (weights_v[k] - wminv) * wscale + 1.0;
    }
  }

  /* Scale to the ParMETIS range. */
  double wscale = 1.0;
  if ((wmaxe - wmine) > 0.0) {
    wscale = (parmetis_maxweight - 1.0) / (wmaxe - wmine);
  }
  for (int k = 0; k < 26 * nr_cells; k++) {
    weights_e[k] = (weights_e[k] - wmine) * wscale + 1.0;
  }

  /* And partition, use both weights or not as requested. */
  if (bothweights)
    pick_parmetis(nodeID, s, nr_nodes, weights_v, weights_e, refine, -1,
                  repartition->celllist);
  else
    pick_parmetis(nodeID, s, nr_nodes, NULL, weights_e, refine, -1,
                  repartition->celllist);

  /* Check that all cells have good values. All nodes have same copy, so just
   * check on one. */
  if (nodeID == 0) {
    for (int k = 0; k < nr_cells; k++)
      if (repartition->celllist[k] < 0 || repartition->celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", repartition->celllist[k], k);
  }

  /* Check that the partition is complete and all nodes have some work. */
  int present[nr_nodes];
  int failed = 0;
  for (int i = 0; i < nr_nodes; i++) present[i] = 0;
  for (int i = 0; i < nr_cells; i++) present[repartition->celllist[i]]++;
  for (int i = 0; i < nr_nodes; i++) {
    if (!present[i]) {
      failed = 1;
      if (nodeID == 0) message("Node %d is not present after repartition", i);
    }
  }

  /* If partition failed continue with the current one, but make this clear. */
  if (failed) {
    if (nodeID == 0)
      message(
          "WARNING: ParMETIS repartition has failed, "
          "continuing with the current partition, "
          "load balance will not be optimal");
    for (int k = 0; k < nr_cells; k++)
      repartition->celllist[k] = cells[k].nodeID;
  }

  /* And apply to our cells */
  split_metis(s, nr_nodes, repartition->celllist);

  /* Clean up. */
  free(inds);
  if (bothweights) free(weights_v);
  free(weights_e);
}
#endif

/**
 * @brief Repartition the space using the given repartition type.
 *
 * Note that at the end of this process all the cells will be re-distributed
 * across the nodes, but the particles themselves will not be.
 *
 * @param reparttype #repartition struct
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
void partition_repartition(struct repartition *reparttype, int nodeID,
                           int nr_nodes, struct space *s, struct task *tasks,
                           int nr_tasks) {

#if defined(WITH_MPI) && defined(HAVE_PARMETIS)

  if (reparttype->type == REPART_PARMETIS_VERTEX_EDGE_COSTS) {
    repart_edge_parmetis(0, 1, reparttype, nodeID, nr_nodes, s, tasks,
                         nr_tasks);

  } else if (reparttype->type == REPART_PARMETIS_EDGE_COSTS) {
    repart_edge_parmetis(0, 0, reparttype, nodeID, nr_nodes, s, tasks,
                         nr_tasks);

  } else if (reparttype->type == REPART_PARMETIS_VERTEX_COSTS_TIMEBINS) {
    repart_edge_parmetis(0, 1, reparttype, nodeID, nr_nodes, s, tasks,
                         nr_tasks);

  } else if (reparttype->type == REPART_NONE) {
    /* Doing nothing. */

  } else {
    error("Impossible repartition type");
  }
#else
  error("SWIFT was not compiled with ParMETIS support.");
#endif
}

/**
 * @brief Initial partition of space cells.
 *
 * Cells are assigned to a node on the basis of various schemes, all of which
 * should attempt to distribute them in geometrically close regions to
 * minimise the movement of particles.
 *
 * Note that the partition type is a suggestion and will be ignored if that
 * scheme fails. In that case we fallback to a vectorised scheme, that is
 * guaranteed to work provided we have more cells than nodes.
 *
 * @param initial_partition the type of partitioning to try.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells.
 */
void partition_initial_partition(struct partition *initial_partition,
                                 int nodeID, int nr_nodes, struct space *s) {

  /* Geometric grid partitioning. */
  if (initial_partition->type == INITPART_GRID) {
    int j, k;
    int ind[3];
    struct cell *c;

    /* If we've got the wrong number of nodes, fail. */
    if (nr_nodes !=
        initial_partition->grid[0] * initial_partition->grid[1] *
            initial_partition->grid[2])
      error("Grid size does not match number of nodes.");

    /* Run through the cells and set their nodeID. */
    // message("s->dim = [%e,%e,%e]", s->dim[0], s->dim[1], s->dim[2]);
    for (k = 0; k < s->nr_cells; k++) {
      c = &s->cells_top[k];
      for (j = 0; j < 3; j++)
        ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
      c->nodeID = ind[0] +
                  initial_partition->grid[0] *
                      (ind[1] + initial_partition->grid[1] * ind[2]);
      // message("cell at [%e,%e,%e]: ind = [%i,%i,%i], nodeID = %i", c->loc[0],
      // c->loc[1], c->loc[2], ind[0], ind[1], ind[2], c->nodeID);
    }

    /* The grid technique can fail, so check for this before proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("Grid initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
      return;
    }

  } else if (initial_partition->type == INITPART_PARMETIS_WEIGHT ||
             initial_partition->type == INITPART_PARMETIS_NOWEIGHT) {
#if defined(WITH_MPI) && defined(HAVE_PARMETIS)
    /* Simple k-way partition selected by ParMETIS using cell particle counts
     * as weights or not. Should be best when starting with a inhomogeneous
     * dist.
     */

    /* Space for particles per cell counts, which will be used as weights or
     * not. */
    double *weights = NULL;
    if (initial_partition->type == INITPART_PARMETIS_WEIGHT) {
      if ((weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate weights buffer.");
      bzero(weights, sizeof(double) * s->nr_cells);

      /* Check each particle and accumulate the counts per cell. */
      accumulate_counts(s, weights);

      /* Get all the counts from all the nodes. */
      if (MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
        error("Failed to allreduce particle cell weights.");
    }

    /* Do the calculation. */
    int *celllist = NULL;
    if ((celllist = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate celllist");
    pick_parmetis(nodeID, s, nr_nodes, weights, NULL, 0, -1, celllist);

    /* And apply to our cells */
    split_metis(s, nr_nodes, celllist);

    /* It's not known if this can fail, but check for this before
     * proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("METIS initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
    }

    if (weights != NULL) free(weights);
    free(celllist);
#else
    error("SWIFT was not compiled with METIS support");
#endif

  } else if (initial_partition->type == INITPART_VECTORIZE) {

#if defined(WITH_MPI)
    /* Vectorised selection, guaranteed to work for samples less than the
     * number of cells, but not very clumpy in the selection of regions. */
    int *samplecells = NULL;
    if ((samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3)) == NULL)
      error("Failed to allocate samplecells");

    if (nodeID == 0) {
      pick_vector(s, nr_nodes, samplecells);
    }

    /* Share the samplecells around all the nodes. */
    int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the partition sample cells.");

    /* And apply to our cells */
    split_vector(s, nr_nodes, samplecells);
    free(samplecells);
#else
    error("SWIFT was not compiled with MPI support");
#endif
  }
}

/**
 * @brief Initialises the partition and re-partition scheme from the parameter
 *        file
 *
 * @param partition The #partition scheme to initialise.
 * @param repartition The #repartition scheme to initialise.
 * @param params The parsed parameter file.
 * @param nr_nodes The number of MPI nodes we are running on.
 */
void partition_init(struct partition *partition,
                    struct repartition *repartition,
                    const struct swift_params *params, int nr_nodes) {

#ifdef WITH_MPI

/* Defaults make use of METIS if available */
#ifdef HAVE_PARMETIS
  const char *default_repart = "costs/costs";
  const char *default_part = "regions";
#else
  const char *default_repart = "none/none";
  const char *default_part = "grid";
#endif

  /* Set a default grid so that grid[0]*grid[1]*grid[2] == nr_nodes. */
  factor(nr_nodes, &partition->grid[0], &partition->grid[1]);
  factor(nr_nodes / partition->grid[1], &partition->grid[0],
         &partition->grid[2]);
  factor(partition->grid[0] * partition->grid[1], &partition->grid[1],
         &partition->grid[0]);

  /* Now let's check what the user wants as an initial domain. */
  char part_type[20];
  parser_get_opt_param_string(params, "DomainDecomposition:initial_type",
                              part_type, default_part);
  switch (part_type[0]) {
    case 'g':
      partition->type = INITPART_GRID;
      break;
    case 'v':
      partition->type = INITPART_VECTORIZE;
      break;
#ifdef HAVE_PARMETIS
    case 'r':
      partition->type = INITPART_PARMETIS_NOWEIGHT;
      break;
    case 'b':
      partition->type = INITPART_PARMETIS_WEIGHT;
      break;
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid', 'regions', 'balanced' or "
          "'vectorized'");
#else
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid' or 'vectorized' when compiled "
          "without ParMETIS.");
#endif
  }

  /* In case of grid, read more parameters */
  if (part_type[0] == 'g') {
    partition->grid[0] = parser_get_opt_param_int(
        params, "DomainDecomposition:initial_grid_x", partition->grid[0]);
    partition->grid[1] = parser_get_opt_param_int(
        params, "DomainDecomposition:initial_grid_y", partition->grid[1]);
    partition->grid[2] = parser_get_opt_param_int(
        params, "DomainDecomposition:initial_grid_z", partition->grid[2]);
  }

  /* Now let's check what the user wants as a repartition strategy */
  parser_get_opt_param_string(params, "DomainDecomposition:repartition_type",
                              part_type, default_repart);

  if (strcmp("none/none", part_type) == 0) {
    repartition->type = REPART_NONE;

#ifdef HAVE_PARMETIS
  } else if (strcmp("costs/costs", part_type) == 0) {
    repartition->type = REPART_PARMETIS_VERTEX_EDGE_COSTS;

  } else if (strcmp("none/costs", part_type) == 0) {
    repartition->type = REPART_PARMETIS_EDGE_COSTS;

  } else if (strcmp("costs/time", part_type) == 0) {
    repartition->type = REPART_PARMETIS_VERTEX_COSTS_TIMEBINS;

  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none/none', 'costs/costs', 'none/costs' "
        "or 'costs/time'");
#else
  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error("Permitted values are: 'none/none' when compiled without ParMETIS.");
#endif
  }

  /* Get the fraction CPU time difference between nodes (<1) or the number
   * of steps between repartitions (>1). */
  repartition->trigger =
      parser_get_opt_param_float(params, "DomainDecomposition:trigger", 0.05f);
  if (repartition->trigger <= 0)
    error("Invalid DomainDecomposition:trigger, must be greater than zero");
  if (repartition->trigger < 2 && repartition->trigger >= 1)
    error(
        "Invalid DomainDecomposition:trigger, must be 2 or greater or less"
        " than 1");

  /* Fraction of particles that should be updated before a repartition
   * based on CPU time is considered. */
  repartition->minfrac =
      parser_get_opt_param_float(params, "DomainDecomposition:minfrac", 0.9f);
  if (repartition->minfrac <= 0 || repartition->minfrac > 1)
    error(
        "Invalid DomainDecomposition:minfrac, must be greater than 0 and less "
        "than equal to 1");

  /* Clear the celllist for use. */
  repartition->ncelllist = 0;
  repartition->celllist = NULL;

#else
  error("SWIFT was not compiled with MPI support");
#endif
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
static int check_complete(struct space *s, int verbose, int nregions) {

  int *present = NULL;
  if ((present = (int *)malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate present array");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    if (s->cells_top[i].nodeID <= nregions)
      present[s->cells_top[i].nodeID]++;
    else
      message("Bad nodeID: s->cells_top[%d].nodeID = %d", i,
              s->cells_top[i].nodeID);
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

/**
 * @brief Partition a space of cells based on another space of cells.
 *
 * The two spaces are expected to be at different cell sizes, so what we'd
 * like to do is assign the second space to geometrically closest nodes
 * of the first, with the effect of minimizing particle movement when
 * rebuilding the second space from the first.
 *
 * Since two spaces cannot exist simultaneously the old space is actually
 * required in a decomposed state. These are the old cells sizes and counts
 * per dimension, along with a list of the old nodeIDs. The old nodeIDs are
 * indexed by the cellid (see cell_getid()), so should be stored that way.
 *
 * On exit the new space cells will have their nodeIDs assigned.
 *
 * @param oldh the cell dimensions of old space.
 * @param oldcdim number of cells per dimension in old space.
 * @param oldnodeIDs the nodeIDs of cells in the old space, indexed by old
 *cellid.
 * @param s the space to be partitioned.
 *
 * @return 1 if the new space contains nodeIDs from all nodes, 0 otherwise.
 */
int partition_space_to_space(double *oldh, double *oldcdim, int *oldnodeIDs,
                             struct space *s) {

  /* Loop over all the new cells. */
  int nr_nodes = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {

        /* Scale indices to old cell space. */
        const int ii = rint(i * s->iwidth[0] * oldh[0]);
        const int jj = rint(j * s->iwidth[1] * oldh[1]);
        const int kk = rint(k * s->iwidth[2] * oldh[2]);

        const int cid = cell_getid(s->cdim, i, j, k);
        const int oldcid = cell_getid(oldcdim, ii, jj, kk);
        s->cells_top[cid].nodeID = oldnodeIDs[oldcid];

        if (oldnodeIDs[oldcid] > nr_nodes) nr_nodes = oldnodeIDs[oldcid];
      }
    }
  }

  /* Check we have all nodeIDs present in the resample. */
  return check_complete(s, 1, nr_nodes + 1);
}

/**
 * @brief save the nodeIDs of the current top-level cells by adding them to a
 *             repartition struct. Used when restarting application.
 *
 * @param s the space with the top-level cells.
 * @param reparttype struct to update with the a list of nodeIDs.
 *
 */
void partition_store_celllist(struct space *s, struct repartition *reparttype) {
  if (reparttype->ncelllist != s->nr_cells) {
    free(reparttype->celllist);
    if ((reparttype->celllist = (int *)malloc(sizeof(int) * s->nr_cells)) ==
        NULL)
      error("Failed to allocate celllist");
    reparttype->ncelllist = s->nr_cells;
  }
  for (int i = 0; i < s->nr_cells; i++) {
    reparttype->celllist[i] = s->cells_top[i].nodeID;
  }
}

/**
 * @brief restore the saved list of nodeIDs by applying them to the
 *        top-level cells of a space. Used when restarting application.
 *
 * @param s the space with the top-level cells.
 * @param reparttype struct with the list of nodeIDs saved,
 *
 */
void partition_restore_celllist(struct space *s,
                                struct repartition *reparttype) {
  if (reparttype->ncelllist > 0) {
    if (reparttype->ncelllist == s->nr_cells) {
      for (int i = 0; i < s->nr_cells; i++) {
        s->cells_top[i].nodeID = reparttype->celllist[i];
      }
      if (!check_complete(s, 1, s->e->nr_nodes)) {
        error("Not all ranks are present in the restored partition");
      }
    } else {
      error(
          "Cannot apply the saved partition celllist as the "
          "number of top-level cells (%d) is different to the "
          "saved number (%d)",
          s->nr_cells, reparttype->ncelllist);
    }
  }
}

/**
 * @brief Write a repartition struct to the given FILE as a stream of bytes.
 *
 * @param reparttype the struct
 * @param stream the file stream
 */
void partition_struct_dump(struct repartition *reparttype, FILE *stream) {
  restart_write_blocks(reparttype, sizeof(struct repartition), 1, stream,
                       "repartition", "repartition params");

  /* Also save the celllist, if we have one. */
  if (reparttype->ncelllist > 0)
    restart_write_blocks(reparttype->celllist,
                         sizeof(int) * reparttype->ncelllist, 1, stream,
                         "celllist", "repartition celllist");
}

/**
 * @brief Restore a repartition struct from the given FILE as a stream of
 * bytes.
 *
 * @param reparttype the struct
 * @param stream the file stream
 */
void partition_struct_restore(struct repartition *reparttype, FILE *stream) {
  restart_read_blocks(reparttype, sizeof(struct repartition), 1, stream, NULL,
                      "repartition params");

  /* Also restore the celllist, if we have one. */
  if (reparttype->ncelllist > 0) {
    if ((reparttype->celllist = malloc(sizeof(int) * reparttype->ncelllist)) ==
        NULL)
      error("Failed to allocate celllist");
    restart_read_blocks(reparttype->celllist,
                        sizeof(int) * reparttype->ncelllist, 1, stream, NULL,
                        "repartition celllist");
  }
}
