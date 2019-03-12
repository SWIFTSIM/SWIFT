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
 *  Currently supported partitioning types: grid, vectorise and METIS/ParMETIS.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/* Include int min and max values. Define these limits in C++ as well. */
#define __STDC_LIMIT_MACROS
#include <stdint.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "debug.h"
#include "engine.h"
#include "error.h"
#include "partition.h"
#include "restart.h"
#include "space.h"
#include "tools.h"

/* Simple descriptions of initial partition types for reports. */
const char *initial_partition_name[] = {
    "axis aligned grids of cells", "vectorized point associated cells",
    "memory balanced, using particle weighted cells",
    "similar sized regions, using unweighted cells"};

/* Simple descriptions of repartition types for reports. */
const char *repartition_name[] = {
    "none", "edge and vertex task cost weights", "task cost edge weights",
    "memory balanced, using particle vertex weights",
    "vertex task costs and edge delta timebin weights"};

/* Local functions, if needed. */
static int check_complete(struct space *s, int verbose, int nregions);

/*
 * Repartition fixed costs per type/subtype. These are determined from the
 * statistics output produced when running with task debugging enabled.
 */
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
static double repartition_costs[task_type_count][task_subtype_count];
#endif
#if defined(WITH_MPI)
static int repart_init_fixed_costs(void);
#endif

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

  /* METIS/ParMETIS support (optional)
   * =================================
   *
   * METIS/ParMETIS partitions using a multi-level k-way scheme. We support
   * using this in a unweighted scheme, which works well and seems to be
   * guaranteed, and a weighted by the number of particles scheme.
   *
   * Repartitioning is based on ParMETIS and uses weights determined from the
   * estimated costs that a cells tasks will take or the relative time bins of
   * the cells next updates.
   */

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Fill the adjncy array defining the graph of cells in a space.
 *
 * See the ParMETIS and METIS manuals if you want to understand this
 * format. The cell graph consists of all nodes as vertices with edges as the
 * connections to all neighbours, so we have 26 per vertex. Note you will
 * also need an xadj array, for METIS that would be:
 *
 *   xadj[0] = 0;
 *   for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
 *
 * but each rank needs a different xadj when using ParMETIS.
 *
 * @param s the space of cells.
 * @param adjncy the adjncy array to fill, must be of size 26 * the number of
 *               cells in the space.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 */
static void graph_init(struct space *s, idx_t *adjncy, idx_t *xadj) {

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

  /* If given set METIS xadj. */
  if (xadj != NULL) {
    xadj[0] = 0;
    for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
  }
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
struct counts_mapper_data {
  double *counts;
  size_t size;
  struct space *s;
};

/* Generic function for accumulating sized counts for TYPE parts. Note uses
 * local memory to reduce contention, the amount of memory required is
 * precalculated by an additional loop determining the range of cell IDs. */
#define ACCUMULATE_SIZES_MAPPER(TYPE)                                          \
  accumulate_sizes_mapper_##TYPE(void *map_data, int num_elements,             \
                                 void *extra_data) {                           \
    struct TYPE *parts = (struct TYPE *)map_data;                              \
    struct counts_mapper_data *mydata =                                        \
        (struct counts_mapper_data *)extra_data;                               \
    double size = mydata->size;                                                \
    int *cdim = mydata->s->cdim;                                               \
    double iwidth[3] = {mydata->s->iwidth[0], mydata->s->iwidth[1],            \
                        mydata->s->iwidth[2]};                                 \
    double dim[3] = {mydata->s->dim[0], mydata->s->dim[1], mydata->s->dim[2]}; \
    double *lcounts = NULL;                                                    \
    int lcid = mydata->s->nr_cells;                                            \
    int ucid = 0;                                                              \
    for (int k = 0; k < num_elements; k++) {                                   \
      for (int j = 0; j < 3; j++) {                                            \
        if (parts[k].x[j] < 0.0)                                               \
          parts[k].x[j] += dim[j];                                             \
        else if (parts[k].x[j] >= dim[j])                                      \
          parts[k].x[j] -= dim[j];                                             \
      }                                                                        \
      const int cid =                                                          \
          cell_getid(cdim, parts[k].x[0] * iwidth[0],                          \
                     parts[k].x[1] * iwidth[1], parts[k].x[2] * iwidth[2]);    \
      if (cid > ucid) ucid = cid;                                              \
      if (cid < lcid) lcid = cid;                                              \
    }                                                                          \
    int nused = ucid - lcid + 1;                                               \
    if ((lcounts = (double *)calloc(sizeof(double), nused)) == NULL)           \
      error("Failed to allocate counts thread-specific buffer");               \
    for (int k = 0; k < num_elements; k++) {                                   \
      const int cid =                                                          \
          cell_getid(cdim, parts[k].x[0] * iwidth[0],                          \
                     parts[k].x[1] * iwidth[1], parts[k].x[2] * iwidth[2]);    \
      lcounts[cid - lcid] += size;                                             \
    }                                                                          \
    for (int k = 0; k < nused; k++)                                            \
      atomic_add_d(&mydata->counts[k + lcid], lcounts[k]);                     \
    free(lcounts);                                                             \
  }

/**
 * @brief Accumulate the sized counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * part version.
 */
static void ACCUMULATE_SIZES_MAPPER(part);

/**
 * @brief Accumulate the sized counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * gpart version.
 */
static void ACCUMULATE_SIZES_MAPPER(gpart);

/**
 * @brief Accumulate the sized counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * spart version.
 */
static void ACCUMULATE_SIZES_MAPPER(spart);

/**
 * @brief Accumulate total memory size in particles per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of bytes in particles per cell. Should be
 *               allocated as size s->nr_cells.
 */
static void accumulate_sizes(struct space *s, double *counts) {

  bzero(counts, sizeof(double) * s->nr_cells);

  struct counts_mapper_data mapper_data;
  mapper_data.counts = counts;
  mapper_data.s = s;

  double hsize = (double)sizeof(struct part);
  if (s->nr_parts > 0) {
    mapper_data.size = hsize;
    threadpool_map(&s->e->threadpool, accumulate_sizes_mapper_part, s->parts,
                   s->nr_parts, sizeof(struct part), space_splitsize,
                   &mapper_data);
  }

  double gsize = (double)sizeof(struct gpart);
  if (s->nr_gparts > 0) {
    mapper_data.size = gsize;
    threadpool_map(&s->e->threadpool, accumulate_sizes_mapper_gpart, s->gparts,
                   s->nr_gparts, sizeof(struct gpart), space_splitsize,
                   &mapper_data);
  }

  double ssize = (double)sizeof(struct spart);
  if (s->nr_sparts > 0) {
    mapper_data.size = ssize;
    threadpool_map(&s->e->threadpool, accumulate_sizes_mapper_spart, s->sparts,
                   s->nr_sparts, sizeof(struct spart), space_splitsize,
                   &mapper_data);
  }

  /* Keep the sum of particles across all ranks in the range of IDX_MAX. */
  if ((s->e->total_nr_parts * hsize + s->e->total_nr_gparts * gsize +
       s->e->total_nr_sparts * ssize) > (double)IDX_MAX) {
    double vscale =
        (double)(IDX_MAX - 1000) /
        (double)(s->e->total_nr_parts * hsize + s->e->total_nr_gparts * gsize +
                 s->e->total_nr_sparts * ssize);
    for (int k = 0; k < s->nr_cells; k++) counts[k] *= vscale;
  }
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Apply METIS cell-list partitioning to a cell structure.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
static void split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) s->cells_top[i].nodeID = celllist[i];

  /* To check or visualise the partition dump all the cells. */
  /*dumpCellRanks("metis_partition", s->cells_top, s->nr_cells);*/
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/* qsort support. */
struct indexval {
  int index;
  int count;
  int old_val;
  int new_val;
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
    ivs[index].old_val = oldlist[k];
    ivs[index].new_val = newlist[k];
  }
  qsort(ivs, indmax, sizeof(struct indexval), indexvalcmp);

  /* Go through the ivs using the largest counts first, these are the
   * regions with the most cells in common, old partition to new. If not
   * returning the permutation, avoid the associated work. */
  int *oldmap = NULL;
  int *newmap = NULL;
  oldmap = permlist; /* Reuse this */
  if ((newmap = (int *)malloc(sizeof(int) * nregions)) == NULL)
    error("Failed to allocate newmap array");

  for (int k = 0; k < nregions; k++) {
    oldmap[k] = -1;
    newmap[k] = -1;
  }

  for (int k = 0; k < indmax; k++) {

    /* Stop when all regions with common cells have been considered. */
    if (ivs[k].count == 0) break;

    /* Store old and new IDs, if not already used. */
    if (newmap[ivs[k].new_val] == -1 && oldmap[ivs[k].old_val] == -1) {
      newmap[ivs[k].new_val] = ivs[k].old_val;
      oldmap[ivs[k].old_val] = ivs[k].new_val;
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
 * @param adaptive whether to use an adaptive reparitition of an existing
 *        partition or simple refinement. Adaptive repartition is controlled
 *        by the itr parameter.
 * @param itr the ratio of inter-process communication time to data
 *            redistribution time. Used to weight repartitioning edge cuts
 *            when refine and adaptive are true.
 * @param celllist on exit this contains the ids of the selected regions,
 *        size of number of cells. If refine is 1, then this should contain
 *        the old partition on entry.
 */
static void pick_parmetis(int nodeID, struct space *s, int nregions,
                          double *vertexw, double *edgew, int refine,
                          int adaptive, float itr, int *celllist) {

  int res;
  MPI_Comm comm;
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

  /* Prepare MPI requests for the asynchronous communications */
  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 5 * nregions)) ==
      NULL)
    error("Failed to allocate MPI request list.");
  for (int k = 0; k < 5 * nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  MPI_Status *stats;
  if ((stats = (MPI_Status *)malloc(sizeof(MPI_Status) * 5 * nregions)) == NULL)
    error("Failed to allocate MPI status list.");

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
    graph_init(s, full_adjncy, NULL);

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

    /* Dump graphs to disk files for testing. ParMETIS xadj isn't right for
     * a dump, so make a serial-like version. */
    /*{
      idx_t *tmp_xadj =
          (idx_t *)malloc(sizeof(idx_t) * (ncells + nregions + 1));
      tmp_xadj[0] = 0;
      for (int k = 0; k < ncells; k++) tmp_xadj[k + 1] = tmp_xadj[k] + 26;
      dumpMETISGraph("parmetis_graph", ncells, 1, tmp_xadj, full_adjncy,
                     full_weights_v, NULL, full_weights_e);
      free(tmp_xadj);
      }*/

    /* Send ranges to the other ranks and keep our own. */
    for (int rank = 0, j1 = 0, j2 = 0, j3 = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      if (refine)
        for (int i = 0; i < nvt; i++) full_regionid[j3 + i] = celllist[j3 + i];

      if (rank == 0) {
        memcpy(xadj, &full_xadj[j1], sizeof(idx_t) * (nvt + 1));
        memcpy(adjncy, &full_adjncy[j2], sizeof(idx_t) * nvt * 26);
        if (weights_e != NULL)
          memcpy(weights_e, &full_weights_e[j2], sizeof(idx_t) * nvt * 26);
        if (weights_v != NULL)
          memcpy(weights_v, &full_weights_v[j3], sizeof(idx_t) * nvt);
        if (refine) memcpy(regionid, full_regionid, sizeof(idx_t) * nvt);

      } else {
        res = MPI_Isend(&full_xadj[j1], nvt + 1, IDX_T, rank, 0, comm,
                        &reqs[5 * rank + 0]);
        if (res == MPI_SUCCESS)
          res = MPI_Isend(&full_adjncy[j2], nvt * 26, IDX_T, rank, 1, comm,
                          &reqs[5 * rank + 1]);
        if (res == MPI_SUCCESS && weights_e != NULL)
          res = MPI_Isend(&full_weights_e[j2], nvt * 26, IDX_T, rank, 2, comm,
                          &reqs[5 * rank + 2]);
        if (res == MPI_SUCCESS && weights_v != NULL)
          res = MPI_Isend(&full_weights_v[j3], nvt, IDX_T, rank, 3, comm,
                          &reqs[5 * rank + 3]);
        if (refine && res == MPI_SUCCESS)
          res = MPI_Isend(&full_regionid[j3], nvt, IDX_T, rank, 4, comm,
                          &reqs[5 * rank + 4]);
        if (res != MPI_SUCCESS) mpi_error(res, "Failed to send graph data");
      }
      j1 += nvt + 1;
      j2 += nvt * 26;
      j3 += nvt;
    }

    /* Wait for all sends to complete. */
    int result;
    if ((result = MPI_Waitall(5 * nregions, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5 * nregions; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &result);
        message("send request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall sending repartition data.");
    }

    /* Clean up. */
    if (weights_v != NULL) free(full_weights_v);
    if (weights_e != NULL) free(full_weights_e);
    free(full_xadj);
    free(full_adjncy);
    if (refine) free(full_regionid);

  } else {

    /* Receive stuff from rank 0. */
    res = MPI_Irecv(xadj, nverts + 1, IDX_T, 0, 0, comm, &reqs[0]);
    if (res == MPI_SUCCESS)
      res = MPI_Irecv(adjncy, nverts * 26, IDX_T, 0, 1, comm, &reqs[1]);
    if (res == MPI_SUCCESS && weights_e != NULL)
      res = MPI_Irecv(weights_e, nverts * 26, IDX_T, 0, 2, comm, &reqs[2]);
    if (res == MPI_SUCCESS && weights_v != NULL)
      res = MPI_Irecv(weights_v, nverts, IDX_T, 0, 3, comm, &reqs[3]);
    if (refine && res == MPI_SUCCESS)
      res += MPI_Irecv((void *)regionid, nverts, IDX_T, 0, 4, comm, &reqs[4]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive graph data");

    /* Wait for all recvs to complete. */
    int result;
    if ((result = MPI_Waitall(5, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &result);
        message("recv request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall receiving repartition data.");
    }
  }

  /* Set up the tpwgts array. This is just 1/nregions. */
  real_t *tpwgts;
  if ((tpwgts = (real_t *)malloc(sizeof(real_t) * nregions)) == NULL)
    error("Failed to allocate tpwgts array");
  for (int i = 0; i < nregions; i++) tpwgts[i] = 1.0 / (real_t)nregions;

  /* Common parameters. */
  idx_t options[4];
  options[0] = 1;
  options[1] = 0;

  idx_t edgecut;
  idx_t ncon = 1;
  idx_t nparts = nregions;
  idx_t numflag = 0;
  idx_t wgtflag = 0;
  if (edgew != NULL) wgtflag += 1;
  if (vertexw != NULL) wgtflag += 2;

  real_t ubvec[1];
  ubvec[0] = 1.001;

  if (refine) {
    /* Refine an existing partition, uncouple as we do not have the cells
     * present on their expected ranks. */
    options[3] = PARMETIS_PSR_UNCOUPLED;

    /* Seed for randoms. */
    options[2] = clocks_random_seed();

    /* Choice is whether to use an adaptive repartition or a simple
     * refinement. */
    if (adaptive) {

      /* Balance between cuts and movement. */
      real_t itr_real_t = itr;
      if (ParMETIS_V3_AdaptiveRepart(
              vtxdist, xadj, adjncy, weights_v, NULL, weights_e, &wgtflag,
              &numflag, &ncon, &nparts, tpwgts, ubvec, &itr_real_t, options,
              &edgecut, regionid, &comm) != METIS_OK)
        error("Call to ParMETIS_V3_AdaptiveRepart failed.");
    } else {
      if (ParMETIS_V3_RefineKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                                 &wgtflag, &numflag, &ncon, &nparts, tpwgts,
                                 ubvec, options, &edgecut, regionid,
                                 &comm) != METIS_OK)
        error("Call to ParMETIS_V3_RefineKway failed.");
    }
  } else {

    /* Create a new partition. Use a number of guesses as that is similar to
     * the way that serial METIS works (serial METIS usually gives the best
     * quality partitions). */
    idx_t best_edgecut = 0;
    idx_t *best_regionid = NULL;
    if ((best_regionid = (idx_t *)malloc(sizeof(idx_t) * (nverts + 1))) == NULL)
      error("Failed to allocate best_regionid array");

    for (int i = 0; i < 10; i++) {
      options[2] = clocks_random_seed();

      if (ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, weights_v, weights_e,
                               &wgtflag, &numflag, &ncon, &nparts, tpwgts,
                               ubvec, options, &edgecut, regionid,
                               &comm) != METIS_OK)
        error("Call to ParMETIS_V3_PartKway failed.");

      if (i == 0 || (best_edgecut > edgecut)) {
        best_edgecut = edgecut;
        memcpy(best_regionid, regionid, sizeof(idx_t) * (nverts + 1));
      }
    }

    /* Keep the best edgecut. */
    memcpy(regionid, best_regionid, sizeof(idx_t) * (nverts + 1));
    free(best_regionid);
  }

  /* Need to gather all the regionid arrays from the ranks. */
  for (int k = 0; k < nregions; k++) reqs[k] = MPI_REQUEST_NULL;

  if (nodeID != 0) {

    /* Send our regions to node 0. */
    res = MPI_Isend(regionid, vtxdist[nodeID + 1] - vtxdist[nodeID], IDX_T, 0,
                    1, comm, &reqs[0]);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to send new regionids");

    /* Wait for send to complete. */
    int err;
    if ((err = MPI_Wait(reqs, stats)) != MPI_SUCCESS) {
      mpi_error(err, "Failed during wait sending regionids.");
    }

  } else {

    /* Node 0 */
    idx_t *remoteids = NULL;
    if ((remoteids = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate remoteids buffer");

    int nvt = vtxdist[1] - vtxdist[0];
    memcpy(remoteids, regionid, sizeof(idx_t) * nvt);

    /* Receive from other ranks. */
    for (int rank = 1, j = nvt; rank < nregions; rank++) {
      nvt = vtxdist[rank + 1] - vtxdist[rank];
      res = MPI_Irecv((void *)&remoteids[j], nvt, IDX_T, rank, 1, comm,
                      &reqs[rank]);
      if (res != MPI_SUCCESS) mpi_error(res, "Failed to receive new regionids");
      j += nvt;
    }

    int err;
    if ((err = MPI_Waitall(nregions, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 5; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &err);
        message("recv request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall receiving regionid data.");
    }

    /* Copy: idx_t -> int. */
    int *newcelllist = NULL;
    if ((newcelllist = (int *)malloc(sizeof(int) * ncells)) == NULL)
      error("Failed to allocate new celllist");
    for (int k = 0; k < ncells; k++) newcelllist[k] = remoteids[k];
    free(remoteids);

    /* Check that the region ids are all good. */
    int bad = 0;
    for (int k = 0; k < ncells; k++) {
      if (newcelllist[k] < 0 || newcelllist[k] >= nregions) {
        message("Got bad nodeID %" PRIDX " for cell %i.", newcelllist[k], k);
        bad++;
      }
    }
    if (bad) error("Bad node IDs located");

    /* Now check the similarity to the old partition and permute if necessary.
     * Checks show that refinement can return a permutation of the partition,
     * we need to check that and correct as necessary. */
    int permute = 1;
    if (!refine) {

      /* No old partition was given, so we need to construct the existing
       * partition from the cells, if one existed. */
      int nsum = 0;
      for (int i = 0; i < s->nr_cells; i++) {
        celllist[i] = s->cells_top[i].nodeID;
        nsum += celllist[i];
      }

      /* If no previous partition then all nodeIDs will be set to 0. */
      if (nsum == 0) permute = 0;
    }

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
    free(newcelllist);
  }

  /* And everyone gets a copy. */
  res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");

  /* Clean up. */
  free(reqs);
  free(stats);
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(vtxdist);
  free(tpwgts);
  free(xadj);
  free(adjncy);
  free(regionid);
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the given edge and
 * vertex weights. If no weights are given then an unweighted partition is
 * performed.
 *
 * @param nodeID the rank of our node.
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights. Need to be in the range of idx_t.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array. Need to be in the range of
 *        idx_t.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 */
static void pick_metis(int nodeID, struct space *s, int nregions,
                       double *vertexw, double *edgew, int *celllist) {

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* Only one node needs to calculate this. */
  if (nodeID == 0) {

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

    /* Define the cell graph. */
    graph_init(s, adjncy, xadj);

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          weights_v[k] = vertexw[k];
        } else {
          weights_v[k] = 1;
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
        if (weights_v[k] < 1) {
          message("Used vertex weight  out of range: %" PRIDX, weights_v[k]);
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
          weights_e[k] = edgew[k];
        } else {
          weights_e[k] = 1;
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
        if (weights_e[k] < 1) {
          message("Used edge weight out of range: %" PRIDX, weights_e[k]);
          failed++;
        }
      }
      if (failed > 0) error("%d edge weights are out of range", failed);
#endif
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
    /*dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy, weights_v,
      NULL, weights_e);*/

    if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v, NULL,
                            weights_e, &idx_nregions, NULL, NULL, options,
                            &objval, regionid) != METIS_OK)
      error("Call to METIS_PartGraphKway failed.");

    /* Check that the regionids are ok. */
    for (int k = 0; k < ncells; k++) {
      if (regionid[k] < 0 || regionid[k] >= nregions)
        error("Got bad nodeID %" PRIDX " for cell %i.", regionid[k], k);

      /* And keep. */
      celllist[k] = regionid[k];
    }

    /* Clean up. */
    if (weights_v != NULL) free(weights_v);
    if (weights_e != NULL) free(weights_e);
    free(xadj);
    free(adjncy);
    free(regionid);
  }

  /* Calculations all done, now everyone gets a copy. */
  int res = MPI_Bcast(celllist, ncells, MPI_INT, 0, MPI_COMM_WORLD);
  if (res != MPI_SUCCESS) mpi_error(res, "Failed to broadcast new celllist");
}
#endif

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

/* Helper struct for partition_gather weights. */
struct weights_mapper_data {
  double *weights_e;
  double *weights_v;
  idx_t *inds;
  int eweights;
  int nodeID;
  int timebins;
  int vweights;
  int nr_cells;
  int use_ticks;
  struct cell *cells;
};

#ifdef SWIFT_DEBUG_CHECKS
static void check_weights(struct task *tasks, int nr_tasks,
                          struct weights_mapper_data *weights_data,
                          double *weights_v, double *weights_e);
#endif

/**
 * @brief Threadpool mapper function to gather cell edge and vertex weights
 *        from the associated tasks.
 *
 * @param map_data part of the data to process in this mapper.
 * @param num_elements the number of data elements to process.
 * @param extra_data additional data for the mapper context.
 */
static void partition_gather_weights(void *map_data, int num_elements,
                                     void *extra_data) {

  struct task *tasks = (struct task *)map_data;
  struct weights_mapper_data *mydata = (struct weights_mapper_data *)extra_data;

  double *weights_e = mydata->weights_e;
  double *weights_v = mydata->weights_v;
  idx_t *inds = mydata->inds;
  int eweights = mydata->eweights;
  int nodeID = mydata->nodeID;
  int nr_cells = mydata->nr_cells;
  int timebins = mydata->timebins;
  int vweights = mydata->vweights;
  int use_ticks = mydata->use_ticks;

  struct cell *cells = mydata->cells;

  /* Loop over the tasks... */
  for (int i = 0; i < num_elements; i++) {
    struct task *t = &tasks[i];

    /* Skip un-interesting tasks. */
    if (t->type == task_type_send || t->type == task_type_recv ||
        t->type == task_type_logger || t->implicit || t->ci == NULL)
      continue;

    /* Get weight for this task. Either based on fixed costs or task timings. */
    double w = 0.0;
    if (use_ticks) {
      w = (double)t->toc - (double)t->tic;
    } else {
      w = repartition_costs[t->type][t->subtype];
    }
    if (w <= 0.0) continue;

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
        t->type == task_type_end_hydro_force ||
        t->type == task_type_end_grav_force || t->type == task_type_cooling ||
        t->type == task_type_star_formation || t->type == task_type_timestep ||
        t->type == task_type_init_grav || t->type == task_type_grav_down ||
        t->type == task_type_grav_long_range) {

      /* Particle updates add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);
    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub_self && cj == NULL &&
              ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);

    }

    /* Pair? */
    else if (t->type == task_type_pair || (t->type == task_type_sub_pair)) {
      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (vweights) atomic_add_d(&weights_v[cid], w);

      }

      /* Distinct cells. */
      else {
        /* Index of the jth cell. */
        int cjd = cj - cells;

        /* Local cells add weight to vertices. */
        if (vweights && ci->nodeID == nodeID) {
          atomic_add_d(&weights_v[cid], 0.5 * w);
          if (cj->nodeID == nodeID) atomic_add_d(&weights_v[cjd], 0.5 * w);
        }

        if (eweights) {

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
              int dti = num_time_bins - get_time_bin(ci->hydro.ti_end_min);
              int dtj = num_time_bins - get_time_bin(cj->hydro.ti_end_min);
              double dt = (double)(1 << dti) + (double)(1 << dtj);
              atomic_add_d(&weights_e[ik], dt);
              atomic_add_d(&weights_e[jk], dt);

            } else {

              /* Add weights from task costs to the edge. */
              atomic_add_d(&weights_e[ik], w);
              atomic_add_d(&weights_e[jk], w);
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Repartition the cells amongst the nodes using weights of
 *        various kinds.
 *
 * @param vweights whether vertex weights will be used.
 * @param eweights whether weights will be used.
 * @param timebins use timebins as the edge weights.
 * @param repartition the partition struct of the local engine.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
static void repart_edge_metis(int vweights, int eweights, int timebins,
                              struct repartition *repartition, int nodeID,
                              int nr_nodes, struct space *s, struct task *tasks,
                              int nr_tasks) {

  /* Create weight arrays using task ticks for vertices and edges (edges
   * assume the same graph structure as used in the part_ calls). */
  int nr_cells = s->nr_cells;
  struct cell *cells = s->cells_top;

  /* Allocate and fill the adjncy indexing array defining the graph of
   * cells. */
  idx_t *inds;
  if ((inds = (idx_t *)malloc(sizeof(idx_t) * 26 * nr_cells)) == NULL)
    error("Failed to allocate the inds array");
  graph_init(s, inds, NULL);

  /* Allocate and init weights. */
  double *weights_v = NULL;
  double *weights_e = NULL;
  if (vweights) {
    if ((weights_v = (double *)malloc(sizeof(double) * nr_cells)) == NULL)
      error("Failed to allocate vertex weights arrays.");
    bzero(weights_v, sizeof(double) * nr_cells);
  }
  if (eweights) {
    if ((weights_e = (double *)malloc(sizeof(double) * 26 * nr_cells)) == NULL)
      error("Failed to allocate edge weights arrays.");
    bzero(weights_e, sizeof(double) * 26 * nr_cells);
  }

  /* Gather weights. */
  struct weights_mapper_data weights_data;

  weights_data.cells = cells;
  weights_data.eweights = eweights;
  weights_data.inds = inds;
  weights_data.nodeID = nodeID;
  weights_data.nr_cells = nr_cells;
  weights_data.timebins = timebins;
  weights_data.vweights = vweights;
  weights_data.weights_e = weights_e;
  weights_data.weights_v = weights_v;
  weights_data.use_ticks = repartition->use_ticks;

  ticks tic = getticks();

  threadpool_map(&s->e->threadpool, partition_gather_weights, tasks, nr_tasks,
                 sizeof(struct task), 0, &weights_data);
  if (s->e->verbose)
    message("weight mapper took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  check_weights(tasks, nr_tasks, &weights_data, weights_v, weights_e);
#endif

  /* Merge the weights arrays across all nodes. */
  int res;
  if (vweights) {
    res = MPI_Allreduce(MPI_IN_PLACE, weights_v, nr_cells, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to allreduce vertex weights.");
  }

  if (eweights) {
    res = MPI_Allreduce(MPI_IN_PLACE, weights_e, 26 * nr_cells, MPI_DOUBLE,
                        MPI_SUM, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to allreduce edge weights.");
  }

    /* Allocate cell list for the partition. If not already done. */
#ifdef HAVE_PARMETIS
  int refine = 1;
#endif
  if (repartition->ncelllist != nr_cells) {
#ifdef HAVE_PARMETIS
    refine = 0;
#endif
    free(repartition->celllist);
    repartition->ncelllist = 0;
    if ((repartition->celllist = (int *)malloc(sizeof(int) * nr_cells)) == NULL)
      error("Failed to allocate celllist");
    repartition->ncelllist = nr_cells;
  }

  /* We need to rescale the sum of the weights so that the sums of the two
   * types of weights are less than IDX_MAX, that is the range of idx_t.  */
  double vsum = 0.0;
  if (vweights)
    for (int k = 0; k < nr_cells; k++) vsum += weights_v[k];
  double esum = 0.0;
  if (eweights)
    for (int k = 0; k < 26 * nr_cells; k++) esum += weights_e[k];

  /* Do the scaling, if needed, keeping both weights in proportion. */
  double vscale = 1.0;
  double escale = 1.0;
  if (vweights && eweights) {
    if (vsum > esum) {
      if (vsum > (double)IDX_MAX) {
        vscale = (double)(IDX_MAX - 1000) / vsum;
        escale = vscale;
      }
    } else {
      if (esum > (double)IDX_MAX) {
        escale = (double)(IDX_MAX - 1000) / esum;
        vscale = escale;
      }
    }
  } else if (vweights) {
    if (vsum > (double)IDX_MAX) {
      vscale = (double)(IDX_MAX - 1000) / vsum;
    }
  } else if (eweights) {
    if (esum > (double)IDX_MAX) {
      escale = (double)(IDX_MAX - 1000) / esum;
    }
  }

  if (vweights && vscale != 1.0) {
    vsum = 0.0;
    for (int k = 0; k < nr_cells; k++) {
      weights_v[k] *= vscale;
      vsum += weights_v[k];
    }
    vscale = 1.0;
  }
  if (eweights && escale != 1.0) {
    esum = 0.0;
    for (int k = 0; k < 26 * nr_cells; k++) {
      weights_e[k] *= escale;
      esum += weights_e[k];
    }
    escale = 1.0;
  }

  /* Balance edges and vertices when the edge weights are timebins, as these
   * have no reason to have equivalent scales, we use an equipartition. */
  if (timebins && eweights) {

    /* Make sums the same. */
    if (vsum > esum) {
      escale = vsum / esum;
      for (int k = 0; k < 26 * nr_cells; k++) weights_e[k] *= escale;
    } else {
      vscale = esum / vsum;
      for (int k = 0; k < nr_cells; k++) weights_v[k] *= vscale;
    }
  }

    /* And repartition/ partition, using both weights or not as requested. */
#ifdef HAVE_PARMETIS
  if (repartition->usemetis) {
    pick_metis(nodeID, s, nr_nodes, weights_v, weights_e,
               repartition->celllist);
  } else {
    pick_parmetis(nodeID, s, nr_nodes, weights_v, weights_e, refine,
                  repartition->adaptive, repartition->itr,
                  repartition->celllist);
  }
#else
  pick_metis(nodeID, s, nr_nodes, weights_v, weights_e, repartition->celllist);
#endif

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
          "WARNING: repartition has failed, continuing with the current"
          " partition, load balance will not be optimal");
    for (int k = 0; k < nr_cells; k++)
      repartition->celllist[k] = cells[k].nodeID;
  }

  /* And apply to our cells */
  split_metis(s, nr_nodes, repartition->celllist);

  /* Clean up. */
  free(inds);
  if (vweights) free(weights_v);
  if (eweights) free(weights_e);
}

/**
 * @brief Repartition the cells amongst the nodes using weights based on
 *        the memory use of particles in the cells.
 *
 * @param repartition the partition struct of the local engine.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 */
static void repart_memory_metis(struct repartition *repartition, int nodeID,
                                int nr_nodes, struct space *s) {

  /* Space for counts of particle memory use per cell. */
  double *weights = NULL;
  if ((weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
    error("Failed to allocate cell weights buffer.");
  bzero(weights, sizeof(double) * s->nr_cells);

  /* Check each particle and accumulate the sizes per cell. */
  accumulate_sizes(s, weights);

  /* Get all the counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce particle cell weights.");

    /* Allocate cell list for the partition. If not already done. */
#ifdef HAVE_PARMETIS
  int refine = 1;
#endif
  if (repartition->ncelllist != s->nr_cells) {
#ifdef HAVE_PARMETIS
    refine = 0;
#endif
    free(repartition->celllist);
    repartition->ncelllist = 0;
    if ((repartition->celllist = (int *)malloc(sizeof(int) * s->nr_cells)) ==
        NULL)
      error("Failed to allocate celllist");
    repartition->ncelllist = s->nr_cells;
  }

  /* We need to rescale the sum of the weights so that the sum is
   * less than IDX_MAX, that is the range of idx_t. */
  double sum = 0.0;
  for (int k = 0; k < s->nr_cells; k++) sum += weights[k];
  if (sum > (double)IDX_MAX) {
    double scale = (double)(IDX_MAX - 1000) / sum;
    for (int k = 0; k < s->nr_cells; k++) weights[k] *= scale;
  }

    /* And repartition. */
#ifdef HAVE_PARMETIS
  if (repartition->usemetis) {
    pick_metis(nodeID, s, nr_nodes, weights, NULL, repartition->celllist);
  } else {
    pick_parmetis(nodeID, s, nr_nodes, weights, NULL, refine,
                  repartition->adaptive, repartition->itr,
                  repartition->celllist);
  }
#else
  pick_metis(nodeID, s, nr_nodes, weights, NULL, repartition->celllist);
#endif

  /* Check that all cells have good values. All nodes have same copy, so just
   * check on one. */
  if (nodeID == 0) {
    for (int k = 0; k < s->nr_cells; k++)
      if (repartition->celllist[k] < 0 || repartition->celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", repartition->celllist[k], k);
  }

  /* Check that the partition is complete and all nodes have some cells. */
  int present[nr_nodes];
  int failed = 0;
  for (int i = 0; i < nr_nodes; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) present[repartition->celllist[i]]++;
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
          "WARNING: repartition has failed, continuing with the current"
          " partition, load balance will not be optimal");
    for (int k = 0; k < s->nr_cells; k++)
      repartition->celllist[k] = s->cells_top[k].nodeID;
  }

  /* And apply to our cells */
  split_metis(s, nr_nodes, repartition->celllist);
}
#endif /* WITH_MPI && (HAVE_METIS || HAVE_PARMETIS) */

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

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

  ticks tic = getticks();

  if (reparttype->type == REPART_METIS_VERTEX_EDGE_COSTS) {
    repart_edge_metis(1, 1, 0, reparttype, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype->type == REPART_METIS_EDGE_COSTS) {
    repart_edge_metis(0, 1, 0, reparttype, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COSTS_TIMEBINS) {
    repart_edge_metis(1, 1, 1, reparttype, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COUNTS) {
    repart_memory_metis(reparttype, nodeID, nr_nodes, s);

  } else if (reparttype->type == REPART_NONE) {
    /* Doing nothing. */

  } else {
    error("Impossible repartition type");
  }

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with METIS or ParMETIS support.");
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
  ticks tic = getticks();

  /* Geometric grid partitioning. */
  if (initial_partition->type == INITPART_GRID) {
    int j, k;
    int ind[3];
    struct cell *c;

    /* If we've got the wrong number of nodes, fail. */
    if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                        initial_partition->grid[2])
      error("Grid size does not match number of nodes.");

    /* Run through the cells and set their nodeID. */
    // message("s->dim = [%e,%e,%e]", s->dim[0], s->dim[1], s->dim[2]);
    for (k = 0; k < s->nr_cells; k++) {
      c = &s->cells_top[k];
      for (j = 0; j < 3; j++)
        ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
      c->nodeID = ind[0] + initial_partition->grid[0] *
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

  } else if (initial_partition->type == INITPART_METIS_WEIGHT ||
             initial_partition->type == INITPART_METIS_NOWEIGHT) {
#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
    /* Simple k-way partition selected by METIS using cell particle
     * counts as weights or not. Should be best when starting with a
     * inhomogeneous dist.
     */

    /* Space for particles sizes per cell, which will be used as weights. */
    double *weights = NULL;
    if (initial_partition->type == INITPART_METIS_WEIGHT) {
      if ((weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate weights buffer.");
      bzero(weights, sizeof(double) * s->nr_cells);

      /* Check each particle and accumilate the sizes per cell. */
      accumulate_sizes(s, weights);

      /* Get all the counts from all the nodes. */
      if (MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
        error("Failed to allreduce particle cell weights.");
    }

    /* Do the calculation. */
    int *celllist = NULL;
    if ((celllist = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
      error("Failed to allocate celllist");
#ifdef HAVE_PARMETIS
    if (initial_partition->usemetis) {
      pick_metis(nodeID, s, nr_nodes, weights, NULL, celllist);
    } else {
      pick_parmetis(nodeID, s, nr_nodes, weights, NULL, 0, 0, 0.0f, celllist);
    }
#else
    pick_metis(nodeID, s, nr_nodes, weights, NULL, celllist);
#endif

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
    error("SWIFT was not compiled with METIS or ParMETIS support");
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

  if (s->e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Initialises the partition and re-partition scheme from the parameter
 *        file.
 *
 * @param partition The #partition scheme to initialise.
 * @param repartition The #repartition scheme to initialise.
 * @param params The parsed parameter file.
 * @param nr_nodes The number of MPI nodes we are running on.
 */
void partition_init(struct partition *partition,
                    struct repartition *repartition,
                    struct swift_params *params, int nr_nodes) {

#ifdef WITH_MPI

/* Defaults make use of METIS if available */
#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  const char *default_repart = "fullcosts";
  const char *default_part = "memory";
#else
  const char *default_repart = "none";
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
#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
    case 'r':
      partition->type = INITPART_METIS_NOWEIGHT;
      break;
    case 'm':
      partition->type = INITPART_METIS_WEIGHT;
      break;
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid', 'region', 'memory' or "
          "'vectorized'");
#else
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid' or 'vectorized' when compiled "
          "without METIS or ParMETIS.");
#endif
  }

  /* In case of grid, read more parameters */
  if (part_type[0] == 'g') {
    parser_get_opt_param_int_array(params, "DomainDecomposition:initial_grid",
                                   3, partition->grid);
  }

  /* Now let's check what the user wants as a repartition strategy */
  parser_get_opt_param_string(params, "DomainDecomposition:repartition_type",
                              part_type, default_repart);

  if (strcmp("none", part_type) == 0) {
    repartition->type = REPART_NONE;

#if defined(HAVE_METIS) || defined(HAVE_PARMETIS)
  } else if (strcmp("fullcosts", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_EDGE_COSTS;

  } else if (strcmp("edgecosts", part_type) == 0) {
    repartition->type = REPART_METIS_EDGE_COSTS;

  } else if (strcmp("memory", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COUNTS;

  } else if (strcmp("timecosts", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COSTS_TIMEBINS;

  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none', 'fullcosts', 'edgecosts' "
        "'memory' or 'timecosts'");
#else
  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none' when compiled without "
        "METIS or ParMETIS.");
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
   * based on CPU time is considered, needs to be high. */
  repartition->minfrac =
      parser_get_opt_param_float(params, "DomainDecomposition:minfrac", 0.95f);
  if (repartition->minfrac <= 0.5 || repartition->minfrac > 1)
    error(
        "Invalid DomainDecomposition:minfrac, must be greater than 0.5 "
        "and less than equal to 1");

  /* Use METIS or ParMETIS when ParMETIS is also available. */
  repartition->usemetis =
      parser_get_opt_param_int(params, "DomainDecomposition:usemetis", 0);
  partition->usemetis = repartition->usemetis;

  /* Use adaptive or simple refinement when repartitioning. */
  repartition->adaptive =
      parser_get_opt_param_int(params, "DomainDecomposition:adaptive", 1);

  /* Ratio of interprocess communication time to data redistribution time. */
  repartition->itr =
      parser_get_opt_param_float(params, "DomainDecomposition:itr", 100.0f);

  /* Clear the celllist for use. */
  repartition->ncelllist = 0;
  repartition->celllist = NULL;

  /* Do we have fixed costs available? These can be used to force
   * repartitioning at any time. Not required if not repartitioning.*/
  repartition->use_fixed_costs = parser_get_opt_param_int(
      params, "DomainDecomposition:use_fixed_costs", 0);
  if (repartition->type == REPART_NONE) repartition->use_fixed_costs = 0;

  /* Check if this is true or required and initialise them. */
  if (repartition->use_fixed_costs || repartition->trigger > 1) {
    if (!repart_init_fixed_costs()) {
      if (repartition->trigger <= 1) {
        if (engine_rank == 0)
          message(
              "WARNING: fixed cost repartitioning was requested but is"
              " not available.");
        repartition->use_fixed_costs = 0;
      } else {
        error(
            "Forced fixed cost repartitioning was requested but is"
            " not available.");
      }
    }
  }

#else
  error("SWIFT was not compiled with MPI support");
#endif
}

#ifdef WITH_MPI
/**
 * @brief Set the fixed costs for repartition using METIS.
 *
 *  These are determined using a run with the -y flag on which produces
 *  a statistical analysis that is condensed into a .h file for inclusion.
 *
 *  If the default include file is used then no fixed costs are set and this
 *  function will return 0.
 */
static int repart_init_fixed_costs(void) {

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
  /* Set the default fixed cost. */
  for (int j = 0; j < task_type_count; j++) {
    for (int k = 0; k < task_subtype_count; k++) {
      repartition_costs[j][k] = 1.0;
    }
  }

#include <partition_fixed_costs.h>
  return HAVE_FIXED_COSTS;
#endif

  return 0;
}
#endif /* WITH_MPI */

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

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Check that the threadpool version of the weights construction is
 *        correct by comparing to the old serial code.
 *
 * @param tasks the list of tasks
 * @param nr_tasks number of tasks
 * @param mydata additional values as passed to threadpool
 * @param ref_weights_v vertex weights to check
 * @param ref_weights_e edge weights to check
 */
static void check_weights(struct task *tasks, int nr_tasks,
                          struct weights_mapper_data *mydata,
                          double *ref_weights_v, double *ref_weights_e) {

  idx_t *inds = mydata->inds;
  int eweights = mydata->eweights;
  int nodeID = mydata->nodeID;
  int nr_cells = mydata->nr_cells;
  int timebins = mydata->timebins;
  int vweights = mydata->vweights;
  int use_ticks = mydata->use_ticks;

  struct cell *cells = mydata->cells;

  /* Allocate and init weights. */
  double *weights_v = NULL;
  double *weights_e = NULL;
  if (vweights) {
    if ((weights_v = (double *)malloc(sizeof(double) * nr_cells)) == NULL)
      error("Failed to allocate vertex weights arrays.");
    bzero(weights_v, sizeof(double) * nr_cells);
  }
  if (eweights) {
    if ((weights_e = (double *)malloc(sizeof(double) * 26 * nr_cells)) == NULL)
      error("Failed to allocate edge weights arrays.");
    bzero(weights_e, sizeof(double) * 26 * nr_cells);
  }

  /* Loop over the tasks... */
  for (int j = 0; j < nr_tasks; j++) {

    /* Get a pointer to the kth task. */
    struct task *t = &tasks[j];

    /* Skip un-interesting tasks. */
    if (t->type == task_type_send || t->type == task_type_recv ||
        t->type == task_type_logger || t->implicit || t->ci == NULL)
      continue;

    /* Get weight for this task. Either based on fixed costs or task timings. */
    double w = 0.0;
    if (use_ticks) {
      w = (double)t->toc - (double)t->tic;
    } else {
      w = repartition_costs[t->type][t->subtype];
    }
    if (w <= 0.0) continue;

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
        t->type == task_type_end_hydro_force ||
        t->type == task_type_end_grav_force || t->type == task_type_cooling ||
        t->type == task_type_star_formation || t->type == task_type_timestep ||
        t->type == task_type_init_grav || t->type == task_type_grav_down ||
        t->type == task_type_grav_long_range) {

      /* Particle updates add only to vertex weight. */
      if (vweights) weights_v[cid] += w;
    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub_self && cj == NULL &&
              ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (vweights) weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair || (t->type == task_type_sub_pair)) {
      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (vweights) weights_v[cid] += w;

      }

      /* Distinct cells. */
      else {
        /* Index of the jth cell. */
        int cjd = cj - cells;

        /* Local cells add weight to vertices. */
        if (vweights && ci->nodeID == nodeID) {
          weights_v[cid] += 0.5 * w;
          if (cj->nodeID == nodeID) weights_v[cjd] += 0.5 * w;
        }

        if (eweights) {

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
              int dti = num_time_bins - get_time_bin(ci->hydro.ti_end_min);
              int dtj = num_time_bins - get_time_bin(cj->hydro.ti_end_min);
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
  }

  /* Now do the comparisons. */
  double refsum = 0.0;
  double sum = 0.0;
  for (int k = 0; k < nr_cells; k++) {
    refsum += ref_weights_v[k];
    sum += weights_v[k];
  }
  if (fabs(sum - refsum) > 1.0) {
    error("vertex partition weights are not consistent (%f!=%f)", sum, refsum);
  } else {
    refsum = 0.0;
    sum = 0.0;
    for (int k = 0; k < 26 * nr_cells; k++) {
      refsum += ref_weights_e[k];
      sum += weights_e[k];
    }
    if (fabs(sum - refsum) > 1.0) {
      error("edge partition weights are not consistent (%f!=%f)", sum, refsum);
    }
  }
  message("partition weights checked successfully");
}
#endif
#endif

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
    if ((reparttype->celllist =
             (int *)malloc(sizeof(int) * reparttype->ncelllist)) == NULL)
      error("Failed to allocate celllist");
    restart_read_blocks(reparttype->celllist,
                        sizeof(int) * reparttype->ncelllist, 1, stream, NULL,
                        "repartition celllist");
  }
}
