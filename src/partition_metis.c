/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Will J. Roper (w.roper@sussex.ac.uk)
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
#include <config.h>

/* Standard headers. */
#include <math.h>
#include <stdlib.h>
#include <strings.h>

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
#include "engine.h"
#include "error.h"
#include "partition.h"
#include "space.h"
#include "threadpool.h"
#include "tools.h"

/* METIS/ParMETIS support (optional)
 * =================================
 *
 * METIS/ParMETIS partitions using a multi-level k-way scheme. We support
 * using this in a unweighted scheme, which works well and seems to be
 * guaranteed, and a weighted by the number of particles scheme.
 *
 * Repartitioning is based on ParMETIS (but can also use METIS) and uses weights
 * determined from the estimated costs that a cells tasks will take or the
 * relative time bins of the cells next updates.
 */

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/**
 * @brief Fill the adjncy array defining the graph of cells in a space.
 *
 * See the ParMETIS and METIS manuals if you want to understand this
 * format. The cell graph consists of all nodes as vertices with edges as the
 * connections to all neighbours, so we have 26 per vertex for periodic
 * boundary, fewer than 26 on the space edges when non-periodic. Note you will
 * also need an xadj array, for METIS that would be:
 *
 *   xadj[0] = 0;
 *   for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
 *
 * but each rank needs a different xadj when using ParMETIS (each segment
 * should be rezeroed).
 *
 * @param s the space of cells.
 * @param periodic whether to assume a periodic space (fixed 26 edges).
 * @param weights_e the edge weights for the cells, if used. On input
 *                  assumed to be ordered with a fixed 26 edges per cell, so
 *                  will need reordering for non-periodic spaces.
 * @param adjncy the adjncy array to fill, must be of size 26 * the number of
 *               cells in the space.
 * @param nadjcny number of adjncy elements used, can be less if not periodic.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 * @param nxadj the number of xadj element used.
 */
void graph_init(struct space *s, int periodic, idx_t *weights_e, idx_t *adjncy,
                int *nadjcny, idx_t *xadj, int *nxadj) {

  /* Loop over all cells in the space. */
  *nadjcny = 0;
  if (periodic) {
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
    *nadjcny = cid * 26;

    /* If given set METIS xadj. */
    if (xadj != NULL) {
      xadj[0] = 0;
      for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
      *nxadj = s->nr_cells;
    }

  } else {

    /* Non periodic. */
    int ind = 0;
    int cid = 0;
    if (xadj != NULL) xadj[0] = 0;

    /* May need to reorder weights, shuffle in place as moving to left. */
    int shuffle = 0;
    if (weights_e != NULL) shuffle = 1;

    for (int l = 0; l < s->cdim[0]; l++) {
      for (int m = 0; m < s->cdim[1]; m++) {
        for (int n = 0; n < s->cdim[2]; n++) {

          /* Visit all neighbours of this cell. */
          int p = 0;
          for (int i = -1; i <= 1; i++) {
            int ii = l + i;
            if (ii >= 0 && ii < s->cdim[0]) {
              for (int j = -1; j <= 1; j++) {
                int jj = m + j;
                if (jj >= 0 && jj < s->cdim[1]) {
                  for (int k = -1; k <= 1; k++) {
                    int kk = n + k;
                    if (kk >= 0 && kk < s->cdim[2]) {

                      /* If not self, record id of neighbour. */
                      if (i || j || k) {
                        adjncy[ind] = cell_getid(s->cdim, ii, jj, kk);

                        if (shuffle) {
                          /* Keep this weight, need index for periodic
                           * version for input weights... */
                          int oldp = ((i + 1) * 9 + (j + 1) * 3 + (k + 1));
                          oldp = oldp - (oldp / 14);
                          weights_e[ind] = weights_e[cid * 26 + oldp];
                        }

                        ind++;
                        p++;
                      }
                    }
                  }
                }
              }
            }
          }

          /* Keep xadj in sync. */
          if (xadj != NULL) {
            xadj[cid + 1] = xadj[cid] + p;
          }
          cid++;
        }
      }
    }
    *nadjcny = ind;
    *nxadj = cid;
  }
}

/* Data structure for accumulating counts in the mappers. */
struct counts_mapper_data {
  double *counts;
  size_t size;
  struct space *s;
};

/* Generic function for accumulating sized counts for TYPE parts. Note uses
 * local memory to reduce contention, the amount of memory required is
 * precalculated by an additional loop determining the range of cell IDs. */
#define ACCUMULATE_SIZES_MAPPER(TYPE)                                          \
  partition_accumulate_sizes_mapper_##TYPE(void *map_data, int num_elements,   \
                                           void *extra_data) {                 \
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
    if ((lcounts = (double *)calloc(nused, sizeof(double))) == NULL)           \
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
void ACCUMULATE_SIZES_MAPPER(part);

/**
 * @brief Accumulate the sized counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * gpart version.
 */
void ACCUMULATE_SIZES_MAPPER(gpart);

/**
 * @brief Accumulate the sized counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * spart version.
 */
void ACCUMULATE_SIZES_MAPPER(spart);

/* qsort support. */
static int ptrcmp(const void *p1, const void *p2) {
  const double *v1 = *(const double **)p1;
  const double *v2 = *(const double **)p2;
  return (*v1) - (*v2);
}

/**
 * @brief Accumulate total memory size in particles per cell.
 *
 * @param s the space containing the cells.
 * @param verbose whether to report any clipped cell counts.
 * @param counts the number of bytes in particles per cell. Should be
 *               allocated as size s->nr_cells.
 */
void accumulate_sizes(struct space *s, int verbose, double *counts) {

  bzero(counts, sizeof(double) * s->nr_cells);

  struct counts_mapper_data mapper_data;
  mapper_data.s = s;
  double gsize = 0.0;
  double *gcounts = NULL;
  double hsize = 0.0;
  double ssize = 0.0;

  if (s->nr_gparts > 0) {
    /* Self-gravity gets more efficient with density (see gitlab issue #640)
     * so we end up with too much weight in cells with larger numbers of
     * gparts, to suppress this we fix a upper weight limit based on a
     * percentile clip to on the numbers of cells. Should be relatively
     * harmless when not really needed. */
    if ((gcounts = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
      error("Failed to allocate gcounts buffer.");
    bzero(gcounts, sizeof(double) * s->nr_cells);
    gsize = (double)sizeof(struct gpart);

    mapper_data.counts = gcounts;
    mapper_data.size = gsize;
    threadpool_map(&s->e->threadpool, partition_accumulate_sizes_mapper_gpart,
                   s->gparts, s->nr_gparts, sizeof(struct gpart),
                   space_splitsize, &mapper_data);

    /* Get all the counts from all the nodes. */
    if (MPI_Allreduce(MPI_IN_PLACE, gcounts, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to allreduce particle cell gpart weights.");

    /* Now we need to sort... */
    double **ptrs = NULL;
    if ((ptrs = (double **)malloc(sizeof(double *) * s->nr_cells)) == NULL)
      error("Failed to allocate pointers buffer.");
    for (int k = 0; k < s->nr_cells; k++) {
      ptrs[k] = &gcounts[k];
    }

    /* Sort pointers, not counts... */
    qsort(ptrs, s->nr_cells, sizeof(double *), ptrcmp);

    /* Percentile cut keeps 99.8% of cells and clips above. */
    int cut = ceil(s->nr_cells * 0.998);
    if (cut == s->nr_cells) cut = s->nr_cells - 1;

    /* And clip. */
    int nadj = 0;
    double clip = *ptrs[cut];
    for (int k = cut + 1; k < s->nr_cells; k++) {
      *ptrs[k] = clip;
      nadj++;
    }
    if (verbose) message("clipped gravity part counts of %d cells", nadj);
    free(ptrs);
  }

  /* Other particle types are assumed to correlate with processing time. */
  if (s->nr_parts > 0) {
    mapper_data.counts = counts;
    hsize = (double)sizeof(struct part);
    mapper_data.size = hsize;
    threadpool_map(&s->e->threadpool, partition_accumulate_sizes_mapper_part,
                   s->parts, s->nr_parts, sizeof(struct part), space_splitsize,
                   &mapper_data);
  }

  if (s->nr_sparts > 0) {
    ssize = (double)sizeof(struct spart);
    mapper_data.size = ssize;
    threadpool_map(&s->e->threadpool, partition_accumulate_sizes_mapper_spart,
                   s->sparts, s->nr_sparts, sizeof(struct spart),
                   space_splitsize, &mapper_data);
  }

  /* Merge the counts arrays across all nodes, if needed. Doesn't include any
   * gparts. */
  if (s->nr_parts > 0 || s->nr_sparts > 0) {
    if (MPI_Allreduce(MPI_IN_PLACE, counts, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to allreduce particle cell weights.");
  }

  /* And merge with gravity counts. */
  double sum = 0.0;
  if (s->nr_gparts > 0) {
    for (int k = 0; k < s->nr_cells; k++) {
      counts[k] += gcounts[k];
      sum += counts[k];
    }
    free(gcounts);
  } else {
    for (int k = 0; k < s->nr_cells; k++) {
      sum += counts[k];
    }
  }

  /* Keep the sum of particles across all ranks in the range of IDX_MAX. */
  if (sum > (double)(IDX_MAX - 10000)) {
    double vscale = (double)(IDX_MAX - 10000) / sum;
    for (int k = 0; k < s->nr_cells; k++) counts[k] *= vscale;
  }
}

/**
 * @brief Make edge weights from the accumulated particle sizes per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of bytes in particles per cell.
 * @param edges weights for the edges of these regions. Should be 26 * counts.
 */
void sizes_to_edges(struct space *s, double *counts, double *edges) {

  bzero(edges, sizeof(double) * s->nr_cells * 26);

  for (int l = 0; l < s->nr_cells; l++) {
    int p = 0;
    for (int i = -1; i <= 1; i++) {
      int isid = ((i < 0) ? 0 : ((i > 0) ? 2 : 1));
      for (int j = -1; j <= 1; j++) {
        int jsid = isid * 3 + ((j < 0) ? 0 : ((j > 0) ? 2 : 1));
        for (int k = -1; k <= 1; k++) {
          int ksid = jsid * 3 + ((k < 0) ? 0 : ((k > 0) ? 2 : 1));

          /* If not self, we work out the sort indices to get the expected
           * fractional weight and add that. Scale to keep sum less than
           * counts and a bit of tuning... */
          if (i || j || k) {
            edges[l * 26 + p] = counts[l] * sid_scale[sortlistID[ksid]] / 26.0;
            p++;
          }
        }
      }
    }
  }
}

/**
 * @brief Apply METIS cell-list partitioning to a cell structure.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
void split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) s->cells_top[i].nodeID = celllist[i];

  /* To check or visualise the partition dump all the cells. */
  /*if (engine_rank == 0) dumpCellRanks("metis_partition", s->cells_top,
                                      s->nr_cells);*/
}

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

#ifdef HAVE_PARMETIS
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
void pick_parmetis(int nodeID, struct space *s, int nregions, double *vertexw,
                   double *edgew, int refine, int adaptive, float itr,
                   int *celllist) {

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
      error("Failed to allocate full xadj buffer.");
    idx_t *std_xadj = NULL;
    if ((std_xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
      error("Failed to allocate std xadj buffer.");
    idx_t *full_adjncy = NULL;
    if ((full_adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
      error("Failed to allocate full adjncy array.");
    idx_t *full_weights_v = NULL;
    if (weights_v != NULL)
      if ((full_weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate full vertex weights array");
    idx_t *full_weights_e = NULL;
    if (weights_e != NULL)
      if ((full_weights_e = (idx_t *)malloc(26 * sizeof(idx_t) * ncells)) ==
          NULL)
        error("Failed to allocate full edge weights array");

    idx_t *full_regionid = NULL;
    if (refine) {
      if ((full_regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
        error("Failed to allocate full regionid array");
    }

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          full_weights_v[k] = vertexw[k];
        } else {
          full_weights_v[k] = 0;
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
        if (full_weights_v[k] < 0) {
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
      for (int k = 0; k < ncells * 26; k++) {
        if (edgew[k] > 1) {
          full_weights_e[k] = edgew[k];
        } else {
          full_weights_e[k] = 1;
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check weights are all in range. */
      int failed = 0;
      for (int k = 0; k < ncells * 26; k++) {

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

    /* Define the cell graph. Keeping the edge weights association. */
    int nadjcny = 0;
    int nxadj = 0;
    graph_init(s, s->periodic, full_weights_e, full_adjncy, &nadjcny, std_xadj,
               &nxadj);

    /* Dump graphs to disk files for testing. */
    /*dumpMETISGraph("parmetis_graph", ncells, 1, std_xadj, full_adjncy,
                   full_weights_v, NULL, full_weights_e); */

    /* xadj is set for each rank, different to serial version in that each
     * rank starts with 0, so we need to re-offset. */
    for (int rank = 0, i = 0, j = 0; rank < nregions; rank++) {

      /* Number of vertices for this rank. */
      int nvt = vtxdist[rank + 1] - vtxdist[rank];

      /* Each xadj section starts at 0 and terminates like a complete one. */
      int offset = std_xadj[j];
      for (int k = 0; k < nvt; k++) {
        full_xadj[i] = std_xadj[j] - offset;
        j++;
        i++;
      }
      full_xadj[i] = std_xadj[j] - offset;
      i++;
    }

    /* Send ranges to the other ranks and keep our own. */
    for (int rank = 0, j1 = 0, j2 = 0, j3 = 0; rank < nregions; rank++) {
      int nvt = vtxdist[rank + 1] - vtxdist[rank];
      int nedge = std_xadj[vtxdist[rank + 1]] - std_xadj[vtxdist[rank]];

      if (refine)
        for (int i = 0; i < nvt; i++) full_regionid[j3 + i] = celllist[j3 + i];

      if (rank == 0) {
        memcpy(xadj, &full_xadj[j1], sizeof(idx_t) * (nvt + 1));
        memcpy(adjncy, &full_adjncy[j2], sizeof(idx_t) * nedge);
        if (weights_e != NULL)
          memcpy(weights_e, &full_weights_e[j2], sizeof(idx_t) * nedge);
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

      /* Note we send 26 edges, but only increment by the correct number. */
      j2 += nedge;
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
    free(std_xadj);
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
        message("Got bad nodeID %d for cell %i.", newcelllist[k], k);
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
void pick_metis(int nodeID, struct space *s, int nregions, double *vertexw,
                double *edgew, int *celllist) {

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

    /* Allocate adjacency and weights arrays . */
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

    /* Init the vertex weights array. */
    if (vertexw != NULL) {
      for (int k = 0; k < ncells; k++) {
        if (vertexw[k] > 1) {
          weights_v[k] = vertexw[k];
        } else {
          weights_v[k] = 0;
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
        if (weights_v[k] < 0) {
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

    /* Define the cell graph. Keeping the edge weights association. */
    int nadjcny = 0;
    int nxadj = 0;
    graph_init(s, s->periodic, weights_e, adjncy, &nadjcny, xadj, &nxadj);

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
