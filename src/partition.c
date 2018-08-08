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
 *  Currently supported partitioning types: grid, vectorise and METIS.
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
/* METIS headers only used when MPI is also available. */
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

/* Maximum weight used for METIS. */
#define metis_maxweight 10000.0f

/* Simple descriptions of initial partition types for reports. */
const char *initial_partition_name[] = {
    "gridded cells", "vectorized point associated cells",
    "METIS particle weighted cells", "METIS unweighted cells"};

/* Simple descriptions of repartition types for reports. */
const char *repartition_name[] = {
    "no",
    "METIS edge and vertex task cost weights",
    "METIS particle count vertex weights",
    "METIS task cost edge weights",
    "METIS particle count vertex and task cost edge weights",
    "METIS vertex task costs and edge delta timebin weights",
    "METIS particle count vertex and edge delta timebin weights",
    "METIS edge delta timebin weights",
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

  /* METIS support
   * =============
   *
   * METIS partitions using a multi-level k-way scheme. We support using this in
   * a unweighted scheme, which works well and seems to be guaranteed, and a
   * weighted by the number of particles scheme. Note METIS is optional.
   *
   * Repartitioning is based on METIS and uses weights determined from the times
   * that cell tasks have taken. These weight the graph edges and vertices, or
   * just the edges, with vertex weights from the particle counts or none.
   */

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Fill the METIS xadj and adjncy arrays defining the graph of cells
 *        in a space.
 *
 * See the METIS manual if you want to understand this format. The cell graph
 * consists of all nodes as vertices with edges as the connections to all
 * neighbours, so we have 26 per vertex.
 *
 * @param s the space of cells.
 * @param adjncy the METIS adjncy array to fill, must be of size
 *               26 * the number of cells in the space.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 */
static void graph_init_metis(struct space *s, idx_t *adjncy, idx_t *xadj) {

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

  /* If given set xadj. */
  if (xadj != NULL) {
    xadj[0] = 0;
    for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
  }
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Accumulate the counts of particles per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of particles per cell. Should be
 *               allocated as size s->nr_parts.
 */
static void accumulate_counts(struct space *s, double *counts) {

  struct part *parts = s->parts;
  int *cdim = s->cdim;
  double iwidth[3] = {s->iwidth[0], s->iwidth[1], s->iwidth[2]};
  double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  bzero(counts, sizeof(double) * s->nr_cells);

  for (size_t k = 0; k < s->nr_parts; k++) {
    for (int j = 0; j < 3; j++) {
      if (parts[k].x[j] < 0.0)
        parts[k].x[j] += dim[j];
      else if (parts[k].x[j] >= dim[j])
        parts[k].x[j] -= dim[j];
    }
    const int cid =
        cell_getid(cdim, parts[k].x[0] * iwidth[0], parts[k].x[1] * iwidth[1],
                   parts[k].x[2] * iwidth[2]);
    counts[cid]++;
  }
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
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
static void split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) s->cells_top[i].nodeID = celllist[i];

  /* To check or visualise the partition dump all the cells. */
  /* dumpCellRanks("metis_partition", s->cells_top, s->nr_cells); */
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)

/* qsort support. */
struct indexval {
  int index;
  int count;
};
static int indexvalcmp(const void *p1, const void *p2) {
  const struct indexval *iv1 = (const struct indexval *)p1;
  const struct indexval *iv2 = (const struct indexval *)p2;
  return iv2->count - iv1->count;
}

/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the
 * given edge and vertex weights. If no weights are given then an
 * unweighted partition is performed.
 *
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
static void pick_metis(struct space *s, int nregions, double *vertexw,
                       double *edgew, int *celllist) {

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

  /* Define the cell graph. */
  graph_init_metis(s, adjncy, xadj);

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
  /*dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy,
   *               weights_v, NULL, weights_e);
   */
  if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v, NULL,
                          weights_e, &idx_nregions, NULL, NULL, options,
                          &objval, regionid) != METIS_OK)
    error("Call to METIS_PartGraphKway failed.");

  /* Check that the regionids are ok. */
  for (int k = 0; k < ncells; k++)
    if (regionid[k] < 0 || regionid[k] >= nregions)
      error("Got bad nodeID %" PRIDX " for cell %i.", regionid[k], k);

  /* We want a solution in which the current regions of the space are
   * preserved when possible, to avoid unneccesary particle movement.
   * So create a 2d-array of cells counts that are common to all pairs
   * of old and new ranks. Each element of the array has a cell count and
   * an unique index so we can sort into decreasing counts. */
  int indmax = nregions * nregions;
  struct indexval *ivs =
      (struct indexval *)malloc(sizeof(struct indexval) * indmax);
  bzero(ivs, sizeof(struct indexval) * indmax);
  for (int k = 0; k < ncells; k++) {
    int index = regionid[k] + nregions * s->cells_top[k].nodeID;
    ivs[index].count++;
    ivs[index].index = index;
  }
  qsort(ivs, indmax, sizeof(struct indexval), indexvalcmp);

  /* Go through the ivs using the largest counts first, these are the
   * regions with the most cells in common, old partition to new. */
  int *oldmap = (int *)malloc(sizeof(int) * nregions);
  int *newmap = (int *)malloc(sizeof(int) * nregions);
  for (int k = 0; k < nregions; k++) {
    oldmap[k] = -1;
    newmap[k] = -1;
  }
  for (int k = 0; k < indmax; k++) {

    /* Stop when all regions with common cells have been considered. */
    if (ivs[k].count == 0) break;

    /* Store old and new IDs, if not already used. */
    int oldregion = ivs[k].index / nregions;
    int newregion = ivs[k].index - oldregion * nregions;
    if (newmap[newregion] == -1 && oldmap[oldregion] == -1) {
      newmap[newregion] = oldregion;
      oldmap[oldregion] = newregion;
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

  /* Set the cell list to the region index. */
  for (int k = 0; k < ncells; k++) {
    celllist[k] = newmap[regionid[k]];
  }

  /* Clean up. */
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(ivs);
  free(oldmap);
  free(newmap);
  free(xadj);
  free(adjncy);
  free(regionid);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Repartition the cells amongst the nodes using task costs
 *        as edge weights and vertex weights also from task costs
 *        or particle cells counts.
 *
 * @param partweights whether particle counts will be used as vertex weights.
 * @param bothweights whether vertex and edge weights will be used, otherwise
 *                    only edge weights will be used.
 * @param timebins use timebins as edge weights.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
static void repart_edge_metis(int partweights, int bothweights, int timebins,
                              int nodeID, int nr_nodes, struct space *s,
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
  graph_init_metis(s, inds, NULL);

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

  /* Generate task weights for vertices. */
  int taskvweights = (bothweights && !partweights);

  /* Loop over the tasks... */
  for (int j = 0; j < nr_tasks; j++) {
    /* Get a pointer to the kth task. */
    struct task *t = &tasks[j];

    /* Skip un-interesting tasks. */
    if (t->cost == 0.f) continue;

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
      if (taskvweights) weights_v[cid] += w;
    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub_self && cj == NULL &&
              ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (taskvweights) weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair || (t->type == task_type_sub_pair)) {
      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (taskvweights) weights_v[cid] += w;

      }

      /* Distinct cells. */
      else {
        /* Index of the jth cell. */
        int cjd = cj - cells;

        /* Local cells add weight to vertices. */
        if (taskvweights && ci->nodeID == nodeID) {
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

  /* Re-calculate the vertices if using particle counts. */
  if (partweights && bothweights) accumulate_counts(s, weights_v);

  /* Merge the weights arrays across all nodes. */
  int res;
  if (bothweights) {
    if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_v, weights_v,
                          nr_cells, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD)) !=
        MPI_SUCCESS)
      mpi_error(res, "Failed to allreduce vertex weights.");
  }

  if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_e, weights_e,
                        26 * nr_cells, MPI_DOUBLE, MPI_SUM, 0,
                        MPI_COMM_WORLD)) != MPI_SUCCESS)
    mpi_error(res, "Failed to allreduce edge weights.");

  /* Allocate cell list for the partition. */
  int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
  if (celllist == NULL) error("Failed to allocate celllist");

  /* As of here, only one node needs to compute the partition. */
  if (nodeID == 0) {

    /* We need to rescale the weights into the range of an integer for METIS
     * (really range of idx_t). Also we would like the range of vertex and
     * edges weights to be similar so they balance. */
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

      /* Scale to the METIS range. */
      double wscale = 1.0;
      if ((wmaxv - wminv) > 0.0) {
        wscale = (metis_maxweight - 1.0) / (wmaxv - wminv);
      }
      for (int k = 0; k < nr_cells; k++) {
        weights_v[k] = (weights_v[k] - wminv) * wscale + 1.0;
      }
    }

    /* Scale to the METIS range. */
    double wscale = 1.0;
    if ((wmaxe - wmine) > 0.0) {
      wscale = (metis_maxweight - 1.0) / (wmaxe - wmine);
    }
    for (int k = 0; k < 26 * nr_cells; k++) {
      weights_e[k] = (weights_e[k] - wmine) * wscale + 1.0;
    }

    /* And partition, use both weights or not as requested. */
    if (bothweights)
      pick_metis(s, nr_nodes, weights_v, weights_e, celllist);
    else
      pick_metis(s, nr_nodes, NULL, weights_e, celllist);

    /* Check that all cells have good values. */
    for (int k = 0; k < nr_cells; k++)
      if (celllist[k] < 0 || celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", celllist[k], k);

    /* Check that the partition is complete and all nodes have some work. */
    int present[nr_nodes];
    int failed = 0;
    for (int i = 0; i < nr_nodes; i++) present[i] = 0;
    for (int i = 0; i < nr_cells; i++) present[celllist[i]]++;
    for (int i = 0; i < nr_nodes; i++) {
      if (!present[i]) {
        failed = 1;
        message("Node %d is not present after repartition", i);
      }
    }

    /* If partition failed continue with the current one, but make this
     * clear. */
    if (failed) {
      message(
          "WARNING: METIS repartition has failed, continuing with "
          "the current partition, load balance will not be optimal");
      for (int k = 0; k < nr_cells; k++) celllist[k] = cells[k].nodeID;
    }
  }

  /* Distribute the celllist partition and apply. */
  if ((res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the cell list");

  /* And apply to our cells */
  split_metis(s, nr_nodes, celllist);

  /* Clean up. */
  free(inds);
  if (bothweights) free(weights_v);
  free(weights_e);
  free(celllist);
}
#endif

/**
 * @brief Repartition the cells amongst the nodes using vertex weights
 *
 * @param s The space containing the local cells.
 * @param nodeID our MPI node id.
 * @param nr_nodes number of MPI nodes.
 */
#if defined(WITH_MPI) && defined(HAVE_METIS)
static void repart_vertex_metis(struct space *s, int nodeID, int nr_nodes) {

  /* Use particle counts as vertex weights. */
  /* Space for particles per cell counts, which will be used as weights. */
  double *weights = NULL;
  if ((weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
    error("Failed to allocate weights buffer.");

  /* Check each particle and accumulate the counts per cell. */
  accumulate_counts(s, weights);

  /* Get all the counts from all the nodes. */
  int res;
  if ((res = MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_DOUBLE,
                           MPI_SUM, MPI_COMM_WORLD)) != MPI_SUCCESS)
    mpi_error(res, "Failed to allreduce particle cell weights.");

  /* Main node does the partition calculation. */
  int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
  if (celllist == NULL) error("Failed to allocate celllist");

  if (nodeID == 0) pick_metis(s, nr_nodes, weights, NULL, celllist);

  /* Distribute the celllist partition and apply. */
  if ((res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the cell list");

  /* And apply to our cells */
  split_metis(s, nr_nodes, celllist);

  free(weights);
  free(celllist);
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

#if defined(WITH_MPI) && defined(HAVE_METIS)

  if (reparttype->type == REPART_METIS_VERTEX_COSTS_EDGE_COSTS) {
    repart_edge_metis(0, 1, 0, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_EDGE_COSTS) {
    repart_edge_metis(0, 0, 0, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COUNTS_EDGE_COSTS) {
    repart_edge_metis(1, 1, 0, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COSTS_EDGE_TIMEBINS) {
    repart_edge_metis(0, 1, 1, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COUNTS_EDGE_TIMEBINS) {
    repart_edge_metis(1, 1, 1, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_EDGE_TIMEBINS) {
    repart_edge_metis(0, 0, 1, nodeID, nr_nodes, s, tasks, nr_tasks);

  } else if (reparttype->type == REPART_METIS_VERTEX_COUNTS) {
    repart_vertex_metis(s, nodeID, nr_nodes);

  } else if (reparttype->type == REPART_NONE) {
    /* Doing nothing. */

  } else {
    error("Impossible repartition type");
  }
#else
  error("SWIFT was not compiled with METIS support.");
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
#if defined(WITH_MPI) && defined(HAVE_METIS)
    /* Simple k-way partition selected by METIS using cell particle counts as
     * weights or not. Should be best when starting with a inhomogeneous dist.
     */

    /* Space for particles per cell counts, which will be used as weights or
     * not. */
    double *weights = NULL;
    if (initial_partition->type == INITPART_METIS_WEIGHT) {
      if ((weights = (double *)malloc(sizeof(double) * s->nr_cells)) == NULL)
        error("Failed to allocate weights buffer.");
      bzero(weights, sizeof(double) * s->nr_cells);

      /* Check each particle and accumilate the counts per cell. */
      accumulate_counts(s, weights);

      /* Get all the counts from all the nodes. */
      if (MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
        error("Failed to allreduce particle cell weights.");
    }

    /* Main node does the partition calculation. */
    int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
    if (celllist == NULL) error("Failed to allocate celllist");
    if (nodeID == 0) pick_metis(s, nr_nodes, weights, NULL, celllist);

    /* Distribute the celllist partition and apply. */
    int res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to bcast the cell list");

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
    int *samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3);
    if (samplecells == NULL) error("Failed to allocate samplecells");

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
                    struct swift_params *params, int nr_nodes) {

#ifdef WITH_MPI

/* Defaults make use of METIS if available */
#ifdef HAVE_METIS
  const char *default_repart = "costs/costs";
  const char *default_part = "simple_metis";
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
#ifdef HAVE_METIS
    case 's':
      partition->type = INITPART_METIS_NOWEIGHT;
      break;
    case 'w':
      partition->type = INITPART_METIS_WEIGHT;
      break;
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid', 'simple_metis', 'weighted_metis'"
          " or 'vectorized'");
#else
    default:
      message("Invalid choice of initial partition type '%s'.", part_type);
      error(
          "Permitted values are: 'grid' or 'vectorized' when compiled "
          "without METIS.");
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

  if (strcmp("none/none", part_type) == 0) {
    repartition->type = REPART_NONE;

#ifdef HAVE_METIS
  } else if (strcmp("costs/costs", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COSTS_EDGE_COSTS;

  } else if (strcmp("counts/none", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COUNTS;

  } else if (strcmp("none/costs", part_type) == 0) {
    repartition->type = REPART_METIS_EDGE_COSTS;

  } else if (strcmp("counts/costs", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COUNTS_EDGE_COSTS;

  } else if (strcmp("costs/time", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COSTS_EDGE_TIMEBINS;

  } else if (strcmp("counts/time", part_type) == 0) {
    repartition->type = REPART_METIS_VERTEX_COUNTS_EDGE_TIMEBINS;

  } else if (strcmp("none/time", part_type) == 0) {
    repartition->type = REPART_METIS_EDGE_TIMEBINS;
  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error(
        "Permitted values are: 'none/none', 'costs/costs',"
        "'counts/none', 'none/costs', 'counts/costs', "
        "'costs/time', 'counts/time' or 'none/time'");
#else
  } else {
    message("Invalid choice of re-partition type '%s'.", part_type);
    error("Permitted values are: 'none/none' when compiled without METIS.");
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

  int *present = (int *)malloc(sizeof(int) * nregions);
  if (present == NULL) error("Failed to allocate present array");

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
  if (reparttype->celllist != NULL) free(reparttype->celllist);
  reparttype->celllist = (int *)malloc(sizeof(int) * s->nr_cells);
  reparttype->ncelllist = s->nr_cells;
  if (reparttype->celllist == NULL) error("Failed to allocate celllist");

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
          "Cannot apply the saved partition celllist as the number of"
          "top-level cells (%d) is different to the saved number (%d)",
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
    reparttype->celllist = (int *)malloc(sizeof(int) * reparttype->ncelllist);
    if (reparttype->celllist == NULL) error("Failed to allocate celllist");
    restart_read_blocks(reparttype->celllist,
                        sizeof(int) * reparttype->ncelllist, 1, stream, NULL,
                        "repartition celllist");
  }
}
