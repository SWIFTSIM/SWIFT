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
#include "threadpool.h"
#include "tools.h"

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

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
  int nadjcny = 0;
  int nxadj = 0;
  graph_init(s, 1 /* periodic */, NULL /* no edge weights */, inds, &nadjcny,
             NULL /* no xadj needed */, &nxadj);

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
                 sizeof(struct task), threadpool_auto_chunk_size,
                 &weights_data);
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
        vscale = (double)(IDX_MAX - 10000) / vsum;
        escale = vscale;
      }
    } else {
      if (esum > (double)IDX_MAX) {
        escale = (double)(IDX_MAX - 10000) / esum;
        vscale = escale;
      }
    }
  } else if (vweights) {
    if (vsum > (double)IDX_MAX) {
      vscale = (double)(IDX_MAX - 10000) / vsum;
    }
  } else if (eweights) {
    if (esum > (double)IDX_MAX) {
      escale = (double)(IDX_MAX - 10000) / esum;
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

  /* Check each particle and accumulate the sizes per cell. */
  accumulate_sizes(s, s->e->verbose, weights);

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
