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
#include "threadpool.h"
#include "tools.h"

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))

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
        t->type == task_type_csds || t->implicit || t->ci == NULL)
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
    for (ci = t->ci; ci->parent != NULL; ci = ci->parent) {
      /* Nothing to do here. */
    }

    if (t->cj != NULL) {
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent) {
        /* Nothing to do here. */
      }
    } else {
      cj = NULL;
    }

    /* Get the cell IDs. */
    int cid = ci - cells;

    /* Different weights for different tasks. */
    if (t->type == task_type_init_grav || t->type == task_type_ghost ||
        t->type == task_type_extra_ghost || t->type == task_type_drift_part ||
        t->type == task_type_drift_spart || t->type == task_type_drift_sink ||
        t->type == task_type_drift_bpart || t->type == task_type_drift_gpart ||
        t->type == task_type_end_hydro_force || t->type == task_type_kick1 ||
        t->type == task_type_kick2 || t->type == task_type_timestep ||
        t->type == task_type_timestep_limiter ||
        t->type == task_type_timestep_sync ||
        t->type == task_type_grav_long_range || t->type == task_type_grav_mm ||
        t->type == task_type_grav_down || t->type == task_type_end_grav_force ||
        t->type == task_type_cooling || t->type == task_type_star_formation ||
        t->type == task_type_star_formation_sink ||
        t->type == task_type_stars_ghost ||
        t->type == task_type_bh_density_ghost ||
        t->type == task_type_bh_swallow_ghost2 ||
        t->type == task_type_sink_density_ghost ||
        t->type == task_type_sink_ghost2 ||
        t->type == task_type_neutrino_weight ||
        t->type == task_type_sink_formation || t->type == task_type_rt_ghost1 ||
        t->type == task_type_rt_ghost2 || t->type == task_type_rt_tchem) {

      /* Particle updates add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);
    }

    /* Self interaction? */
    else if (t->type == task_type_self && ci->nodeID == nodeID) {
      /* Self interactions add only to vertex weight. */
      if (vweights) atomic_add_d(&weights_v[cid], w);
    }

    /* Pair? */
    else if (t->type == task_type_pair) {

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
  int nadjcny = 0;
  int nxadj = 0;
  partition_graph_init(s, 1 /* periodic */, NULL /* no edge weights */, inds,
                       &nadjcny, NULL /* no xadj needed */, &nxadj);

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
  partition_check_weights(tasks, nr_tasks, &weights_data, weights_v, weights_e);
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
    partition_pick_metis(nodeID, s, nr_nodes, weights_v, weights_e,
                         repartition->celllist);
  } else {
    partition_pick_parmetis(nodeID, s, nr_nodes, weights_v, weights_e, refine,
                            repartition->adaptive, repartition->itr,
                            repartition->celllist);
  }
#else
  partition_pick_metis(nodeID, s, nr_nodes, weights_v, weights_e,
                       repartition->celllist);
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
  partition_split_metis(s, nr_nodes, repartition->celllist);

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
  partition_accumulate_sizes(s, s->e->verbose, weights);

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
    partition_pick_metis(nodeID, s, nr_nodes, weights, NULL,
                         repartition->celllist);
  } else {
    partition_pick_parmetis(nodeID, s, nr_nodes, weights, NULL, refine,
                            repartition->adaptive, repartition->itr,
                            repartition->celllist);
  }
#else
  partition_pick_metis(nodeID, s, nr_nodes, weights, NULL,
                       repartition->celllist);
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
  partition_split_metis(s, nr_nodes, repartition->celllist);
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
