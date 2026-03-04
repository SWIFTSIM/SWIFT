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
void check_weights(struct task *tasks, int nr_tasks,
                   struct weights_mapper_data *mydata, double *ref_weights_v,
                   double *ref_weights_e) {

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
      /* Nothing to do here */
    }
    if (t->cj != NULL) {
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent) {
        /* Nothing to do here */
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
      if (vweights) weights_v[cid] += w;
    }

    /* Self interaction? */
    else if (t->type == task_type_self && ci->nodeID == nodeID) {
      /* Self interactions add only to vertex weight. */
      if (vweights) weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair) {
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
  if (vweights) {
    if (engine_rank == 0) message("checking vertex weight consistency");
    if (ref_weights_v == NULL)
      error("vertex partition weights are inconsistent");
    for (int k = 0; k < nr_cells; k++) {
      refsum += ref_weights_v[k];
      sum += weights_v[k];
    }
    if (fabs(sum - refsum) > 1.0) {
      error("vertex partition weights are not consistent (%f!=%f)", sum,
            refsum);
    }
  }
  if (eweights) {
    if (engine_rank == 0) message("checking edge weight consistency");
    refsum = 0.0;
    sum = 0.0;
    if (ref_weights_e == NULL) error("edge partition weights are inconsistent");
    for (int k = 0; k < 26 * nr_cells; k++) {
      refsum += ref_weights_e[k];
      sum += weights_e[k];
    }
    if (fabs(sum - refsum) > 1.0) {
      error("edge partition weights are not consistent (%f!=%f)", sum, refsum);
    }
  }
  if (engine_rank == 0) message("partition weights checked successfully");
}
#endif
#endif
