/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS headers only used when MPI is also available. */
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "timers.h"

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/** The rank of the engine as a global variable (for messages). */
int engine_rank;

/**
 * @brief Link a density/force task to a cell.
 *
 * @param e The #engine.
 * @param l The #link.
 * @param t The #task.
 *
 * @return The new #link pointer.
 */

struct link *engine_addlink(struct engine *e, struct link *l, struct task *t) {

  const int ind = atomic_inc(&e->nr_links);
  if (ind >= e->size_links) {
    error("Link table overflow.");
  }
  struct link *res = &e->links[ind];
  res->next = l;
  res->t = t;
  return res;
}

/**
 * @brief Generate the ghost and kick tasks for a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param super The super #cell.
 */

void engine_mkghosts(struct engine *e, struct cell *c, struct cell *super) {

  int k;
  struct scheduler *s = &e->sched;

  //  message("in here");

  /* Am I the super-cell? */
  if (super == NULL && c->nr_tasks > 0) {

    /* Remember me. */
    super = c;

    // message("Adding tasks");

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Generate the ghost task. */
      c->ghost = scheduler_addtask(s, task_type_ghost, task_subtype_none, 0, 0,
                                   c, NULL, 0);
      /* Add the drift task. */
      c->drift = scheduler_addtask(s, task_type_drift, task_subtype_none, 0, 0,
                                   c, NULL, 0);
      /* Add the init task. */
      c->init = scheduler_addtask(s, task_type_init, task_subtype_none, 0, 0, c,
                                  NULL, 0);
      /* Add the kick task. */
      c->kick = scheduler_addtask(s, task_type_kick, task_subtype_none, 0, 0, c,
                                  NULL, 0);
    }
  }

  /* Set the super-cell. */
  c->super = super;

  /* Recurse. */
  if (c->split)
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) engine_mkghosts(e, c->progeny[k], super);
}

/**
 * @brief Redistribute the particles amongst the nodes according
 *      to their cell's node IDs.
 *
 * @param e The #engine.
 */

void engine_redistribute(struct engine *e) {

#ifdef WITH_MPI

  int nr_nodes = e->nr_nodes, nodeID = e->nodeID;
  struct space *s = e->s;
  int my_cells = 0;
  int *cdim = s->cdim;
  struct cell *cells = s->cells;
  int nr_cells = s->nr_cells;

  /* Start by sorting the particles according to their nodes and
     getting the counts. The counts array is indexed as
     count[from * nr_nodes + to]. */
  int *counts, *dest;
  double ih[3], dim[3];
  ih[0] = s->ih[0];
  ih[1] = s->ih[1];
  ih[2] = s->ih[2];
  dim[0] = s->dim[0];
  dim[1] = s->dim[1];
  dim[2] = s->dim[2];
  if ((counts = (int *)malloc(sizeof(int) *nr_nodes *nr_nodes)) == NULL ||
      (dest = (int *)malloc(sizeof(int) * s->nr_parts)) == NULL)
    error("Failed to allocate count and dest buffers.");
  bzero(counts, sizeof(int) * nr_nodes * nr_nodes);
  struct part *parts = s->parts;
  for (int k = 0; k < s->nr_parts; k++) {
    for (int j = 0; j < 3; j++) {
      if (parts[k].x[j] < 0.0)
        parts[k].x[j] += dim[j];
      else if (parts[k].x[j] >= dim[j])
        parts[k].x[j] -= dim[j];
    }
    const int cid = cell_getid(cdim, parts[k].x[0] * ih[0],
                               parts[k].x[1] * ih[1], parts[k].x[2] * ih[2]);
    dest[k] = cells[cid].nodeID;
    counts[nodeID * nr_nodes + dest[k]] += 1;
  }
  parts_sort(s->parts, s->xparts, dest, s->nr_parts, 0, nr_nodes - 1);

  /* Get all the counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, counts, nr_nodes * nr_nodes, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce particle transfer counts.");

  /* Get the new number of parts for this node, be generous in allocating. */
  int nr_parts = 0;
  for (int k = 0; k < nr_nodes; k++) nr_parts += counts[k * nr_nodes + nodeID];
  struct part *parts_new;
  struct xpart *xparts_new, *xparts = s->xparts;
  if (posix_memalign((void **)&parts_new, part_align,
                     sizeof(struct part) * nr_parts * 1.2) != 0 ||
      posix_memalign((void **)&xparts_new, part_align,
                     sizeof(struct xpart) * nr_parts * 1.2) != 0)
    error("Failed to allocate new part data.");

  /* Emit the sends and recvs for the particle data. */
  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 4 * nr_nodes)) ==
      NULL)
    error("Failed to allocate MPI request list.");
  for (int k = 0; k < 4 * nr_nodes; k++) reqs[k] = MPI_REQUEST_NULL;
  for (int offset_send = 0, offset_recv = 0, k = 0; k < nr_nodes; k++) {
    int ind_send = nodeID * nr_nodes + k;
    int ind_recv = k * nr_nodes + nodeID;
    if (counts[ind_send] > 0) {
      if (k == nodeID) {
        memcpy(&parts_new[offset_recv], &s->parts[offset_send],
               sizeof(struct part) * counts[ind_recv]);
        memcpy(&xparts_new[offset_recv], &s->xparts[offset_send],
               sizeof(struct xpart) * counts[ind_recv]);
        offset_send += counts[ind_send];
        offset_recv += counts[ind_recv];
      } else {
        if (MPI_Isend(&s->parts[offset_send],
                      sizeof(struct part) * counts[ind_send], MPI_BYTE, k,
                      2 * ind_send + 0, MPI_COMM_WORLD,
                      &reqs[4 * k]) != MPI_SUCCESS)
          error("Failed to isend parts to node %i.", k);
        if (MPI_Isend(&s->xparts[offset_send],
                      sizeof(struct xpart) * counts[ind_send], MPI_BYTE, k,
                      2 * ind_send + 1, MPI_COMM_WORLD,
                      &reqs[4 * k + 1]) != MPI_SUCCESS)
          error("Failed to isend xparts to node %i.", k);
        offset_send += counts[ind_send];
      }
    }
    if (k != nodeID && counts[ind_recv] > 0) {
      if (MPI_Irecv(&parts_new[offset_recv],
                    sizeof(struct part) * counts[ind_recv], MPI_BYTE, k,
                    2 * ind_recv + 0, MPI_COMM_WORLD,
                    &reqs[4 * k + 2]) != MPI_SUCCESS)
        error("Failed to emit irecv of parts from node %i.", k);
      if (MPI_Irecv(&xparts_new[offset_recv],
                    sizeof(struct xpart) * counts[ind_recv], MPI_BYTE, k,
                    2 * ind_recv + 1, MPI_COMM_WORLD,
                    &reqs[4 * k + 3]) != MPI_SUCCESS)
        error("Failed to emit irecv of parts from node %i.", k);
      offset_recv += counts[ind_recv];
    }
  }

  /* Wait for all the sends and recvs to tumble in. */
  MPI_Status stats[4 * nr_nodes];
  int res;
  if ((res = MPI_Waitall(4 * nr_nodes, reqs, stats)) != MPI_SUCCESS) {
    for (int k = 0; k < 4 * nr_nodes; k++) {
      char buff[MPI_MAX_ERROR_STRING];
      int res;
      MPI_Error_string(stats[k].MPI_ERROR, buff, &res);
      message("request %i has error '%s'.", k, buff);
    }
    error("Failed during waitall for part data.");
  }

  /* Verify that all parts are in the right place. */
  /* for ( k = 0 ; k < nr_parts ; k++ ) {
      cid = cell_getid( cdim , parts_new[k].x[0]*ih[0] , parts_new[k].x[1]*ih[1]
     , parts_new[k].x[2]*ih[2] );
      if ( cells[ cid ].nodeID != nodeID )
          error( "Received particle (%i) that does not belong here (nodeID=%i)."
     , k , cells[ cid ].nodeID );
      } */

  /* Set the new part data, free the old. */
  free(parts);
  free(xparts);
  s->parts = parts_new;
  s->xparts = xparts_new;
  s->nr_parts = nr_parts;
  s->size_parts = 1.2 * nr_parts;

  /* Be verbose about what just happened. */
  for (int k = 0; k < nr_cells; k++)
    if (cells[k].nodeID == nodeID) my_cells += 1;
  message("node %i now has %i parts in %i cells.", nodeID, nr_parts, my_cells);

  /* Clean up other stuff. */
  free(reqs);
  free(counts);
  free(dest);

#else
  error("SWIFT was not compiled with MPI and METIS support.");
#endif
}

/**
 * @brief Repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */

void engine_repartition(struct engine *e) {

#if defined(WITH_MPI) && defined(HAVE_METIS)

  int i, j, k, l, cid, cjd, ii, jj, kk, res;
  idx_t *inds, *nodeIDs;
  idx_t *weights_v = NULL, *weights_e = NULL;
  struct space *s = e->s;
  int nr_cells = s->nr_cells, my_cells = 0;
  struct cell *cells = s->cells;
  int ind[3], *cdim = s->cdim;
  struct task *t, *tasks = e->sched.tasks;
  struct cell *ci, *cj;
  int nr_nodes = e->nr_nodes, nodeID = e->nodeID;
  float wscale = 1e-3, vscale = 1e-3, wscale_buff;
  idx_t wtot = 0;
  const idx_t wmax = 1e9 / e->nr_nodes;

  /* Clear the repartition flag. */
  e->forcerepart = 0;

  /* Allocate the inds and weights. */
  if ((inds = (idx_t *)malloc(sizeof(idx_t) * 26 *nr_cells)) == NULL ||
      (weights_v = (idx_t *)malloc(sizeof(idx_t) *nr_cells)) == NULL ||
      (weights_e = (idx_t *)malloc(sizeof(idx_t) * 26 *nr_cells)) == NULL ||
      (nodeIDs = (idx_t *)malloc(sizeof(idx_t) * nr_cells)) == NULL)
    error("Failed to allocate inds and weights arrays.");

  /* Fill the inds array. */
  for (cid = 0; cid < nr_cells; cid++) {
    ind[0] = cells[cid].loc[0] / s->cells[cid].h[0] + 0.5;
    ind[1] = cells[cid].loc[1] / s->cells[cid].h[1] + 0.5;
    ind[2] = cells[cid].loc[2] / s->cells[cid].h[2] + 0.5;
    l = 0;
    for (i = -1; i <= 1; i++) {
      ii = ind[0] + i;
      if (ii < 0)
        ii += cdim[0];
      else if (ii >= cdim[0])
        ii -= cdim[0];
      for (j = -1; j <= 1; j++) {
        jj = ind[1] + j;
        if (jj < 0)
          jj += cdim[1];
        else if (jj >= cdim[1])
          jj -= cdim[1];
        for (k = -1; k <= 1; k++) {
          kk = ind[2] + k;
          if (kk < 0)
            kk += cdim[2];
          else if (kk >= cdim[2])
            kk -= cdim[2];
          if (i || j || k) {
            inds[cid * 26 + l] = cell_getid(cdim, ii, jj, kk);
            l += 1;
          }
        }
      }
    }
  }

  /* Init the weights arrays. */
  bzero(weights_e, sizeof(idx_t) * 26 * nr_cells);
  bzero(weights_v, sizeof(idx_t) * nr_cells);

  /* Loop over the tasks... */
  for (j = 0; j < e->sched.nr_tasks; j++) {

    /* Get a pointer to the kth task. */
    t = &tasks[j];

    /* Skip un-interesting tasks. */
    if (t->type != task_type_self && t->type != task_type_pair &&
        t->type != task_type_sub && t->type != task_type_ghost &&
        t->type != task_type_drift && t->type != task_type_kick &&
        t->type != task_type_init)
      continue;

    /* Get the task weight. */
    idx_t w = (t->toc - t->tic) * wscale;
    if (w < 0) error("Bad task weight (%" SCIDX ").", w);

    /* Do we need to re-scale? */
    wtot += w;
    while (wtot > wmax) {
      wscale /= 2;
      wtot /= 2;
      w /= 2;
      for (k = 0; k < 26 * nr_cells; k++) weights_e[k] *= 0.5;
      for (k = 0; k < nr_cells; k++) weights_v[k] *= 0.5;
    }

    /* Get the top-level cells involved. */
    for (ci = t->ci; ci->parent != NULL; ci = ci->parent)
      ;
    if (t->cj != NULL)
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent)
        ;
    else
      cj = NULL;

    /* Get the cell IDs. */
    cid = ci - cells;

    /* Different weights for different tasks. */
    if (t->type == task_type_ghost || t->type == task_type_drift ||
        t->type == task_type_kick || t->type == task_type_drift) {

      /* Particle updates add only to vertex weight. */
      weights_v[cid] += w;

    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub && cj == NULL && ci->nodeID == nodeID)) {

      /* Self interactions add only to vertex weight. */
      weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair ||
             (t->type == task_type_sub && cj != NULL)) {

      /* In-cell pair? */
      if (ci == cj) {

        /* Add weight to vertex for ci. */
        weights_v[cid] += w;

      }

      /* Distinct cells with local ci? */
      else if (ci->nodeID == nodeID) {

        /* Index of the jth cell. */
        cjd = cj - cells;

        /* Add half of weight to each cell. */
        if (ci->nodeID == nodeID) weights_v[cid] += 0.5 * w;
        if (cj->nodeID == nodeID) weights_v[cjd] += 0.5 * w;

        /* Add Weight to edge. */
        for (k = 26 * cid; inds[k] != cjd; k++)
          ;
        weights_e[k] += w;
        for (k = 26 * cjd; inds[k] != cid; k++)
          ;
        weights_e[k] += w;
      }
    }
  }

  /* Get the minimum scaling and re-scale if necessary. */
  if ((res = MPI_Allreduce(&wscale, &wscale_buff, 1, MPI_FLOAT, MPI_MIN,
                           MPI_COMM_WORLD)) != MPI_SUCCESS) {
    char buff[MPI_MAX_ERROR_STRING];
    MPI_Error_string(res, buff, &i);
    error("Failed to allreduce the weight scales (%s).", buff);
  }
  if (wscale_buff != wscale) {
    float scale = wscale_buff / wscale;
    for (k = 0; k < 26 * nr_cells; k++) weights_e[k] *= scale;
    for (k = 0; k < nr_cells; k++) weights_v[k] *= scale;
  }

/* Merge the weights arrays across all nodes. */
#if IDXTYPEWIDTH == 32
  if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_v, weights_v,
                        nr_cells, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS) {
#else
  if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_v, weights_v,
                        nr_cells, MPI_LONG_LONG_INT, MPI_SUM, 0,
                        MPI_COMM_WORLD)) != MPI_SUCCESS) {
#endif
    char buff[MPI_MAX_ERROR_STRING];
    MPI_Error_string(res, buff, &i);
    error("Failed to allreduce vertex weights (%s).", buff);
  }
#if IDXTYPEWIDTH == 32
  if (MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_e, weights_e,
                 26 * nr_cells, MPI_INT, MPI_SUM, 0,
                 MPI_COMM_WORLD) != MPI_SUCCESS)
#else
  if (MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_e, weights_e,
                 26 * nr_cells, MPI_LONG_LONG_INT, MPI_SUM, 0,
                 MPI_COMM_WORLD) != MPI_SUCCESS)
#endif
    error("Failed to allreduce edge weights.");

  /* As of here, only one node needs to compute the partition. */
  if (nodeID == 0) {

    /* Check that the edge weights are fully symmetric. */
    /* for ( cid = 0 ; cid < nr_cells ; cid++ )
        for ( k = 0 ; k < 26 ; k++ ) {
            cjd = inds[ cid*26 + k ];
            for ( j = 26*cjd ; inds[j] != cid ; j++ );
            if ( weights_e[ cid*26+k ] != weights_e[ j ] )
                error( "Unsymmetric edge weights detected (%i vs %i)." ,
       weights_e[ cid*26+k ] , weights_e[ j ] );
            } */
    /* int w_min = weights_e[0], w_max = weights_e[0], w_tot = weights_e[0];
    for ( k = 1 ; k < 26*nr_cells ; k++ ) {
        w_tot += weights_e[k];
        if ( weights_e[k] < w_min )
            w_min = weights_e[k];
        else if ( weights_e[k] > w_max )
            w_max = weights_e[k];
        }
    message( "edge weights in [ %i , %i ], tot=%i." , w_min , w_max , w_tot );
    w_min = weights_e[0], w_max = weights_e[0]; w_tot = weights_v[0];
    for ( k = 1 ; k < nr_cells ; k++ ) {
        w_tot += weights_v[k];
        if ( weights_v[k] < w_min )
            w_min = weights_v[k];
        else if ( weights_v[k] > w_max )
            w_max = weights_v[k];
        }
    message( "vertex weights in [ %i , %i ], tot=%i." , w_min , w_max , w_tot );
    */

    /* Make sure there are no zero weights. */
    for (k = 0; k < 26 * nr_cells; k++)
      if (weights_e[k] == 0) weights_e[k] = 1;
    for (k = 0; k < nr_cells; k++)
      if ((weights_v[k] *= vscale) == 0) weights_v[k] = 1;

    /* Allocate and fill the connection array. */
    idx_t *offsets;
    if ((offsets = (idx_t *)malloc(sizeof(idx_t) * (nr_cells + 1))) == NULL)
      error("Failed to allocate offsets buffer.");
    offsets[0] = 0;
    for (k = 0; k < nr_cells; k++) offsets[k + 1] = offsets[k] + 26;

    /* Set the METIS options. +1 to keep the GCC sanitizer happy. */
    idx_t options[METIS_NOPTIONS + 1];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_NCUTS] = 10;
    options[METIS_OPTION_NITER] = 20;
    // options[ METIS_OPTION_UFACTOR ] = 1;

    /* Set the initial partition, although this is probably ignored. */
    for (k = 0; k < nr_cells; k++) nodeIDs[k] = cells[k].nodeID;

    /* Call METIS. */
    idx_t one = 1, idx_nr_cells = nr_cells, idx_nr_nodes = nr_nodes;
    idx_t objval;
    if (METIS_PartGraphRecursive(&idx_nr_cells, &one, offsets, inds, weights_v,
                                 NULL, weights_e, &idx_nr_nodes, NULL, NULL,
                                 options, &objval, nodeIDs) != METIS_OK)
      error("Call to METIS_PartGraphKway failed.");

    /* Dump the 3d array of cell IDs. */
    /* printf( "engine_repartition: nodeIDs = reshape( [" );
    for ( i = 0 ; i < cdim[0]*cdim[1]*cdim[2] ; i++ )
        printf( "%i " , (int)nodeIDs[ i ] );
    printf("] ,%i,%i,%i);\n",cdim[0],cdim[1],cdim[2]); */
  }

/* Broadcast the result of the partition. */
#if IDXTYPEWIDTH == 32
  if (MPI_Bcast(nodeIDs, nr_cells, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to bcast the node IDs.");
#else
  if (MPI_Bcast(nodeIDs, nr_cells, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to bcast the node IDs.");
#endif

  /* Set the cell nodeIDs and clear any non-local parts. */
  for (k = 0; k < nr_cells; k++) {
    cells[k].nodeID = nodeIDs[k];
    if (nodeIDs[k] == nodeID) my_cells += 1;
  }

  /* Clean up. */
  free(inds);
  free(weights_v);
  free(weights_e);
  free(nodeIDs);

  /* Now comes the tricky part: Exchange particles between all nodes.
     This is done in two steps, first allreducing a matrix of
     how many particles go from where to where, then re-allocating
     the parts array, and emitting the sends and receives.
     Finally, the space, tasks, and proxies need to be rebuilt. */

  /* Redistribute the particles between the nodes. */
  engine_redistribute(e);

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Tell the engine it should re-build whenever possible */
  e->forcerebuild = 1;

#else
  error("SWIFT was not compiled with MPI and METIS support.");
#endif
}

/* /\** */
/*  * @brief Add up/down gravity tasks to a cell hierarchy. */
/*  * */
/*  * @param e The #engine. */
/*  * @param c The #cell */
/*  * @param up The upward gravity #task. */
/*  * @param down The downward gravity #task. */
/*  *\/ */

/* void engine_addtasks_grav(struct engine *e, struct cell *c, struct task *up,
 */
/*                           struct task *down) { */

/*   /\* Link the tasks to this cell. *\/ */
/*   c->grav_up = up; */
/*   c->grav_down = down; */

/*   /\* Recurse? *\/ */
/*   if (c->split) */
/*     for (int k = 0; k < 8; k++) */
/*       if (c->progeny[k] != NULL) */
/*         engine_addtasks_grav(e, c->progeny[k], up, down); */
/* } */

/**
 * @brief Add send tasks to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj The receiving #cell
 */

void engine_addtasks_send(struct engine *e, struct cell *ci, struct cell *cj) {

#ifdef WITH_MPI
  int k;
  struct link *l = NULL;
  struct scheduler *s = &e->sched;

  /* Check if any of the density tasks are for the target node. */
  for (l = ci->density; l != NULL; l = l->next)
    if (l->t->ci->nodeID == cj->nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == cj->nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks. */
    struct task *t_xv =
        scheduler_addtask(&e->sched, task_type_send, task_subtype_none,
                          2 * ci->tag, 0, ci, cj, 0);
    struct task *t_rho =
        scheduler_addtask(&e->sched, task_type_send, task_subtype_none,
                          2 * ci->tag + 1, 0, ci, cj, 0);

    /* The send_rho task depends on the cell's ghost task. */
    scheduler_addunlock(s, ci->super->ghost, t_rho);

    /* The send_rho task should unlock the super-cell's kick task. */
    scheduler_addunlock(s, t_rho, ci->super->kick);

    /* The send_xv task should unlock the super-cell's ghost task. */
    scheduler_addunlock(s, t_xv, ci->super->ghost);

  }

  /* Recurse? */
  else if (ci->split)
    for (k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) engine_addtasks_send(e, ci->progeny[k], cj);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param t_xv The recv_xv #task, if it has already been created.
 * @param t_rho The recv_rho #task, if it has already been created.
 */

void engine_addtasks_recv(struct engine *e, struct cell *c, struct task *t_xv,
                          struct task *t_rho) {

#ifdef WITH_MPI
  int k;
  struct scheduler *s = &e->sched;

  /* Do we need to construct a recv task? */
  if (t_xv == NULL && c->nr_density > 0) {

    /* Create the tasks. */
    t_xv = c->recv_xv =
        scheduler_addtask(&e->sched, task_type_recv, task_subtype_none,
                          2 * c->tag, 0, c, NULL, 0);
    t_rho = c->recv_rho =
        scheduler_addtask(&e->sched, task_type_recv, task_subtype_none,
                          2 * c->tag + 1, 0, c, NULL, 0);
  }

  /* Add dependencies. */
  for (struct link *l = c->density; l != NULL; l = l->next) {
    scheduler_addunlock(s, t_xv, l->t);
    scheduler_addunlock(s, l->t, t_rho);
  }
  for (struct link *l = c->force; l != NULL; l = l->next)
    scheduler_addunlock(s, t_rho, l->t);
  if (c->sorts != NULL) scheduler_addunlock(s, t_xv, c->sorts);

  /* Recurse? */
  if (c->split)
    for (k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv(e, c->progeny[k], t_xv, t_rho);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchange cell structures with other nodes.
 *
 * @param e The #engine.
 */

void engine_exchange_cells(struct engine *e) {

#ifdef WITH_MPI

  int j, k, pid, count = 0;
  struct pcell *pcells;
  struct space *s = e->s;
  struct cell *cells = s->cells;
  int nr_cells = s->nr_cells;
  int nr_proxies = e->nr_proxies;
  int offset[nr_cells];
  MPI_Request reqs_in[engine_maxproxies];
  MPI_Request reqs_out[engine_maxproxies];
  MPI_Status status;

  /* Run through the cells and get the size of the ones that will be sent off.
   */
  for (k = 0; k < nr_cells; k++) {
    offset[k] = count;
    if (cells[k].sendto)
      count += (cells[k].pcell_size = cell_getsize(&cells[k]));
  }

  /* Allocate the pcells. */
  if ((pcells = (struct pcell *)malloc(sizeof(struct pcell) * count)) == NULL)
    error("Failed to allocate pcell buffer.");

  /* Pack the cells. */
  cell_next_tag = 0;
  for (k = 0; k < nr_cells; k++)
    if (cells[k].sendto) {
      cell_pack(&cells[k], &pcells[offset[k]]);
      cells[k].pcell = &pcells[offset[k]];
    }

  /* Launch the proxies. */
  for (k = 0; k < nr_proxies; k++) {
    proxy_cells_exch1(&e->proxies[k]);
    reqs_in[k] = e->proxies[k].req_cells_count_in;
    reqs_out[k] = e->proxies[k].req_cells_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (k = 0; k < nr_proxies; k++) {
    if (MPI_Waitany(nr_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "request from proxy %i has arrived." , pid );
    proxy_cells_exch2(&e->proxies[pid]);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(nr_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Set the requests for the cells. */
  for (k = 0; k < nr_proxies; k++) {
    reqs_in[k] = e->proxies[k].req_cells_in;
    reqs_out[k] = e->proxies[k].req_cells_out;
  }

  /* Wait for each pcell array to come in from the proxies. */
  for (k = 0; k < nr_proxies; k++) {
    if (MPI_Waitany(nr_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "cell data from proxy %i has arrived." , pid );
    for (count = 0, j = 0; j < e->proxies[pid].nr_cells_in; j++)
      count += cell_unpack(&e->proxies[pid].pcells_in[count],
                           e->proxies[pid].cells_in[j], e->s);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(nr_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Count the number of particles we need to import and re-allocate
     the buffer if needed. */
  for (count = 0, k = 0; k < nr_proxies; k++)
    for (j = 0; j < e->proxies[k].nr_cells_in; j++)
      count += e->proxies[k].cells_in[j]->count;
  if (count > s->size_parts_foreign) {
    if (s->parts_foreign != NULL) free(s->parts_foreign);
    s->size_parts_foreign = 1.1 * count;
    if (posix_memalign((void **)&s->parts_foreign, part_align,
                       sizeof(struct part) * s->size_parts_foreign) != 0)
      error("Failed to allocate foreign part data.");
  }

  /* Unpack the cells and link to the particle data. */
  struct part *parts = s->parts_foreign;
  for (k = 0; k < nr_proxies; k++) {
    for (count = 0, j = 0; j < e->proxies[k].nr_cells_in; j++) {
      count += cell_link(e->proxies[k].cells_in[j], parts);
      parts = &parts[e->proxies[k].cells_in[j]->count];
    }
  }
  s->nr_parts_foreign = parts - s->parts_foreign;

  /* Is the parts buffer large enough? */
  if (s->nr_parts_foreign > s->size_parts_foreign)
    error("Foreign parts buffer too small.");

  /* Free the pcell buffer. */
  free(pcells);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchange straying parts with other nodes.
 *
 * @param e The #engine.
 * @param offset The index in the parts array as of which the foreign parts
 *reside.
 * @param ind The ID of the foreign #cell.
 * @param N The number of stray parts.
 *
 * @return The number of arrived parts copied to parts and xparts.
 */

int engine_exchange_strays(struct engine *e, int offset, int *ind, int N) {

#ifdef WITH_MPI

  int k, pid, count = 0, nr_in = 0, nr_out = 0;
  MPI_Request reqs_in[2 * engine_maxproxies];
  MPI_Request reqs_out[2 * engine_maxproxies];
  MPI_Status status;
  struct proxy *p;
  struct space *s = e->s;

  /* Re-set the proxies. */
  for (k = 0; k < e->nr_proxies; k++) e->proxies[k].nr_parts_out = 0;

  /* Put the parts into the corresponding proxies. */
  for (k = 0; k < N; k++) {
    int node_id = e->s->cells[ind[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    pid = e->proxy_ind[node_id];
    if (pid < 0)
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%llu, x=[%e,%e,%e].",
          node_id, s->parts[offset + k].id, s->parts[offset + k].x[0],
          s->parts[offset + k].x[1], s->parts[offset + k].x[2]);
    proxy_parts_load(&e->proxies[pid], &s->parts[offset + k],
                     &s->xparts[offset + k], 1);
  }

  /* Launch the proxies. */
  for (k = 0; k < e->nr_proxies; k++) {
    proxy_parts_exch1(&e->proxies[k]);
    reqs_in[k] = e->proxies[k].req_parts_count_in;
    reqs_out[k] = e->proxies[k].req_parts_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (k = 0; k < e->nr_proxies; k++) {
    if (MPI_Waitany(e->nr_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "request from proxy %i has arrived." , pid );
    proxy_parts_exch2(&e->proxies[pid]);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(e->nr_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Count the total number of incoming particles and make sure we have
     enough space to accommodate them. */
  int count_in = 0;
  for (k = 0; k < e->nr_proxies; k++) count_in += e->proxies[k].nr_parts_in;
  message("sent out %i particles, got %i back.", N, count_in);
  if (offset + count_in > s->size_parts) {
    s->size_parts = (offset + count_in) * 1.05;
    struct part *parts_new;
    struct xpart *xparts_new;
    if (posix_memalign((void **)&parts_new, part_align,
                       sizeof(struct part) * s->size_parts) != 0 ||
        posix_memalign((void **)&xparts_new, part_align,
                       sizeof(struct xpart) * s->size_parts) != 0)
      error("Failed to allocate new part data.");
    memcpy(parts_new, s->parts, sizeof(struct part) * offset);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * offset);
    free(s->parts);
    free(s->xparts);
    s->parts = parts_new;
    s->xparts = xparts_new;
  }

  /* Collect the requests for the particle data from the proxies. */
  for (k = 0; k < e->nr_proxies; k++) {
    if (e->proxies[k].nr_parts_in > 0) {
      reqs_in[2 * k] = e->proxies[k].req_parts_in;
      reqs_in[2 * k + 1] = e->proxies[k].req_xparts_in;
      nr_in += 1;
    } else
      reqs_in[2 * k] = reqs_in[2 * k + 1] = MPI_REQUEST_NULL;
    if (e->proxies[k].nr_parts_out > 0) {
      reqs_out[2 * k] = e->proxies[k].req_parts_out;
      reqs_out[2 * k + 1] = e->proxies[k].req_xparts_out;
      nr_out += 1;
    } else
      reqs_out[2 * k] = reqs_out[2 * k + 1] = MPI_REQUEST_NULL;
  }

  /* Wait for each part array to come in and collect the new
     parts from the proxies. */
  for (k = 0; k < 2 * (nr_in + nr_out); k++) {
    int err;
    if ((err = MPI_Waitany(2 * e->nr_proxies, reqs_in, &pid, &status)) !=
        MPI_SUCCESS) {
      char buff[MPI_MAX_ERROR_STRING];
      int res;
      MPI_Error_string(err, buff, &res);
      error("MPI_Waitany failed (%s).", buff);
    }
    if (pid == MPI_UNDEFINED) break;
    // message( "request from proxy %i has arrived." , pid );
    if (reqs_in[pid & ~1] == MPI_REQUEST_NULL &&
        reqs_in[pid | 1] == MPI_REQUEST_NULL) {
      p = &e->proxies[pid >> 1];
      memcpy(&s->parts[offset + count], p->parts_in,
             sizeof(struct part) * p->nr_parts_in);
      memcpy(&s->xparts[offset + count], p->xparts_in,
             sizeof(struct xpart) * p->nr_parts_in);
      /* for (int k = offset; k < offset + count; k++)
         message(
            "received particle %lli, x=[%.3e %.3e %.3e], h=%.3e, from node %i.",
            s->parts[k].id, s->parts[k].x[0], s->parts[k].x[1],
            s->parts[k].x[2], s->parts[k].h, p->nodeID); */
      count += p->nr_parts_in;
    }
  }

  /* Wait for all the sends to have finished too. */
  if (nr_out > 0)
    if (MPI_Waitall(2 * e->nr_proxies, reqs_out, MPI_STATUSES_IGNORE) !=
        MPI_SUCCESS)
      error("MPI_Waitall on sends failed.");

  /* Return the number of harvested parts. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Fill the #space's task list.
 *
 * @param e The #engine we are working with.
 */

void engine_maketasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  struct cell *cells = s->cells;
  int nr_cells = s->nr_cells;
  int nodeID = e->nodeID;
  int i, j, k, ii, jj, kk, iii, jjj, kkk, cid, cjd, sid;
  int *cdim = s->cdim;
  struct task *t, *t2;
  struct cell *ci, *cj;

  /* Re-set the scheduler. */
  scheduler_reset(sched, s->tot_cells * engine_maxtaskspercell);

  /* Run through the highest level of cells and add pairs. */
  for (i = 0; i < cdim[0]; i++)
    for (j = 0; j < cdim[1]; j++)
      for (k = 0; k < cdim[2]; k++) {
        cid = cell_getid(cdim, i, j, k);
        if (cells[cid].count == 0) continue;
        ci = &cells[cid];
        if (ci->count == 0) continue;
        if (ci->nodeID == nodeID)
          scheduler_addtask(sched, task_type_self, task_subtype_density, 0, 0,
                            ci, NULL, 0);
        for (ii = -1; ii < 2; ii++) {
          iii = i + ii;
          if (!s->periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (jj = -1; jj < 2; jj++) {
            jjj = j + jj;
            if (!s->periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (kk = -1; kk < 2; kk++) {
              kkk = k + kk;
              if (!s->periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];
              cjd = cell_getid(cdim, iii, jjj, kkk);
              cj = &cells[cjd];
              if (cid >= cjd || cj->count == 0 ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;
              sid = sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
              scheduler_addtask(sched, task_type_pair, task_subtype_density,
                                sid, 0, ci, cj, 1);
            }
          }
        }
      }

  /* /\* Add the gravity mm tasks. *\/ */
  /* for (i = 0; i < nr_cells; i++) */
  /*   if (cells[i].gcount > 0) { */
  /*     scheduler_addtask(sched, task_type_grav_mm, task_subtype_none, -1, 0,
   */
  /*                       &cells[i], NULL, 0); */
  /*     for (j = i + 1; j < nr_cells; j++) */
  /*       if (cells[j].gcount > 0) */
  /*         scheduler_addtask(sched, task_type_grav_mm, task_subtype_none, -1,
   * 0, */
  /*                           &cells[i], &cells[j], 0); */
  /*   } */

  /* Split the tasks. */
  scheduler_splittasks(sched);

  /* Allocate the list of cell-task links. The maximum number of links
     is the number of cells (s->tot_cells) times the number of neighbours (27)
     times the number of interaction types (2, density and force). */
  free(e->links);
  e->size_links = s->tot_cells * 27 * 2;
  if ((e->links = malloc(sizeof(struct link) * e->size_links)) == NULL)
    error("Failed to allocate cell-task links.");
  e->nr_links = 0;

  space_link_cleanup(s);

  /* /\* Add the gravity up/down tasks at the top-level cells and push them
   * down. *\/ */
  /* for (k = 0; k < nr_cells; k++) */
  /*   if (cells[k].nodeID == nodeID && cells[k].gcount > 0) { */

  /*     /\* Create tasks at top level. *\/ */
  /*     struct task *up = */
  /*         scheduler_addtask(sched, task_type_grav_up, task_subtype_none, 0,
   * 0, */
  /*                           &cells[k], NULL, 0); */
  /*     struct task *down = */
  /*         scheduler_addtask(sched, task_type_grav_down, task_subtype_none, 0,
   * 0, */
  /*                           &cells[k], NULL, 0); */

  /*     /\* Push tasks down the cell hierarchy. *\/ */
  /*     engine_addtasks_grav(e, &cells[k], up, down); */
  /*   } */

  message("nb tasks: %d", sched->nr_tasks);

  /* Count the number of tasks associated with each cell and
     store the density tasks in each cell, and make each sort
     depend on the sorts of its progeny. */
  for (k = 0; k < sched->nr_tasks; k++) {

    /* Get the current task. */
    t = &sched->tasks[k];
    if (t->skip) continue;

    /* Link sort tasks together. */
    if (t->type == task_type_sort && t->ci->split)
      for (j = 0; j < 8; j++)
        if (t->ci->progeny[j] != NULL && t->ci->progeny[j]->sorts != NULL) {
          t->ci->progeny[j]->sorts->skip = 0;
          scheduler_addunlock(sched, t->ci->progeny[j]->sorts, t);
        }

    /* Link density tasks to cells. */
    if (t->type == task_type_self) {
      atomic_inc(&t->ci->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->link_density = engine_addlink(e, t->ci->link_density, t);
        atomic_inc(&t->ci->nr_link_density);
      }
    } else if (t->type == task_type_pair) {
      atomic_inc(&t->ci->nr_tasks);
      atomic_inc(&t->cj->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->link_density = engine_addlink(e, t->ci->link_density, t);
        atomic_inc(&t->ci->nr_link_density);
        t->cj->link_density = engine_addlink(e, t->cj->link_density, t);
        atomic_inc(&t->cj->nr_link_density);
      }
    } else if (t->type == task_type_sub) {
      atomic_inc(&t->ci->nr_tasks);
      if (t->cj != NULL) atomic_inc(&t->cj->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->link_density = engine_addlink(e, t->ci->link_density, t);
        atomic_inc(&t->ci->nr_link_density);
        if (t->cj != NULL) {
          t->cj->link_density = engine_addlink(e, t->cj->link_density, t);
          atomic_inc(&t->cj->nr_link_density);
        }
      }
    }

    /* /\* Link gravity multipole tasks to the up/down tasks. *\/ */
    /* if (t->type == task_type_grav_mm || */
    /*     (t->type == task_type_sub && t->subtype == task_subtype_grav)) { */
    /*   atomic_inc(&t->ci->nr_tasks); */
    /*   scheduler_addunlock(sched, t->ci->grav_up, t); */
    /*   scheduler_addunlock(sched, t, t->ci->grav_down); */
    /*   if (t->cj != NULL && t->ci->grav_up != t->cj->grav_up) { */
    /*     scheduler_addunlock(sched, t->cj->grav_up, t); */
    /*     scheduler_addunlock(sched, t, t->cj->grav_down); */
    /*   } */
    /* } */
  }

  /* Append a ghost task to each cell, and add kick tasks to the
     super cells. */
  for (k = 0; k < nr_cells; k++) engine_mkghosts(e, &cells[k], NULL);

  /* Run through the tasks and make force tasks for each density task.
     Each force task depends on the cell ghosts and unlocks the kick task
     of its super-cell. */
  kk = sched->nr_tasks;
  for (k = 0; k < kk; k++) {

    /* Get a pointer to the task. */
    t = &sched->tasks[k];

    /* Skip? */
    if (t->skip) continue;

    /* Self-interaction? */
    if (t->type == task_type_self && t->subtype == task_subtype_density) {
      scheduler_addunlock(sched, t->ci->super->init, t);
      scheduler_addunlock(sched, t, t->ci->super->ghost);
      t2 = scheduler_addtask(sched, task_type_self, task_subtype_force, 0, 0,
                             t->ci, NULL, 0);
      scheduler_addunlock(sched, t->ci->super->ghost, t2);
      scheduler_addunlock(sched, t2, t->ci->super->kick);
      t->ci->link_force = engine_addlink(e, t->ci->link_force, t2);
      atomic_inc(&t->ci->nr_link_force);
    }

    /* Otherwise, pair interaction? */
    else if (t->type == task_type_pair && t->subtype == task_subtype_density) {
      t2 = scheduler_addtask(sched, task_type_pair, task_subtype_force, 0, 0,
                             t->ci, t->cj, 0);
      if (t->ci->nodeID == nodeID) {
        scheduler_addunlock(sched, t->ci->super->init, t);
        scheduler_addunlock(sched, t, t->ci->super->ghost);
        scheduler_addunlock(sched, t->ci->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->ci->super->kick);
      }
      if (t->cj->nodeID == nodeID && t->ci->super != t->cj->super) {
        scheduler_addunlock(sched, t->cj->super->init, t);
        scheduler_addunlock(sched, t, t->cj->super->ghost);
        scheduler_addunlock(sched, t->cj->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->cj->super->kick);
      }
      t->ci->link_force = engine_addlink(e, t->ci->link_force, t2);
      atomic_inc(&t->ci->nr_link_force);
      t->cj->link_force = engine_addlink(e, t->cj->link_force, t2);
      atomic_inc(&t->cj->nr_link_force);
    }

    /* Otherwise, sub interaction? */
    else if (t->type == task_type_sub && t->subtype == task_subtype_density) {
      t2 = scheduler_addtask(sched, task_type_sub, task_subtype_force, t->flags,
                             0, t->ci, t->cj, 0);
      if (t->ci->nodeID == nodeID) {
        scheduler_addunlock(sched, t->ci->super->init, t);
        scheduler_addunlock(sched, t, t->ci->super->ghost);
        scheduler_addunlock(sched, t->ci->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->ci->super->kick);
      }
      if (t->cj != NULL && t->cj->nodeID == nodeID &&
          t->ci->super != t->cj->super) {
        scheduler_addunlock(sched, t->cj->super->init, t);
        scheduler_addunlock(sched, t, t->cj->super->ghost);
        scheduler_addunlock(sched, t->cj->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->cj->super->kick);
      }
      t->ci->link_force = engine_addlink(e, t->ci->link_force, t2);
      atomic_inc(&t->ci->nr_link_force);
      if (t->cj != NULL) {
        t->cj->link_force = engine_addlink(e, t->cj->link_force, t2);
        atomic_inc(&t->cj->nr_link_force);
      }
    }

    /* /\* Kick tasks should rely on the grav_down tasks of their cell. *\/ */
    /* else if (t->type == task_type_kick && t->ci->grav_down != NULL) */
    /*   scheduler_addunlock(sched, t->ci->grav_down, t); */
  }

/* Add the communication tasks if MPI is being used. */
#ifdef WITH_MPI

  /* Loop over the proxies. */
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    struct proxy *p = &e->proxies[pid];

    /* Loop through the proxy's incoming cells and add the
       recv tasks. */
    for (k = 0; k < p->nr_cells_in; k++)
      engine_addtasks_recv(e, p->cells_in[k], NULL, NULL);

    /* Loop through the proxy's outgoing cells and add the
       send tasks. */
    for (k = 0; k < p->nr_cells_out; k++)
      engine_addtasks_send(e, p->cells_out[k], p->cells_in[0]);
  }

#endif

  /* Rank the tasks. */
  scheduler_ranktasks(sched);

  /* Weight the tasks. */
  scheduler_reweight(sched);

  /* Set the tasks age. */
  e->tasks_age = 0;
}

/**
 * @brief Mark tasks to be skipped and set the sort flags accordingly.
 *
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */

int engine_marktasks(struct engine *e) {

  struct scheduler *s = &e->sched;
  int k, nr_tasks = s->nr_tasks, *ind = s->tasks_ind;
  struct task *t, *tasks = s->tasks;
  // float t_end = e->time;
  struct cell *ci, *cj;
  // ticks tic = getticks();

  /* /\* Much less to do here if we're on a fixed time-step. *\/ */
  /* if (!(e->policy & engine_policy_multistep)) { */

  /*   /\* Run through the tasks and mark as skip or not. *\/ */
  /*   for (k = 0; k < nr_tasks; k++) { */

  /*     /\* Get a handle on the kth task. *\/ */
  /*     t = &tasks[ind[k]]; */

  /*     /\* Pair? *\/ */
  /*     if (t->type == task_type_pair || */
  /*         (t->type == task_type_sub && t->cj != NULL)) { */

  /*       /\* Local pointers. *\/ */
  /*       ci = t->ci; */
  /*       cj = t->cj; */

  /*       /\* Too much particle movement? *\/ */
  /*       if (t->tight && */
  /*           (fmaxf(ci->h_max, cj->h_max) + ci->dx_max + cj->dx_max > cj->dmin
   * || */
  /*            ci->dx_max > space_maxreldx * ci->h_max || */
  /*            cj->dx_max > space_maxreldx * cj->h_max)) */
  /*         return 1; */

  /*     } */

  /*     /\* Sort? *\/ */
  /*     else if (t->type == task_type_sort) { */

  /*       /\* If all the sorts have been done, make this task implicit. *\/ */
  /*       if (!(t->flags & (t->flags ^ t->ci->sorted))) t->implicit = 1; */
  /*     } */
  /*   } */

  /* } else { */

  /* Run through the tasks and mark as skip or not. */
  for (k = 0; k < nr_tasks; k++) {

    /* Get a handle on the kth task. */
    t = &tasks[ind[k]];

    /* Sort-task? Note that due to the task ranking, the sorts
       will all come before the pairs. */
    if (t->type == task_type_sort) {

      /* Re-set the flags. */
      t->flags = 0;
      t->skip = 1;

    }

    /* Single-cell task? */
    else if (t->type == task_type_self || t->type == task_type_ghost ||
             (t->type == task_type_sub && t->cj == NULL)) {

      /* Set this task's skip. */
      // t->skip = (t->ci->t_end_min >= t_end);

    }

    /* Pair? */
    else if (t->type == task_type_pair ||
             (t->type == task_type_sub && t->cj != NULL)) {

      /* Local pointers. */
      ci = t->ci;
      cj = t->cj;

      /* Set this task's skip. */
      // t->skip = (ci->t_end_min >= t_end && cj->t_end_min >= t_end);

      /* Too much particle movement? */
      if (t->tight &&
          (fmaxf(ci->h_max, cj->h_max) + ci->dx_max + cj->dx_max > cj->dmin ||
           ci->dx_max > space_maxreldx * ci->h_max ||
           cj->dx_max > space_maxreldx * cj->h_max))
        return 1;

      /* Set the sort flags. */
      if (!t->skip && t->type == task_type_pair) {
        if (!(ci->sorted & (1 << t->flags))) {
          ci->sorts->flags |= (1 << t->flags);
          ci->sorts->skip = 0;
        }
        if (!(cj->sorted & (1 << t->flags))) {
          cj->sorts->flags |= (1 << t->flags);
          cj->sorts->skip = 0;
        }
      }

    }

    /* Kick? */
    else if (t->type == task_type_kick)
      t->skip = 0;

    /* Drift? */
    else if (t->type == task_type_drift)
      t->skip = 0;

    /* Init? */
    else if (t->type == task_type_init)
      t->skip = 0;

    /* None? */
    else if (t->type == task_type_none)
      t->skip = 1;
  }
  //}

  // message( "took %.3f ms." , (double)(getticks() - tic)/CPU_TPS*1000 );

  /* All is well... */
  return 0;
}

/**
 * @brief Prints the number of tasks in the engine
 *
 * @param e The #engine.
 */

void engine_print(struct engine *e) {

  int k;
  struct scheduler *sched = &e->sched;

  /* Count and print the number of each task type. */
  int counts[task_type_count + 1];
  for (k = 0; k <= task_type_count; k++) counts[k] = 0;
  for (k = 0; k < sched->nr_tasks; k++)
    if (!sched->tasks[k].skip)
      counts[(int)sched->tasks[k].type] += 1;
    else
      counts[task_type_count] += 1;
#ifdef WITH_MPI
  printf("[%03i] engine_print: task counts are [ %s=%i", e->nodeID,
         taskID_names[0], counts[0]);
#else
  printf("engine_print: task counts are [ %s=%i", taskID_names[0], counts[0]);
#endif
  for (k = 1; k < task_type_count; k++)
    printf(" %s=%i", taskID_names[k], counts[k]);
  printf(" skipped=%i ]\n", counts[task_type_count]);
  fflush(stdout);
  message("nr_parts = %i.", e->s->nr_parts);
}

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 */

void engine_rebuild(struct engine *e) {

  message("REBUILD !!!");

  /* Clear the forcerebuild flag, whatever it was. */
  e->forcerebuild = 0;

  /* Re-build the space. */
  // tic = getticks();
  space_rebuild(e->s, 0.0, e->nodeID == 0);
// message( "space_rebuild took %.3f ms." , (double)(getticks() -
// tic)/CPU_TPS*1000 );

/* If in parallel, exchange the cell structure. */
#ifdef WITH_MPI
  // tic = getticks();
  engine_exchange_cells(e);
// message( "engine_exchange_cells took %.3f ms." , (double)(getticks() -
// tic)/CPU_TPS*1000 );
#endif

  /* Re-build the tasks. */
  // tic = getticks();
  engine_maketasks(e);
  // message( "engine_maketasks took %.3f ms." , (double)(getticks() -
  // tic)/CPU_TPS*1000 );

  /* Run through the tasks and mark as skip or not. */
  // tic = getticks();
  if (engine_marktasks(e))
    error("engine_marktasks failed after space_rebuild.");
  // message( "engine_marktasks took %.3f ms." , (double)(getticks() -
  // tic)/CPU_TPS*1000 );

  /* Print the status of the system */
  engine_print(e);
}

/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 */

void engine_prepare(struct engine *e) {

  int rebuild;

  TIMER_TIC

  /* Run through the tasks and mark as skip or not. */
  // tic = getticks();
  rebuild = (e->forcerebuild || engine_marktasks(e));
// message( "space_marktasks took %.3f ms." , (double)(getticks() -
// tic)/CPU_TPS*1000 );

/* Collect the values of rebuild from all nodes. */
#ifdef WITH_MPI
  // tic = getticks();
  int buff;
  if (MPI_Allreduce(&rebuild, &buff, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to aggregate the rebuild flag across nodes.");
  rebuild = buff;
// message( "rebuild allreduce took %.3f ms." , (double)(getticks() -
// tic)/CPU_TPS*1000 );
#endif
  e->tic_step = getticks();

  /* Did this not go through? */
  if (rebuild) {
    // tic = getticks();
    engine_rebuild(e);
    // message( "engine_rebuild took %.3f ms." , (double)(getticks() -
    // tic)/CPU_TPS*1000 );
  }

  /* Re-rank the tasks every now and then. */
  if (e->tasks_age % engine_tasksreweight == 1) {
    // tic = getticks();
    scheduler_reweight(&e->sched);
    // message( "scheduler_reweight took %.3f ms." , (double)(getticks() -
    // tic)/CPU_TPS*1000 );
  }
  e->tasks_age += 1;

  TIMER_TOC(timer_prepare);
}

/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 * @param tid The thread ID
 */

void engine_barrier(struct engine *e, int tid) {

  /* First, get the barrier mutex. */
  if (pthread_mutex_lock(&e->barrier_mutex) != 0)
    error("Failed to get barrier mutex.");

  /* This thread is no longer running. */
  e->barrier_running -= 1;

  /* If all threads are in, send a signal... */
  if (e->barrier_running == 0)
    if (pthread_cond_broadcast(&e->barrier_cond) != 0)
      error("Failed to broadcast barrier full condition.");

  /* Wait for the barrier to open. */
  while (e->barrier_launch == 0 || tid >= e->barrier_launchcount)
    if (pthread_cond_wait(&e->barrier_cond, &e->barrier_mutex) != 0)
      error("Error waiting for barrier to close.");

  /* This thread has been launched. */
  e->barrier_running += 1;
  e->barrier_launch -= 1;

  /* If I'm the last one out, signal the condition again. */
  if (e->barrier_launch == 0)
    if (pthread_cond_broadcast(&e->barrier_cond) != 0)
      error("Failed to broadcast empty barrier condition.");

  /* Last but not least, release the mutex. */
  if (pthread_mutex_unlock(&e->barrier_mutex) != 0)
    error("Failed to get unlock the barrier mutex.");
}

/**
 * @brief Mapping function to collect the data from the second kick.
 */

void engine_collect_kick(struct cell *c) {

  int k, updated = 0;
  float t_end_min = FLT_MAX, t_end_max = 0.0f;
  double ekin = 0.0, epot = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f}, ang[3] = {0.0f, 0.0f, 0.0f};
  struct cell *cp;

  /* If I am a super-cell, return immediately. */
  if (c->kick != NULL || c->count == 0) return;

  /* If this cell is not split, I'm in trouble. */
  if (!c->split) error("Cell has no super-cell.");

  /* Collect the values from the progeny. */
  for (k = 0; k < 8; k++)
    if ((cp = c->progeny[k]) != NULL) {
      engine_collect_kick(cp);
      t_end_min = fminf(t_end_min, cp->t_end_min);
      t_end_max = fmaxf(t_end_max, cp->t_end_max);
      updated += cp->updated;
      ekin += cp->ekin;
      epot += cp->epot;
      mom[0] += cp->mom[0];
      mom[1] += cp->mom[1];
      mom[2] += cp->mom[2];
      ang[0] += cp->ang[0];
      ang[1] += cp->ang[1];
      ang[2] += cp->ang[2];
    }

  /* Store the collected values in the cell. */
  c->t_end_min = t_end_min;
  c->t_end_max = t_end_max;
  c->updated = updated;
  c->ekin = ekin;
  c->epot = epot;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];
}

/**
 * @brief Compute the force on a single particle brute-force.
 */

// void engine_single_density ( double *dim , long long int pid , struct part
// *__restrict__ parts , int N , int periodic ) {
//
//     int i, k;
//     double r2, dx[3];
//     float fdx[3], ih;
//     struct part p;
//
//     /* Find "our" part. */
//     for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
//     if ( k == N )
//         error( "Part not found." );
//     p = parts[k];
//
//     /* Clear accumulators. */
//     ih = 1.0f / p.h;
//     p.rho = 0.0f; p.rho_dh = 0.0f;
//     p.density.wcount = 0.0f; p.density.wcount_dh = 0.0f;
//     p.density.div_v = 0.0;
//     for ( k=0 ; k < 3 ; k++)
//         p.density.curl_v[k] = 0.0;
//
//     /* Loop over all particle pairs (force). */
//     for ( k = 0 ; k < N ; k++ ) {
//         if ( parts[k].id == p.id )
//             continue;
//         for ( i = 0 ; i < 3 ; i++ ) {
//             dx[i] = p.x[i] - parts[k].x[i];
//             if ( periodic ) {
//                 if ( dx[i] < -dim[i]/2 )
//                     dx[i] += dim[i];
//                 else if ( dx[i] > dim[i]/2 )
//                     dx[i] -= dim[i];
//                 }
//             fdx[i] = dx[i];
//             }
//         r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
//         if ( r2 < p.h*p.h*kernel_gamma2 ) {
//             runner_iact_nonsym_density( r2 , fdx , p.h , parts[k].h , &p ,
// &parts[k] );
//             }
//         }
//
//     /* Dump the result. */
//     p.rho = ih * ih * ih * ( p.rho + p.mass*kernel_root );
//     p.rho_dh = p.rho_dh * ih * ih * ih * ih;
//     p.density.wcount = ( p.density.wcount + kernel_root ) * ( 4.0f / 3.0 *
// M_PI * kernel_gamma3 );
//     message( "part %lli (h=%e) has wcount=%e, rho=%e, rho_dh=%e." , p.id ,
// p.h , p.density.wcount , p.rho , p.rho_dh );
//     fflush(stdout);
//
//     }

// void engine_single_force ( double *dim , long long int pid , struct part
// *__restrict__ parts , int N , int periodic ) {
//
//     int i, k;
//     double r2, dx[3];
//     float fdx[3];
//     struct part p;
//
//     /* Find "our" part. */
//     for ( k = 0 ; k < N && parts[k].id != pid ; k++ );
//     if ( k == N )
//         error( "Part not found." );
//     p = parts[k];
//
//     /* Clear accumulators. */
//     p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
//     p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
//
//     /* Loop over all particle pairs (force). */
//     for ( k = 0 ; k < N ; k++ ) {
//     // for ( k = N-1 ; k >= 0 ; k-- ) {
//         if ( parts[k].id == p.id )
//             continue;
//         for ( i = 0 ; i < 3 ; i++ ) {
//             dx[i] = p.x[i] - parts[k].x[i];
//             if ( periodic ) {
//                 if ( dx[i] < -dim[i]/2 )
//                     dx[i] += dim[i];
//                 else if ( dx[i] > dim[i]/2 )
//                     dx[i] -= dim[i];
//                 }
//             fdx[i] = dx[i];
//             }
//         r2 = fdx[0]*fdx[0] + fdx[1]*fdx[1] + fdx[2]*fdx[2];
//         if ( r2 < p.h*p.h*kernel_gamma2 || r2 <
// parts[k].h*parts[k].h*kernel_gamma2 ) {
//             p.a[0] = 0.0f; p.a[1] = 0.0f; p.a[2] = 0.0f;
//             p.force.u_dt = 0.0f; p.force.h_dt = 0.0f; p.force.v_sig = 0.0f;
//             runner_iact_nonsym_force( r2 , fdx , p.h , parts[k].h , &p ,
// &parts[k] );
//             double dvdr = ( (p.v[0]-parts[k].v[0])*fdx[0] +
// (p.v[1]-parts[k].v[1])*fdx[1] + (p.v[2]-parts[k].v[2])*fdx[2] ) / sqrt(r2);
//             message( "part %lli and %lli interact (r=%.3e,dvdr=%.3e) with
// a=[%.3e,%.3e,%.3e], dudt=%.3e." ,
//                 p.id , parts[k].id , sqrt(r2) , dvdr , p.a[0] , p.a[1],
// p.a[2] , p.force.u_dt );
//             }
//         }
//
//     /* Dump the result. */
//     // message( "part %lli (h=%e) has a=[%.3e,%.3e,%.3e], udt=%e." , p.id ,
// p.h , p.a[0] , p.a[1] , p.a[2] , p.force.u_dt );
//     fflush(stdout);
//
//     }

/**
 * @brief Launch the runners.
 *
 * @param e The #engine.
 * @param nr_runners The number of #runner to let loose.
 * @param mask The task mask to launch.
 */

void engine_launch(struct engine *e, int nr_runners, unsigned int mask) {

  /* Prepare the scheduler. */
  atomic_inc(&e->sched.waiting);

  /* Cry havoc and let loose the dogs of war. */
  e->barrier_launch = nr_runners;
  e->barrier_launchcount = nr_runners;
  if (pthread_cond_broadcast(&e->barrier_cond) != 0)
    error("Failed to broadcast barrier open condition.");

  /* Print out what we do */
  printf("\n");
  task_print_mask(mask);
  printf("\n");
  fflush(stdout);

  engine_print(e);

  /* Load the tasks. */
  pthread_mutex_unlock(&e->barrier_mutex);
  scheduler_start(&e->sched, mask);
  pthread_mutex_lock(&e->barrier_mutex);

  /* Remove the safeguard. */
  pthread_mutex_lock(&e->sched.sleep_mutex);
  atomic_dec(&e->sched.waiting);
  pthread_cond_broadcast(&e->sched.sleep_cond);
  pthread_mutex_unlock(&e->sched.sleep_mutex);

  message("Before barrier");
  fflush(stdout);

  /* Sit back and wait for the runners to come home. */
  while (e->barrier_launch || e->barrier_running)
    if (pthread_cond_wait(&e->barrier_cond, &e->barrier_mutex) != 0)
      error("Error while waiting for barrier.");

  message("After barrier");
  fflush(stdout);
}

/* void hassorted(struct cell *c) { */

/*   if (c->sorted) error("Suprious sorted flags."); */

/*   if (c->split) */
/*     for (int k = 0; k < 8; k++) */
/*       if (c->progeny[k] != NULL) hassorted(c->progeny[k]); */
/* } */

/**
 * @brief Initialises the particles and set them in a state ready to move
 *forward in time.
 *
 * @param e The #engine
 */
void engine_init_particles(struct engine *e) {

  struct space *s = e->s;

  // engine_repartition(e);

  engine_prepare(e);

  engine_print(e);

  // engine_maketasks(e);

  engine_print(e);

  engine_marktasks(e);

  engine_print(e);

  fflush(stdout);
  message("Engine prepared");

  /* Make sure all particles are ready to go */
  void initParts(struct part * p, struct xpart * xp, struct cell * c) {
    p->t_begin = 0.;
    p->t_end = 0.;
    p->rho = -1.;
    xp->v_full[0] = p->v[0];
    xp->v_full[1] = p->v[1];
    xp->v_full[2] = p->v[2];
    c->t_end_min = 0.;
  }

  message("Initialising particles");
  space_map_parts_xparts(s, initParts);

  /* Now everybody should have sensible smoothing length */
  void printParts(struct part * p, struct xpart * xp, struct cell * c) {
    if (p->id == 1000)
      message("id=%lld h=%f rho=%f t_begin=%f t_end=%f", p->id, p->h, p->rho,
              p->t_begin, p->t_end);
  }
  // space_map_parts_xparts(s, printParts);

  void printCells(struct part * p, struct xpart * xp, struct cell * c) {
    if (c->super != NULL && 0)
      message(
          "c->t_end_min=%f c->t_end_max=%f c->super=%p sort=%p ghost=%p "
          "kick=%p",
          c->t_end_min, c->t_end_max, c->super, c->sorts, c->ghost, c->kick);
  }
  // space_map_parts_xparts(s, printCells);

  /* Now do a density calculation */
  TIMER_TIC;
  engine_launch(e, e->nr_threads,
                (1 << task_type_sort) | (1 << task_type_self) |
                    (1 << task_type_pair) | (1 << task_type_sub) |
                    (1 << task_type_init) | (1 << task_type_ghost) |
                    (1 << task_type_send) | (1 << task_type_recv) |
                    (1 << task_type_link));

  TIMER_TOC(timer_runners);

  // space_map_parts_xparts(s, printParts);

  printf("\n\n");

  /* Ready to go */
  e->step = -1;
}

/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 */
void engine_step(struct engine *e) {

  int k;
  float t_end_min = FLT_MAX, t_end_max = 0.f;
  double epot = 0.0, ekin = 0.0;
  float mom[3] = {0.0, 0.0, 0.0};
  float ang[3] = {0.0, 0.0, 0.0};
  int count = 0;
  struct cell *c;
  struct space *s = e->s;

  TIMER_TIC2;

  /* Collect the cell data. */
  for (k = 0; k < s->nr_cells; k++)
    if (s->cells[k].nodeID == e->nodeID) {
      c = &s->cells[k];
      engine_collect_kick(c);
      t_end_min = fminf(t_end_min, c->t_end_min);
      t_end_max = fmaxf(t_end_max, c->t_end_max);
      ekin += c->ekin;
      epot += c->epot;
      count += c->updated;
      mom[0] += c->mom[0];
      mom[1] += c->mom[1];
      mom[2] += c->mom[2];
      ang[0] += c->ang[0];
      ang[1] += c->ang[1];
      ang[2] += c->ang[2];
    }

/* Aggregate the data from the different nodes. */
#ifdef WITH_MPI
  double in[3], out[3];
  out[0] = t_end_min;
  if (MPI_Allreduce(out, in, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to aggregate t_end_min.");
  t_end_min = in[0];
  out[0] = t_end_max;
  if (MPI_Allreduce(out, in, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to aggregate t_end_max.");
  t_end_max = in[0];
  out[0] = count;
  out[1] = ekin;
  out[2] = epot;
  if (MPI_Allreduce(out, in, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to aggregate energies.");
  count = in[0];
  ekin = in[1];
  epot = in[2];
#endif

  message("t_end_min=%f t_end_max=%f", t_end_min, t_end_max);

  /* Move forward in time */
  e->timeOld = e->time;
  e->time = t_end_min;
  e->step += 1;
  message("Step: %d e->time=%f", e->step, e->time);

  /* Drift everybody */
  engine_launch(e, e->nr_threads, 1 << task_type_drift);

  /* Re-distribute the particles amongst the nodes? */
  if (e->forcerepart) engine_repartition(e);

  /* Prepare the space. */
  engine_prepare(e);

  engine_maketasks(e);

  // engine_marktasks(e);

  engine_print(e);
  message("Go !");
  fflush(stdout);

  /* Send off the runners. */
  TIMER_TIC;
  engine_launch(e, e->nr_threads,
                (1 << task_type_sort) | (1 << task_type_self) |
                    (1 << task_type_pair) | (1 << task_type_sub) |
                    (1 << task_type_init) | (1 << task_type_ghost) |
                    (1 << task_type_kick) | (1 << task_type_send) |
                    (1 << task_type_recv) | (1 << task_type_link));

  TIMER_TOC(timer_runners);

  TIMER_TOC2(timer_step);
}

/**
 * @brief Create and fill the proxies.
 *
 * @param e The #engine.
 */

void engine_makeproxies(struct engine *e) {

  int i, j, k, ii, jj, kk;
  int cid, cjd, pid, ind[3], *cdim = e->s->cdim;
  struct space *s = e->s;
  struct cell *cells = s->cells;
  struct proxy *proxies = e->proxies;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* The following loop is super-clunky, but it's necessary
     to ensure that the order of the send and recv cells in
     the proxies is identical for all nodes! */

  /* Loop over each cell in the space. */
  for (ind[0] = 0; ind[0] < cdim[0]; ind[0]++)
    for (ind[1] = 0; ind[1] < cdim[1]; ind[1]++)
      for (ind[2] = 0; ind[2] < cdim[2]; ind[2]++) {

        /* Get the cell ID. */
        cid = cell_getid(cdim, ind[0], ind[1], ind[2]);

        /* Loop over all its neighbours (periodic). */
        for (i = -1; i <= 1; i++) {
          ii = ind[0] + i;
          if (ii >= cdim[0])
            ii -= cdim[0];
          else if (ii < 0)
            ii += cdim[0];
          for (j = -1; j <= 1; j++) {
            jj = ind[1] + j;
            if (jj >= cdim[1])
              jj -= cdim[1];
            else if (jj < 0)
              jj += cdim[1];
            for (k = -1; k <= 1; k++) {
              kk = ind[2] + k;
              if (kk >= cdim[2])
                kk -= cdim[2];
              else if (kk < 0)
                kk += cdim[2];

              /* Get the cell ID. */
              cjd = cell_getid(cdim, ii, jj, kk);

              /* Add to proxies? */
              if (cells[cid].nodeID == e->nodeID &&
                  cells[cjd].nodeID != e->nodeID) {
                pid = e->proxy_ind[cells[cjd].nodeID];
                if (pid < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cjd].nodeID);
                  e->proxy_ind[cells[cjd].nodeID] = e->nr_proxies;
                  pid = e->nr_proxies;
                  e->nr_proxies += 1;
                }
                proxy_addcell_in(&proxies[pid], &cells[cjd]);
                proxy_addcell_out(&proxies[pid], &cells[cid]);
                cells[cid].sendto |= (1ULL << pid);
              }

              if (cells[cjd].nodeID == e->nodeID &&
                  cells[cid].nodeID != e->nodeID) {
                pid = e->proxy_ind[cells[cid].nodeID];
                if (pid < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cid].nodeID);
                  e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
                  pid = e->nr_proxies;
                  e->nr_proxies += 1;
                }
                proxy_addcell_in(&proxies[pid], &cells[cid]);
                proxy_addcell_out(&proxies[pid], &cells[cjd]);
                cells[cjd].sendto |= (1ULL << pid);
              }
            }
          }
        }
      }
}

/**
 * @brief Split the underlying space according to the given grid.
 *
 * @param e The #engine.
 * @param grid The grid.
 */

void engine_split(struct engine *e, int *grid) {

  int j, k;
  int ind[3];
  struct space *s = e->s;
  struct cell *c;

  /* If we've got the wrong number of nodes, fail. */
  if (e->nr_nodes != grid[0] * grid[1] * grid[2])
    error("Grid size does not match number of nodes.");

  /* Run through the cells and set their nodeID. */
  // message("s->dim = [%e,%e,%e]", s->dim[0], s->dim[1], s->dim[2]);
  for (k = 0; k < s->nr_cells; k++) {
    c = &s->cells[k];
    for (j = 0; j < 3; j++) ind[j] = c->loc[j] / s->dim[j] * grid[j];
    c->nodeID = ind[0] + grid[0] * (ind[1] + grid[1] * ind[2]);
    // message("cell at [%e,%e,%e]: ind = [%i,%i,%i], nodeID = %i", c->loc[0],
    // c->loc[1], c->loc[2], ind[0], ind[1], ind[2], c->nodeID);
  }

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Re-allocate the local parts. */
  if (e->nodeID == 0)
    message("Re-allocating parts array from %i to %i.", s->size_parts,
            (int)(s->nr_parts * 1.2));
  s->size_parts = s->nr_parts * 1.2;
  struct part *parts_new;
  struct xpart *xparts_new;
  if (posix_memalign((void **)&parts_new, part_align,
                     sizeof(struct part) * s->size_parts) != 0 ||
      posix_memalign((void **)&xparts_new, part_align,
                     sizeof(struct xpart) * s->size_parts) != 0)
    error("Failed to allocate new part data.");
  memcpy(parts_new, s->parts, sizeof(struct part) * s->nr_parts);
  memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->nr_parts);
  free(s->parts);
  free(s->xparts);
  s->parts = parts_new;
  s->xparts = xparts_new;
}

/**
 * @brief init an engine with the given number of threads, queues, and
 *      the given policy.
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param dt The initial time step to use.
 * @param nr_threads The number of threads to spawn.
 * @param nr_queues The number of task queues to create.
 * @param nr_nodes The number of MPI ranks.
 * @param nodeID The MPI rank of this node.
 * @param policy The queuing policy to use.
 * @param timeBegin Time at the begininning of the simulation.
 * @param timeEnd Time at the end of the simulation.
 */

void engine_init(struct engine *e, struct space *s, float dt, int nr_threads,
                 int nr_queues, int nr_nodes, int nodeID, int policy,
                 float timeBegin, float timeEnd) {

  int k;
#if defined(HAVE_SETAFFINITY)
  int nr_cores = sysconf(_SC_NPROCESSORS_ONLN);
  int i, j, cpuid[nr_cores];
  cpu_set_t cpuset;
  if (policy & engine_policy_cputight) {
    for (k = 0; k < nr_cores; k++) cpuid[k] = k;
  } else {
    /*  Get next highest power of 2. */
    int maxint = 1;
    while (maxint < nr_cores) maxint *= 2;

    cpuid[0] = 0;
    k = 1;
    for (i = 1; i < maxint; i *= 2)
      for (j = maxint / i / 2; j < maxint; j += maxint / i)
        if (j < nr_cores && j != 0) cpuid[k++] = j;

    if (nodeID == 0) {
#ifdef WITH_MPI
      message("engine_init: cpu map is [ ");
#else
      printf("[%03i] engine_init: cpu map is [ ", nodeID);
#endif
      for (i = 0; i < nr_cores; i++) printf("%i ", cpuid[i]);
      printf("].\n");
    }
  }
#endif

  /* Store the values. */
  e->s = s;
  e->nr_threads = nr_threads;
  e->policy = policy;
  e->step = 0;
  e->nullstep = 0;
  e->time = 0.0;
  e->nr_nodes = nr_nodes;
  e->nodeID = nodeID;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->forcerebuild = 1;
  e->forcerepart = 0;
  e->links = NULL;
  e->nr_links = 0;
  e->timeBegin = timeBegin;
  e->timeEnd = timeEnd;
  engine_rank = nodeID;

  /* Make the space link back to the engine. */
  s->e = e;

  /* Are we doing stuff in parallel? */
  if (nr_nodes > 1) {
#ifndef HAVE_MPI
    error("SWIFT was not compiled with MPI support.");
#else
    e->policy |= engine_policy_mpi;
    if ((e->proxies = (struct proxy *)malloc(sizeof(struct proxy) *
                                             engine_maxproxies)) == NULL)
      error("Failed to allocate memory for proxies.");
    bzero(e->proxies, sizeof(struct proxy) * engine_maxproxies);
    e->nr_proxies = 0;
#endif
  }

  /* First of all, init the barrier and lock it. */
  if (pthread_mutex_init(&e->barrier_mutex, NULL) != 0)
    error("Failed to initialize barrier mutex.");
  if (pthread_cond_init(&e->barrier_cond, NULL) != 0)
    error("Failed to initialize barrier condition variable.");
  if (pthread_mutex_lock(&e->barrier_mutex) != 0)
    error("Failed to lock barrier mutex.");
  e->barrier_running = 0;
  e->barrier_launch = 0;
  e->barrier_launchcount = 0;

  /* Init the scheduler. */
  scheduler_init(&e->sched, e->s, nr_queues, scheduler_flag_steal, e->nodeID);
  s->nr_queues = nr_queues;

  /* Allocate and init the threads. */
  if ((e->runners =
           (struct runner *)malloc(sizeof(struct runner) * nr_threads)) == NULL)
    error("Failed to allocate threads array.");
  for (k = 0; k < nr_threads; k++) {
    e->runners[k].id = k;
    e->runners[k].e = e;
    e->barrier_running += 1;
    if (pthread_create(&e->runners[k].thread, NULL, &runner_main,
                       &e->runners[k]) != 0)
      error("Failed to create runner thread.");
    if (e->policy & engine_policy_setaffinity) {
#if defined(HAVE_SETAFFINITY)

      /* Set a reasonable queue ID. */
      e->runners[k].cpuid = cpuid[k % nr_cores];
      if (nr_queues < nr_threads)
        e->runners[k].qid = cpuid[k % nr_cores] * nr_queues / nr_cores;
      else
        e->runners[k].qid = k;

      /* Set the cpu mask to zero | e->id. */
      CPU_ZERO(&cpuset);
      CPU_SET(cpuid[k % nr_cores], &cpuset);

      /* Apply this mask to the runner's pthread. */
      if (pthread_setaffinity_np(e->runners[k].thread, sizeof(cpu_set_t),
                                 &cpuset) != 0)
        error("Failed to set thread affinity.");

#else
      error("SWIFT was not compiled with affinity enabled.");
#endif
    } else {
      e->runners[k].cpuid = k;
      e->runners[k].qid = k * nr_queues / nr_threads;
    }
    // message( "runner %i on cpuid=%i with qid=%i." , e->runners[k].id ,
    // e->runners[k].cpuid , e->runners[k].qid );
  }

  /* Wait for the runner threads to be in place. */
  while (e->barrier_running || e->barrier_launch)
    if (pthread_cond_wait(&e->barrier_cond, &e->barrier_mutex) != 0)
      error("Error while waiting for runner threads to get in place.");
}
