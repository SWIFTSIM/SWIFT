/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <stdbool.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef HAVE_LIBNUMA
#include <numa.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "atomic.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "hydro.h"
#include "minmax.h"
#include "part.h"
#include "partition.h"
#include "timers.h"

const char *engine_policy_names[12] = {
    "none",          "rand",   "steal",        "keep",
    "block",         "fix_dt", "cpu_tight",    "mpi",
    "numa_affinity", "hydro",  "self_gravity", "external_gravity"};

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
 * @brief Generate the ghosts all the O(Npart) tasks for a hierarchy of cells.
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * @param e The #engine.
 * @param c The #cell.
 * @param super The super #cell.
 */

void engine_make_ghost_tasks(struct engine *e, struct cell *c,
                             struct cell *super) {

  struct scheduler *s = &e->sched;

  /* Am I the super-cell? */
  if (super == NULL && (c->count > 0 || c->gcount > 0)) {

    /* Remember me. */
    super = c;

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      if(c->count > 0) {
        /* Generate the ghost task. */
        c->ghost = scheduler_addtask(s, task_type_ghost, task_subtype_none, 0, 0,
                                     c, NULL, 0);
        /* Add the init task. */
        c->init = scheduler_addtask(s, task_type_init, task_subtype_none, 0, 0, c,
                                    NULL, 0);
      }

      /* Add the drift task. */
      c->drift = scheduler_addtask(s, task_type_drift, task_subtype_none, 0, 0,
                                   c, NULL, 0);

      /* Add the kick task. */
      c->kick = scheduler_addtask(s, task_type_kick, task_subtype_none, 0, 0, c,
                                  NULL, 0);

      if(c->gcount > 0) {
        /* Add the gravity tasks */
        c->grav_external = scheduler_addtask(s, task_type_grav_external, task_subtype_none, 0, 0,
                                             c, NULL, 0);

        /* Enforce gravity calculated before kick  */
        scheduler_addunlock(s, c->grav_external, c->kick);
      }
    }
  }

  /* Set the super-cell. */
  c->super = super;

  /* Recurse. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_make_ghost_tasks(e, c->progeny[k], super);
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
  ticks tic = getticks();

  /* Start by sorting the particles according to their nodes and
     getting the counts. The counts array is indexed as
     count[from * nr_nodes + to]. */
  int *counts;
  size_t *dest;
  double ih[3], dim[3];
  ih[0] = s->ih[0];
  ih[1] = s->ih[1];
  ih[2] = s->ih[2];
  dim[0] = s->dim[0];
  dim[1] = s->dim[1];
  dim[2] = s->dim[2];
  if ((counts = (int *)malloc(sizeof(int) *nr_nodes *nr_nodes)) == NULL ||
      (dest = (size_t *)malloc(sizeof(size_t) * s->nr_parts)) == NULL)
    error("Failed to allocate count and dest buffers.");
  bzero(counts, sizeof(int) * nr_nodes * nr_nodes);
  struct part *parts = s->parts;
  for (size_t k = 0; k < s->nr_parts; k++) {
    for (int j = 0; j < 3; j++) {
      if (parts[k].x[j] < 0.0)
        parts[k].x[j] += dim[j];
      else if (parts[k].x[j] >= dim[j])
        parts[k].x[j] -= dim[j];
    }
    const int cid = cell_getid(cdim, parts[k].x[0] * ih[0],
                               parts[k].x[1] * ih[1], parts[k].x[2] * ih[2]);
    /* if (cid < 0 || cid >= s->nr_cells)
       error("Bad cell id %i for part %i at [%.3e,%.3e,%.3e].",
             cid, k, parts[k].x[0], parts[k].x[1], parts[k].x[2]); */
    dest[k] = cells[cid].nodeID;
    counts[nodeID * nr_nodes + dest[k]] += 1;
  }
  space_parts_sort(s, dest, s->nr_parts, 0, nr_nodes - 1, e->verbose);

  /* Get all the counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, counts, nr_nodes * nr_nodes, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce particle transfer counts.");

  /* Get the new number of parts for this node, be generous in allocating. */
  size_t nr_parts = 0;
  for (int k = 0; k < nr_nodes; k++) nr_parts += counts[k * nr_nodes + nodeID];
  struct part *parts_new = NULL;
  struct xpart *xparts_new = NULL, *xparts = s->xparts;
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
  for (size_t offset_send = 0, offset_recv = 0, k = 0; k < nr_nodes; k++) {
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
        if (MPI_Isend(&s->parts[offset_send], counts[ind_send],
                      e->part_mpi_type, k, 2 * ind_send + 0, MPI_COMM_WORLD,
                      &reqs[4 * k]) != MPI_SUCCESS)
          error("Failed to isend parts to node %zi.", k);
        if (MPI_Isend(&s->xparts[offset_send], counts[ind_send],
                      e->xpart_mpi_type, k, 2 * ind_send + 1, MPI_COMM_WORLD,
                      &reqs[4 * k + 1]) != MPI_SUCCESS)
          error("Failed to isend xparts to node %zi.", k);
        offset_send += counts[ind_send];
      }
    }
    if (k != nodeID && counts[ind_recv] > 0) {
      if (MPI_Irecv(&parts_new[offset_recv], counts[ind_recv], e->part_mpi_type,
                    k, 2 * ind_recv + 0, MPI_COMM_WORLD,
                    &reqs[4 * k + 2]) != MPI_SUCCESS)
        error("Failed to emit irecv of parts from node %zi.", k);
      if (MPI_Irecv(&xparts_new[offset_recv], counts[ind_recv],
                    e->xpart_mpi_type, k, 2 * ind_recv + 1, MPI_COMM_WORLD,
                    &reqs[4 * k + 3]) != MPI_SUCCESS)
        error("Failed to emit irecv of parts from node %zi.", k);
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
  /* for ( int k = 0 ; k < nr_parts ; k++ ) {
      int cid = cell_getid( cdim , parts_new[k].x[0]*ih[0], parts_new[k].x[1]*ih[1], parts_new[k].x[2]*ih[2] );
      if ( cells[ cid ].nodeID != nodeID )
          error( "Received particle (%i) that does not belong here (nodeID=%i).", k , cells[ cid ].nodeID );
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
  if (e->verbose)
    message("node %i now has %zi parts in %i cells.", nodeID, nr_parts,
            my_cells);

  /* Clean up other stuff. */
  free(reqs);
  free(counts);
  free(dest);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */

void engine_repartition(struct engine *e) {

#if defined(WITH_MPI) && defined(HAVE_METIS)

  ticks tic = getticks();

  /* Clear the repartition flag. */
  enum repartition_type reparttype = e->forcerepart;
  e->forcerepart = REPART_NONE;

  /* Nothing to do if only using a single node. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (e->nr_nodes == 1) return;

  /* Do the repartitioning. */
  partition_repartition(reparttype, e->nodeID, e->nr_nodes, e->s,
                        e->sched.tasks, e->sched.nr_tasks);

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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI and METIS support.");
#endif
}

/**
 * @brief Add up/down gravity tasks to a cell hierarchy.
 *
 * @param e The #engine.
 * @param c The #cell
 * @param up The upward gravity #task.
 * @param down The downward gravity #task.
 */

void engine_addtasks_grav(struct engine *e, struct cell *c, struct task *up,
                          struct task *down) {

  /* Link the tasks to this cell. */
  c->grav_up = up;
  c->grav_down = down;

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_grav(e, c->progeny[k], up, down);
}

/**
 * @brief Add send tasks to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj The receiving #cell
 */

void engine_addtasks_send(struct engine *e, struct cell *ci, struct cell *cj) {

#ifdef WITH_MPI
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
    struct task *t_xv = scheduler_addtask(s, task_type_send, task_subtype_none,
                                          2 * ci->tag, 0, ci, cj, 0);
    struct task *t_rho = scheduler_addtask(s, task_type_send, task_subtype_none,
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
    for (int k = 0; k < 8; k++)
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
  struct scheduler *s = &e->sched;

  /* Do we need to construct a recv task? */
  if (t_xv == NULL && c->nr_density > 0) {

    /* Create the tasks. */
    t_xv = c->recv_xv = scheduler_addtask(s, task_type_recv, task_subtype_none,
                                          2 * c->tag, 0, c, NULL, 0);
    t_rho = c->recv_rho = scheduler_addtask(
        s, task_type_recv, task_subtype_none, 2 * c->tag + 1, 0, c, NULL, 0);
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
    for (int k = 0; k < 8; k++)
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

  struct space *s = e->s;
  struct cell *cells = s->cells;
  const int nr_cells = s->nr_cells;
  const int nr_proxies = e->nr_proxies;
  int offset[nr_cells];
  MPI_Request reqs_in[engine_maxproxies];
  MPI_Request reqs_out[engine_maxproxies];
  MPI_Status status;
  ticks tic = getticks();

  /* Run through the cells and get the size of the ones that will be sent off.
   */
  int count_out = 0;
  for (int k = 0; k < nr_cells; k++) {
    offset[k] = count_out;
    if (cells[k].sendto)
      count_out += (cells[k].pcell_size = cell_getsize(&cells[k]));
  }

  /* Allocate the pcells. */
  struct pcell *pcells;
  if ((pcells = (struct pcell *)malloc(sizeof(struct pcell) * count_out)) ==
      NULL)
    error("Failed to allocate pcell buffer.");

  /* Pack the cells. */
  cell_next_tag = 0;
  for (int k = 0; k < nr_cells; k++)
    if (cells[k].sendto) {
      cell_pack(&cells[k], &pcells[offset[k]]);
      cells[k].pcell = &pcells[offset[k]];
    }

  /* Launch the proxies. */
  for (int k = 0; k < nr_proxies; k++) {
    proxy_cells_exch1(&e->proxies[k]);
    reqs_in[k] = e->proxies[k].req_cells_count_in;
    reqs_out[k] = e->proxies[k].req_cells_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (int k = 0; k < nr_proxies; k++) {
    int pid = MPI_UNDEFINED;
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
  for (int k = 0; k < nr_proxies; k++) {
    reqs_in[k] = e->proxies[k].req_cells_in;
    reqs_out[k] = e->proxies[k].req_cells_out;
  }

  /* Wait for each pcell array to come in from the proxies. */
  for (int k = 0; k < nr_proxies; k++) {
    int pid = MPI_UNDEFINED;
    if (MPI_Waitany(nr_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "cell data from proxy %i has arrived." , pid );
    for (int count = 0, j = 0; j < e->proxies[pid].nr_cells_in; j++)
      count += cell_unpack(&e->proxies[pid].pcells_in[count],
                           e->proxies[pid].cells_in[j], e->s);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(nr_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Count the number of particles we need to import and re-allocate
     the buffer if needed. */
  int count_in = 0;
  for (int k = 0; k < nr_proxies; k++)
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++)
      count_in += e->proxies[k].cells_in[j]->count;
  if (count_in > s->size_parts_foreign) {
    if (s->parts_foreign != NULL) free(s->parts_foreign);
    s->size_parts_foreign = 1.1 * count_in;
    if (posix_memalign((void **)&s->parts_foreign, part_align,
                       sizeof(struct part) * s->size_parts_foreign) != 0)
      error("Failed to allocate foreign part data.");
  }

  /* Unpack the cells and link to the particle data. */
  struct part *parts = s->parts_foreign;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {
      cell_link(e->proxies[k].cells_in[j], parts);
      parts = &parts[e->proxies[k].cells_in[j]->count];
    }
  }
  s->nr_parts_foreign = parts - s->parts_foreign;

  /* Is the parts buffer large enough? */
  if (s->nr_parts_foreign > s->size_parts_foreign)
    error("Foreign parts buffer too small.");

  /* Free the pcell buffer. */
  free(pcells);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

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

int engine_exchange_strays(struct engine *e, int offset, size_t *ind,
                           size_t N) {

#ifdef WITH_MPI

  struct space *s = e->s;
  ticks tic = getticks();

  /* Re-set the proxies. */
  for (int k = 0; k < e->nr_proxies; k++) e->proxies[k].nr_parts_out = 0;

  /* Put the parts into the corresponding proxies. */
  for (size_t k = 0; k < N; k++) {
    const int node_id = e->s->cells[ind[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
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
  MPI_Request reqs_in[2 * engine_maxproxies];
  MPI_Request reqs_out[2 * engine_maxproxies];
  for (int k = 0; k < e->nr_proxies; k++) {
    proxy_parts_exch1(&e->proxies[k]);
    reqs_in[k] = e->proxies[k].req_parts_count_in;
    reqs_out[k] = e->proxies[k].req_parts_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (int k = 0; k < e->nr_proxies; k++) {
    int pid = MPI_UNDEFINED;
    if (MPI_Waitany(e->nr_proxies, reqs_in, &pid, MPI_STATUS_IGNORE) !=
            MPI_SUCCESS ||
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
  size_t count_in = 0;
  for (int k = 0; k < e->nr_proxies; k++) count_in += e->proxies[k].nr_parts_in;
  if (e->verbose) message("sent out %zi particles, got %zi back.", N, count_in);
  if (offset + count_in > s->size_parts) {
    s->size_parts = (offset + count_in) * 1.05;
    struct part *parts_new = NULL;
    struct xpart *xparts_new = NULL;
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
  int nr_in = 0, nr_out = 0;
  for (int k = 0; k < e->nr_proxies; k++) {
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
  size_t count = 0;
  for (int k = 0; k < 2 * (nr_in + nr_out); k++) {
    int err = 0, pid = MPI_UNDEFINED;
    if ((err = MPI_Waitany(2 * e->nr_proxies, reqs_in, &pid,
                           MPI_STATUS_IGNORE)) != MPI_SUCCESS) {
      char buff[MPI_MAX_ERROR_STRING];
      int res;
      MPI_Error_string(err, buff, &res);
      error("MPI_Waitany failed (%s).", buff);
    }
    if (pid == MPI_UNDEFINED) break;
    // message( "request from proxy %i has arrived." , pid );
    if (reqs_in[pid & ~1] == MPI_REQUEST_NULL &&
        reqs_in[pid | 1] == MPI_REQUEST_NULL) {
      struct proxy *p = &e->proxies[pid >> 1];
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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Return the number of harvested parts. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Constructs the top-level pair tasks for the first hydro loop over
 *neighbours
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * @param e The #engine.
 */
void engine_make_hydroloop_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->cdim;
  struct cell *cells = s->cells;

  /* Run through the highest level of cells and add pairs. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell */
        const int cid = cell_getid(cdim, i, j, k);
        struct cell *ci = &cells[cid];

        /* Skip cells without hydro particles */
        if (ci->count == 0) continue;

        /* If the cells is local build a self-interaction */
        if (ci->nodeID == nodeID)
          scheduler_addtask(sched, task_type_self, task_subtype_density, 0, 0,
                            ci, NULL, 0);

        /* Now loop over all the neighbours of this cell */
        for (int ii = -1; ii < 2; ii++) {
          int iii = i + ii;
          if (!s->periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -1; jj < 2; jj++) {
            int jjj = j + jj;
            if (!s->periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -1; kk < 2; kk++) {
              int kkk = k + kk;
              if (!s->periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];

              /* Get the neighbouring cell */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);
              struct cell *cj = &cells[cjd];

              /* Is that neighbour local and does it have particles ? */
              if (cid >= cjd || cj->count == 0 ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              /* Construct the pair task */
              const int sid =
                  sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
              scheduler_addtask(sched, task_type_pair, task_subtype_density,
                                sid, 0, ci, cj, 1);
            }
          }
        }
      }
    }
  }
}

/**
 * @brief Counts the tasks associated with one cell and constructs the links
 *
 * For each hydrodynamic task, construct the links with the corresponding cell.
 * Similarly, construct the dependencies for all the sorting tasks.
 *
 * @param e The #engine.
 */
void engine_count_and_link_tasks(struct engine *e) {

  struct scheduler *sched = &e->sched;
  const int nr_tasks = sched->nr_tasks;

  for (int k = 0; k < nr_tasks; k++) {

    /* Get the current task. */
    struct task *t = &sched->tasks[k];
    if (t->skip) continue;

    /* Link sort tasks together. */
    if (t->type == task_type_sort && t->ci->split)
      for (int j = 0; j < 8; j++)
        if (t->ci->progeny[j] != NULL && t->ci->progeny[j]->sorts != NULL) {
          t->ci->progeny[j]->sorts->skip = 0;
          scheduler_addunlock(sched, t->ci->progeny[j]->sorts, t);
        }

    /* Link density tasks to cells. */
    if (t->type == task_type_self) {
      atomic_inc(&t->ci->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->density = engine_addlink(e, t->ci->density, t);
        atomic_inc(&t->ci->nr_density);
      }
    } else if (t->type == task_type_pair) {
      atomic_inc(&t->ci->nr_tasks);
      atomic_inc(&t->cj->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->density = engine_addlink(e, t->ci->density, t);
        atomic_inc(&t->ci->nr_density);
        t->cj->density = engine_addlink(e, t->cj->density, t);
        atomic_inc(&t->cj->nr_density);
      }
    } else if (t->type == task_type_sub) {
      atomic_inc(&t->ci->nr_tasks);
      if (t->cj != NULL) atomic_inc(&t->cj->nr_tasks);
      if (t->subtype == task_subtype_density) {
        t->ci->density = engine_addlink(e, t->ci->density, t);
        atomic_inc(&t->ci->nr_density);
        if (t->cj != NULL) {
          t->cj->density = engine_addlink(e, t->cj->density, t);
          atomic_inc(&t->cj->nr_density);
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
}

/**
 * @brief Duplicates the first hydro loop and construct all the
 * dependencies for the hydro part
 *
 * This is done by looping over all the previously constructed tasks
 * and adding another task involving the same cells but this time
 * corresponding to the second hydro loop over neighbours.
 * With all the relevant tasks for a given cell available, we construct
 * all the dependencies for that cell.
 *
 * @param e The #engine.
 */
void engine_make_extra_hydroloop_tasks(struct engine *e) {

  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int nr_tasks = sched->nr_tasks;

  for (int k = 0; k < nr_tasks; k++) {

    /* Get a pointer to the task. */
    struct task *t = &sched->tasks[k];

    /* Skip? */
    if (t->skip) continue;

    /* Self-interaction? */
    if (t->type == task_type_self && t->subtype == task_subtype_density) {

      /* Start by constructing the task for the second hydro loop */
      struct task *t2 = scheduler_addtask(
          sched, task_type_self, task_subtype_force, 0, 0, t->ci, NULL, 0);

      /* Add the link between the new loop and the cell */
      t->ci->force = engine_addlink(e, t->ci->force, t2);
      atomic_inc(&t->ci->nr_force);

      /* Now, build all the dependencies for the hydro */
      /* init --> t (density loop) --> ghost --> t2 (force loop) --> kick */
      scheduler_addunlock(sched, t->ci->super->init, t);
      scheduler_addunlock(sched, t, t->ci->super->ghost);
      scheduler_addunlock(sched, t->ci->super->ghost, t2);
      scheduler_addunlock(sched, t2, t->ci->super->kick);
    }

    /* Otherwise, pair interaction? */
    else if (t->type == task_type_pair && t->subtype == task_subtype_density) {

      /* Start by constructing the task for the second hydro loop */
      struct task *t2 = scheduler_addtask(
          sched, task_type_pair, task_subtype_force, 0, 0, t->ci, t->cj, 0);

      /* Add the link between the new loop and both cells */
      t->ci->force = engine_addlink(e, t->ci->force, t2);
      atomic_inc(&t->ci->nr_force);
      t->cj->force = engine_addlink(e, t->cj->force, t2);
      atomic_inc(&t->cj->nr_force);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super-cells */
      /* init --> t (density loop) --> ghost --> t2 (force loop) --> kick */
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
    }

    /* Otherwise, sub interaction? */
    else if (t->type == task_type_sub && t->subtype == task_subtype_density) {

      /* Start by constructing the task for the second hydro loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_sub, task_subtype_force, t->flags,
                            0, t->ci, t->cj, 0);

      /* Add the link between the new loop and both cells */
      t->ci->force = engine_addlink(e, t->ci->force, t2);
      atomic_inc(&t->ci->nr_force);
      if (t->cj != NULL) {
        t->cj->force = engine_addlink(e, t->cj->force, t2);
        atomic_inc(&t->cj->nr_force);
      }

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super-cells */
      /* init --> t (density loop) --> ghost --> t2 (force loop) --> kick */
      if (t->ci->nodeID == nodeID) {
        scheduler_addunlock(sched, t, t->ci->super->ghost);
        scheduler_addunlock(sched, t->ci->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->ci->super->kick);
      }
      if (t->cj != NULL && t->cj->nodeID == nodeID &&
          t->ci->super != t->cj->super) {
        scheduler_addunlock(sched, t, t->cj->super->ghost);
        scheduler_addunlock(sched, t->cj->super->ghost, t2);
        scheduler_addunlock(sched, t2, t->cj->super->kick);
      }
    }

    /* /\* Kick tasks should rely on the grav_down tasks of their cell. *\/ */
    /* else if (t->type == task_type_kick && t->ci->grav_down != NULL) */
    /*   scheduler_addunlock(sched, t->ci->grav_down, t); */
  }
}

/**
 * @brief Constructs the top-level pair tasks for the gravity M-M interactions
 *
 * Correct implementation is still lacking here.
 *
 * @param e The #engine.
 */
void engine_make_gravityinteraction_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nr_cells = s->nr_cells;
  struct cell *cells = s->cells;

  /* Loop over all cells. */
  for (int i = 0; i < nr_cells; i++) {

    /* If it has gravity particles, add a self-task */
    if (cells[i].gcount > 0) {
      scheduler_addtask(sched, task_type_grav_mm, task_subtype_none, -1, 0,
                        &cells[i], NULL, 0);

      /* Loop over all remainding cells */
      for (int j = i + 1; j < nr_cells; j++) {

        /* If that other cell has gravity parts, add a pair interaction */
        if (cells[j].gcount > 0) {
          scheduler_addtask(sched, task_type_grav_mm, task_subtype_none, -1, 0,
                            &cells[i], &cells[j], 0);
        }
      }
    }
  }
}

/**
 * @brief Constructs the gravity tasks building the multipoles and propagating
 *them to the children
 *
 * Correct implementation is still lacking here.
 *
 * @param e The #engine.
 */
void engine_make_gravityrecursive_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int nr_cells = s->nr_cells;
  struct cell *cells = s->cells;

  for (int k = 0; k < nr_cells; k++) {

    /* Only do this for local cells containing gravity particles */
    if (cells[k].nodeID == nodeID && cells[k].gcount > 0) {

      /* Create tasks at top level. */
      struct task *up =
          scheduler_addtask(sched, task_type_grav_up, task_subtype_none, 0, 0,
                            &cells[k], NULL, 0);
      struct task *down =
          scheduler_addtask(sched, task_type_grav_down, task_subtype_none, 0, 0,
                            &cells[k], NULL, 0);

      /* Push tasks down the cell hierarchy. */
      engine_addtasks_grav(e, &cells[k], up, down);
    }
  }
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
  const int nr_cells = s->nr_cells;
  const ticks tic = getticks();

  /* Re-set the scheduler. */
  scheduler_reset(sched, s->tot_cells * engine_maxtaskspercell);

  /* Add the space sorting tasks. */
  for (int i = 0; i < e->nr_threads; i++)
    scheduler_addtask(sched, task_type_psort, task_subtype_none, i, 0, NULL,
                      NULL, 0);

  /* Construct the firt hydro loop over neighbours */
  engine_make_hydroloop_tasks(e);

  /* Add the gravity mm tasks. */
  if ((e->policy & engine_policy_self_gravity) == engine_policy_self_gravity)
    engine_make_gravityinteraction_tasks(e);

  /* Split the tasks. */
  scheduler_splittasks(sched);

  /* Allocate the list of cell-task links. The maximum number of links
     is the number of cells (s->tot_cells) times the number of neighbours (27)
     times the number of interaction types (2, density and force). */
  if (e->links != NULL) free(e->links);
  e->size_links = s->tot_cells * 27 * 2;
  if ((e->links = malloc(sizeof(struct link) * e->size_links)) == NULL)
    error("Failed to allocate cell-task links.");
  e->nr_links = 0;

  /* Add the gravity up/down tasks at the top-level cells and push them down. */
  if ((e->policy & engine_policy_self_gravity) == engine_policy_self_gravity)
    engine_make_gravityrecursive_tasks(e);

  /* Count the number of tasks associated with each cell and
     store the density tasks in each cell, and make each sort
     depend on the sorts of its progeny. */
  engine_count_and_link_tasks(e);

  /* Append a ghost task to each cell, and add kick tasks to the
     super cells. */
  for (int k = 0; k < nr_cells; k++)
    engine_make_ghost_tasks(e, &cells[k], NULL);

  /* Run through the tasks and make force tasks for each density task.
     Each force task depends on the cell ghosts and unlocks the kick task
     of its super-cell. */
  engine_make_extra_hydroloop_tasks(e);

  /* Add the communication tasks if MPI is being used. */
  if ((e->policy & engine_policy_mpi) == engine_policy_mpi) {

    /* Loop over the proxies. */
    for (int pid = 0; pid < e->nr_proxies; pid++) {

      /* Get a handle on the proxy. */
      struct proxy *p = &e->proxies[pid];

      /* Loop through the proxy's incoming cells and add the
         recv tasks. */
      for (int k = 0; k < p->nr_cells_in; k++)
        engine_addtasks_recv(e, p->cells_in[k], NULL, NULL);

      /* Loop through the proxy's outgoing cells and add the
         send tasks. */
      for (int k = 0; k < p->nr_cells_out; k++)
        engine_addtasks_send(e, p->cells_out[k], p->cells_in[0]);
    }
  }

  /* Set the unlocks per task. */
  scheduler_set_unlocks(sched);

  /* Rank the tasks. */
  scheduler_ranktasks(sched);

  /* Weight the tasks. */
  scheduler_reweight(sched);

  /* Set the tasks age. */
  e->tasks_age = 0;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Mark tasks to be skipped and set the sort flags accordingly.
 *
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */

int engine_marktasks(struct engine *e) {

  struct scheduler *s = &e->sched;
  const int ti_end = e->ti_current;
  const int nr_tasks = s->nr_tasks;
  const int *const ind = s->tasks_ind;
  struct task *tasks = s->tasks;
  const ticks tic = getticks();

  /* Much less to do here if we're on a fixed time-step. */
  if ((e->policy & engine_policy_fixdt) == engine_policy_fixdt) {

    /* Run through the tasks and mark as skip or not. */
    for (int k = 0; k < nr_tasks; k++) {

      /* Get a handle on the kth task. */
      struct task *t = &tasks[ind[k]];

      /* Pair? */
      if (t->type == task_type_pair ||
          (t->type == task_type_sub && t->cj != NULL)) {

        /* Local pointers. */
        const struct cell *ci = t->ci;
        const struct cell *cj = t->cj;

        /* Too much particle movement? */
        if (t->tight &&
            (fmaxf(ci->h_max, cj->h_max) + ci->dx_max + cj->dx_max > cj->dmin ||
             ci->dx_max > space_maxreldx * ci->h_max ||
             cj->dx_max > space_maxreldx * cj->h_max))
          return 1;

      }

      /* Sort? */
      else if (t->type == task_type_sort) {

        /* If all the sorts have been done, make this task implicit. */
        if (!(t->flags & (t->flags ^ t->ci->sorted))) t->implicit = 1;
      }
    }

    /* Multiple-timestep case */
  } else {

    /* Run through the tasks and mark as skip or not. */
    for (int k = 0; k < nr_tasks; k++) {

      /* Get a handle on the kth task. */
      struct task *t = &tasks[ind[k]];

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
        t->skip = (t->ci->ti_end_min > ti_end);
      }

      /* Pair? */
      else if (t->type == task_type_pair ||
               (t->type == task_type_sub && t->cj != NULL)) {

        /* Local pointers. */
        const struct cell *ci = t->ci;
        const struct cell *cj = t->cj;

        /* Set this task's skip. */
        t->skip = (ci->ti_end_min > ti_end && cj->ti_end_min > ti_end);

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
      else if (t->type == task_type_kick) {
        t->skip = (t->ci->ti_end_min > ti_end);
        t->ci->updated = 0;
      }

      /* Drift? */
      else if (t->type == task_type_drift)
        t->skip = 0;

      /* Init? */
      else if (t->type == task_type_init) {
        /* Set this task's skip. */
        t->skip = (t->ci->ti_end_min > ti_end);
      }

      /* None? */
      else if (t->type == task_type_none)
        t->skip = 1;
    }
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* All is well... */
  return 0;
}

/**
 * @brief Prints the number of tasks in the engine
 *
 * @param e The #engine.
 */

void engine_print_task_counts(struct engine *e) {

  struct scheduler *sched = &e->sched;

  /* Count and print the number of each task type. */
  int counts[task_type_count + 1];
  for (int k = 0; k <= task_type_count; k++) counts[k] = 0;
  for (int k = 0; k < sched->nr_tasks; k++)
    if (!sched->tasks[k].skip)
      counts[(int)sched->tasks[k].type] += 1;
    else
      counts[task_type_count] += 1;
#ifdef WITH_MPI
  printf("[%04i] %s engine_print_task_counts: task counts are [ %s=%i",
         e->nodeID, clocks_get_timesincestart(), taskID_names[0], counts[0]);
#else
  printf("%s engine_print_task_counts: task counts are [ %s=%i",
         clocks_get_timesincestart(), taskID_names[0], counts[0]);
#endif
  for (int k = 1; k < task_type_count; k++)
    printf(" %s=%i", taskID_names[k], counts[k]);
  printf(" skipped=%i ]\n", counts[task_type_count]);
  fflush(stdout);
  message("nr_parts = %zi.", e->s->nr_parts);
}

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 */

void engine_rebuild(struct engine *e) {

  const ticks tic = getticks();

  /* Clear the forcerebuild flag, whatever it was. */
  e->forcerebuild = 0;

  /* Re-build the space. */
  space_rebuild(e->s, 0.0, e->verbose);

/* If in parallel, exchange the cell structure. */
#ifdef WITH_MPI
  engine_exchange_cells(e);
#endif

  /* Re-build the tasks. */
  engine_maketasks(e);

  /* Run through the tasks and mark as skip or not. */
  if (engine_marktasks(e))
    error("engine_marktasks failed after space_rebuild.");

  /* Print the status of the system */
  engine_print_task_counts(e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 */

void engine_prepare(struct engine *e) {

  TIMER_TIC;

  /* Run through the tasks and mark as skip or not. */
  int rebuild = (e->forcerebuild || engine_marktasks(e));

/* Collect the values of rebuild from all nodes. */
#ifdef WITH_MPI
  int buff = 0;
  if (MPI_Allreduce(&rebuild, &buff, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    error("Failed to aggregate the rebuild flag across nodes.");
  rebuild = buff;
#endif
  e->tic_step = getticks();

  /* Did this not go through? */
  if (rebuild) {
    engine_rebuild(e);
  }

  /* Re-rank the tasks every now and then. */
  if (e->tasks_age % engine_tasksreweight == 1) {
    scheduler_reweight(&e->sched);
  }
  e->tasks_age += 1;

  TIMER_TOC(timer_prepare);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
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
 * @brief Mapping function to collect the data from the kick.
 */

void engine_collect_kick(struct cell *c) {

  /* Skip super-cells (Their values are already set) */
  if (c->kick != NULL) return;

  /* Counters for the different quantities. */
  int updated = 0;
  double e_kin = 0.0, e_int = 0.0, e_pot = 0.0;
  float mom[3] = {0.0f, 0.0f, 0.0f}, ang[3] = {0.0f, 0.0f, 0.0f};
  int ti_end_min = max_nr_timesteps, ti_end_max = 0;

  /* Only do something is the cell is non-empty */
  if (c->count != 0 || c->gcount != 0) {

    /* If this cell is not split, I'm in trouble. */
    if (!c->split) error("Cell has no super-cell.");

    /* Collect the values from the progeny. */
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      if (cp != NULL) {

        /* Recurse */
        engine_collect_kick(cp);

        /* And update */
        ti_end_min = min(ti_end_min, cp->ti_end_min);
        ti_end_max = max(ti_end_max, cp->ti_end_max);
        updated += cp->updated;
        e_kin += cp->e_kin;
        e_int += cp->e_int;
        e_pot += cp->e_pot;
        mom[0] += cp->mom[0];
        mom[1] += cp->mom[1];
        mom[2] += cp->mom[2];
        ang[0] += cp->ang[0];
        ang[1] += cp->ang[1];
        ang[2] += cp->ang[2];
      }
    }
  }

  /* Store the collected values in the cell. */
  c->ti_end_min = ti_end_min;
  c->ti_end_max = ti_end_max;
  c->updated = updated;
  c->e_kin = e_kin;
  c->e_int = e_int;
  c->e_pot = e_pot;
  c->mom[0] = mom[0];
  c->mom[1] = mom[1];
  c->mom[2] = mom[2];
  c->ang[0] = ang[0];
  c->ang[1] = ang[1];
  c->ang[2] = ang[2];
}

/**
 * @brief Launch the runners.
 *
 * @param e The #engine.
 * @param nr_runners The number of #runner to let loose.
 * @param mask The task mask to launch.
 * @param submask The sub-task mask to launch.
 */

void engine_launch(struct engine *e, int nr_runners, unsigned int mask,
                   unsigned int submask) {

  /* Prepare the scheduler. */
  atomic_inc(&e->sched.waiting);

  /* Cry havoc and let loose the dogs of war. */
  e->barrier_launch = nr_runners;
  e->barrier_launchcount = nr_runners;
  if (pthread_cond_broadcast(&e->barrier_cond) != 0)
    error("Failed to broadcast barrier open condition.");

  /* Load the tasks. */
  pthread_mutex_unlock(&e->barrier_mutex);
  scheduler_start(&e->sched, mask, submask);
  pthread_mutex_lock(&e->barrier_mutex);

  /* Remove the safeguard. */
  pthread_mutex_lock(&e->sched.sleep_mutex);
  atomic_dec(&e->sched.waiting);
  pthread_cond_broadcast(&e->sched.sleep_cond);
  pthread_mutex_unlock(&e->sched.sleep_mutex);

  /* Sit back and wait for the runners to come home. */
  while (e->barrier_launch || e->barrier_running)
    if (pthread_cond_wait(&e->barrier_cond, &e->barrier_mutex) != 0)
      error("Error while waiting for barrier.");
}

/**
 * @brief Initialises the particles and set them in a state ready to move
 *forward in time.
 *
 * @param e The #engine
 */
void engine_init_particles(struct engine *e) {

  struct space *s = e->s;

  if (e->nodeID == 0) message("Initialising particles");

  /* Make sure all particles are ready to go */
  /* i.e. clean-up any stupid state in the ICs */
  space_map_cells_pre(s, 1, cell_init_parts, NULL);

  engine_prepare(e);

  engine_marktasks(e);

  /* Build the masks corresponding to the policy */
  unsigned int mask = 0;
  unsigned int submask = 0;

  /* We always have sort tasks */
  mask |= 1 << task_type_sort;

  /* Add the tasks corresponding to hydro operations to the masks */
  if ((e->policy & engine_policy_hydro) == engine_policy_hydro) {

    mask |= 1 << task_type_init;
    mask |= 1 << task_type_self;
    mask |= 1 << task_type_pair;
    mask |= 1 << task_type_sub;
    mask |= 1 << task_type_ghost;

    submask |= 1 << task_subtype_density;
  }

  /* Add the tasks corresponding to self-gravity to the masks */
  if ((e->policy & engine_policy_self_gravity) == engine_policy_self_gravity) {

    /* Nothing here for now */
  }

  /* Add the tasks corresponding to external gravity to the masks */
  if ((e->policy & engine_policy_external_gravity) ==
      engine_policy_external_gravity) {
    printf("%s: JR: Excellent lets add the external gravity tasks here.....\n", __FUNCTION__);
    mask |= 1 << task_type_grav_external;
  }

  /* Add MPI tasks if need be */
  if ((e->policy & engine_policy_mpi) == engine_policy_mpi) {

    mask |= 1 << task_type_send;
    mask |= 1 << task_type_recv;
  }

  /* Now, launch the calculation */
  TIMER_TIC;
  engine_launch(e, e->nr_threads, mask, submask);
  TIMER_TOC(timer_runners);

  /* Apply some conversions (e.g. internal energy -> entropy) */
  space_map_cells_pre(s, 1, cell_convert_hydro, NULL);

  // printParticle(e->s->parts, e->s->xparts,1000, e->s->nr_parts);
  // printParticle(e->s->parts, e->s->xparts,515050, e->s->nr_parts);

  /* Ready to go */
  e->step = -1;
}

/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 */
void engine_step(struct engine *e) {

  int updates = 0;
  int ti_end_min = max_nr_timesteps, ti_end_max = 0;
  double e_pot = 0.0, e_int = 0.0, e_kin = 0.0;
  float mom[3] = {0.0, 0.0, 0.0};
  float ang[3] = {0.0, 0.0, 0.0};
  struct space *s = e->s;

  TIMER_TIC2;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  /* Collect the cell data. */
  for (int k = 0; k < s->nr_cells; k++)
    if (s->cells[k].nodeID == e->nodeID) {
      struct cell *c = &s->cells[k];

      /* Recurse */
      engine_collect_kick(c);

      /* And aggregate */
      ti_end_min = min(ti_end_min, c->ti_end_min);
      ti_end_max = max(ti_end_max, c->ti_end_max);
      e_kin += c->e_kin;
      e_int += c->e_int;
      e_pot += c->e_pot;
      updates += c->updated;
      mom[0] += c->mom[0];
      mom[1] += c->mom[1];
      mom[2] += c->mom[2];
      ang[0] += c->ang[0];
      ang[1] += c->ang[1];
      ang[2] += c->ang[2];
    }

/* Aggregate the data from the different nodes. */
#ifdef WITH_MPI
  {
    int in_i[1], out_i[1];
    in_i[0] = 0;
    out_i[0] = ti_end_min;
    if (MPI_Allreduce(out_i, in_i, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate t_end_min.");
    ti_end_min = in_i[0];
    out_i[0] = ti_end_max;
    if (MPI_Allreduce(out_i, in_i, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate t_end_max.");
    ti_end_max = in_i[0];
  }
  {
    double in_d[4], out_d[4];
    out_d[0] = updates;
    out_d[1] = e_kin;
    out_d[2] = e_int;
    out_d[3] = e_pot;
    if (MPI_Allreduce(out_d, in_d, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) !=
        MPI_SUCCESS)
      error("Failed to aggregate energies.");
    updates = in_d[0];
    e_kin = in_d[1];
    e_int = in_d[2];
    e_pot = in_d[3];
  }
#endif

  // message("\nDRIFT\n");

  /* Move forward in time */
  e->ti_old = e->ti_current;
  e->ti_current = ti_end_min;
  e->step += 1;
  e->time = e->ti_current * e->timeBase + e->timeBegin;
  e->timeOld = e->ti_old * e->timeBase + e->timeBegin;
  e->timeStep = (e->ti_current - e->ti_old) * e->timeBase;

  /* Drift everybody */
  engine_launch(e, e->nr_threads, 1 << task_type_drift, 0);

  // printParticle(e->s->parts, e->s->xparts, 1000, e->s->nr_parts);
  // printParticle(e->s->parts, e->s->xparts, 515050, e->s->nr_parts);

  // if(e->step == 2)   exit(0);

  if (e->nodeID == 0) {

    /* Print some information to the screen */
    printf("%d %e %e %d %.3f\n", e->step, e->time, e->timeStep, updates,
           e->wallclock_time);
    fflush(stdout);

    /* Write some energy statistics */
    fprintf(e->file_stats, "%d %f %f %f %f %f %f %f %f %f %f %f\n", e->step,
            e->time, e_kin, e_int, e_pot, e_kin + e_int + e_pot, mom[0], mom[1],
            mom[2], ang[0], ang[1], ang[2]);
    fflush(e->file_stats);
  }

  // message("\nACCELERATION AND KICK\n");

  /* Re-distribute the particles amongst the nodes? */
  if (e->forcerepart != REPART_NONE) engine_repartition(e);

  /* Prepare the space. */
  engine_prepare(e);

  /* Build the masks corresponding to the policy */
  unsigned int mask = 0, submask = 0;

  /* We always have sort tasks and kick tasks */
  mask |= 1 << task_type_sort;
  mask |= 1 << task_type_kick;

  /* Add the tasks corresponding to hydro operations to the masks */
  if ((e->policy & engine_policy_hydro) == engine_policy_hydro) {

    mask |= 1 << task_type_init;
    mask |= 1 << task_type_self;
    mask |= 1 << task_type_pair;
    mask |= 1 << task_type_sub;
    mask |= 1 << task_type_ghost;

    submask |= 1 << task_subtype_density;
    submask |= 1 << task_subtype_force;
  }

  /* Add the tasks corresponding to self-gravity to the masks */
  if ((e->policy & engine_policy_self_gravity) == engine_policy_self_gravity) {

    /* Nothing here for now */
  }

  /* Add the tasks corresponding to self-gravity to the masks */
  if ((e->policy & engine_policy_external_gravity) ==
      engine_policy_external_gravity) {
    mask |= 1 << task_type_grav_external;
  }

  /* Add MPI tasks if need be */
  if ((e->policy & engine_policy_mpi) == engine_policy_mpi) {

    mask |= 1 << task_type_send;
    mask |= 1 << task_type_recv;
  }

  /* Send off the runners. */
  TIMER_TIC;
  engine_launch(e, e->nr_threads, mask, submask);
  TIMER_TOC(timer_runners);

  TIMER_TOC2(timer_step);

  clocks_gettime(&time2);

  e->wallclock_time = (float)clocks_diff(&time1, &time2);
  // printParticle(e->s->parts, e->s->xparts,1000, e->s->nr_parts);
  // printParticle(e->s->parts, e->s->xparts,515050, e->s->nr_parts);
}

/**
 * @brief Returns 1 if the simulation has reached its end point, 0 otherwise
 */
int engine_is_done(struct engine *e) {
  return !(e->ti_current < max_nr_timesteps);
}

/**
 * @brief Create and fill the proxies.
 *
 * @param e The #engine.
 */

void engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  const int *cdim = e->s->cdim;
  const struct space *s = e->s;
  struct cell *cells = s->cells;
  struct proxy *proxies = e->proxies;
  ticks tic = getticks();

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* The following loop is super-clunky, but it's necessary
     to ensure that the order of the send and recv cells in
     the proxies is identical for all nodes! */

  /* Loop over each cell in the space. */
  int ind[3];
  for (ind[0] = 0; ind[0] < cdim[0]; ind[0]++)
    for (ind[1] = 0; ind[1] < cdim[1]; ind[1]++)
      for (ind[2] = 0; ind[2] < cdim[2]; ind[2]++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, ind[0], ind[1], ind[2]);

        /* Loop over all its neighbours (periodic). */
        for (int i = -1; i <= 1; i++) {
          int ii = ind[0] + i;
          if (ii >= cdim[0])
            ii -= cdim[0];
          else if (ii < 0)
            ii += cdim[0];
          for (int j = -1; j <= 1; j++) {
            int jj = ind[1] + j;
            if (jj >= cdim[1])
              jj -= cdim[1];
            else if (jj < 0)
              jj += cdim[1];
            for (int k = -1; k <= 1; k++) {
              int kk = ind[2] + k;
              if (kk >= cdim[2])
                kk -= cdim[2];
              else if (kk < 0)
                kk += cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, ii, jj, kk);

              /* Add to proxies? */
              if (cells[cid].nodeID == e->nodeID &&
                  cells[cjd].nodeID != e->nodeID) {
                int pid = e->proxy_ind[cells[cjd].nodeID];
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
                int pid = e->proxy_ind[cells[cid].nodeID];
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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Split the underlying space into regions and assign to separate nodes.
 *
 * @param e The #engine.
 * @param initial_partition structure defining the cell partition technique
 */

void engine_split(struct engine *e, struct partition *initial_partition) {

#ifdef WITH_MPI
  struct space *s = e->s;

  /* Do the initial partition of the cells. */
  partition_initial_partition(initial_partition, e->nodeID, e->nr_nodes, s);

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Re-allocate the local parts. */
  if (e->nodeID == 0)
    message("Re-allocating parts array from %zi to %zi.", s->size_parts,
            (size_t)(s->nr_parts * 1.2));
  s->size_parts = s->nr_parts * 1.2;
  struct part *parts_new = NULL;
  struct xpart *xparts_new = NULL;
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
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
static bool hyperthreads_present(void) {
#ifdef __linux__
  FILE *f =
      fopen("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list", "r");

  int c;
  while ((c = fgetc(f)) != EOF && c != ',')
    ;
  fclose(f);

  return c == ',';
#else
  return true;  // just guess
#endif
}
#endif

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
 * @param dt_min Minimal allowed timestep (unsed with fixdt policy)
 * @param dt_max Maximal allowed timestep
 * @param verbose Is this #engine talkative ?
 */

void engine_init(struct engine *e, struct space *s, float dt, int nr_threads,
                 int nr_queues, int nr_nodes, int nodeID, int policy,
                 float timeBegin, float timeEnd, float dt_min, float dt_max,
                 int verbose) {

  /* Store the values. */
  e->s = s;
  e->nr_threads = nr_threads;
  e->policy = policy;
  e->step = 0;
  e->nullstep = 0;
  e->nr_nodes = nr_nodes;
  e->nodeID = nodeID;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->forcerebuild = 1;
  e->forcerepart = REPART_NONE;
  e->links = NULL;
  e->nr_links = 0;
  e->timeBegin = timeBegin;
  e->timeEnd = timeEnd;
  e->timeOld = timeBegin;
  e->time = timeBegin;
  e->ti_old = 0;
  e->ti_current = 0;
  e->timeStep = 0.;
  e->dt_min = dt_min;
  e->dt_max = dt_max;
  e->file_stats = NULL;
  e->verbose = verbose;
  e->count_step = 0;
  e->wallclock_time = 0.f;
  engine_rank = nodeID;

  /* Make the space link back to the engine. */
  s->e = e;

#if defined(HAVE_SETAFFINITY)
  const int nr_cores = sysconf(_SC_NPROCESSORS_ONLN);
  int cpuid[nr_cores];
  cpu_set_t cpuset;
  if ((policy & engine_policy_cputight) == engine_policy_cputight) {
    for (int k = 0; k < nr_cores; k++) cpuid[k] = k;
  } else {
    /*  Get next highest power of 2. */
    int maxint = 1;
    while (maxint < nr_cores) maxint *= 2;

    cpuid[0] = 0;
    int k = 1;
    for (int i = 1; i < maxint; i *= 2)
      for (int j = maxint / i / 2; j < maxint; j += maxint / i)
        if (j < nr_cores && j != 0) cpuid[k++] = j;

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
    /* Ascending NUMA distance. Bubblesort(!) for stable equidistant CPUs. */
    if (numa_available() >= 0) {
      if (nodeID == 0) message("prefer NUMA-local CPUs");

      const int home = numa_node_of_cpu(sched_getcpu());
      const int half = nr_cores / 2;
      const bool swap_hyperthreads = hyperthreads_present();
      bool done = false;
      if (swap_hyperthreads && nodeID == 0)
        message("prefer physical cores to hyperthreads");

      while (!done) {
        done = true;
        for (int i = 1; i < nr_cores; i++) {
          const int node_a = numa_node_of_cpu(cpuid[i - 1]);
          const int node_b = numa_node_of_cpu(cpuid[i]);

          /* Avoid using local hyperthreads over unused remote physical cores.
           * Assume two hyperthreads, and that cpuid >= half partitions them.
           */
          const int thread_a = swap_hyperthreads && cpuid[i - 1] >= half;
          const int thread_b = swap_hyperthreads && cpuid[i] >= half;

          bool swap = thread_a > thread_b;
          if (thread_a == thread_b)
            swap = numa_distance(home, node_a) > numa_distance(home, node_b);

          if (swap) {
            const int t = cpuid[i - 1];
            cpuid[i - 1] = cpuid[i];
            cpuid[i] = t;
            done = false;
          }
        }
      }
    }
#endif

    if (nodeID == 0) {
#ifdef WITH_MPI
      printf("[%04i] %s engine_init: cpu map is [ ", nodeID,
             clocks_get_timesincestart());
#else
      printf("%s engine_init: cpu map is [ ", clocks_get_timesincestart());
#endif
      for (int i = 0; i < nr_cores; i++) printf("%i ", cpuid[i]);
      printf("].\n");
    }
  }
#endif

  /* Are we doing stuff in parallel? */
  if (nr_nodes > 1) {
#ifndef WITH_MPI
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

  /* Open some files */
  if (e->nodeID == 0) {
    e->file_stats = fopen("energy.txt", "w");
    fprintf(e->file_stats,
            "# Step Time E_kin E_int E_pot E_tot "
            "p_x p_y p_z ang_x ang_y ang_z\n");
  }

  /* Print policy */
  engine_print_policy(e);

  /* Print information about the hydro scheme */
  if (e->nodeID == 0) message("Hydrodynamic scheme: %s", SPH_IMPLEMENTATION);

  /* Check we have sensible time bounds */
  if (timeBegin >= timeEnd)
    error(
        "Final simulation time (t_end = %e) must be larger than the start time "
        "(t_beg = %e)",
        timeEnd, timeBegin);

  /* Check we have sensible time step bounds */
  if (e->dt_min > e->dt_max)
    error(
        "Minimal time step size must be smaller than maximal time step size ");

  /* Deal with timestep */
  e->timeBase = (timeEnd - timeBegin) / max_nr_timesteps;
  e->ti_current = 0;

  /* Fixed time-step case */
  if ((e->policy & engine_policy_fixdt) == engine_policy_fixdt) {
    e->dt_min = e->dt_max;

    /* Find timestep on the timeline */
    int dti_timeline = max_nr_timesteps;
    while (e->dt_min < dti_timeline * e->timeBase) dti_timeline /= 2;

    e->dt_min = e->dt_max = dti_timeline * e->timeBase;

    if (e->nodeID == 0) message("Timestep set to %e", e->dt_max);
  } else {

    if (e->nodeID == 0) {
      message("Absolute minimal timestep size: %e", e->timeBase);

      float dt_min = timeEnd - timeBegin;
      while (dt_min > e->dt_min) dt_min /= 2.f;

      message("Minimal timestep size (on time-line): %e", dt_min);

      float dt_max = timeEnd - timeBegin;
      while (dt_max > e->dt_max) dt_max /= 2.f;

      message("Maximal timestep size (on time-line): %e", dt_max);
    }
  }

  if (e->dt_min < e->timeBase && e->nodeID == 0)
    error(
        "Minimal time-step size smaller than the absolute possible minimum "
        "dt=%e",
        e->timeBase);

  if (e->dt_max > (e->timeEnd - e->timeBegin) && e->nodeID == 0)
    error("Maximal time-step size larger than the simulation run time t=%e",
          e->timeEnd - e->timeBegin);

/* Construct types for MPI communications */
#ifdef WITH_MPI
  part_create_mpi_type(&e->part_mpi_type);
  xpart_create_mpi_type(&e->xpart_mpi_type);
#endif

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

  /* Init the scheduler with enough tasks for the initial sorting tasks. */
  int nr_tasks = 2 * s->tot_cells + e->nr_threads;
  scheduler_init(&e->sched, e->s, nr_tasks, nr_queues, scheduler_flag_steal,
                 e->nodeID);
  s->nr_queues = nr_queues;

  /* Create the sorting tasks. */
  for (int i = 0; i < e->nr_threads; i++)
    scheduler_addtask(&e->sched, task_type_psort, task_subtype_none, i, 0, NULL,
                      NULL, 0);

  scheduler_ranktasks(&e->sched);

  /* Allocate and init the threads. */
  if ((e->runners =
           (struct runner *)malloc(sizeof(struct runner) * nr_threads)) == NULL)
    error("Failed to allocate threads array.");
  for (int k = 0; k < nr_threads; k++) {
    e->runners[k].id = k;
    e->runners[k].e = e;
    e->barrier_running += 1;
    if (pthread_create(&e->runners[k].thread, NULL, &runner_main,
                       &e->runners[k]) != 0)
      error("Failed to create runner thread.");
    if ((e->policy & engine_policy_setaffinity) == engine_policy_setaffinity) {
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

/**
 * @brief Prints the current policy of an engine
 *
 * @param e The engine to print information about
 */
void engine_print_policy(struct engine *e) {

#ifdef WITH_MPI
  if (e->nodeID == 0) {
    printf("[0000] %s engine_policy: engine policies are [ ",
           clocks_get_timesincestart());
    for (int k = 1; k < 32; k++)
      if (e->policy & (1 << k)) printf(" %s ", engine_policy_names[k + 1]);
    printf(" ]\n");
    fflush(stdout);
  }
#else
  printf("%s engine_policy: engine policies are [ ",
         clocks_get_timesincestart());
  for (int k = 1; k < 32; k++)
    if (e->policy & (1 << k)) printf(" %s ", engine_policy_names[k + 1]);
  printf(" ]\n");
  fflush(stdout);
#endif
}
