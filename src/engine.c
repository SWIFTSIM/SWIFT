/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef HAVE_LIBNUMA
#include <numa.h>
#endif

/* Load the profiler header, if needed. */
#ifdef WITH_PROFILER
#include <gperftools/profiler.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "cell.h"
#include "chemistry.h"
#include "clocks.h"
#include "cooling.h"
#include "cosmology.h"
#include "cycle.h"
#include "debug.h"
#include "error.h"
#include "gravity.h"
#include "gravity_cache.h"
#include "hydro.h"
#include "map.h"
#include "memswap.h"
#include "minmax.h"
#include "outputlist.h"
#include "parallel_io.h"
#include "part.h"
#include "partition.h"
#include "profiler.h"
#include "proxy.h"
#include "restart.h"
#include "runner.h"
#include "serial_io.h"
#include "single_io.h"
#include "sort_part.h"
#include "sourceterms.h"
#include "statistics.h"
#include "timers.h"
#include "tools.h"
#include "units.h"
#include "velociraptor_interface.h"
#include "version.h"

/* Particle cache size. */
#define CACHE_SIZE 512

const char *engine_policy_names[] = {"none",
                                     "rand",
                                     "steal",
                                     "keep",
                                     "block",
                                     "cpu tight",
                                     "mpi",
                                     "numa affinity",
                                     "hydro",
                                     "self gravity",
                                     "external gravity",
                                     "cosmological integration",
                                     "drift everything",
                                     "reconstruct multi-poles",
                                     "cooling",
                                     "sourceterms",
                                     "stars",
                                     "structure finding"};

/** The rank of the engine as a global variable (for messages). */
int engine_rank;

/**
 * @brief Data collected from the cells at the end of a time-step
 */
struct end_of_step_data {

  size_t updates, g_updates, s_updates;
  integertime_t ti_hydro_end_min, ti_hydro_end_max, ti_hydro_beg_max;
  integertime_t ti_gravity_end_min, ti_gravity_end_max, ti_gravity_beg_max;
  struct engine *e;
};

/**
 * @brief Link a density/force task to a cell.
 *
 * @param e The #engine.
 * @param l A pointer to the #link, will be modified atomically.
 * @param t The #task.
 *
 * @return The new #link pointer.
 */
void engine_addlink(struct engine *e, struct link **l, struct task *t) {

  /* Get the next free link. */
  const size_t ind = atomic_inc(&e->nr_links);
  if (ind >= e->size_links) {
    error("Link table overflow.");
  }
  struct link *res = &e->links[ind];

  /* Set it atomically. */
  res->t = t;
  res->next = atomic_swap(l, res);
}

/**
 * @brief Recursively add non-implicit ghost tasks to a cell hierarchy.
 */
void engine_add_ghosts(struct engine *e, struct cell *c, struct task *ghost_in,
                       struct task *ghost_out) {

  /* If we have reached the leaf OR have to few particles to play with*/
  if (!c->split || c->count < engine_max_parts_per_ghost) {

    /* Add the ghost task and its dependencies */
    struct scheduler *s = &e->sched;
    c->ghost =
        scheduler_addtask(s, task_type_ghost, task_subtype_none, 0, 0, c, NULL);
    scheduler_addunlock(s, ghost_in, c->ghost);
    scheduler_addunlock(s, c->ghost, ghost_out);
  } else {
    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_add_ghosts(e, c->progeny[k], ghost_in, ghost_out);
  }
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- timestep version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_common(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;
  const int is_with_cooling = (e->policy & engine_policy_cooling);

  /* Are we in a super-cell ? */
  if (c->super == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the two half kicks */
      c->kick1 = scheduler_addtask(s, task_type_kick1, task_subtype_none, 0, 0,
                                   c, NULL);

      c->kick2 = scheduler_addtask(s, task_type_kick2, task_subtype_none, 0, 0,
                                   c, NULL);

      /* Add the time-step calculation task and its dependency */
      c->timestep = scheduler_addtask(s, task_type_timestep, task_subtype_none,
                                      0, 0, c, NULL);

      /* Add the task finishing the force calculation */
      c->end_force = scheduler_addtask(s, task_type_end_force,
                                       task_subtype_none, 0, 0, c, NULL);

      if (!is_with_cooling) scheduler_addunlock(s, c->end_force, c->kick2);
      scheduler_addunlock(s, c->kick2, c->timestep);
      scheduler_addunlock(s, c->timestep, c->kick1);
    }

  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_common(e, c->progeny[k]);
  }
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- hydro version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_hydro(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;
  const int is_with_cooling = (e->policy & engine_policy_cooling);
  const int is_with_sourceterms = (e->policy & engine_policy_sourceterms);

  /* Are we in a super-cell ? */
  if (c->super_hydro == c) {

    /* Add the sort task. */
    c->sorts =
        scheduler_addtask(s, task_type_sort, task_subtype_none, 0, 0, c, NULL);

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      /* Add the drift task. */
      c->drift_part = scheduler_addtask(s, task_type_drift_part,
                                        task_subtype_none, 0, 0, c, NULL);

      /* Generate the ghost tasks. */
      c->ghost_in =
          scheduler_addtask(s, task_type_ghost_in, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      c->ghost_out =
          scheduler_addtask(s, task_type_ghost_out, task_subtype_none, 0,
                            /* implicit = */ 1, c, NULL);
      engine_add_ghosts(e, c, c->ghost_in, c->ghost_out);

#ifdef EXTRA_HYDRO_LOOP
      /* Generate the extra ghost task. */
      c->extra_ghost = scheduler_addtask(s, task_type_extra_ghost,
                                         task_subtype_none, 0, 0, c, NULL);
#endif

      /* Cooling task */
      if (is_with_cooling) {
        c->cooling = scheduler_addtask(s, task_type_cooling, task_subtype_none,
                                       0, 0, c, NULL);

        scheduler_addunlock(s, c->super->end_force, c->cooling);
        scheduler_addunlock(s, c->cooling, c->super->kick2);
      }

      /* add source terms */
      if (is_with_sourceterms) {
        c->sourceterms = scheduler_addtask(s, task_type_sourceterms,
                                           task_subtype_none, 0, 0, c, NULL);
      }
    }

  } else { /* We are above the super-cell so need to go deeper */

    /* Recurse. */
    if (c->split)
      for (int k = 0; k < 8; k++)
        if (c->progeny[k] != NULL)
          engine_make_hierarchical_tasks_hydro(e, c->progeny[k]);
  }
}

/**
 * @brief Generate the hydro hierarchical tasks for a hierarchy of cells -
 * i.e. all the O(Npart) tasks -- gravity version
 *
 * Tasks are only created here. The dependencies will be added later on.
 *
 * Note that there is no need to recurse below the super-cell. Note also
 * that we only add tasks if the relevant particles are present in the cell.
 *
 * @param e The #engine.
 * @param c The #cell.
 */
void engine_make_hierarchical_tasks_gravity(struct engine *e, struct cell *c) {

  struct scheduler *s = &e->sched;
  const int periodic = e->s->periodic;
  const int is_self_gravity = (e->policy & engine_policy_self_gravity);

  /* Are we in a super-cell ? */
  if (c->super_gravity == c) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      c->drift_gpart = scheduler_addtask(s, task_type_drift_gpart,
                                         task_subtype_none, 0, 0, c, NULL);

      if (is_self_gravity) {

        /* Initialisation of the multipoles */
        c->init_grav = scheduler_addtask(s, task_type_init_grav,
                                         task_subtype_none, 0, 0, c, NULL);

        /* Gravity non-neighbouring pm calculations */
        c->grav_long_range = scheduler_addtask(
            s, task_type_grav_long_range, task_subtype_none, 0, 0, c, NULL);

        /* Gravity recursive down-pass */
        c->grav_down = scheduler_addtask(s, task_type_grav_down,
                                         task_subtype_none, 0, 0, c, NULL);

        /* Implicit tasks for the up and down passes */
        c->init_grav_out = scheduler_addtask(s, task_type_init_grav_out,
                                             task_subtype_none, 0, 1, c, NULL);
        c->grav_down_in = scheduler_addtask(s, task_type_grav_down_in,
                                            task_subtype_none, 0, 1, c, NULL);

        /* Gravity mesh force propagation */
        if (periodic)
          c->grav_mesh = scheduler_addtask(s, task_type_grav_mesh,
                                           task_subtype_none, 0, 0, c, NULL);

        if (periodic) scheduler_addunlock(s, c->drift_gpart, c->grav_mesh);
        if (periodic) scheduler_addunlock(s, c->grav_mesh, c->grav_down);
        scheduler_addunlock(s, c->init_grav, c->grav_long_range);
        scheduler_addunlock(s, c->grav_long_range, c->grav_down);
        scheduler_addunlock(s, c->grav_down, c->super->end_force);

        /* Link in the implicit tasks */
        scheduler_addunlock(s, c->init_grav, c->init_grav_out);
        scheduler_addunlock(s, c->grav_down_in, c->grav_down);
      }
    }
  }

  /* We are below the super-cell but not below the maximal splitting depth */
  else if (c->super_gravity != NULL && c->depth <= space_subdepth_grav) {

    /* Local tasks only... */
    if (c->nodeID == e->nodeID) {

      if (is_self_gravity) {

        c->init_grav_out = scheduler_addtask(s, task_type_init_grav_out,
                                             task_subtype_none, 0, 1, c, NULL);

        c->grav_down_in = scheduler_addtask(s, task_type_grav_down_in,
                                            task_subtype_none, 0, 1, c, NULL);

        scheduler_addunlock(s, c->parent->init_grav_out, c->init_grav_out);
        scheduler_addunlock(s, c->grav_down_in, c->parent->grav_down_in);
      }
    }
  }

  /* Recurse but not below the maximal splitting depth */
  if (c->split && c->depth <= space_subdepth_grav)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_make_hierarchical_tasks_gravity(e, c->progeny[k]);
}

void engine_make_hierarchical_tasks_mapper(void *map_data, int num_elements,
                                           void *extra_data) {
  struct engine *e = (struct engine *)extra_data;
  const int is_with_hydro = (e->policy & engine_policy_hydro);
  const int is_with_self_gravity = (e->policy & engine_policy_self_gravity);
  const int is_with_external_gravity =
      (e->policy & engine_policy_external_gravity);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &((struct cell *)map_data)[ind];
    /* Make the common tasks (time integration) */
    engine_make_hierarchical_tasks_common(e, c);
    /* Add the hydro stuff */
    if (is_with_hydro) engine_make_hierarchical_tasks_hydro(e, c);
    /* And the gravity stuff */
    if (is_with_self_gravity || is_with_external_gravity)
      engine_make_hierarchical_tasks_gravity(e, c);
  }
}

#ifdef WITH_MPI
/**
 * Do the exchange of one type of particles with all the other nodes.
 *
 * @param counts 2D array with the counts of particles to exchange with
 *               each other node.
 * @param parts the particle data to exchange
 * @param new_nr_parts the number of particles this node will have after all
 *                     exchanges have completed.
 * @param sizeofparts sizeof the particle struct.
 * @param alignsize the memory alignment required for this particle type.
 * @param mpi_type the MPI_Datatype for these particles.
 * @param nr_nodes the number of nodes to exchange with.
 * @param nodeID the id of this node.
 *
 * @result new particle data constructed from all the exchanges with the
 *         given alignment.
 */
static void *engine_do_redistribute(int *counts, char *parts,
                                    size_t new_nr_parts, size_t sizeofparts,
                                    size_t alignsize, MPI_Datatype mpi_type,
                                    int nr_nodes, int nodeID) {

  /* Allocate a new particle array with some extra margin */
  char *parts_new = NULL;
  if (posix_memalign(
          (void **)&parts_new, alignsize,
          sizeofparts * new_nr_parts * engine_redistribute_alloc_margin) != 0)
    error("Failed to allocate new particle data.");

  /* Prepare MPI requests for the asynchronous communications */
  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 2 * nr_nodes)) ==
      NULL)
    error("Failed to allocate MPI request list.");

  /* Only send and receive only "chunk" particles per request. So we need to
   * loop as many times as necessary here. Make 2Gb/sizeofparts so we only
   * send 2Gb packets. */
  const int chunk = INT_MAX / sizeofparts;
  int sent = 0;
  int recvd = 0;

  int activenodes = 1;
  while (activenodes) {

    for (int k = 0; k < 2 * nr_nodes; k++) reqs[k] = MPI_REQUEST_NULL;

    /* Emit the sends and recvs for the data. */
    size_t offset_send = sent;
    size_t offset_recv = recvd;
    activenodes = 0;

    for (int k = 0; k < nr_nodes; k++) {

      /* Indices in the count arrays of the node of interest */
      const int ind_send = nodeID * nr_nodes + k;
      const int ind_recv = k * nr_nodes + nodeID;

      /* Are we sending any data this loop? */
      int sending = counts[ind_send] - sent;
      if (sending > 0) {
        activenodes++;
        if (sending > chunk) sending = chunk;

        /* If the send and receive is local then just copy. */
        if (k == nodeID) {
          int receiving = counts[ind_recv] - recvd;
          if (receiving > chunk) receiving = chunk;
          memcpy(&parts_new[offset_recv * sizeofparts],
                 &parts[offset_send * sizeofparts], sizeofparts * receiving);
        } else {
          /* Otherwise send it. */
          int res =
              MPI_Isend(&parts[offset_send * sizeofparts], sending, mpi_type, k,
                        ind_send, MPI_COMM_WORLD, &reqs[2 * k + 0]);
          if (res != MPI_SUCCESS)
            mpi_error(res, "Failed to isend parts to node %i.", k);
        }
      }

      /* If we're sending to this node, then move past it to next. */
      if (counts[ind_send] > 0) offset_send += counts[ind_send];

      /* Are we receiving any data from this node? Note already done if coming
       * from this node. */
      if (k != nodeID) {
        int receiving = counts[ind_recv] - recvd;
        if (receiving > 0) {
          activenodes++;
          if (receiving > chunk) receiving = chunk;
          int res = MPI_Irecv(&parts_new[offset_recv * sizeofparts], receiving,
                              mpi_type, k, ind_recv, MPI_COMM_WORLD,
                              &reqs[2 * k + 1]);
          if (res != MPI_SUCCESS)
            mpi_error(res, "Failed to emit irecv of parts from node %i.", k);
        }
      }

      /* If we're receiving from this node, then move past it to next. */
      if (counts[ind_recv] > 0) offset_recv += counts[ind_recv];
    }

    /* Wait for all the sends and recvs to tumble in. */
    MPI_Status stats[2 * nr_nodes];
    int res;
    if ((res = MPI_Waitall(2 * nr_nodes, reqs, stats)) != MPI_SUCCESS) {
      for (int k = 0; k < 2 * nr_nodes; k++) {
        char buff[MPI_MAX_ERROR_STRING];
        MPI_Error_string(stats[k].MPI_ERROR, buff, &res);
        message("request from source %i, tag %i has error '%s'.",
                stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
      }
      error("Failed during waitall for part data.");
    }

    /* Move to next chunks. */
    sent += chunk;
    recvd += chunk;
  }

  /* Free temps. */
  free(reqs);

  /* And return new memory. */
  return parts_new;
}
#endif

#ifdef WITH_MPI /* redist_mapper */

/* Support for engine_redistribute threadpool dest mappers. */
struct redist_mapper_data {
  int *counts;
  int *dest;
  int nodeID;
  int nr_nodes;
  struct cell *cells;
  struct space *s;
  void *base;
};

/* Generic function for accumulating counts for TYPE parts. Note
 * we use a local counts array to avoid the atomic_add in the parts
 * loop. */
#define ENGINE_REDISTRIBUTE_DEST_MAPPER(TYPE)                              \
  engine_redistribute_dest_mapper_##TYPE(void *map_data, int num_elements, \
                                         void *extra_data) {               \
    struct TYPE *parts = (struct TYPE *)map_data;                          \
    struct redist_mapper_data *mydata =                                    \
        (struct redist_mapper_data *)extra_data;                           \
    struct space *s = mydata->s;                                           \
    int *dest =                                                            \
        mydata->dest + (ptrdiff_t)(parts - (struct TYPE *)mydata->base);   \
    int *lcounts = NULL;                                                   \
    if ((lcounts = (int *)calloc(                                          \
             sizeof(int), mydata->nr_nodes * mydata->nr_nodes)) == NULL)   \
      error("Failed to allocate counts thread-specific buffer");           \
    for (int k = 0; k < num_elements; k++) {                               \
      for (int j = 0; j < 3; j++) {                                        \
        if (parts[k].x[j] < 0.0)                                           \
          parts[k].x[j] += s->dim[j];                                      \
        else if (parts[k].x[j] >= s->dim[j])                               \
          parts[k].x[j] -= s->dim[j];                                      \
      }                                                                    \
      const int cid = cell_getid(s->cdim, parts[k].x[0] * s->iwidth[0],    \
                                 parts[k].x[1] * s->iwidth[1],             \
                                 parts[k].x[2] * s->iwidth[2]);            \
      dest[k] = s->cells_top[cid].nodeID;                                  \
      size_t ind = mydata->nodeID * mydata->nr_nodes + dest[k];            \
      lcounts[ind] += 1;                                                   \
    }                                                                      \
    for (int k = 0; k < (mydata->nr_nodes * mydata->nr_nodes); k++)        \
      atomic_add(&mydata->counts[k], lcounts[k]);                          \
    free(lcounts);                                                         \
  }

/**
 * @brief Accumulate the counts of particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * part version.
 */
static void ENGINE_REDISTRIBUTE_DEST_MAPPER(part);

/**
 * @brief Accumulate the counts of star particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * spart version.
 */
static void ENGINE_REDISTRIBUTE_DEST_MAPPER(spart);

/**
 * @brief Accumulate the counts of gravity particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * gpart version.
 */
static void ENGINE_REDISTRIBUTE_DEST_MAPPER(gpart);

#endif /* redist_mapper_data */

#ifdef WITH_MPI /* savelink_mapper_data */

/* Support for saving the linkage between gparts and parts/sparts. */
struct savelink_mapper_data {
  int nr_nodes;
  int *counts;
  void *parts;
  int nodeID;
};

/**
 * @brief Save the offset of each gravity partner of a part or spart.
 *
 * The offset is from the start of the sorted particles to be sent to a node.
 * This is possible as parts without gravity partners have a positive id.
 * These offsets are used to restore the pointers on the receiving node.
 *
 * CHECKS should be eliminated as dead code when optimizing.
 */
#define ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(TYPE, CHECKS)                      \
  engine_redistribute_savelink_mapper_##TYPE(void *map_data, int num_elements, \
                                             void *extra_data) {               \
    int *nodes = (int *)map_data;                                              \
    struct savelink_mapper_data *mydata =                                      \
        (struct savelink_mapper_data *)extra_data;                             \
    int nodeID = mydata->nodeID;                                               \
    int nr_nodes = mydata->nr_nodes;                                           \
    int *counts = mydata->counts;                                              \
    struct TYPE *parts = (struct TYPE *)mydata->parts;                         \
                                                                               \
    for (int j = 0; j < num_elements; j++) {                                   \
      int node = nodes[j];                                                     \
      int count = 0;                                                           \
      size_t offset = 0;                                                       \
      for (int i = 0; i < node; i++) offset += counts[nodeID * nr_nodes + i];  \
                                                                               \
      for (int k = 0; k < counts[nodeID * nr_nodes + node]; k++) {             \
        if (parts[k + offset].gpart != NULL) {                                 \
          if (CHECKS)                                                          \
            if (parts[k].gpart->id_or_neg_offset > 0)                          \
              error("Trying to link a partnerless " #TYPE "!");                \
          parts[k + offset].gpart->id_or_neg_offset = -count;                  \
          count++;                                                             \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  }

/**
 * @brief Save position of part-gpart links.
 * Threadpool helper for accumulating the counts of particles per cell.
 */
#ifdef SWIFT_DEBUG_CHECKS
static void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(part, 1);
#else
static void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(part, 0);
#endif

/**
 * @brief Save position of spart-gpart links.
 * Threadpool helper for accumulating the counts of particles per cell.
 */
#ifdef SWIFT_DEBUG_CHECKS
static void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(spart, 1);
#else
static void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(spart, 0);
#endif

#endif /* savelink_mapper_data */

#ifdef WITH_MPI /* relink_mapper_data */

/* Support for relinking parts, gparts and sparts after moving between nodes. */
struct relink_mapper_data {
  int nodeID;
  int nr_nodes;
  int *counts;
  int *s_counts;
  int *g_counts;
  struct space *s;
};

/**
 * @brief Restore the part/gpart and spart/gpart links for a list of nodes.
 *
 * @param map_data address of nodes to process.
 * @param num_elements the number nodes to process.
 * @param extra_data additional data defining the context (a
 * relink_mapper_data).
 */
static void engine_redistribute_relink_mapper(void *map_data, int num_elements,
                                              void *extra_data) {

  int *nodes = (int *)map_data;
  struct relink_mapper_data *mydata = (struct relink_mapper_data *)extra_data;

  int nodeID = mydata->nodeID;
  int nr_nodes = mydata->nr_nodes;
  int *counts = mydata->counts;
  int *g_counts = mydata->g_counts;
  int *s_counts = mydata->s_counts;
  struct space *s = mydata->s;

  for (int i = 0; i < num_elements; i++) {

    int node = nodes[i];

    /* Get offsets to correct parts of the counts arrays for this node. */
    size_t offset_parts = 0;
    size_t offset_gparts = 0;
    size_t offset_sparts = 0;
    for (int n = 0; n < node; n++) {
      int ind_recv = n * nr_nodes + nodeID;
      offset_parts += counts[ind_recv];
      offset_gparts += g_counts[ind_recv];
      offset_sparts += s_counts[ind_recv];
    }

    /* Number of gparts sent from this node. */
    int ind_recv = node * nr_nodes + nodeID;
    const size_t count_gparts = g_counts[ind_recv];

    /* Loop over the gparts received from this node */
    for (size_t k = offset_gparts; k < offset_gparts + count_gparts; k++) {

      /* Does this gpart have a gas partner ? */
      if (s->gparts[k].type == swift_type_gas) {

        const ptrdiff_t partner_index =
            offset_parts - s->gparts[k].id_or_neg_offset;

        /* Re-link */
        s->gparts[k].id_or_neg_offset = -partner_index;
        s->parts[partner_index].gpart = &s->gparts[k];
      }

      /* Does this gpart have a star partner ? */
      else if (s->gparts[k].type == swift_type_star) {

        const ptrdiff_t partner_index =
            offset_sparts - s->gparts[k].id_or_neg_offset;

        /* Re-link */
        s->gparts[k].id_or_neg_offset = -partner_index;
        s->sparts[partner_index].gpart = &s->gparts[k];
      }
    }
  }
}

#endif /* relink_mapper_data */

/**
 * @brief Redistribute the particles amongst the nodes according
 *      to their cell's node IDs.
 *
 * The strategy here is as follows:
 * 1) Each node counts the number of particles it has to send to each other
 * node.
 * 2) The number of particles of each type is then exchanged.
 * 3) The particles to send are placed in a temporary buffer in which the
 * part-gpart links are preserved.
 * 4) Each node allocates enough space for the new particles.
 * 5) (Asynchronous) communications are issued to transfer the data.
 *
 *
 * @param e The #engine.
 */
void engine_redistribute(struct engine *e) {

#ifdef WITH_MPI

  const int nr_nodes = e->nr_nodes;
  const int nodeID = e->nodeID;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;
  struct part *parts = s->parts;
  struct gpart *gparts = s->gparts;
  struct spart *sparts = s->sparts;
  ticks tic = getticks();

  /* Allocate temporary arrays to store the counts of particles to be sent
   * and the destination of each particle */
  int *counts;
  if ((counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate counts temporary buffer.");

  int *dest;
  if ((dest = (int *)malloc(sizeof(int) * s->nr_parts)) == NULL)
    error("Failed to allocate dest temporary buffer.");

  /* Simple index of node IDs, used for mappers over nodes. */
  int *nodes = NULL;
  if ((nodes = (int *)malloc(sizeof(int) * nr_nodes)) == NULL)
    error("Failed to allocate nodes temporary buffer.");
  for (int k = 0; k < nr_nodes; k++) nodes[k] = k;

  /* Get destination of each particle */
  struct redist_mapper_data redist_data;
  redist_data.s = s;
  redist_data.nodeID = nodeID;
  redist_data.nr_nodes = nr_nodes;

  redist_data.counts = counts;
  redist_data.dest = dest;
  redist_data.base = (void *)parts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_part, parts,
                 s->nr_parts, sizeof(struct part), 0, &redist_data);

  /* Sort the particles according to their cell index. */
  if (s->nr_parts > 0)
    space_parts_sort(s->parts, s->xparts, dest, &counts[nodeID * nr_nodes],
                     nr_nodes, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < s->nr_parts; k++) {
    const struct part *p = &s->parts[k];

    /* New cell index */
    const int new_cid =
        cell_getid(s->cdim, p->x[0] * s->iwidth[0], p->x[1] * s->iwidth[1],
                   p->x[2] * s->iwidth[2]);

    /* New cell of this part */
    const struct cell *c = &s->cells_top[new_cid];
    const int new_node = c->nodeID;

    if (dest[k] != new_node)
      error("part's new node index not matching sorted index.");

    if (p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] < c->loc[1] || p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] < c->loc[2] || p->x[2] > c->loc[2] + c->width[2])
      error("part not sorted into the right top-level cell!");
  }
#endif

  /* We will need to re-link the gpart partners of parts, so save their
   * relative positions in the sent lists. */
  if (s->nr_parts > 0 && s->nr_gparts > 0) {

    struct savelink_mapper_data savelink_data;
    savelink_data.nr_nodes = nr_nodes;
    savelink_data.counts = counts;
    savelink_data.parts = (void *)parts;
    savelink_data.nodeID = nodeID;
    threadpool_map(&e->threadpool, engine_redistribute_savelink_mapper_part,
                   nodes, nr_nodes, sizeof(int), 0, &savelink_data);
  }
  free(dest);

  /* Get destination of each s-particle */
  int *s_counts;
  if ((s_counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate s_counts temporary buffer.");

  int *s_dest;
  if ((s_dest = (int *)malloc(sizeof(int) * s->nr_sparts)) == NULL)
    error("Failed to allocate s_dest temporary buffer.");

  redist_data.counts = s_counts;
  redist_data.dest = s_dest;
  redist_data.base = (void *)sparts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_spart, sparts,
                 s->nr_sparts, sizeof(struct spart), 0, &redist_data);

  /* Sort the particles according to their cell index. */
  if (s->nr_sparts > 0)
    space_sparts_sort(s->sparts, s_dest, &s_counts[nodeID * nr_nodes], nr_nodes,
                      0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < s->nr_sparts; k++) {
    const struct spart *sp = &s->sparts[k];

    /* New cell index */
    const int new_cid =
        cell_getid(s->cdim, sp->x[0] * s->iwidth[0], sp->x[1] * s->iwidth[1],
                   sp->x[2] * s->iwidth[2]);

    /* New cell of this spart */
    const struct cell *c = &s->cells_top[new_cid];
    const int new_node = c->nodeID;

    if (s_dest[k] != new_node)
      error("spart's new node index not matching sorted index.");

    if (sp->x[0] < c->loc[0] || sp->x[0] > c->loc[0] + c->width[0] ||
        sp->x[1] < c->loc[1] || sp->x[1] > c->loc[1] + c->width[1] ||
        sp->x[2] < c->loc[2] || sp->x[2] > c->loc[2] + c->width[2])
      error("spart not sorted into the right top-level cell!");
  }
#endif

  /* We need to re-link the gpart partners of sparts. */
  if (s->nr_sparts > 0) {

    struct savelink_mapper_data savelink_data;
    savelink_data.nr_nodes = nr_nodes;
    savelink_data.counts = s_counts;
    savelink_data.parts = (void *)sparts;
    savelink_data.nodeID = nodeID;
    threadpool_map(&e->threadpool, engine_redistribute_savelink_mapper_spart,
                   nodes, nr_nodes, sizeof(int), 0, &savelink_data);
  }
  free(s_dest);

  /* Get destination of each g-particle */
  int *g_counts;
  if ((g_counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate g_gcount temporary buffer.");

  int *g_dest;
  if ((g_dest = (int *)malloc(sizeof(int) * s->nr_gparts)) == NULL)
    error("Failed to allocate g_dest temporary buffer.");

  redist_data.counts = g_counts;
  redist_data.dest = g_dest;
  redist_data.base = (void *)gparts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_gpart, gparts,
                 s->nr_gparts, sizeof(struct gpart), 0, &redist_data);

  /* Sort the gparticles according to their cell index. */
  if (s->nr_gparts > 0)
    space_gparts_sort(s->gparts, s->parts, s->sparts, g_dest,
                      &g_counts[nodeID * nr_nodes], nr_nodes);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < s->nr_gparts; k++) {
    const struct gpart *gp = &s->gparts[k];

    /* New cell index */
    const int new_cid =
        cell_getid(s->cdim, gp->x[0] * s->iwidth[0], gp->x[1] * s->iwidth[1],
                   gp->x[2] * s->iwidth[2]);

    /* New cell of this gpart */
    const struct cell *c = &s->cells_top[new_cid];
    const int new_node = c->nodeID;

    if (g_dest[k] != new_node)
      error("gpart's new node index not matching sorted index (%d != %d).",
            g_dest[k], new_node);

    if (gp->x[0] < c->loc[0] || gp->x[0] > c->loc[0] + c->width[0] ||
        gp->x[1] < c->loc[1] || gp->x[1] > c->loc[1] + c->width[1] ||
        gp->x[2] < c->loc[2] || gp->x[2] > c->loc[2] + c->width[2])
      error("gpart not sorted into the right top-level cell!");
  }
#endif

  free(g_dest);

  /* Get all the counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, counts, nr_nodes * nr_nodes, MPI_INT, MPI_SUM,
                    MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce particle transfer counts.");

  /* Get all the s_counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, g_counts, nr_nodes * nr_nodes, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce gparticle transfer counts.");

  /* Get all the g_counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, s_counts, nr_nodes * nr_nodes, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce sparticle transfer counts.");

  /* Report how many particles will be moved. */
  if (e->verbose) {
    if (e->nodeID == 0) {
      size_t total = 0, g_total = 0, s_total = 0;
      size_t unmoved = 0, g_unmoved = 0, s_unmoved = 0;
      for (int p = 0, r = 0; p < nr_nodes; p++) {
        for (int n = 0; n < nr_nodes; n++) {
          total += counts[r];
          g_total += g_counts[r];
          s_total += s_counts[r];
          if (p == n) {
            unmoved += counts[r];
            g_unmoved += g_counts[r];
            s_unmoved += s_counts[r];
          }
          r++;
        }
      }
      if (total > 0)
        message("%ld of %ld (%.2f%%) of particles moved", total - unmoved,
                total, 100.0 * (double)(total - unmoved) / (double)total);
      if (g_total > 0)
        message("%ld of %ld (%.2f%%) of g-particles moved", g_total - g_unmoved,
                g_total,
                100.0 * (double)(g_total - g_unmoved) / (double)g_total);
      if (s_total > 0)
        message("%ld of %ld (%.2f%%) of s-particles moved", s_total - s_unmoved,
                s_total,
                100.0 * (double)(s_total - s_unmoved) / (double)s_total);
    }
  }

  /* Now each node knows how many parts, sparts and gparts will be transferred
   * to every other node.
   * Get the new numbers of particles for this node. */
  size_t nr_parts = 0, nr_gparts = 0, nr_sparts = 0;
  for (int k = 0; k < nr_nodes; k++) nr_parts += counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_gparts += g_counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_sparts += s_counts[k * nr_nodes + nodeID];

  /* Now exchange the particles, type by type to keep the memory required
   * under control. */

  /* SPH particles. */
  void *new_parts = engine_do_redistribute(counts, (char *)s->parts, nr_parts,
                                           sizeof(struct part), part_align,
                                           part_mpi_type, nr_nodes, nodeID);
  free(s->parts);
  s->parts = (struct part *)new_parts;
  s->nr_parts = nr_parts;
  s->size_parts = engine_redistribute_alloc_margin * nr_parts;

  /* Extra SPH particle properties. */
  new_parts = engine_do_redistribute(counts, (char *)s->xparts, nr_parts,
                                     sizeof(struct xpart), xpart_align,
                                     xpart_mpi_type, nr_nodes, nodeID);
  free(s->xparts);
  s->xparts = (struct xpart *)new_parts;

  /* Gravity particles. */
  new_parts = engine_do_redistribute(g_counts, (char *)s->gparts, nr_gparts,
                                     sizeof(struct gpart), gpart_align,
                                     gpart_mpi_type, nr_nodes, nodeID);
  free(s->gparts);
  s->gparts = (struct gpart *)new_parts;
  s->nr_gparts = nr_gparts;
  s->size_gparts = engine_redistribute_alloc_margin * nr_gparts;

  /* Star particles. */
  new_parts = engine_do_redistribute(s_counts, (char *)s->sparts, nr_sparts,
                                     sizeof(struct spart), spart_align,
                                     spart_mpi_type, nr_nodes, nodeID);
  free(s->sparts);
  s->sparts = (struct spart *)new_parts;
  s->nr_sparts = nr_sparts;
  s->size_sparts = engine_redistribute_alloc_margin * nr_sparts;

  /* All particles have now arrived. Time for some final operations on the
     stuff we just received */

  /* Restore the part<->gpart and spart<->gpart links.
   * Generate indices and counts for threadpool tasks. Note we process a node
   * at a time. */
  struct relink_mapper_data relink_data;
  relink_data.s = s;
  relink_data.counts = counts;
  relink_data.g_counts = g_counts;
  relink_data.s_counts = s_counts;
  relink_data.nodeID = nodeID;
  relink_data.nr_nodes = nr_nodes;

  threadpool_map(&e->threadpool, engine_redistribute_relink_mapper, nodes,
                 nr_nodes, sizeof(int), 1, &relink_data);
  free(nodes);

  /* Clean up the counts now we are done. */
  free(counts);
  free(g_counts);
  free(s_counts);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that all parts are in the right place. */
  for (size_t k = 0; k < nr_parts; k++) {
    const int cid = cell_getid(s->cdim, s->parts[k].x[0] * s->iwidth[0],
                               s->parts[k].x[1] * s->iwidth[1],
                               s->parts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received particle (%zu) that does not belong here (nodeID=%i).", k,
            cells[cid].nodeID);
  }
  for (size_t k = 0; k < nr_gparts; k++) {
    const int cid = cell_getid(s->cdim, s->gparts[k].x[0] * s->iwidth[0],
                               s->gparts[k].x[1] * s->iwidth[1],
                               s->gparts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received g-particle (%zu) that does not belong here (nodeID=%i).",
            k, cells[cid].nodeID);
  }
  for (size_t k = 0; k < nr_sparts; k++) {
    const int cid = cell_getid(s->cdim, s->sparts[k].x[0] * s->iwidth[0],
                               s->sparts[k].x[1] * s->iwidth[1],
                               s->sparts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received s-particle (%zu) that does not belong here (nodeID=%i).",
            k, cells[cid].nodeID);
  }

  /* Verify that the links are correct */
  part_verify_links(s->parts, s->gparts, s->sparts, nr_parts, nr_gparts,
                    nr_sparts, e->verbose);
#endif

  /* Be verbose about what just happened. */
  if (e->verbose) {
    int my_cells = 0;
    for (int k = 0; k < nr_cells; k++)
      if (cells[k].nodeID == nodeID) my_cells += 1;
    message("node %i now has %zu parts, %zu sparts and %zu gparts in %i cells.",
            nodeID, nr_parts, nr_sparts, nr_gparts, my_cells);
  }

  /* Flag that a redistribute has taken place */
  e->step_props |= engine_step_prop_redistribute;

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

#ifdef SWIFT_DEBUG_CHECKS
  /* Be verbose about this. */
  if (e->nodeID == 0 || e->verbose) message("repartitioning space");
  fflush(stdout);

  /* Check that all cells have been drifted to the current time */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);
#endif

  /* Clear the repartition flag. */
  e->forcerepart = 0;

  /* Nothing to do if only using a single node. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (e->nr_nodes == 1) return;

  /* Do the repartitioning. */
  partition_repartition(e->reparttype, e->nodeID, e->nr_nodes, e->s,
                        e->sched.tasks, e->sched.nr_tasks);

  /* Partitioning requires copies of the particles, so we need to reduce the
   * memory in use to the minimum, we can free the sorting indices and the
   * tasks as these will be regenerated at the next rebuild. */

  /* Sorting indices. */
  if (e->s->cells_top != NULL) space_free_cells(e->s);

  /* Task arrays. */
  scheduler_free_tasks(&e->sched);

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

  /* Flag that a repartition has taken place */
  e->step_props |= engine_step_prop_repartition;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  if (e->reparttype->type != REPART_NONE)
    error("SWIFT was not compiled with MPI and METIS support.");

  /* Clear the repartition flag. */
  e->forcerepart = 0;
#endif
}

/**
 * @brief Decide whether trigger a repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */
void engine_repartition_trigger(struct engine *e) {

#ifdef WITH_MPI

  /* Do nothing if there have not been enough steps since the last
   * repartition, don't want to repeat this too often or immediately after
   * a repartition step. Also nothing to do when requested. */
  if (e->step - e->last_repartition >= 2 &&
      e->reparttype->type != REPART_NONE) {

    /* Old style if trigger is >1 or this is the second step (want an early
     * repartition following the initial repartition). */
    if (e->reparttype->trigger > 1 || e->step == 2) {
      if (e->reparttype->trigger > 1) {
        if ((e->step % (int)e->reparttype->trigger) == 0) e->forcerepart = 1;
      } else {
        e->forcerepart = 1;
      }

    } else {

      /* Use cputimes from ranks to estimate the imbalance. */
      /* First check if we are going to skip this stage anyway, if so do that
       * now. If is only worth checking the CPU loads when we have processed a
       * significant number of all particles. */
      if ((e->updates > 1 &&
           e->updates >= e->total_nr_parts * e->reparttype->minfrac) ||
          (e->g_updates > 1 &&
           e->g_updates >= e->total_nr_gparts * e->reparttype->minfrac)) {

        /* Get CPU time used since the last call to this function. */
        double elapsed_cputime =
            clocks_get_cputime_used() - e->cputime_last_step;

        /* Gather the elapsed CPU times from all ranks for the last step. */
        double elapsed_cputimes[e->nr_nodes];
        MPI_Gather(&elapsed_cputime, 1, MPI_DOUBLE, elapsed_cputimes, 1,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (e->nodeID == 0) {

          /* Get the range and mean of cputimes. */
          double mintime = elapsed_cputimes[0];
          double maxtime = elapsed_cputimes[0];
          double sum = elapsed_cputimes[0];
          for (int k = 1; k < e->nr_nodes; k++) {
            if (elapsed_cputimes[k] > maxtime) maxtime = elapsed_cputimes[k];
            if (elapsed_cputimes[k] < mintime) mintime = elapsed_cputimes[k];
            sum += elapsed_cputimes[k];
          }
          double mean = sum / (double)e->nr_nodes;

          /* Are we out of balance? */
          if (((maxtime - mintime) / mean) > e->reparttype->trigger) {
            if (e->verbose)
              message("trigger fraction %.3f exceeds %.3f will repartition",
                      (maxtime - mintime) / mintime, e->reparttype->trigger);
            e->forcerepart = 1;
          }
        }

        /* All nodes do this together. */
        MPI_Bcast(&e->forcerepart, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
    }

    /* Remember we did this. */
    if (e->forcerepart) e->last_repartition = e->step;
  }

  /* We always reset CPU time for next check, unless it will not be used. */
  if (e->reparttype->type != REPART_NONE)
    e->cputime_last_step = clocks_get_cputime_used();
#endif
}

/**
 * @brief Add send tasks for the hydro pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_xv The send_xv #task, if it has already been created.
 * @param t_rho The send_rho #task, if it has already been created.
 * @param t_gradient The send_gradient #task, if already created.
 */
void engine_addtasks_send_hydro(struct engine *e, struct cell *ci,
                                struct cell *cj, struct task *t_xv,
                                struct task *t_rho, struct task *t_gradient) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Check if any of the density tasks are for the target node. */
  for (l = ci->density; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_xv == NULL) {

      t_xv = scheduler_addtask(s, task_type_send, task_subtype_xv,
                               6 * ci->tag + 0, 0, ci, cj);
      t_rho = scheduler_addtask(s, task_type_send, task_subtype_rho,
                                6 * ci->tag + 1, 0, ci, cj);
#ifdef EXTRA_HYDRO_LOOP
      t_gradient = scheduler_addtask(s, task_type_send, task_subtype_gradient,
                                     6 * ci->tag + 3, 0, ci, cj);
#endif

#ifdef EXTRA_HYDRO_LOOP

      scheduler_addunlock(s, t_gradient, ci->super->kick2);

      scheduler_addunlock(s, ci->super_hydro->extra_ghost, t_gradient);

      /* The send_rho task should unlock the super_hydro-cell's extra_ghost
       * task. */
      scheduler_addunlock(s, t_rho, ci->super_hydro->extra_ghost);

      /* The send_rho task depends on the cell's ghost task. */
      scheduler_addunlock(s, ci->super_hydro->ghost_out, t_rho);

      /* The send_xv task should unlock the super_hydro-cell's ghost task. */
      scheduler_addunlock(s, t_xv, ci->super_hydro->ghost_in);

#else
      /* The send_rho task should unlock the super_hydro-cell's kick task. */
      scheduler_addunlock(s, t_rho, ci->super->end_force);

      /* The send_rho task depends on the cell's ghost task. */
      scheduler_addunlock(s, ci->super_hydro->ghost_out, t_rho);

      /* The send_xv task should unlock the super_hydro-cell's ghost task. */
      scheduler_addunlock(s, t_xv, ci->super_hydro->ghost_in);

#endif

      /* Drift before you send */
      scheduler_addunlock(s, ci->super_hydro->drift_part, t_xv);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->send_xv, t_xv);
    engine_addlink(e, &ci->send_rho, t_rho);
#ifdef EXTRA_HYDRO_LOOP
    engine_addlink(e, &ci->send_gradient, t_gradient);
#endif
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_hydro(e, ci->progeny[k], cj, t_xv, t_rho,
                                   t_gradient);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_grav The send_grav #task, if it has already been created.
 */
void engine_addtasks_send_gravity(struct engine *e, struct cell *ci,
                                  struct cell *cj, struct task *t_grav) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Check if any of the gravity tasks are for the target node. */
  for (l = ci->grav; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* If so, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_grav == NULL) {

      t_grav = scheduler_addtask(s, task_type_send, task_subtype_gpart,
                                 6 * ci->tag + 4, 0, ci, cj);

      /* The sends should unlock the down pass. */
      scheduler_addunlock(s, t_grav, ci->super_gravity->grav_down);

      /* Drift before you send */
      scheduler_addunlock(s, ci->super_gravity->drift_gpart, t_grav);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->send_grav, t_grav);
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_gravity(e, ci->progeny[k], cj, t_grav);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add send tasks for the time-step to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param ci The sending #cell.
 * @param cj Dummy cell containing the nodeID of the receiving node.
 * @param t_ti The send_ti #task, if it has already been created.
 */
void engine_addtasks_send_timestep(struct engine *e, struct cell *ci,
                                   struct cell *cj, struct task *t_ti) {

#ifdef WITH_MPI
  struct link *l = NULL;
  struct scheduler *s = &e->sched;
  const int nodeID = cj->nodeID;

  /* Check if any of the gravity tasks are for the target node. */
  for (l = ci->grav; l != NULL; l = l->next)
    if (l->t->ci->nodeID == nodeID ||
        (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
      break;

  /* Check whether instead any of the hydro tasks are for the target node. */
  if (l == NULL)
    for (l = ci->density; l != NULL; l = l->next)
      if (l->t->ci->nodeID == nodeID ||
          (l->t->cj != NULL && l->t->cj->nodeID == nodeID))
        break;

  /* If found anything, attach send tasks. */
  if (l != NULL) {

    /* Create the tasks and their dependencies? */
    if (t_ti == NULL) {

      t_ti = scheduler_addtask(s, task_type_send, task_subtype_tend,
                               6 * ci->tag + 2, 0, ci, cj);

      /* The super-cell's timestep task should unlock the send_ti task. */
      scheduler_addunlock(s, ci->super->timestep, t_ti);
    }

    /* Add them to the local cell. */
    engine_addlink(e, &ci->send_ti, t_ti);
  }

  /* Recurse? */
  if (ci->split)
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL)
        engine_addtasks_send_timestep(e, ci->progeny[k], cj, t_ti);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for hydro pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_xv The recv_xv #task, if it has already been created.
 * @param t_rho The recv_rho #task, if it has already been created.
 * @param t_gradient The recv_gradient #task, if it has already been created.
 */
void engine_addtasks_recv_hydro(struct engine *e, struct cell *c,
                                struct task *t_xv, struct task *t_rho,
                                struct task *t_gradient) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Have we reached a level where there are any hydro tasks ? */
  if (t_xv == NULL && c->density != NULL) {

    /* Create the tasks. */
    t_xv = scheduler_addtask(s, task_type_recv, task_subtype_xv, 6 * c->tag + 0,
                             0, c, NULL);
    t_rho = scheduler_addtask(s, task_type_recv, task_subtype_rho,
                              6 * c->tag + 1, 0, c, NULL);
#ifdef EXTRA_HYDRO_LOOP
    t_gradient = scheduler_addtask(s, task_type_recv, task_subtype_gradient,
                                   6 * c->tag + 3, 0, c, NULL);
#endif
  }

  c->recv_xv = t_xv;
  c->recv_rho = t_rho;
  c->recv_gradient = t_gradient;

  /* Add dependencies. */
  if (c->sorts != NULL) scheduler_addunlock(s, t_xv, c->sorts);

  for (struct link *l = c->density; l != NULL; l = l->next) {
    scheduler_addunlock(s, t_xv, l->t);
    scheduler_addunlock(s, l->t, t_rho);
  }
#ifdef EXTRA_HYDRO_LOOP
  for (struct link *l = c->gradient; l != NULL; l = l->next) {
    scheduler_addunlock(s, t_rho, l->t);
    scheduler_addunlock(s, l->t, t_gradient);
  }
  for (struct link *l = c->force; l != NULL; l = l->next)
    scheduler_addunlock(s, t_gradient, l->t);
#else
  for (struct link *l = c->force; l != NULL; l = l->next)
    scheduler_addunlock(s, t_rho, l->t);
#endif

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_hydro(e, c->progeny[k], t_xv, t_rho, t_gradient);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_grav The recv_gpart #task, if it has already been created.
 */
void engine_addtasks_recv_gravity(struct engine *e, struct cell *c,
                                  struct task *t_grav) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Have we reached a level where there are any gravity tasks ? */
  if (t_grav == NULL && c->grav != NULL) {

    /* Create the tasks. */
    t_grav = scheduler_addtask(s, task_type_recv, task_subtype_gpart,
                               6 * c->tag + 4, 0, c, NULL);
  }

  c->recv_grav = t_grav;

  for (struct link *l = c->grav; l != NULL; l = l->next)
    scheduler_addunlock(s, t_grav, l->t);

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_gravity(e, c->progeny[k], t_grav);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add recv tasks for gravity pairs to a hierarchy of cells.
 *
 * @param e The #engine.
 * @param c The foreign #cell.
 * @param t_ti The recv_ti #task, if already been created.
 */
void engine_addtasks_recv_timestep(struct engine *e, struct cell *c,
                                   struct task *t_ti) {

#ifdef WITH_MPI
  struct scheduler *s = &e->sched;

  /* Have we reached a level where there are any self/pair tasks ? */
  if (t_ti == NULL && (c->grav != NULL || c->density != NULL)) {

    t_ti = scheduler_addtask(s, task_type_recv, task_subtype_tend,
                             6 * c->tag + 2, 0, c, NULL);
  }

  c->recv_ti = t_ti;

  for (struct link *l = c->grav; l != NULL; l = l->next)
    scheduler_addunlock(s, l->t, t_ti);

  for (struct link *l = c->force; l != NULL; l = l->next)
    scheduler_addunlock(s, l->t, t_ti);

  /* Recurse? */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        engine_addtasks_recv_timestep(e, c->progeny[k], t_ti);

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
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;
  const int nr_proxies = e->nr_proxies;
  int offset[nr_cells];
  MPI_Request reqs_in[engine_maxproxies];
  MPI_Request reqs_out[engine_maxproxies];
  MPI_Status status;
  const ticks tic = getticks();

  /* Run through the cells and get the size of the ones that will be sent off.
   */
  int count_out = 0;
  for (int k = 0; k < nr_cells; k++) {
    offset[k] = count_out;
    if (cells[k].sendto)
      count_out += (cells[k].pcell_size = cell_getsize(&cells[k]));
  }

  /* Allocate the pcells. */
  struct pcell *pcells = NULL;
  if (posix_memalign((void **)&pcells, SWIFT_CACHE_ALIGNMENT,
                     sizeof(struct pcell) * count_out) != 0)
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
  size_t count_parts_in = 0, count_gparts_in = 0, count_sparts_in = 0;
  for (int k = 0; k < nr_proxies; k++)
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {
      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro)
        count_parts_in += e->proxies[k].cells_in[j]->count;
      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity)
        count_gparts_in += e->proxies[k].cells_in[j]->gcount;
      count_sparts_in += e->proxies[k].cells_in[j]->scount;
    }
  if (count_parts_in > s->size_parts_foreign) {
    if (s->parts_foreign != NULL) free(s->parts_foreign);
    s->size_parts_foreign = 1.1 * count_parts_in;
    if (posix_memalign((void **)&s->parts_foreign, part_align,
                       sizeof(struct part) * s->size_parts_foreign) != 0)
      error("Failed to allocate foreign part data.");
  }
  if (count_gparts_in > s->size_gparts_foreign) {
    if (s->gparts_foreign != NULL) free(s->gparts_foreign);
    s->size_gparts_foreign = 1.1 * count_gparts_in;
    if (posix_memalign((void **)&s->gparts_foreign, gpart_align,
                       sizeof(struct gpart) * s->size_gparts_foreign) != 0)
      error("Failed to allocate foreign gpart data.");
  }
  if (count_sparts_in > s->size_sparts_foreign) {
    if (s->sparts_foreign != NULL) free(s->sparts_foreign);
    s->size_sparts_foreign = 1.1 * count_sparts_in;
    if (posix_memalign((void **)&s->sparts_foreign, spart_align,
                       sizeof(struct spart) * s->size_sparts_foreign) != 0)
      error("Failed to allocate foreign spart data.");
  }

  /* Unpack the cells and link to the particle data. */
  struct part *parts = s->parts_foreign;
  struct gpart *gparts = s->gparts_foreign;
  struct spart *sparts = s->sparts_foreign;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro) {
        cell_link_parts(e->proxies[k].cells_in[j], parts);
        parts = &parts[e->proxies[k].cells_in[j]->count];
      }

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity) {
        cell_link_gparts(e->proxies[k].cells_in[j], gparts);
        gparts = &gparts[e->proxies[k].cells_in[j]->gcount];
      }

      cell_link_sparts(e->proxies[k].cells_in[j], sparts);
      sparts = &sparts[e->proxies[k].cells_in[j]->scount];
    }
  }
  s->nr_parts_foreign = parts - s->parts_foreign;
  s->nr_gparts_foreign = gparts - s->gparts_foreign;
  s->nr_sparts_foreign = sparts - s->sparts_foreign;

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
 * @brief Exchange straying particles with other nodes.
 *
 * @param e The #engine.
 * @param offset_parts The index in the parts array as of which the foreign
 *        parts reside.
 * @param ind_part The foreign #cell ID of each part.
 * @param Npart The number of stray parts, contains the number of parts received
 *        on return.
 * @param offset_gparts The index in the gparts array as of which the foreign
 *        parts reside.
 * @param ind_gpart The foreign #cell ID of each gpart.
 * @param Ngpart The number of stray gparts, contains the number of gparts
 *        received on return.
 * @param offset_sparts The index in the sparts array as of which the foreign
 *        parts reside.
 * @param ind_spart The foreign #cell ID of each spart.
 * @param Nspart The number of stray sparts, contains the number of sparts
 *        received on return.
 *
 * Note that this function does not mess-up the linkage between parts and
 * gparts, i.e. the received particles have correct linkeage.
 */
void engine_exchange_strays(struct engine *e, size_t offset_parts,
                            int *ind_part, size_t *Npart, size_t offset_gparts,
                            int *ind_gpart, size_t *Ngpart,
                            size_t offset_sparts, int *ind_spart,
                            size_t *Nspart) {

#ifdef WITH_MPI

  struct space *s = e->s;
  ticks tic = getticks();

  /* Re-set the proxies. */
  for (int k = 0; k < e->nr_proxies; k++) {
    e->proxies[k].nr_parts_out = 0;
    e->proxies[k].nr_gparts_out = 0;
    e->proxies[k].nr_sparts_out = 0;
  }

  /* Put the parts into the corresponding proxies. */
  for (size_t k = 0; k < *Npart; k++) {
    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_part[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->parts[offset_parts + k].id,
          s->parts[offset_parts + k].x[0], s->parts[offset_parts + k].x[1],
          s->parts[offset_parts + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the part. */
    if (s->parts[offset_parts + k].gpart != NULL) {
      s->parts[offset_parts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_parts_out;
    }

    /* Load the part and xpart into the proxy. */
    proxy_parts_load(&e->proxies[pid], &s->parts[offset_parts + k],
                     &s->xparts[offset_parts + k], 1);
  }

  /* Put the sparts into the corresponding proxies. */
  for (size_t k = 0; k < *Nspart; k++) {
    const int node_id = e->s->cells_top[ind_spart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0)
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->sparts[offset_sparts + k].id,
          s->sparts[offset_sparts + k].x[0], s->sparts[offset_sparts + k].x[1],
          s->sparts[offset_sparts + k].x[2]);

    /* Re-link the associated gpart with the buffer offset of the spart. */
    if (s->sparts[offset_sparts + k].gpart != NULL) {
      s->sparts[offset_sparts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_sparts_out;
    }

    /* Load the spart into the proxy */
    proxy_sparts_load(&e->proxies[pid], &s->sparts[offset_sparts + k], 1);
  }

  /* Put the gparts into the corresponding proxies. */
  for (size_t k = 0; k < *Ngpart; k++) {
    const int node_id = e->s->cells_top[ind_gpart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0)
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lli, x=[%e,%e,%e].",
          node_id, s->gparts[offset_gparts + k].id_or_neg_offset,
          s->gparts[offset_gparts + k].x[0], s->gparts[offset_gparts + k].x[1],
          s->gparts[offset_gparts + k].x[2]);

    /* Load the gpart into the proxy */
    proxy_gparts_load(&e->proxies[pid], &s->gparts[offset_gparts + k], 1);
  }

  /* Launch the proxies. */
  MPI_Request reqs_in[4 * engine_maxproxies];
  MPI_Request reqs_out[4 * engine_maxproxies];
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
  int count_parts_in = 0;
  int count_gparts_in = 0;
  int count_sparts_in = 0;
  for (int k = 0; k < e->nr_proxies; k++) {
    count_parts_in += e->proxies[k].nr_parts_in;
    count_gparts_in += e->proxies[k].nr_gparts_in;
    count_sparts_in += e->proxies[k].nr_sparts_in;
  }
  if (e->verbose) {
    message("sent out %zu/%zu/%zu parts/gparts/sparts, got %i/%i/%i back.",
            *Npart, *Ngpart, *Nspart, count_parts_in, count_gparts_in,
            count_sparts_in);
  }

  /* Reallocate the particle arrays if necessary */
  if (offset_parts + count_parts_in > s->size_parts) {
    message("re-allocating parts array.");
    s->size_parts = (offset_parts + count_parts_in) * engine_parts_size_grow;
    struct part *parts_new = NULL;
    struct xpart *xparts_new = NULL;
    if (posix_memalign((void **)&parts_new, part_align,
                       sizeof(struct part) * s->size_parts) != 0 ||
        posix_memalign((void **)&xparts_new, xpart_align,
                       sizeof(struct xpart) * s->size_parts) != 0)
      error("Failed to allocate new part data.");
    memcpy(parts_new, s->parts, sizeof(struct part) * offset_parts);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * offset_parts);
    free(s->parts);
    free(s->xparts);
    s->parts = parts_new;
    s->xparts = xparts_new;
    for (size_t k = 0; k < offset_parts; k++) {
      if (s->parts[k].gpart != NULL) {
        s->parts[k].gpart->id_or_neg_offset = -k;
      }
    }
  }
  if (offset_sparts + count_sparts_in > s->size_sparts) {
    message("re-allocating sparts array.");
    s->size_sparts = (offset_sparts + count_sparts_in) * engine_parts_size_grow;
    struct spart *sparts_new = NULL;
    if (posix_memalign((void **)&sparts_new, spart_align,
                       sizeof(struct spart) * s->size_sparts) != 0)
      error("Failed to allocate new spart data.");
    memcpy(sparts_new, s->sparts, sizeof(struct spart) * offset_sparts);
    free(s->sparts);
    s->sparts = sparts_new;
    for (size_t k = 0; k < offset_sparts; k++) {
      if (s->sparts[k].gpart != NULL) {
        s->sparts[k].gpart->id_or_neg_offset = -k;
      }
    }
  }
  if (offset_gparts + count_gparts_in > s->size_gparts) {
    message("re-allocating gparts array.");
    s->size_gparts = (offset_gparts + count_gparts_in) * engine_parts_size_grow;
    struct gpart *gparts_new = NULL;
    if (posix_memalign((void **)&gparts_new, gpart_align,
                       sizeof(struct gpart) * s->size_gparts) != 0)
      error("Failed to allocate new gpart data.");
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * offset_gparts);
    free(s->gparts);
    s->gparts = gparts_new;

    for (size_t k = 0; k < offset_gparts; k++) {
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_star) {
        s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }
    }
  }

  /* Collect the requests for the particle data from the proxies. */
  int nr_in = 0, nr_out = 0;
  for (int k = 0; k < e->nr_proxies; k++) {
    if (e->proxies[k].nr_parts_in > 0) {
      reqs_in[4 * k] = e->proxies[k].req_parts_in;
      reqs_in[4 * k + 1] = e->proxies[k].req_xparts_in;
      nr_in += 2;
    } else {
      reqs_in[4 * k] = reqs_in[4 * k + 1] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_gparts_in > 0) {
      reqs_in[4 * k + 2] = e->proxies[k].req_gparts_in;
      nr_in += 1;
    } else {
      reqs_in[4 * k + 2] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sparts_in > 0) {
      reqs_in[4 * k + 3] = e->proxies[k].req_sparts_in;
      nr_in += 1;
    } else {
      reqs_in[4 * k + 3] = MPI_REQUEST_NULL;
    }

    if (e->proxies[k].nr_parts_out > 0) {
      reqs_out[4 * k] = e->proxies[k].req_parts_out;
      reqs_out[4 * k + 1] = e->proxies[k].req_xparts_out;
      nr_out += 2;
    } else {
      reqs_out[4 * k] = reqs_out[4 * k + 1] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_gparts_out > 0) {
      reqs_out[4 * k + 2] = e->proxies[k].req_gparts_out;
      nr_out += 1;
    } else {
      reqs_out[4 * k + 2] = MPI_REQUEST_NULL;
    }
    if (e->proxies[k].nr_sparts_out > 0) {
      reqs_out[4 * k + 3] = e->proxies[k].req_sparts_out;
      nr_out += 1;
    } else {
      reqs_out[4 * k + 3] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for each part array to come in and collect the new
     parts from the proxies. */
  int count_parts = 0, count_gparts = 0, count_sparts = 0;
  for (int k = 0; k < nr_in; k++) {
    int err, pid;
    if ((err = MPI_Waitany(4 * e->nr_proxies, reqs_in, &pid,
                           MPI_STATUS_IGNORE)) != MPI_SUCCESS) {
      char buff[MPI_MAX_ERROR_STRING];
      int res;
      MPI_Error_string(err, buff, &res);
      error("MPI_Waitany failed (%s).", buff);
    }
    if (pid == MPI_UNDEFINED) break;
    // message( "request from proxy %i has arrived." , pid / 4 );
    pid = 4 * (pid / 4);

    /* If all the requests for a given proxy have arrived... */
    if (reqs_in[pid + 0] == MPI_REQUEST_NULL &&
        reqs_in[pid + 1] == MPI_REQUEST_NULL &&
        reqs_in[pid + 2] == MPI_REQUEST_NULL &&
        reqs_in[pid + 3] == MPI_REQUEST_NULL) {
      /* Copy the particle data to the part/xpart/gpart arrays. */
      struct proxy *prox = &e->proxies[pid / 4];
      memcpy(&s->parts[offset_parts + count_parts], prox->parts_in,
             sizeof(struct part) * prox->nr_parts_in);
      memcpy(&s->xparts[offset_parts + count_parts], prox->xparts_in,
             sizeof(struct xpart) * prox->nr_parts_in);
      memcpy(&s->gparts[offset_gparts + count_gparts], prox->gparts_in,
             sizeof(struct gpart) * prox->nr_gparts_in);
      memcpy(&s->sparts[offset_sparts + count_sparts], prox->sparts_in,
             sizeof(struct spart) * prox->nr_sparts_in);
      /* for (int k = offset; k < offset + count; k++)
         message(
            "received particle %lli, x=[%.3e %.3e %.3e], h=%.3e, from node %i.",
            s->parts[k].id, s->parts[k].x[0], s->parts[k].x[1],
            s->parts[k].x[2], s->parts[k].h, p->nodeID); */

      /* Re-link the gparts. */
      for (int kk = 0; kk < prox->nr_gparts_in; kk++) {
        struct gpart *gp = &s->gparts[offset_gparts + count_gparts + kk];

        if (gp->type == swift_type_gas) {
          struct part *p =
              &s->parts[offset_parts + count_parts - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->parts - p;
          p->gpart = gp;
        } else if (gp->type == swift_type_star) {
          struct spart *sp =
              &s->sparts[offset_sparts + count_sparts - gp->id_or_neg_offset];
          gp->id_or_neg_offset = s->sparts - sp;
          sp->gpart = gp;
        }
      }

      /* Advance the counters. */
      count_parts += prox->nr_parts_in;
      count_gparts += prox->nr_gparts_in;
      count_sparts += prox->nr_sparts_in;
    }
  }

  /* Wait for all the sends to have finished too. */
  if (nr_out > 0)
    if (MPI_Waitall(4 * e->nr_proxies, reqs_out, MPI_STATUSES_IGNORE) !=
        MPI_SUCCESS)
      error("MPI_Waitall on sends failed.");

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* Return the number of harvested parts. */
  *Npart = count_parts;
  *Ngpart = count_gparts;
  *Nspart = count_sparts;

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchanges the top-level multipoles between all the nodes
 * such that every node has a multipole for each top-level cell.
 *
 * @param e The #engine.
 */
void engine_exchange_top_multipoles(struct engine *e) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < e->s->nr_cells; ++i) {
    const struct gravity_tensors *m = &e->s->multipoles_top[i];
    if (e->s->cells_top[i].nodeID == engine_rank) {
      if (m->m_pole.M_000 > 0.) {
        if (m->CoM[0] < 0. || m->CoM[0] > e->s->dim[0])
          error("Invalid multipole position in X");
        if (m->CoM[1] < 0. || m->CoM[1] > e->s->dim[1])
          error("Invalid multipole position in Y");
        if (m->CoM[2] < 0. || m->CoM[2] > e->s->dim[2])
          error("Invalid multipole position in Z");
      }
    } else {
      if (m->m_pole.M_000 != 0.) error("Non-zero mass for foreign m-pole");
      if (m->CoM[0] != 0.) error("Non-zero position in X for foreign m-pole");
      if (m->CoM[1] != 0.) error("Non-zero position in Y for foreign m-pole");
      if (m->CoM[2] != 0.) error("Non-zero position in Z for foreign m-pole");
      if (m->m_pole.num_gpart != 0)
        error("Non-zero gpart count in foreign m-pole");
    }
  }
#endif

  /* Each node (space) has constructed its own top-level multipoles.
   * We now need to make sure every other node has a copy of everything.
   *
   * WARNING: Adult stuff ahead: don't do this at home!
   *
   * Since all nodes have their top-level multi-poles computed
   * and all foreign ones set to 0 (all bytes), we can gather all the m-poles
   * by doing a bit-wise OR reduction across all the nodes directly in
   * place inside the multi-poles_top array.
   * This only works if the foreign m-poles on every nodes are zeroed and no
   * multi-pole is present on more than one node (two things guaranteed by the
   * domain decomposition).
   */
  MPI_Allreduce(MPI_IN_PLACE, e->s->multipoles_top,
                e->s->nr_cells * sizeof(struct gravity_tensors), MPI_BYTE,
                MPI_BOR, MPI_COMM_WORLD);

#ifdef SWIFT_DEBUG_CHECKS
  long long counter = 0;

  /* Let's check that what we received makes sense */
  for (int i = 0; i < e->s->nr_cells; ++i) {
    const struct gravity_tensors *m = &e->s->multipoles_top[i];
    counter += m->m_pole.num_gpart;
    if (m->m_pole.M_000 > 0.) {
      if (m->CoM[0] < 0. || m->CoM[0] > e->s->dim[0])
        error("Invalid multipole position in X");
      if (m->CoM[1] < 0. || m->CoM[1] > e->s->dim[1])
        error("Invalid multipole position in Y");
      if (m->CoM[2] < 0. || m->CoM[2] > e->s->dim[2])
        error("Invalid multipole position in Z");
    }
  }
  if (counter != e->total_nr_gparts)
    error("Total particles in multipoles inconsistent with engine");
#endif

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

void engine_exchange_proxy_multipoles(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();

  /* Start by counting the number of cells to send and receive */
  int count_send = 0;
  int count_recv = 0;
  int count_send_requests = 0;
  int count_recv_requests = 0;

  /* Loop over the proxies. */
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    /* Now collect the number of requests associated */
    count_recv_requests += p->nr_cells_in;
    count_send_requests += p->nr_cells_out;

    /* And the actual number of things we are going to ship */
    for (int k = 0; k < p->nr_cells_in; k++)
      count_recv += p->cells_in[k]->pcell_size;

    for (int k = 0; k < p->nr_cells_out; k++)
      count_send += p->cells_out[k]->pcell_size;
  }

  /* Allocate the buffers for the packed data */
  struct gravity_tensors *buffer_send = NULL;
  if (posix_memalign((void **)&buffer_send, SWIFT_CACHE_ALIGNMENT,
                     count_send * sizeof(struct gravity_tensors)) != 0)
    error("Unable to allocate memory for multipole transactions");

  struct gravity_tensors *buffer_recv = NULL;
  if (posix_memalign((void **)&buffer_recv, SWIFT_CACHE_ALIGNMENT,
                     count_recv * sizeof(struct gravity_tensors)) != 0)
    error("Unable to allocate memory for multipole transactions");

  /* Also allocate the MPI requests */
  const int count_requests = count_send_requests + count_recv_requests;
  MPI_Request *requests =
      (MPI_Request *)malloc(sizeof(MPI_Request) * count_requests);
  if (requests == NULL) error("Unable to allocate memory for MPI requests");

  int this_request = 0;
  int this_recv = 0;
  int this_send = 0;

  /* Loop over the proxies to issue the receives. */
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    for (int k = 0; k < p->nr_cells_in; k++) {

      const int num_elements = p->cells_in[k]->pcell_size;

      /* Receive everything */
      MPI_Irecv(&buffer_recv[this_recv],
                num_elements * sizeof(struct gravity_tensors), MPI_BYTE,
                p->cells_in[k]->nodeID, p->cells_in[k]->tag, MPI_COMM_WORLD,
                &requests[this_request]);

      /* Move to the next slot in the buffers */
      this_recv += num_elements;
      this_request++;
    }

    /* Loop over the proxies to issue the sends. */
    for (int k = 0; k < p->nr_cells_out; k++) {

      /* Number of multipoles in this cell hierarchy */
      const int num_elements = p->cells_out[k]->pcell_size;

      /* Let's pack everything recursively */
      cell_pack_multipoles(p->cells_out[k], &buffer_send[this_send]);

      /* Send everything (note the use of cells_in[0] to get the correct node
       * ID. */
      MPI_Isend(&buffer_send[this_send],
                num_elements * sizeof(struct gravity_tensors), MPI_BYTE,
                p->cells_in[0]->nodeID, p->cells_out[k]->tag, MPI_COMM_WORLD,
                &requests[this_request]);

      /* Move to the next slot in the buffers */
      this_send += num_elements;
      this_request++;
    }
  }

  /* Wait for all the requests to arrive home */
  MPI_Status *stats = (MPI_Status *)malloc(count_requests * sizeof(MPI_Status));
  int res;
  if ((res = MPI_Waitall(count_requests, requests, stats)) != MPI_SUCCESS) {
    for (int k = 0; k < count_requests; ++k) {
      char buff[MPI_MAX_ERROR_STRING];
      MPI_Error_string(stats[k].MPI_ERROR, buff, &res);
      message("request from source %i, tag %i has error '%s'.",
              stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
    }
    error("Failed during waitall for multipole data.");
  }

  /* Let's now unpack the multipoles at the right place */
  this_recv = 0;
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    for (int k = 0; k < p->nr_cells_in; k++) {

      const int num_elements = p->cells_in[k]->pcell_size;

#ifdef SWIFT_DEBUG_CHECKS

      /* Check that the first element (top-level cell's multipole) matches what
       * we received */
      if (p->cells_in[k]->multipole->m_pole.num_gpart !=
          buffer_recv[this_recv].m_pole.num_gpart)
        error("Current: M_000=%e num_gpart=%lld\n New: M_000=%e num_gpart=%lld",
              p->cells_in[k]->multipole->m_pole.M_000,
              p->cells_in[k]->multipole->m_pole.num_gpart,
              buffer_recv[this_recv].m_pole.M_000,
              buffer_recv[this_recv].m_pole.num_gpart);
#endif

      /* Unpack recursively */
      cell_unpack_multipoles(p->cells_in[k], &buffer_recv[this_recv]);

      /* Move to the next slot in the buffers */
      this_recv += num_elements;
    }
  }

  /* Free everything */
  free(stats);
  free(buffer_send);
  free(buffer_recv);
  free(requests);

  /* How much time did this take? */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * and long-range gravity interactions.
 *
 * - All top-cells get a self task.
 * - All pairs within range according to the multipole acceptance
 *   criterion get a pair task.
 */
void engine_make_self_gravity_tasks_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct engine *e = ((struct engine **)extra_data)[0];
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->cells_top;
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_distance = e->mesh->r_cut_max;

  /* Compute how many cells away we need to walk */
  const double distance = 2.5 * cells[0].width[0] / theta_crit;
  int delta = (int)(distance / cells[0].width[0]) + 1;
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
    }
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->gcount == 0) continue;

    /* Is that cell local ? */
    if (ci->nodeID != nodeID) continue;

    /* If the cells is local build a self-interaction */
    scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci, NULL);

    /* Recover the multipole information */
    const struct gravity_tensors *const multi_i = ci->multipole;
    const double CoM_i[3] = {multi_i->CoM[0], multi_i->CoM[1], multi_i->CoM[2]};

#ifdef SWIFT_DEBUG_CHECKS
    if (cell_getid(cdim, i, j, k) != cid)
      error("Incorrect calculation of indices (i,j,k)=(%d,%d,%d) cid=%d", i, j,
            k, cid);

    if (multi_i->r_max != multi_i->r_max_rebuild)
      error(
          "Multipole size not equal ot it's size after rebuild. But we just "
          "rebuilt...");
#endif

    /* Loop over every other cell within (Manhattan) range delta */
    for (int x = -delta_m; x <= delta_p; x++) {
      int ii = i + x;
      if (ii >= cdim[0])
        ii -= cdim[0];
      else if (ii < 0)
        ii += cdim[0];
      for (int y = -delta_m; y <= delta_p; y++) {
        int jj = j + y;
        if (jj >= cdim[1])
          jj -= cdim[1];
        else if (jj < 0)
          jj += cdim[1];
        for (int z = -delta_m; z <= delta_p; z++) {
          int kk = k + z;
          if (kk >= cdim[2])
            kk -= cdim[2];
          else if (kk < 0)
            kk += cdim[2];

          /* Get the cell */
          const int cjd = cell_getid(cdim, ii, jj, kk);
          struct cell *cj = &cells[cjd];

#ifdef SWIFT_DEBUG_CHECKS
          const int iii = cjd / (cdim[1] * cdim[2]);
          const int jjj = (cjd / cdim[2]) % cdim[1];
          const int kkk = cjd % cdim[2];

          if (ii != iii || jj != jjj || kk != kkk)
            error(
                "Incorrect calculation of indices (iii,jjj,kkk)=(%d,%d,%d) "
                "cjd=%d",
                iii, jjj, kkk, cjd);
#endif

          /* Avoid duplicates of local pairs*/
          if (cid <= cjd && cj->nodeID == nodeID) continue;

          /* Skip cells without gravity particles */
          if (cj->gcount == 0) continue;

          /* Recover the multipole information */
          const struct gravity_tensors *const multi_j = cj->multipole;

          /* Get the distance between the CoMs */
          double dx = CoM_i[0] - multi_j->CoM[0];
          double dy = CoM_i[1] - multi_j->CoM[1];
          double dz = CoM_i[2] - multi_j->CoM[2];

          /* Apply BC */
          if (periodic) {
            dx = nearest(dx, dim[0]);
            dy = nearest(dy, dim[1]);
            dz = nearest(dz, dim[2]);
          }
          const double r2 = dx * dx + dy * dy + dz * dz;

          /* Minimal distance between any pair of particles */
          const double min_radius =
              sqrt(r2) - (multi_i->r_max + multi_j->r_max);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (periodic && min_radius > max_distance) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm(ci, cj, e, s)) {

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
                              ci, cj);
          }
        }
      }
    }
  }
}

/**
 * @brief Constructs the top-level tasks for the short-range gravity
 * interactions (master function).
 *
 * - Create the FFT task and the array of gravity ghosts.
 * - Call the mapper function to create the other tasks.
 *
 * @param e The #engine.
 */
void engine_make_self_gravity_tasks(struct engine *e) {

  struct space *s = e->s;
  struct task **ghosts = NULL;

  /* Create the multipole self and pair tasks. */
  void *extra_data[2] = {e, ghosts};
  threadpool_map(&e->threadpool, engine_make_self_gravity_tasks_mapper, NULL,
                 s->nr_cells, 1, 0, extra_data);
}

/**
 * @brief Constructs the top-level tasks for the external gravity.
 *
 * @param e The #engine.
 */
void engine_make_external_gravity_tasks(struct engine *e) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;

  for (int cid = 0; cid < nr_cells; ++cid) {

    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->gcount == 0) continue;

    /* Is that neighbour local ? */
    if (ci->nodeID != nodeID) continue;

    /* If the cell is local, build a self-interaction */
    scheduler_addtask(sched, task_type_self, task_subtype_external_grav, 0, 0,
                      ci, NULL);
  }
}

/**
 * @brief Constructs the top-level pair tasks for the first hydro loop over
 * neighbours
 *
 * Here we construct all the tasks for all possible neighbouring non-empty
 * local cells in the hierarchy. No dependencies are being added thus far.
 * Additional loop over neighbours can later be added by simply duplicating
 * all the tasks created by this function.
 *
 * @param map_data Offset of first two indices disguised as a pointer.
 * @param num_elements Number of cells to traverse.
 * @param extra_data The #engine.
 */
void engine_make_hydroloop_tasks_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  /* Extract the engine pointer. */
  struct engine *e = (struct engine *)extra_data;

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int *cdim = s->cdim;
  struct cell *cells = s->cells_top;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without hydro particles */
    if (ci->count == 0) continue;

    /* If the cells is local build a self-interaction */
    if (ci->nodeID == nodeID)
      scheduler_addtask(sched, task_type_self, task_subtype_density, 0, 0, ci,
                        NULL);

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
          const int sid = sortlistID[(kk + 1) + 3 * ((jj + 1) + 3 * (ii + 1))];
          scheduler_addtask(sched, task_type_pair, task_subtype_density, sid, 0,
                            ci, cj);
        }
      }
    }
  }
}

/**
 * @brief Counts the tasks associated with one cell and constructs the links
 *
 * For each hydrodynamic and gravity task, construct the links with
 * the corresponding cell.  Similarly, construct the dependencies for
 * all the sorting tasks.
 */
void engine_count_and_link_tasks_mapper(void *map_data, int num_elements,
                                        void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *const sched = &e->sched;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &((struct task *)map_data)[ind];

    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

    /* Link sort tasks to all the higher sort task. */
    if (t_type == task_type_sort) {
      for (struct cell *finger = t->ci->parent; finger != NULL;
           finger = finger->parent)
        if (finger->sorts != NULL) scheduler_addunlock(sched, t, finger->sorts);
    }

    /* Link self tasks to cells. */
    else if (t_type == task_type_self) {
      atomic_inc(&ci->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav, t);
      } else if (t_subtype == task_subtype_external_grav) {
        engine_addlink(e, &ci->grav, t);
      }

      /* Link pair tasks to cells. */
    } else if (t_type == task_type_pair) {
      atomic_inc(&ci->nr_tasks);
      atomic_inc(&cj->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->density, t);
        engine_addlink(e, &cj->density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav, t);
        engine_addlink(e, &cj->grav, t);
      }
#ifdef SWIFT_DEBUG_CHECKS
      else if (t_subtype == task_subtype_external_grav) {
        error("Found a pair/external-gravity task...");
      }
#endif

      /* Link sub-self tasks to cells. */
    } else if (t_type == task_type_sub_self) {
      atomic_inc(&ci->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav, t);
      } else if (t_subtype == task_subtype_external_grav) {
        engine_addlink(e, &ci->grav, t);
      }

      /* Link sub-pair tasks to cells. */
    } else if (t_type == task_type_sub_pair) {
      atomic_inc(&ci->nr_tasks);
      atomic_inc(&cj->nr_tasks);

      if (t_subtype == task_subtype_density) {
        engine_addlink(e, &ci->density, t);
        engine_addlink(e, &cj->density, t);
      } else if (t_subtype == task_subtype_grav) {
        engine_addlink(e, &ci->grav, t);
        engine_addlink(e, &cj->grav, t);
      }
#ifdef SWIFT_DEBUG_CHECKS
      else if (t_subtype == task_subtype_external_grav) {
        error("Found a sub-pair/external-gravity task...");
      }
#endif

      /* Note that we do not need to link the M-M tasks */
      /* since we already did so when splitting the gravity */
      /* tasks. */
    }
  }
}

/**
 * @brief Creates all the task dependencies for the gravity
 *
 * @param e The #engine
 */
void engine_link_gravity_tasks(struct engine *e) {

  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int nr_tasks = sched->nr_tasks;

  for (int k = 0; k < nr_tasks; k++) {

    /* Get a pointer to the task. */
    struct task *t = &sched->tasks[k];

    /* Get the cells we act on */
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const enum task_types t_type = t->type;
    const enum task_subtypes t_subtype = t->subtype;

/* Node ID (if running with MPI) */
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Self-interaction for self-gravity? */
    if (t_type == task_type_self && t_subtype == task_subtype_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local self task");
#endif

      /* drift ---+-> gravity --> grav_down */
      /* init  --/    */
      scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
      scheduler_addunlock(sched, ci->init_grav_out, t);
      scheduler_addunlock(sched, t, ci->grav_down_in);
    }

    /* Self-interaction for external gravity ? */
    if (t_type == task_type_self && t_subtype == task_subtype_external_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local self task");
#endif

      /* drift -----> gravity --> end_force */
      scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
      scheduler_addunlock(sched, t, ci->end_force);
    }

    /* Otherwise, pair interaction? */
    else if (t_type == task_type_pair && t_subtype == task_subtype_grav) {

      if (ci_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
        scheduler_addunlock(sched, ci->init_grav_out, t);
        scheduler_addunlock(sched, t, ci->grav_down_in);
      }
      if (cj_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        if (ci->super_gravity != cj->super_gravity) /* Avoid double unlock */
          scheduler_addunlock(sched, cj->super_gravity->drift_gpart, t);
        scheduler_addunlock(sched, cj->init_grav_out, t);
        scheduler_addunlock(sched, t, cj->grav_down_in);
      }
    }

    /* Otherwise, sub-self interaction? */
    else if (t_type == task_type_sub_self && t_subtype == task_subtype_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local sub-self task");
#endif
      /* drift ---+-> gravity --> grav_down */
      /* init  --/    */
      scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
      scheduler_addunlock(sched, ci->init_grav_out, t);
      scheduler_addunlock(sched, t, ci->grav_down_in);
    }

    /* Sub-self-interaction for external gravity ? */
    else if (t_type == task_type_sub_self &&
             t_subtype == task_subtype_external_grav) {

#ifdef SWIFT_DEBUG_CHECKS
      if (ci_nodeID != nodeID) error("Non-local sub-self task");
#endif

      /* drift -----> gravity --> end_force */
      scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
      scheduler_addunlock(sched, t, ci->end_force);
    }

    /* Otherwise, sub-pair interaction? */
    else if (t_type == task_type_sub_pair && t_subtype == task_subtype_grav) {

      if (ci_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        scheduler_addunlock(sched, ci->super_gravity->drift_gpart, t);
        scheduler_addunlock(sched, ci->init_grav_out, t);
        scheduler_addunlock(sched, t, ci->grav_down_in);
      }
      if (cj_nodeID == nodeID) {

        /* drift ---+-> gravity --> grav_down */
        /* init  --/    */
        if (ci->super_gravity != cj->super_gravity) /* Avoid double unlock */
          scheduler_addunlock(sched, cj->super_gravity->drift_gpart, t);
        scheduler_addunlock(sched, cj->init_grav_out, t);
        scheduler_addunlock(sched, t, cj->grav_down_in);
      }
    }

    /* Otherwise M-M interaction? */
    else if (t_type == task_type_grav_mm) {

      if (ci_nodeID == nodeID) {

        /* init -----> gravity --> grav_down */
        scheduler_addunlock(sched, ci->init_grav_out, t);
        scheduler_addunlock(sched, t, ci->grav_down_in);
      }
      if (cj_nodeID == nodeID) {

        /* init -----> gravity --> grav_down */
        scheduler_addunlock(sched, cj->init_grav_out, t);
        scheduler_addunlock(sched, t, cj->grav_down_in);
      }
    }
  }
}

#ifdef EXTRA_HYDRO_LOOP

/**
 * @brief Creates the dependency network for the hydro tasks of a given cell.
 *
 * @param sched The #scheduler.
 * @param density The density task to link.
 * @param gradient The gradient task to link.
 * @param force The force task to link.
 * @param c The cell.
 * @param with_cooling Do we have a cooling task ?
 */
static inline void engine_make_hydro_loops_dependencies(
    struct scheduler *sched, struct task *density, struct task *gradient,
    struct task *force, struct cell *c, int with_cooling) {

  /* density loop --> ghost --> gradient loop --> extra_ghost */
  /* extra_ghost --> force loop  */
  scheduler_addunlock(sched, density, c->super_hydro->ghost_in);
  scheduler_addunlock(sched, c->super_hydro->ghost_out, gradient);
  scheduler_addunlock(sched, gradient, c->super_hydro->extra_ghost);
  scheduler_addunlock(sched, c->super_hydro->extra_ghost, force);
}

#else

/**
 * @brief Creates the dependency network for the hydro tasks of a given cell.
 *
 * @param sched The #scheduler.
 * @param density The density task to link.
 * @param force The force task to link.
 * @param c The cell.
 * @param with_cooling Are we running with cooling switched on ?
 */
static inline void engine_make_hydro_loops_dependencies(struct scheduler *sched,
                                                        struct task *density,
                                                        struct task *force,
                                                        struct cell *c,
                                                        int with_cooling) {
  /* density loop --> ghost --> force loop */
  scheduler_addunlock(sched, density, c->super_hydro->ghost_in);
  scheduler_addunlock(sched, c->super_hydro->ghost_out, force);
}

#endif
/**
 * @brief Duplicates the first hydro loop and construct all the
 * dependencies for the hydro part
 *
 * This is done by looping over all the previously constructed tasks
 * and adding another task involving the same cells but this time
 * corresponding to the second hydro loop over neighbours.
 * With all the relevant tasks for a given cell available, we construct
 * all the dependencies for that cell.
 */
void engine_make_extra_hydroloop_tasks_mapper(void *map_data, int num_elements,
                                              void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int with_cooling = (e->policy & engine_policy_cooling);

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &((struct task *)map_data)[ind];

    /* Sort tasks depend on the drift of the cell. */
    if (t->type == task_type_sort && t->ci->nodeID == engine_rank) {
      scheduler_addunlock(sched, t->ci->super_hydro->drift_part, t);
    }

    /* Self-interaction? */
    else if (t->type == task_type_self && t->subtype == task_subtype_density) {

      /* Make the self-density tasks depend on the drift only. */
      scheduler_addunlock(sched, t->ci->super_hydro->drift_part, t);

#ifdef EXTRA_HYDRO_LOOP
      /* Start by constructing the task for the second  and third hydro loop. */
      struct task *t2 = scheduler_addtask(
          sched, task_type_self, task_subtype_gradient, 0, 0, t->ci, NULL);
      struct task *t3 = scheduler_addtask(
          sched, task_type_self, task_subtype_force, 0, 0, t->ci, NULL);

      /* Add the link between the new loops and the cell */
      engine_addlink(e, &t->ci->gradient, t2);
      engine_addlink(e, &t->ci->force, t3);

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->ci,
                                           with_cooling);
      scheduler_addunlock(sched, t3, t->ci->super->end_force);
#else

      /* Start by constructing the task for the second hydro loop */
      struct task *t2 = scheduler_addtask(
          sched, task_type_self, task_subtype_force, 0, 0, t->ci, NULL);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->force, t2);

      /* Now, build all the dependencies for the hydro */
      engine_make_hydro_loops_dependencies(sched, t, t2, t->ci, with_cooling);
      scheduler_addunlock(sched, t2, t->ci->super->end_force);
#endif
    }

    /* Otherwise, pair interaction? */
    else if (t->type == task_type_pair && t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift and the sorts. */
      if (t->ci->nodeID == engine_rank)
        scheduler_addunlock(sched, t->ci->super_hydro->drift_part, t);
      scheduler_addunlock(sched, t->ci->super_hydro->sorts, t);
      if (t->ci->super_hydro != t->cj->super_hydro) {
        if (t->cj->nodeID == engine_rank)
          scheduler_addunlock(sched, t->cj->super_hydro->drift_part, t);
        scheduler_addunlock(sched, t->cj->super_hydro->sorts, t);
      }

#ifdef EXTRA_HYDRO_LOOP
      /* Start by constructing the task for the second and third hydro loop */
      struct task *t2 = scheduler_addtask(
          sched, task_type_pair, task_subtype_gradient, 0, 0, t->ci, t->cj);
      struct task *t3 = scheduler_addtask(
          sched, task_type_pair, task_subtype_force, 0, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->gradient, t2);
      engine_addlink(e, &t->cj->gradient, t2);
      engine_addlink(e, &t->ci->force, t3);
      engine_addlink(e, &t->cj->force, t3);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->ci,
                                             with_cooling);
        scheduler_addunlock(sched, t3, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super_hydro != t->cj->super_hydro)
          engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->cj,
                                               with_cooling);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t3, t->cj->super->end_force);
      }

#else

      /* Start by constructing the task for the second hydro loop */
      struct task *t2 = scheduler_addtask(
          sched, task_type_pair, task_subtype_force, 0, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->force, t2);
      engine_addlink(e, &t->cj->force, t2);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t->ci, with_cooling);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super_hydro != t->cj->super_hydro)
          engine_make_hydro_loops_dependencies(sched, t, t2, t->cj,
                                               with_cooling);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t2, t->cj->super->end_force);
      }

#endif

    }

    /* Otherwise, sub-self interaction? */
    else if (t->type == task_type_sub_self &&
             t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift and sorts. */
      scheduler_addunlock(sched, t->ci->super_hydro->drift_part, t);
      scheduler_addunlock(sched, t->ci->super_hydro->sorts, t);

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_sub_self, task_subtype_gradient,
                            t->flags, 0, t->ci, t->cj);
      struct task *t3 =
          scheduler_addtask(sched, task_type_sub_self, task_subtype_force,
                            t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->gradient, t2);
      engine_addlink(e, &t->ci->force, t3);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->ci,
                                             with_cooling);
        scheduler_addunlock(sched, t3, t->ci->super->end_force);
      }

#else
      /* Start by constructing the task for the second hydro loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_sub_self, task_subtype_force,
                            t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and the cell */
      engine_addlink(e, &t->ci->force, t2);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t->ci, with_cooling);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
#endif
    }

    /* Otherwise, sub-pair interaction? */
    else if (t->type == task_type_sub_pair &&
             t->subtype == task_subtype_density) {

      /* Make all density tasks depend on the drift. */
      if (t->ci->nodeID == engine_rank)
        scheduler_addunlock(sched, t->ci->super_hydro->drift_part, t);
      scheduler_addunlock(sched, t->ci->super_hydro->sorts, t);
      if (t->ci->super_hydro != t->cj->super_hydro) {
        if (t->cj->nodeID == engine_rank)
          scheduler_addunlock(sched, t->cj->super_hydro->drift_part, t);
        scheduler_addunlock(sched, t->cj->super_hydro->sorts, t);
      }

#ifdef EXTRA_HYDRO_LOOP

      /* Start by constructing the task for the second and third hydro loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_sub_pair, task_subtype_gradient,
                            t->flags, 0, t->ci, t->cj);
      struct task *t3 =
          scheduler_addtask(sched, task_type_sub_pair, task_subtype_force,
                            t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->gradient, t2);
      engine_addlink(e, &t->cj->gradient, t2);
      engine_addlink(e, &t->ci->force, t3);
      engine_addlink(e, &t->cj->force, t3);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->ci,
                                             with_cooling);
        scheduler_addunlock(sched, t3, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super_hydro != t->cj->super_hydro)
          engine_make_hydro_loops_dependencies(sched, t, t2, t3, t->cj,
                                               with_cooling);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t3, t->cj->super->end_force);
      }

#else
      /* Start by constructing the task for the second hydro loop */
      struct task *t2 =
          scheduler_addtask(sched, task_type_sub_pair, task_subtype_force,
                            t->flags, 0, t->ci, t->cj);

      /* Add the link between the new loop and both cells */
      engine_addlink(e, &t->ci->force, t2);
      engine_addlink(e, &t->cj->force, t2);

      /* Now, build all the dependencies for the hydro for the cells */
      /* that are local and are not descendant of the same super_hydro-cells */
      if (t->ci->nodeID == nodeID) {
        engine_make_hydro_loops_dependencies(sched, t, t2, t->ci, with_cooling);
        scheduler_addunlock(sched, t2, t->ci->super->end_force);
      }
      if (t->cj->nodeID == nodeID) {
        if (t->ci->super_hydro != t->cj->super_hydro)
          engine_make_hydro_loops_dependencies(sched, t, t2, t->cj,
                                               with_cooling);
        if (t->ci->super != t->cj->super)
          scheduler_addunlock(sched, t2, t->cj->super->end_force);
      }
#endif
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
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;
  const ticks tic = getticks();

  /* Re-set the scheduler. */
  scheduler_reset(sched, engine_estimate_nr_tasks(e));

  ticks tic2 = getticks();

  /* Construct the firt hydro loop over neighbours */
  if (e->policy & engine_policy_hydro)
    threadpool_map(&e->threadpool, engine_make_hydroloop_tasks_mapper, NULL,
                   s->nr_cells, 1, 0, e);

  if (e->verbose)
    message("Making hydro tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Add the self gravity tasks. */
  if (e->policy & engine_policy_self_gravity) engine_make_self_gravity_tasks(e);

  if (e->verbose)
    message("Making gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Add the external gravity tasks. */
  if (e->policy & engine_policy_external_gravity)
    engine_make_external_gravity_tasks(e);

  if (e->sched.nr_tasks == 0 && (s->nr_gparts > 0 || s->nr_parts > 0))
    error("We have particles but no hydro or gravity tasks were created.");

  /* Free the old list of cell-task links. */
  if (e->links != NULL) free(e->links);
  e->size_links = 0;

/* The maximum number of links is the
 * number of cells (s->tot_cells) times the number of neighbours (26) times
 * the number of interaction types, so 26 * 2 (density, force) pairs
 * and 2 (density, force) self.
 */
#ifdef EXTRA_HYDRO_LOOP
  const size_t hydro_tasks_per_cell = 27 * 3;
#else
  const size_t hydro_tasks_per_cell = 27 * 2;
#endif
  const size_t self_grav_tasks_per_cell = 125;
  const size_t ext_grav_tasks_per_cell = 1;

  if (e->policy & engine_policy_hydro)
    e->size_links += s->tot_cells * hydro_tasks_per_cell;
  if (e->policy & engine_policy_external_gravity)
    e->size_links += s->tot_cells * ext_grav_tasks_per_cell;
  if (e->policy & engine_policy_self_gravity)
    e->size_links += s->tot_cells * self_grav_tasks_per_cell;

  /* Allocate the new link list */
  if ((e->links = (struct link *)malloc(sizeof(struct link) * e->size_links)) ==
      NULL)
    error("Failed to allocate cell-task links.");
  e->nr_links = 0;

  tic2 = getticks();

  /* Split the tasks. */
  scheduler_splittasks(sched);

  if (e->verbose)
    message("Splitting tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that we are not left with invalid tasks */
  for (int i = 0; i < e->sched.nr_tasks; ++i) {
    const struct task *t = &e->sched.tasks[i];
    if (t->ci == NULL && t->cj != NULL && !t->skip) error("Invalid task");
  }
#endif

  tic2 = getticks();

  /* Count the number of tasks associated with each cell and
     store the density tasks in each cell, and make each sort
     depend on the sorts of its progeny. */
  threadpool_map(&e->threadpool, engine_count_and_link_tasks_mapper,
                 sched->tasks, sched->nr_tasks, sizeof(struct task), 0, e);

  if (e->verbose)
    message("Counting and linking tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Now that the self/pair tasks are at the right level, set the super
   * pointers. */
  threadpool_map(&e->threadpool, cell_set_super_mapper, cells, nr_cells,
                 sizeof(struct cell), 0, e);

  if (e->verbose)
    message("Setting super-pointers took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Append hierarchical tasks to each cell. */
  threadpool_map(&e->threadpool, engine_make_hierarchical_tasks_mapper, cells,
                 nr_cells, sizeof(struct cell), 0, e);

  tic2 = getticks();

  /* Run through the tasks and make force tasks for each density task.
     Each force task depends on the cell ghosts and unlocks the kick task
     of its super-cell. */
  if (e->policy & engine_policy_hydro)
    threadpool_map(&e->threadpool, engine_make_extra_hydroloop_tasks_mapper,
                   sched->tasks, sched->nr_tasks, sizeof(struct task), 0, e);

  if (e->verbose)
    message("Making extra hydroloop tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Add the dependencies for the gravity stuff */
  if (e->policy & (engine_policy_self_gravity | engine_policy_external_gravity))
    engine_link_gravity_tasks(e);

  if (e->verbose)
    message("Linking gravity tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

#ifdef WITH_MPI

  /* Add the communication tasks if MPI is being used. */
  if (e->policy & engine_policy_mpi) {

    /* Loop over the proxies. */
    for (int pid = 0; pid < e->nr_proxies; pid++) {

      /* Get a handle on the proxy. */
      struct proxy *p = &e->proxies[pid];

      for (int k = 0; k < p->nr_cells_in; k++)
        engine_addtasks_recv_timestep(e, p->cells_in[k], NULL);

      for (int k = 0; k < p->nr_cells_out; k++)
        engine_addtasks_send_timestep(e, p->cells_out[k], p->cells_in[0], NULL);

      /* Loop through the proxy's incoming cells and add the
         recv tasks for the cells in the proxy that have a hydro connection. */
      if (e->policy & engine_policy_hydro)
        for (int k = 0; k < p->nr_cells_in; k++)
          if (p->cells_in_type[k] & proxy_cell_type_hydro)
            engine_addtasks_recv_hydro(e, p->cells_in[k], NULL, NULL, NULL);

      /* Loop through the proxy's incoming cells and add the
         recv tasks for the cells in the proxy that have a gravity connection.
         */
      if (e->policy & engine_policy_self_gravity)
        for (int k = 0; k < p->nr_cells_in; k++)
          if (p->cells_in_type[k] & proxy_cell_type_gravity)
            engine_addtasks_recv_gravity(e, p->cells_in[k], NULL);

      /* Loop through the proxy's outgoing cells and add the
         send tasks for the cells in the proxy that have a hydro connection. */
      if (e->policy & engine_policy_hydro)
        for (int k = 0; k < p->nr_cells_out; k++)
          if (p->cells_out_type[k] & proxy_cell_type_hydro)
            engine_addtasks_send_hydro(e, p->cells_out[k], p->cells_in[0], NULL,
                                       NULL, NULL);

      /* Loop through the proxy's outgoing cells and add the
         send tasks for the cells in the proxy that have a gravity connection.
         */
      if (e->policy & engine_policy_self_gravity)
        for (int k = 0; k < p->nr_cells_out; k++)
          if (p->cells_out_type[k] & proxy_cell_type_gravity)
            engine_addtasks_send_gravity(e, p->cells_out[k], p->cells_in[0],
                                         NULL);
    }
  }
#endif

  tic2 = getticks();

  /* Set the unlocks per task. */
  scheduler_set_unlocks(sched);

  if (e->verbose)
    message("Setting unlocks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  tic2 = getticks();

  /* Rank the tasks. */
  scheduler_ranktasks(sched);

  if (e->verbose)
    message("Ranking the tasks took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Weight the tasks. */
  scheduler_reweight(sched, e->verbose);

  /* Set the tasks age. */
  e->tasks_age = 0;

  if (e->verbose)
    message("took %.3f %s (including reweight).",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

/**
 * @brief Mark tasks to be un-skipped and set the sort flags accordingly.
 *        Threadpool mapper function.
 *
 * @param map_data pointer to the tasks
 * @param num_elements number of tasks
 * @param extra_data pointer to int that will define if a rebuild is needed.
 */
void engine_marktasks_mapper(void *map_data, int num_elements,
                             void *extra_data) {
  /* Unpack the arguments. */
  struct task *tasks = (struct task *)map_data;
  size_t *rebuild_space = &((size_t *)extra_data)[1];
  struct scheduler *s = (struct scheduler *)(((size_t *)extra_data)[2]);
  struct engine *e = (struct engine *)((size_t *)extra_data)[0];

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Single-cell task? */
    if (t->type == task_type_self || t->type == task_type_sub_self) {

      /* Local pointer. */
      struct cell *ci = t->ci;

      if (ci->nodeID != engine_rank) error("Non-local self task found");

      /* Activate the hydro drift */
      if (t->type == task_type_self && t->subtype == task_subtype_density) {
        if (cell_is_active_hydro(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_drift_part(ci, s);
        }
      }

      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_self &&
               t->subtype == task_subtype_density) {
        if (cell_is_active_hydro(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_hydro_tasks(ci, NULL, s);
        }
      }

      else if (t->type == task_type_self && t->subtype == task_subtype_force) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t->type == task_type_sub_self &&
               t->subtype == task_subtype_force) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

#ifdef EXTRA_HYDRO_LOOP
      else if (t->type == task_type_self &&
               t->subtype == task_subtype_gradient) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }

      else if (t->type == task_type_sub_self &&
               t->subtype == task_subtype_gradient) {
        if (cell_is_active_hydro(ci, e)) scheduler_activate(s, t);
      }
#endif

      /* Activate the gravity drift */
      else if (t->type == task_type_self && t->subtype == task_subtype_grav) {
        if (cell_is_active_gravity(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_subcell_grav_tasks(t->ci, NULL, s);
        }
      }

      /* Activate the gravity drift */
      else if (t->type == task_type_self &&
               t->subtype == task_subtype_external_grav) {
        if (cell_is_active_gravity(ci, e)) {
          scheduler_activate(s, t);
          cell_activate_drift_gpart(t->ci, s);
        }
      }

#ifdef SWIFT_DEBUG_CHECKS
      else {
        error("Invalid task type / sub-type encountered");
      }
#endif
    }

    /* Pair? */
    else if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      /* Local pointers. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;
      const int ci_active_hydro = cell_is_active_hydro(ci, e);
      const int cj_active_hydro = cell_is_active_hydro(cj, e);
      const int ci_active_gravity = cell_is_active_gravity(ci, e);
      const int cj_active_gravity = cell_is_active_gravity(cj, e);

      /* Only activate tasks that involve a local active cell. */
      if ((t->subtype == task_subtype_density ||
           t->subtype == task_subtype_gradient ||
           t->subtype == task_subtype_force) &&
          ((ci_active_hydro && ci->nodeID == engine_rank) ||
           (cj_active_hydro && cj->nodeID == engine_rank))) {

        scheduler_activate(s, t);

        /* Set the correct sorting flags */
        if (t->type == task_type_pair && t->subtype == task_subtype_density) {

          /* Store some values. */
          atomic_or(&ci->requires_sorts, 1 << t->flags);
          atomic_or(&cj->requires_sorts, 1 << t->flags);
          ci->dx_max_sort_old = ci->dx_max_sort;
          cj->dx_max_sort_old = cj->dx_max_sort;

          /* Activate the hydro drift tasks. */
          if (ci->nodeID == engine_rank) cell_activate_drift_part(ci, s);
          if (cj->nodeID == engine_rank) cell_activate_drift_part(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_sorts(ci, t->flags, s);
          cell_activate_sorts(cj, t->flags, s);

        }

        /* Store current values of dx_max and h_max. */
        else if (t->type == task_type_sub_pair &&
                 t->subtype == task_subtype_density) {
          cell_activate_subcell_hydro_tasks(t->ci, t->cj, s);
        }
      }

      if ((t->subtype == task_subtype_grav) &&
          ((ci_active_gravity && ci->nodeID == engine_rank) ||
           (cj_active_gravity && cj->nodeID == engine_rank))) {

        scheduler_activate(s, t);

        if (t->type == task_type_pair && t->subtype == task_subtype_grav) {
          /* Activate the gravity drift */
          cell_activate_subcell_grav_tasks(t->ci, t->cj, s);
        }

        else if (t->type == task_type_sub_pair &&
                 t->subtype == task_subtype_grav) {
          error("Invalid task sub-type encountered");
        }
      }

      /* Only interested in density tasks as of here. */
      if (t->subtype == task_subtype_density) {

        /* Too much particle movement? */
        if (cell_need_rebuild_for_pair(ci, cj)) *rebuild_space = 1;

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci->nodeID != engine_rank) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (cj_active_hydro) {
            scheduler_activate(s, ci->recv_xv);
            if (ci_active_hydro) {
              scheduler_activate(s, ci->recv_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate(s, ci->recv_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_hydro) scheduler_activate(s, ci->recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_hydro) {

            struct link *l =
                scheduler_activate_send(s, cj->send_xv, ci->nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (cj_active_hydro) {
              scheduler_activate_send(s, cj->send_rho, ci->nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, cj->send_gradient, ci->nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_hydro)
            scheduler_activate_send(s, cj->send_ti, ci->nodeID);

        } else if (cj->nodeID != engine_rank) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_hydro) {
            scheduler_activate(s, cj->recv_xv);
            if (cj_active_hydro) {
              scheduler_activate(s, cj->recv_rho);
#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate(s, cj->recv_gradient);
#endif
            }
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_hydro) scheduler_activate(s, cj->recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_hydro) {

            struct link *l =
                scheduler_activate_send(s, ci->send_xv, cj->nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_part(l->t->ci, s);

            /* If the local cell is also active, more stuff will be needed. */
            if (ci_active_hydro) {

              scheduler_activate_send(s, ci->send_rho, cj->nodeID);

#ifdef EXTRA_HYDRO_LOOP
              scheduler_activate_send(s, ci->send_gradient, cj->nodeID);
#endif
            }
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_hydro)
            scheduler_activate_send(s, ci->send_ti, cj->nodeID);
        }
#endif
      }

      /* Only interested in gravity tasks as of here. */
      if (t->subtype == task_subtype_grav) {

#ifdef WITH_MPI
        /* Activate the send/recv tasks. */
        if (ci->nodeID != engine_rank) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (cj_active_gravity) {
            scheduler_activate(s, ci->recv_grav);
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (ci_active_gravity) scheduler_activate(s, ci->recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (ci_active_gravity) {

            struct link *l =
                scheduler_activate_send(s, cj->send_grav, ci->nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (cj_active_gravity)
            scheduler_activate_send(s, cj->send_ti, ci->nodeID);

        } else if (cj->nodeID != engine_rank) {

          /* If the local cell is active, receive data from the foreign cell. */
          if (ci_active_gravity) {
            scheduler_activate(s, cj->recv_grav);
          }

          /* If the foreign cell is active, we want its ti_end values. */
          if (cj_active_gravity) scheduler_activate(s, cj->recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          if (cj_active_gravity) {

            struct link *l =
                scheduler_activate_send(s, ci->send_grav, cj->nodeID);

            /* Drift the cell which will be sent at the level at which it is
               sent, i.e. drift the cell specified in the send task (l->t)
               itself. */
            cell_activate_drift_gpart(l->t->ci, s);
          }

          /* If the local cell is active, send its ti_end values. */
          if (ci_active_gravity)
            scheduler_activate_send(s, ci->send_ti, cj->nodeID);
        }
#endif
      }
    }

    /* End force ? */
    else if (t->type == task_type_end_force) {

      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Kick ? */
    else if (t->type == task_type_kick1 || t->type == task_type_kick2) {

      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Hydro ghost tasks ? */
    else if (t->type == task_type_ghost || t->type == task_type_extra_ghost ||
             t->type == task_type_ghost_in || t->type == task_type_ghost_out) {
      if (cell_is_active_hydro(t->ci, e)) scheduler_activate(s, t);
    }

    /* Gravity stuff ? */
    else if (t->type == task_type_grav_down || t->type == task_type_grav_mesh ||
             t->type == task_type_grav_long_range ||
             t->type == task_type_init_grav ||
             t->type == task_type_init_grav_out ||
             t->type == task_type_grav_down_in) {
      if (cell_is_active_gravity(t->ci, e)) scheduler_activate(s, t);
    }

    else if (t->type == task_type_grav_mm) {

      /* Local pointers. */
      const struct cell *ci = t->ci;
      const struct cell *cj = t->cj;
      const int ci_nodeID = ci->nodeID;
      const int cj_nodeID = cj->nodeID;
      const int ci_active_gravity = cell_is_active_gravity(ci, e);
      const int cj_active_gravity = cell_is_active_gravity(cj, e);

      if ((ci_active_gravity && ci_nodeID == engine_rank) ||
          (cj_active_gravity && cj_nodeID == engine_rank))
        scheduler_activate(s, t);
    }

    /* Time-step? */
    else if (t->type == task_type_timestep) {
      t->ci->updated = 0;
      t->ci->g_updated = 0;
      t->ci->s_updated = 0;
      if (cell_is_active_hydro(t->ci, e) || cell_is_active_gravity(t->ci, e))
        scheduler_activate(s, t);
    }

    /* Subgrid tasks */
    else if (t->type == task_type_cooling || t->type == task_type_sourceterms) {
      if (cell_is_active_hydro(t->ci, e)) scheduler_activate(s, t);
    }
  }
}

/**
 * @brief Mark tasks to be un-skipped and set the sort flags accordingly.
 *
 * @return 1 if the space has to be rebuilt, 0 otherwise.
 */
int engine_marktasks(struct engine *e) {

  struct scheduler *s = &e->sched;
  const ticks tic = getticks();
  int rebuild_space = 0;

  /* Run through the tasks and mark as skip or not. */
  size_t extra_data[3] = {(size_t)e, (size_t)rebuild_space, (size_t)&e->sched};
  threadpool_map(&e->threadpool, engine_marktasks_mapper, s->tasks, s->nr_tasks,
                 sizeof(struct task), 0, extra_data);
  rebuild_space = extra_data[1];

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

  /* All is well... */
  return rebuild_space;
}

/**
 * @brief Prints the number of tasks in the engine
 *
 * @param e The #engine.
 */
void engine_print_task_counts(struct engine *e) {

  const ticks tic = getticks();
  struct scheduler *const sched = &e->sched;
  const int nr_tasks = sched->nr_tasks;
  const struct task *const tasks = sched->tasks;

  /* Count and print the number of each task type. */
  int counts[task_type_count + 1];
  for (int k = 0; k <= task_type_count; k++) counts[k] = 0;
  for (int k = 0; k < nr_tasks; k++) {
    if (tasks[k].skip)
      counts[task_type_count] += 1;
    else
      counts[(int)tasks[k].type] += 1;
  }
  message("Total = %d  (per cell = %d)", nr_tasks,
          (int)ceil((double)nr_tasks / e->s->tot_cells));
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
  message("nr_parts = %zu.", e->s->nr_parts);
  message("nr_gparts = %zu.", e->s->nr_gparts);
  message("nr_sparts = %zu.", e->s->nr_sparts);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief if necessary, estimate the number of tasks required given
 *        the current tasks in use and the numbers of cells.
 *
 * If e->tasks_per_cell is set greater than 0 then that value is used
 * as the estimate of the average number of tasks per cell,
 * otherwise we attempt an estimate.
 *
 * @param e the #engine
 *
 * @return the estimated total number of tasks
 */
int engine_estimate_nr_tasks(struct engine *e) {

  int tasks_per_cell = e->tasks_per_cell;
  if (tasks_per_cell > 0) return e->s->tot_cells * tasks_per_cell;

  /* Our guess differs depending on the types of tasks we are using, but we
   * basically use a formula <n1>*ntopcells + <n2>*(totcells - ntopcells).
   * Where <n1> is the expected maximum tasks per top-level/super cell, and
   * <n2> the expected maximum tasks for all other cells. These should give
   * a safe upper limit.
   */
  int n1 = 0;
  int n2 = 0;
  if (e->policy & engine_policy_hydro) {
    n1 += 37;
    n2 += 2;
#ifdef WITH_MPI
    n1 += 6;
#endif

#ifdef EXTRA_HYDRO_LOOP
    n1 += 15;
#ifdef WITH_MPI
    n1 += 2;
#endif
#endif
  }
  if (e->policy & engine_policy_self_gravity) {
    n1 += 125;
    n2 += 8;
#ifdef WITH_MPI
    n2 += 2;
#endif
  }
  if (e->policy & engine_policy_external_gravity) {
    n1 += 2;
  }
  if (e->policy & engine_policy_cosmology) {
    n1 += 2;
  }
  if (e->policy & engine_policy_cooling) {
    n1 += 2;
  }
  if (e->policy & engine_policy_sourceterms) {
    n1 += 2;
  }
  if (e->policy & engine_policy_stars) {
    n1 += 2;
  }

#ifdef WITH_MPI

  /* We need fewer tasks per rank when using MPI, but we could have
   * imbalances, so we need to work using the locally active cells, not just
   * some equipartition amongst the nodes. Don't want to recurse the whole
   * cell tree, so just make a guess of the maximum possible total cells. */
  int ntop = 0;
  int ncells = 0;
  for (int k = 0; k < e->s->nr_cells; k++) {
    struct cell *c = &e->s->cells_top[k];

    /* Any cells with particles will have tasks (local & foreign). */
    int nparts = c->count + c->gcount + c->scount;
    if (nparts > 0) {
      ntop++;
      ncells++;

      /* Count cell depth until we get below the parts per cell threshold. */
      int depth = 0;
      while (nparts > space_splitsize) {
        depth++;
        nparts /= 8;
        ncells += (1 << (depth * 3));
      }
    }
  }

  /* If no local cells, we are probably still initialising, so just keep
   * room for the top-level. */
  if (ncells == 0) {
    ntop = e->s->nr_cells;
    ncells = ntop;
  }
#else
  int ntop = e->s->nr_cells;
  int ncells = e->s->tot_cells;
#endif

  double ntasks = n1 * ntop + n2 * (ncells - ntop);
  if (ncells > 0) tasks_per_cell = ceil(ntasks / ncells);

  if (tasks_per_cell < 1.0) tasks_per_cell = 1.0;
  if (e->verbose)
    message("tasks per cell estimated as: %d, maximum tasks: %d",
            tasks_per_cell, ncells * tasks_per_cell);

  return ncells * tasks_per_cell;
}

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 * @param clean_smoothing_length_values Are we cleaning up the values of
 * the smoothing lengths before building the tasks ?
 */
void engine_rebuild(struct engine *e, int clean_smoothing_length_values) {

  const ticks tic = getticks();

  /* Clear the forcerebuild flag, whatever it was. */
  e->forcerebuild = 0;
  e->restarting = 0;

  /* Re-build the space. */
  space_rebuild(e->s, e->verbose);

  /* Re-compute the mesh forces */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
    pm_mesh_compute_potential(e->mesh, e->s, e->verbose);

  /* Re-compute the maximal RMS displacement constraint */
  if (e->policy & engine_policy_cosmology)
    engine_recompute_displacement_constraint(e);

#ifdef SWIFT_DEBUG_CHECKS
  part_verify_links(e->s->parts, e->s->gparts, e->s->sparts, e->s->nr_parts,
                    e->s->nr_gparts, e->s->nr_sparts, e->verbose);
#endif

  /* Initial cleaning up session ? */
  if (clean_smoothing_length_values) space_sanitize(e->s);

/* If in parallel, exchange the cell structure, top-level and neighbouring
 * multipoles. */
#ifdef WITH_MPI
  engine_exchange_cells(e);

  if (e->policy & engine_policy_self_gravity) engine_exchange_top_multipoles(e);

  if (e->policy & engine_policy_self_gravity)
    engine_exchange_proxy_multipoles(e);
#endif

  /* Re-build the tasks. */
  engine_maketasks(e);

  /* Make the list of top-level cells that have tasks */
  space_list_cells_with_tasks(e->s);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);
#endif

  /* Run through the tasks and mark as skip or not. */
  if (engine_marktasks(e))
    error("engine_marktasks failed after space_rebuild.");

  /* Print the status of the system */
  if (e->verbose) engine_print_task_counts(e);

  /* Clear the counters of updates since the last rebuild */
  e->updates_since_rebuild = 0;
  e->g_updates_since_rebuild = 0;
  e->s_updates_since_rebuild = 0;

  /* Flag that a rebuild has taken place */
  e->step_props |= engine_step_prop_rebuild;

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

  TIMER_TIC2;
  const ticks tic = getticks();

  int drifted_all = 0;

  /* Unskip active tasks and check for rebuild */
  if (!e->forcerebuild && !e->forcerepart && !e->restarting) engine_unskip(e);

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->forcerebuild, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#endif

  /* Do we need repartitioning ? */
  if (e->forcerepart) {

    /* Let's start by drifting everybody to the current time */
    engine_drift_all(e);
    drifted_all = 1;

    /* And repartition */
    engine_repartition(e);
  }

  /* Do we need rebuilding ? */
  if (e->forcerebuild) {

    /* Let's start by drifting everybody to the current time */
    if (!e->restarting && !drifted_all) engine_drift_all(e);

    /* And rebuild */
    engine_rebuild(e, 0);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (e->forcerepart || e->forcerebuild) {
    /* Check that all cells have been drifted to the current time.
     * That can include cells that have not previously been active on this
     * rank. Skip if haven't got any cells (yet). */
    if (e->s->cells_top != NULL)
      space_check_drift_point(e->s, e->ti_current,
                              e->policy & engine_policy_self_gravity);
  }
#endif

  /* Re-rank the tasks every now and then. XXX this never executes. */
  if (e->tasks_age % engine_tasksreweight == 1) {
    scheduler_reweight(&e->sched, e->verbose);
  }
  e->tasks_age += 1;

  TIMER_TOC2(timer_prepare);

  if (e->verbose)
    message("took %.3f %s (including unskip, rebuild and reweight).",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
}

/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 */
void engine_barrier(struct engine *e) {

  /* Wait at the wait barrier. */
  swift_barrier_wait(&e->wait_barrier);

  /* Wait at the run barrier. */
  swift_barrier_wait(&e->run_barrier);
}

/**
 * @brief Mapping function to collect the data from the kick.
 *
 * @param c A super-cell.
 */
void engine_collect_end_of_step_recurse(struct cell *c) {

/* Skip super-cells (Their values are already set) */
#ifdef WITH_MPI
  if (c->timestep != NULL || c->recv_ti != NULL) return;
#else
  if (c->timestep != NULL) return;
#endif /* WITH_MPI */

  /* Counters for the different quantities. */
  size_t updated = 0, g_updated = 0, s_updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL && (cp->count > 0 || cp->gcount > 0 || cp->scount > 0)) {

      /* Recurse */
      engine_collect_end_of_step_recurse(cp);

      /* And update */
      ti_hydro_end_min = min(ti_hydro_end_min, cp->ti_hydro_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, cp->ti_hydro_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, cp->ti_hydro_beg_max);
      ti_gravity_end_min = min(ti_gravity_end_min, cp->ti_gravity_end_min);
      ti_gravity_end_max = max(ti_gravity_end_max, cp->ti_gravity_end_max);
      ti_gravity_beg_max = max(ti_gravity_beg_max, cp->ti_gravity_beg_max);
      updated += cp->updated;
      g_updated += cp->g_updated;
      s_updated += cp->s_updated;

      /* Collected, so clear for next time. */
      cp->updated = 0;
      cp->g_updated = 0;
      cp->s_updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->ti_hydro_end_min = ti_hydro_end_min;
  c->ti_hydro_end_max = ti_hydro_end_max;
  c->ti_hydro_beg_max = ti_hydro_beg_max;
  c->ti_gravity_end_min = ti_gravity_end_min;
  c->ti_gravity_end_max = ti_gravity_end_max;
  c->ti_gravity_beg_max = ti_gravity_beg_max;
  c->updated = updated;
  c->g_updated = g_updated;
  c->s_updated = s_updated;
}

void engine_collect_end_of_step_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct end_of_step_data *data = (struct end_of_step_data *)extra_data;
  struct engine *e = data->e;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;

  /* Local collectible */
  size_t updates = 0, g_updates = 0, s_updates = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];

    if (c->count > 0 || c->gcount > 0 || c->scount > 0) {

      /* Make the top-cells recurse */
      engine_collect_end_of_step_recurse(c);

      /* And aggregate */
      ti_hydro_end_min = min(ti_hydro_end_min, c->ti_hydro_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, c->ti_hydro_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, c->ti_hydro_beg_max);
      ti_gravity_end_min = min(ti_gravity_end_min, c->ti_gravity_end_min);
      ti_gravity_end_max = max(ti_gravity_end_max, c->ti_gravity_end_max);
      ti_gravity_beg_max = max(ti_gravity_beg_max, c->ti_gravity_beg_max);
      updates += c->updated;
      g_updates += c->g_updated;
      s_updates += c->s_updated;

      /* Collected, so clear for next time. */
      c->updated = 0;
      c->g_updated = 0;
      c->s_updated = 0;
    }
  }

  /* Let's write back to the global data.
   * We use the space lock to garanty single access*/
  if (lock_lock(&s->lock) == 0) {
    data->updates += updates;
    data->g_updates += g_updates;
    data->s_updates += s_updates;
    data->ti_hydro_end_min = min(ti_hydro_end_min, data->ti_hydro_end_min);
    data->ti_hydro_end_max = max(ti_hydro_end_max, data->ti_hydro_end_max);
    data->ti_hydro_beg_max = max(ti_hydro_beg_max, data->ti_hydro_beg_max);
    data->ti_gravity_end_min =
        min(ti_gravity_end_min, data->ti_gravity_end_min);
    data->ti_gravity_end_max =
        max(ti_gravity_end_max, data->ti_gravity_end_max);
    data->ti_gravity_beg_max =
        max(ti_gravity_beg_max, data->ti_gravity_beg_max);
  }
  if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space");
}

/**
 * @brief Collects the next time-step and rebuild flag.
 *
 * The next time-step is determined by making each super-cell recurse to
 * collect the minimal of ti_end and the number of updated particles.  When in
 * MPI mode this routines reduces these across all nodes and also collects the
 * forcerebuild flag -- this is so that we only use a single collective MPI
 * call per step for all these values.
 *
 * Note that the results are stored in e->collect_group1 struct not in the
 * engine fields, unless apply is true. These can be applied field-by-field
 * or all at once using collectgroup1_copy();
 *
 * @param e The #engine.
 * @param apply whether to apply the results to the engine or just keep in the
 *              group1 struct.
 */
void engine_collect_end_of_step(struct engine *e, int apply) {

  const ticks tic = getticks();
  const struct space *s = e->s;
  struct end_of_step_data data;
  data.updates = 0, data.g_updates = 0, data.s_updates = 0;
  data.ti_hydro_end_min = max_nr_timesteps, data.ti_hydro_end_max = 0,
  data.ti_hydro_beg_max = 0;
  data.ti_gravity_end_min = max_nr_timesteps, data.ti_gravity_end_max = 0,
  data.ti_gravity_beg_max = 0;
  data.e = e;

  /* Collect information from the local top-level cells */
  threadpool_map(&e->threadpool, engine_collect_end_of_step_mapper,
                 s->local_cells_top, s->nr_local_cells, sizeof(int), 0, &data);

  /* Store these in the temporary collection group. */
  collectgroup1_init(&e->collect_group1, data.updates, data.g_updates,
                     data.s_updates, data.ti_hydro_end_min,
                     data.ti_hydro_end_max, data.ti_hydro_beg_max,
                     data.ti_gravity_end_min, data.ti_gravity_end_max,
                     data.ti_gravity_beg_max, e->forcerebuild);

/* Aggregate collective data from the different nodes for this step. */
#ifdef WITH_MPI
  collectgroup1_reduce(&e->collect_group1);

#ifdef SWIFT_DEBUG_CHECKS
  {
    /* Check the above using the original MPI calls. */
    integertime_t in_i[2], out_i[2];
    in_i[0] = 0;
    in_i[1] = 0;
    out_i[0] = data.ti_hydro_end_min;
    out_i[1] = data.ti_gravity_end_min;
    if (MPI_Allreduce(out_i, in_i, 2, MPI_LONG_LONG_INT, MPI_MIN,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate ti_end_min.");
    if (in_i[0] != (long long)e->collect_group1.ti_hydro_end_min)
      error("Failed to get same ti_hydro_end_min, is %lld, should be %lld",
            in_i[0], e->collect_group1.ti_hydro_end_min);
    if (in_i[1] != (long long)e->collect_group1.ti_gravity_end_min)
      error("Failed to get same ti_gravity_end_min, is %lld, should be %lld",
            in_i[1], e->collect_group1.ti_gravity_end_min);

    long long in_ll[3], out_ll[3];
    out_ll[0] = data.updates;
    out_ll[1] = data.g_updates;
    out_ll[2] = data.s_updates;
    if (MPI_Allreduce(out_ll, in_ll, 3, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.updates)
      error("Failed to get same updates, is %lld, should be %lld", in_ll[0],
            e->collect_group1.updates);
    if (in_ll[1] != (long long)e->collect_group1.g_updates)
      error("Failed to get same g_updates, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_updates);
    if (in_ll[2] != (long long)e->collect_group1.s_updates)
      error("Failed to get same s_updates, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_updates);

    int buff = 0;
    if (MPI_Allreduce(&e->forcerebuild, &buff, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    if (!!buff != !!e->collect_group1.forcerebuild)
      error(
          "Failed to get same rebuild flag from all nodes, is %d,"
          "should be %d",
          buff, e->collect_group1.forcerebuild);
  }
#endif
#endif

  /* Apply to the engine, if requested. */
  if (apply) collectgroup1_apply(&e->collect_group1, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Print the conserved quantities statistics to a log file
 *
 * @param e The #engine.
 */
void engine_print_stats(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);

  /* Be verbose about this */
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Saving statistics at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Saving statistics at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#else
  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Saving statistics at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Saving statistics at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

  struct statistics stats;
  stats_init(&stats);

  /* Collect the stats on this node */
  stats_collect(e->s, &stats);

/* Aggregate the data from the different nodes. */
#ifdef WITH_MPI
  struct statistics global_stats;
  stats_init(&global_stats);

  if (MPI_Reduce(&stats, &global_stats, 1, statistics_mpi_type,
                 statistics_mpi_reduce_op, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to aggregate stats.");
#else
  struct statistics global_stats = stats;
#endif

  /* Finalize operations */
  stats_finalize(&stats);

  /* Print info */
  if (e->nodeID == 0)
    stats_print_to_file(e->file_stats, &global_stats, e->time);

  /* Flag that we dumped some statistics */
  e->step_props |= engine_step_prop_statistics;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Sets all the force, drift and kick tasks to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_skip_force_and_kick(struct engine *e) {

  struct task *tasks = e->sched.tasks;
  const int nr_tasks = e->sched.nr_tasks;

  for (int i = 0; i < nr_tasks; ++i) {

    struct task *t = &tasks[i];

    /* Skip everything that updates the particles */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_kick1 || t->type == task_type_kick2 ||
        t->type == task_type_timestep || t->subtype == task_subtype_force ||
        t->subtype == task_subtype_grav || t->type == task_type_end_force ||
        t->type == task_type_grav_long_range || t->type == task_type_grav_mm ||
        t->type == task_type_grav_down || t->type == task_type_cooling ||
        t->type == task_type_sourceterms)
      t->skip = 1;
  }

  /* Run through the cells and clear some flags. */
  space_map_cells_pre(e->s, 1, cell_clear_drift_flags, NULL);
}

/**
 * @brief Sets all the drift and first kick tasks to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_skip_drift(struct engine *e) {

  struct task *tasks = e->sched.tasks;
  const int nr_tasks = e->sched.nr_tasks;

  for (int i = 0; i < nr_tasks; ++i) {

    struct task *t = &tasks[i];

    /* Skip everything that moves the particles */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart)
      t->skip = 1;
  }

  /* Run through the cells and clear some flags. */
  space_map_cells_pre(e->s, 1, cell_clear_drift_flags, NULL);
}

/**
 * @brief Launch the runners.
 *
 * @param e The #engine.
 */
void engine_launch(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Re-set all the cell task counters to 0 */
  space_reset_task_counters(e->s);
#endif

  /* Prepare the scheduler. */
  atomic_inc(&e->sched.waiting);

  /* Cry havoc and let loose the dogs of war. */
  swift_barrier_wait(&e->run_barrier);

  /* Load the tasks. */
  scheduler_start(&e->sched);

  /* Remove the safeguard. */
  pthread_mutex_lock(&e->sched.sleep_mutex);
  atomic_dec(&e->sched.waiting);
  pthread_cond_broadcast(&e->sched.sleep_cond);
  pthread_mutex_unlock(&e->sched.sleep_mutex);

  /* Sit back and wait for the runners to come home. */
  swift_barrier_wait(&e->wait_barrier);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Calls the 'first init' function on the particles of all types.
 *
 * @param e The #engine.
 */
void engine_first_init_particles(struct engine *e) {

  const ticks tic = getticks();

  /* Set the particles in a state where they are ready for a run. */
  space_first_init_parts(e->s, e->verbose);
  space_first_init_gparts(e->s, e->verbose);
  space_first_init_sparts(e->s, e->verbose);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Initialises the particles and set them in a state ready to move
 *forward in time.
 *
 * @param e The #engine
 * @param flag_entropy_ICs Did the 'Internal Energy' of the particles actually
 * contain entropy ?
 * @param clean_h_values Are we cleaning up the values of h before building
 * the tasks ?
 */
void engine_init_particles(struct engine *e, int flag_entropy_ICs,
                           int clean_h_values) {

  struct space *s = e->s;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_update(e->gravity_properties, e->cosmology);

  /* Start by setting the particles in a good state */
  if (e->nodeID == 0) message("Setting particles to a valid state...");
  engine_first_init_particles(e);

  if (e->nodeID == 0) message("Computing initial gas densities.");

  /* Construct all cells and tasks to start everything */
  engine_rebuild(e, clean_h_values);

  /* No time integration. We just want the density and ghosts */
  engine_skip_force_and_kick(e);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

  /* Init the particle data (by hand). */
  space_init_parts(s, e->verbose);
  space_init_gparts(s, e->verbose);

  /* Now, launch the calculation */
  TIMER_TIC;
  engine_launch(e);
  TIMER_TOC(timer_runners);

  /* Apply some conversions (e.g. internal energy -> entropy) */
  if (!flag_entropy_ICs) {

    if (e->nodeID == 0) message("Converting internal energy variable.");

    space_convert_quantities(e->s, e->verbose);

    /* Correct what we did (e.g. in PE-SPH, need to recompute rho_bar) */
    if (hydro_need_extra_init_loop) {
      engine_marktasks(e);
      engine_skip_force_and_kick(e);
      engine_launch(e);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].multipole->m_pole.num_gpart;
    if (num_gpart_mpole != e->total_nr_gparts)
      error(
          "Top-level multipoles don't contain the total number of gpart "
          "s->nr_gpart=%lld, "
          "m_poles=%lld",
          e->total_nr_gparts, num_gpart_mpole);
  }
#endif

  /* Now time to get ready for the first time-step */
  if (e->nodeID == 0) message("Running initial fake time-step.");

  /* Prepare all the tasks again for a new round */
  engine_marktasks(e);

  /* No drift this time */
  engine_skip_drift(e);

  /* Init the particle data (by hand). */
  space_init_parts(e->s, e->verbose);
  space_init_gparts(e->s, e->verbose);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Run the brute-force gravity calculation for some gparts */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_compute(e->s, e);
#endif

  if (e->nodeID == 0) scheduler_write_dependencies(&e->sched, e->verbose);

  /* Run the 0th time-step */
  TIMER_TIC2;
  engine_launch(e);
  TIMER_TOC2(timer_runners);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Check the accuracy of the gravity calculation */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_check(e->s, e, 1e-1);
#endif

  /* Recover the (integer) end of the next time-step */
  engine_collect_end_of_step(e, 1);

  /* Check if any particles have the same position. This is not
   * allowed (/0) so we abort.*/
  if (s->nr_parts > 0) {

    /* Sorting should put the same positions next to each other... */
    int failed = 0;
    double *prev_x = s->parts[0].x;
    long long *prev_id = &s->parts[0].id;
    for (size_t k = 1; k < s->nr_parts; k++) {
      if (prev_x[0] == s->parts[k].x[0] && prev_x[1] == s->parts[k].x[1] &&
          prev_x[2] == s->parts[k].x[2]) {
        if (e->verbose)
          message("Two particles occupy location: %f %f %f id=%lld id=%lld",
                  prev_x[0], prev_x[1], prev_x[2], *prev_id, s->parts[k].id);
        failed++;
      }
      prev_x = s->parts[k].x;
      prev_id = &s->parts[k].id;
    }
    if (failed > 0)
      error(
          "Have %d particle pairs with the same locations.\n"
          "Cannot continue",
          failed);
  }

  /* Also check any gparts. This is not supposed to be fatal so only warn. */
  if (s->nr_gparts > 0) {
    int failed = 0;
    double *prev_x = s->gparts[0].x;
    for (size_t k = 1; k < s->nr_gparts; k++) {
      if (prev_x[0] == s->gparts[k].x[0] && prev_x[1] == s->gparts[k].x[1] &&
          prev_x[2] == s->gparts[k].x[2]) {
        if (e->verbose)
          message("Two gparts occupy location: %f %f %f / %f %f %f", prev_x[0],
                  prev_x[1], prev_x[2], s->gparts[k].x[0], s->gparts[k].x[1],
                  s->gparts[k].x[2]);
        failed++;
      }
      prev_x = s->gparts[k].x;
    }
    if (failed > 0)
      message(
          "WARNING: found %d gpart pairs at the same location. "
          "That is not optimal",
          failed);
  }

  /* Check the top-level cell h_max matches the particles as these can be
   * updated in the the ghost tasks (only a problem if the ICs estimates for h
   * are too small). Note this must be followed by a rebuild as sub-cells will
   * not be updated until that is done. */
  if (s->cells_top != NULL && s->nr_parts > 0) {
    for (int i = 0; i < s->nr_cells; i++) {
      struct cell *c = &s->cells_top[i];
      if (c->nodeID == engine_rank && c->count > 0) {
        float part_h_max = c->parts[0].h;
        for (int k = 1; k < c->count; k++) {
          if (c->parts[k].h > part_h_max) part_h_max = c->parts[k].h;
        }
        c->h_max = max(part_h_max, c->h_max);
      }
    }
  }

  clocks_gettime(&time2);

#ifdef SWIFT_DEBUG_CHECKS
  space_check_timesteps(e->s);
  part_verify_links(e->s->parts, e->s->gparts, e->s->sparts, e->s->nr_parts,
                    e->s->nr_gparts, e->s->nr_sparts, e->verbose);
#endif

  /* Ready to go */
  e->step = 0;
  e->forcerebuild = 1;
  e->wallclock_time = (float)clocks_diff(&time1, &time2);

  if (e->verbose) message("took %.3f %s.", e->wallclock_time, clocks_getunit());
}

/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 */
void engine_step(struct engine *e) {

  TIMER_TIC2;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

#ifdef SWIFT_DEBUG_TASKS
  e->tic_step = getticks();
#endif

  if (e->nodeID == 0) {

    /* Print some information to the screen */
    printf(
        "  %6d %14e %14e %10.5f %14e %4d %4d %12lld %12lld %12lld %21.3f %6d\n",
        e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
        e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
        e->s_updates, e->wallclock_time, e->step_props);
    fflush(stdout);

    if (!e->restarting)
      fprintf(e->file_timesteps,
              "  %6d %14e %14e %10.5f %14e %4d %4d %12lld %12lld %12lld %21.3f "
              "%6d\n",
              e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
              e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
              e->s_updates, e->wallclock_time, e->step_props);
    fflush(e->file_timesteps);
  }

  /* We need some cells to exist but not the whole task stuff. */
  if (e->restarting) space_rebuild(e->s, e->verbose);

  /* Move forward in time */
  e->ti_old = e->ti_current;
  e->ti_current = e->ti_end_min;
  e->max_active_bin = get_max_active_bin(e->ti_end_min);
  e->min_active_bin = get_min_active_bin(e->ti_current, e->ti_old);
  e->step += 1;
  e->step_props = engine_step_prop_none;

  /* When restarting, move everyone to the current time. */
  if (e->restarting) engine_drift_all(e);

  /* Get the physical value of the time and time-step size */
  if (e->policy & engine_policy_cosmology) {
    e->time_old = e->time;
    cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
    e->time = e->cosmology->time;
    e->time_step = e->time - e->time_old;
  } else {
    e->time = e->ti_current * e->time_base + e->time_begin;
    e->time_old = e->ti_old * e->time_base + e->time_begin;
    e->time_step = (e->ti_current - e->ti_old) * e->time_base;
  }

  /*****************************************************/
  /* OK, we now know what the next end of time-step is */
  /*****************************************************/

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_update(e->gravity_properties, e->cosmology);

  /* Trigger a tree-rebuild if we passed the frequency threshold */
  if ((e->policy & engine_policy_self_gravity) &&
      ((double)e->g_updates_since_rebuild >
       ((double)e->total_nr_gparts) * e->gravity_properties->rebuild_frequency))
    e->forcerebuild = 1;

  /* Are we drifting everything (a la Gadget/GIZMO) ? */
  if (e->policy & engine_policy_drift_all && !e->forcerebuild)
    engine_drift_all(e);

  /* Are we reconstructing the multipoles or drifting them ?*/
  if ((e->policy & engine_policy_self_gravity) && !e->forcerebuild) {

    if (e->policy & engine_policy_reconstruct_mpoles)
      engine_reconstruct_multipoles(e);
    else
      engine_drift_top_multipoles(e);
  }

#ifdef WITH_MPI
  /* Repartition the space amongst the nodes? */
  engine_repartition_trigger(e);
#endif

  /* Prepare the tasks to be launched, rebuild or repartition if needed. */
  engine_prepare(e);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

    /* Dump local cells and active particle counts. */
    /* dumpCells("cells", 0, 0, 0, 0, e->s, e->nodeID, e->step); */

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].multipole->m_pole.num_gpart;
    if (num_gpart_mpole != e->total_nr_gparts)
      error(
          "Multipoles don't contain the total number of gpart mpoles=%lld "
          "ngparts=%lld",
          num_gpart_mpole, e->total_nr_gparts);
  }
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Run the brute-force gravity calculation for some gparts */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_compute(e->s, e);
#endif

  /* Start all the tasks. */
  TIMER_TIC;
  engine_launch(e);
  TIMER_TOC(timer_runners);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Check the accuracy of the gravity calculation */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_check(e->s, e, 1e-1);
#endif

  /* Collect information about the next time-step */
  engine_collect_end_of_step(e, 0);
  e->forcerebuild = e->collect_group1.forcerebuild;

  /* Now apply all the collected time step updates and particle counts. */
  collectgroup1_apply(&e->collect_group1, e);
  e->updates_since_rebuild += e->collect_group1.updates;
  e->g_updates_since_rebuild += e->collect_group1.g_updates;
  e->s_updates_since_rebuild += e->collect_group1.s_updates;

  /********************************************************/
  /* OK, we are done with the regular stuff. Time for i/o */
  /********************************************************/

  /* Create a restart file if needed. */
  engine_dump_restarts(e, 0, e->restart_onexit && engine_is_done(e));

  /* Save some statistics ? */
  int save_stats = 0;
  if (e->ti_end_min > e->ti_next_stats && e->ti_next_stats > 0) save_stats = 1;

  /* Do we want a snapshot? */
  int dump_snapshot = 0;
  if (e->ti_end_min > e->ti_next_snapshot && e->ti_next_snapshot > 0)
    dump_snapshot = 1;

  /* Do we want to perform structure finding? */
  int run_stf = 0;
  if ((e->policy & engine_policy_structure_finding)) {
    if (e->stf_output_freq_format == STEPS && e->step % e->deltaStepSTF == 0)
      run_stf = 1;
    else if (e->stf_output_freq_format == TIME &&
             e->ti_end_min > e->ti_nextSTF && e->ti_nextSTF > 0)
      run_stf = 1;
  }

  /* Store information before attempting extra dump-related drifts */
  integertime_t ti_current = e->ti_current;
  timebin_t max_active_bin = e->max_active_bin;
  double time = e->time;

  /* Write some form of output */
  if (dump_snapshot && save_stats) {

    /* If both, need to figure out which one occurs first */
    if (e->ti_next_stats == e->ti_next_snapshot) {

      /* Let's fake that we are at the common dump time */
      e->ti_current = e->ti_next_snapshot;
      e->max_active_bin = 0;
      if (!(e->policy & engine_policy_cosmology))
        e->time = e->ti_next_snapshot * e->time_base + e->time_begin;

      /* Drift everyone */
      engine_drift_all(e);

      /* Dump everything */
      engine_print_stats(e);
      engine_dump_snapshot(e);

    } else if (e->ti_next_stats < e->ti_next_snapshot) {

      /* Let's fake that we are at the stats dump time */
      e->ti_current = e->ti_next_stats;
      e->max_active_bin = 0;
      if (!(e->policy & engine_policy_cosmology))
        e->time = e->ti_next_stats * e->time_base + e->time_begin;

      /* Drift everyone */
      engine_drift_all(e);

      /* Dump stats */
      engine_print_stats(e);

      /* Let's fake that we are at the snapshot dump time */
      e->ti_current = e->ti_next_snapshot;
      e->max_active_bin = 0;
      if (!(e->policy & engine_policy_cosmology))
        e->time = e->ti_next_snapshot * e->time_base + e->time_begin;

      /* Drift everyone */
      engine_drift_all(e);

      /* Dump snapshot */
      engine_dump_snapshot(e);

    } else if (e->ti_next_stats > e->ti_next_snapshot) {

      /* Let's fake that we are at the snapshot dump time */
      e->ti_current = e->ti_next_snapshot;
      e->max_active_bin = 0;
      if (!(e->policy & engine_policy_cosmology))
        e->time = e->ti_next_snapshot * e->time_base + e->time_begin;

      /* Drift everyone */
      engine_drift_all(e);

      /* Dump snapshot */
      engine_dump_snapshot(e);

      /* Let's fake that we are at the stats dump time */
      e->ti_current = e->ti_next_stats;
      e->max_active_bin = 0;
      if (!(e->policy & engine_policy_cosmology))
        e->time = e->ti_next_stats * e->time_base + e->time_begin;

      /* Drift everyone */
      engine_drift_all(e);

      /* Dump stats */
      engine_print_stats(e);
    }

    /* Let's compute the time of the next outputs */
    engine_compute_next_snapshot_time(e);
    engine_compute_next_statistics_time(e);

  } else if (dump_snapshot) {

    /* Let's fake that we are at the snapshot dump time */
    e->ti_current = e->ti_next_snapshot;
    e->max_active_bin = 0;
    if (!(e->policy & engine_policy_cosmology))
      e->time = e->ti_next_snapshot * e->time_base + e->time_begin;

    /* Drift everyone */
    engine_drift_all(e);

    /* Dump... */
    engine_dump_snapshot(e);

    /* ... and find the next output time */
    engine_compute_next_snapshot_time(e);
  } else if (save_stats) {

    /* Let's fake that we are at the stats dump time */
    e->ti_current = e->ti_next_stats;
    e->max_active_bin = 0;
    if (!(e->policy & engine_policy_cosmology))
      e->time = e->ti_next_stats * e->time_base + e->time_begin;

    /* Drift everyone */
    engine_drift_all(e);

    /* Dump */
    engine_print_stats(e);

    /* and move on */
    engine_compute_next_statistics_time(e);
  }

  /* Perform structure finding? */
  if (run_stf) {

  // MATTHIEU: Add a drift_all here. And check the order with the order i/o
  // options.

#ifdef HAVE_VELOCIRAPTOR
    velociraptor_init(e);
    velociraptor_invoke(e);

    /* ... and find the next output time */
    if (e->stf_output_freq_format == TIME) engine_compute_next_stf_time(e);
#endif
  }

  /* Restore the information we stored */
  e->ti_current = ti_current;
  e->max_active_bin = max_active_bin;
  e->time = time;

  TIMER_TOC2(timer_step);

  clocks_gettime(&time2);
  e->wallclock_time = (float)clocks_diff(&time1, &time2);

#ifdef SWIFT_DEBUG_TASKS
  /* Time in ticks at the end of this step. */
  e->toc_step = getticks();
#endif
}

/**
 * @brief dump restart files if it is time to do so and dumps are enabled.
 *
 * @param e the engine.
 * @param drifted_all true if a drift_all has just been performed.
 * @param force force a dump, if dumping is enabled.
 */
void engine_dump_restarts(struct engine *e, int drifted_all, int force) {

  if (e->restart_dump) {
    ticks tic = getticks();

    /* Dump when the time has arrived, or we are told to. */
    int dump = ((tic > e->restart_next) || force);

#ifdef WITH_MPI
    /* Synchronize this action from rank 0 (ticks may differ between
     * machines). */
    MPI_Bcast(&dump, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (dump) {

      if (e->nodeID == 0) message("Writing restart files");

      /* Clean out the previous saved files, if found. Do this now as we are
       * MPI synchronized. */
      restart_remove_previous(e->restart_file);

      /* Drift all particles first (may have just been done). */
      if (!drifted_all) engine_drift_all(e);
      restart_write(e, e->restart_file);

      if (e->verbose)
        message("Dumping restart files took %.3f %s",
                clocks_from_ticks(getticks() - tic), clocks_getunit());

      /* Time after which next dump will occur. */
      e->restart_next += e->restart_dt;

      /* Flag that we dumped the restarts */
      e->step_props |= engine_step_prop_restarts;
    }
  }
}

/**
 * @brief Returns 1 if the simulation has reached its end point, 0 otherwise
 */
int engine_is_done(struct engine *e) {
  return !(e->ti_current < max_nr_timesteps);
}

/**
 * @brief Unskip all the tasks that act on active cells at this time.
 *
 * @param e The #engine.
 */
void engine_unskip(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;

#ifdef WITH_PROFILER
  static int count = 0;
  char filename[100];
  sprintf(filename, "/tmp/swift_runner_do_usnkip_mapper_%06i.prof", count++);
  ProfilerStart(filename);
#endif  // WITH_PROFILER

  /* Move the active local cells to the top of the list. */
  int *local_cells = e->s->local_cells_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];
    if ((e->policy & engine_policy_hydro && cell_is_active_hydro(c, e)) ||
        (e->policy &
             (engine_policy_self_gravity | engine_policy_external_gravity) &&
         cell_is_active_gravity(c, e))) {
      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, runner_do_unskip_mapper, local_cells,
                 num_active_cells, sizeof(int), 1, e);

#ifdef WITH_PROFILER
  ProfilerStop();
#endif  // WITH_PROFILER

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Mapper function to drift *all* particle types and multipoles forward
 * in time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_all_mapper(void *map_data, int num_elements,
                                void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &cells[ind];
    if (c != NULL && c->nodeID == e->nodeID) {
      /* Drift all the particles */
      cell_drift_part(c, e, 1);

      /* Drift all the g-particles */
      cell_drift_gpart(c, e, 1);
    }

    /* Drift the multipoles */
    if (e->policy & engine_policy_self_gravity) {
      cell_drift_all_multipoles(c, e);
    }
  }
}

/**
 * @brief Drift *all* particles and multipoles at all levels
 * forward to the current time.
 *
 * @param e The #engine.
 */
void engine_drift_all(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Drifting all to a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Drifting all to t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

  threadpool_map(&e->threadpool, engine_do_drift_all_mapper, e->s->cells_top,
                 e->s->nr_cells, sizeof(struct cell), 0, e);

  /* Synchronize particle positions */
  space_synchronize_particle_positions(e->s);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);
  part_verify_links(e->s->parts, e->s->gparts, e->s->sparts, e->s->nr_parts,
                    e->s->nr_gparts, e->s->nr_sparts, e->verbose);
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Mapper function to drift *all* top-level multipoles forward in
 * time.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void engine_do_drift_top_multipoles_mapper(void *map_data, int num_elements,
                                           void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &cells[ind];
    if (c != NULL) {

      /* Drift the multipole at this level only */
      if (c->ti_old_multipole != e->ti_current) cell_drift_multipole(c, e);
    }
  }
}

/**
 * @brief Drift *all* top-level multipoles forward to the current time.
 *
 * @param e The #engine.
 */
void engine_drift_top_multipoles(struct engine *e) {

  const ticks tic = getticks();

  threadpool_map(&e->threadpool, engine_do_drift_top_multipoles_mapper,
                 e->s->cells_top, e->s->nr_cells, sizeof(struct cell), 0, e);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time. */
  space_check_top_multipoles_drift_point(e->s, e->ti_current);
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void engine_do_reconstruct_multipoles_mapper(void *map_data, int num_elements,
                                             void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &cells[ind];
    if (c != NULL && c->nodeID == e->nodeID) {

      /* Construct the multipoles in this cell hierarchy */
      cell_make_multipoles(c, e->ti_current);
    }
  }
}

/**
 * @brief Reconstruct all the multipoles at all the levels in the tree.
 *
 * @param e The #engine.
 */
void engine_reconstruct_multipoles(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Reconstructing multipoles at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Reconstructing multipoles at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

  threadpool_map(&e->threadpool, engine_do_reconstruct_multipoles_mapper,
                 e->s->cells_top, e->s->nr_cells, sizeof(struct cell), 0, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Create and fill the proxies.
 *
 * @param e The #engine.
 */
void engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  const int nodeID = e->nodeID;
  const struct space *s = e->s;
  const int *cdim = s->cdim;
  const int periodic = s->periodic;

  /* Get some info about the physics */
  const double *dim = s->dim;
  const struct gravity_props *props = e->gravity_properties;
  const double theta_crit2 = props->theta_crit2;
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Let's time this */
  const ticks tic = getticks();

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Compute how many cells away we need to walk */
  int delta = 1; /*hydro case */
  if (with_gravity) {
    const double distance = 2.5 * cells[0].width[0] / props->theta_crit;
    delta = (int)(distance / cells[0].width[0]) + 1;
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message("Looking for proxies up to %d top-level cells away", delta);

  /* Loop over each cell in the space. */
  int ind[3];
  for (ind[0] = 0; ind[0] < cdim[0]; ind[0]++) {
    for (ind[1] = 0; ind[1] < cdim[1]; ind[1]++) {
      for (ind[2] = 0; ind[2] < cdim[2]; ind[2]++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, ind[0], ind[1], ind[2]);

        double CoM_i[3] = {0., 0., 0.};
        double r_max_i = 0.;

        if (with_gravity) {

          /* Get ci's multipole */
          const struct gravity_tensors *multi_i = cells[cid].multipole;
          CoM_i[0] = multi_i->CoM[0];
          CoM_i[1] = multi_i->CoM[1];
          CoM_i[2] = multi_i->CoM[2];
          r_max_i = multi_i->r_max;
        }

        /* Loop over all its neighbours (periodic). */
        for (int i = -delta; i <= delta; i++) {
          int ii = ind[0] + i;
          if (ii >= cdim[0])
            ii -= cdim[0];
          else if (ii < 0)
            ii += cdim[0];
          for (int j = -delta; j <= delta; j++) {
            int jj = ind[1] + j;
            if (jj >= cdim[1])
              jj -= cdim[1];
            else if (jj < 0)
              jj += cdim[1];
            for (int k = -delta; k <= delta; k++) {
              int kk = ind[2] + k;
              if (kk >= cdim[2])
                kk -= cdim[2];
              else if (kk < 0)
                kk += cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, ii, jj, kk);

              /* Early abort (same cell) */
              if (cid == cjd) continue;

              /* Early abort (both same node) */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID == nodeID)
                continue;

              /* Early abort (both foreign node) */
              if (cells[cid].nodeID != nodeID && cells[cjd].nodeID != nodeID)
                continue;

              int proxy_type = 0;

              /* In the hydro case, only care about neighbours */
              if (with_hydro) {

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(ind[0] - ii) <= 1 ||
                      abs(ind[0] - ii - cdim[0]) <= 1 ||
                      abs(ind[0] - ii + cdim[0]) <= 1) &&
                     (abs(ind[1] - jj) <= 1 ||
                      abs(ind[1] - jj - cdim[1]) <= 1 ||
                      abs(ind[1] - jj + cdim[1]) <= 1) &&
                     (abs(ind[2] - kk) <= 1 ||
                      abs(ind[2] - kk - cdim[2]) <= 1 ||
                      abs(ind[2] - kk + cdim[2]) <= 1)))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* Get cj's multipole */
                const struct gravity_tensors *multi_j = cells[cjd].multipole;
                const double CoM_j[3] = {multi_j->CoM[0], multi_j->CoM[1],
                                         multi_j->CoM[2]};
                const double r_max_j = multi_j->r_max;

                /* Let's compute the current distance between the cell pair*/
                double dx = CoM_i[0] - CoM_j[0];
                double dy = CoM_i[1] - CoM_j[1];
                double dz = CoM_i[2] - CoM_j[2];

                /* Apply BC */
                if (periodic) {
                  dx = nearest(dx, dim[0]);
                  dy = nearest(dy, dim[1]);
                  dz = nearest(dz, dim[2]);
                }
                const double r2 = dx * dx + dy * dy + dz * dz;

                /* Are we too close for M2L? */
                if (!gravity_M2L_accept(r_max_i, r_max_j, theta_crit2, r2))
                  proxy_type |= (int)proxy_cell_type_gravity;
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int pid = e->proxy_ind[cells[cjd].nodeID];
                if (pid < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cjd].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cjd].nodeID] = e->nr_proxies;
                  pid = e->nr_proxies;
                  e->nr_proxies += 1;
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[pid], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[pid], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                cells[cid].sendto |= (1ULL << pid);
              }

              /* Same for the symmetric case? */
              if (cells[cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int pid = e->proxy_ind[cells[cid].nodeID];
                if (pid < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cid].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
                  pid = e->nr_proxies;
                  e->nr_proxies += 1;
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[pid], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[pid], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cells[cjd].sendto |= (1ULL << pid);
              }
            }
          }
        }
      }
    }
  }

  /* Be clear about the time */
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
  if (e->verbose)
    message("Re-allocating parts array from %zu to %zu.", s->size_parts,
            (size_t)(s->nr_parts * 1.2));
  s->size_parts = s->nr_parts * 1.2;
  struct part *parts_new = NULL;
  struct xpart *xparts_new = NULL;
  if (posix_memalign((void **)&parts_new, part_align,
                     sizeof(struct part) * s->size_parts) != 0 ||
      posix_memalign((void **)&xparts_new, xpart_align,
                     sizeof(struct xpart) * s->size_parts) != 0)
    error("Failed to allocate new part data.");
  if (s->nr_parts > 0) {
    memcpy(parts_new, s->parts, sizeof(struct part) * s->nr_parts);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->nr_parts);
  }
  free(s->parts);
  free(s->xparts);
  s->parts = parts_new;
  s->xparts = xparts_new;

  /* Re-link the gparts to their parts. */
  if (s->nr_parts > 0 && s->nr_gparts > 0)
    part_relink_gparts_to_parts(s->parts, s->nr_parts, 0);

  /* Re-allocate the local sparts. */
  if (e->verbose)
    message("Re-allocating sparts array from %zu to %zu.", s->size_sparts,
            (size_t)(s->nr_sparts * 1.2));
  s->size_sparts = s->nr_sparts * 1.2;
  struct spart *sparts_new = NULL;
  if (posix_memalign((void **)&sparts_new, spart_align,
                     sizeof(struct spart) * s->size_sparts) != 0)
    error("Failed to allocate new spart data.");
  if (s->nr_sparts > 0)
    memcpy(sparts_new, s->sparts, sizeof(struct spart) * s->nr_sparts);
  free(s->sparts);
  s->sparts = sparts_new;

  /* Re-link the gparts to their sparts. */
  if (s->nr_sparts > 0 && s->nr_gparts > 0)
    part_relink_gparts_to_sparts(s->sparts, s->nr_sparts, 0);

  /* Re-allocate the local gparts. */
  if (e->verbose)
    message("Re-allocating gparts array from %zu to %zu.", s->size_gparts,
            (size_t)(s->nr_gparts * 1.2));
  s->size_gparts = s->nr_gparts * 1.2;
  struct gpart *gparts_new = NULL;
  if (posix_memalign((void **)&gparts_new, gpart_align,
                     sizeof(struct gpart) * s->size_gparts) != 0)
    error("Failed to allocate new gpart data.");
  if (s->nr_gparts > 0)
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * s->nr_gparts);
  free(s->gparts);
  s->gparts = gparts_new;

  /* Re-link the parts. */
  if (s->nr_parts > 0 && s->nr_gparts > 0)
    part_relink_parts_to_gparts(s->gparts, s->nr_gparts, s->parts);

  /* Re-link the sparts. */
  if (s->nr_sparts > 0 && s->nr_gparts > 0)
    part_relink_sparts_to_gparts(s->gparts, s->nr_gparts, s->sparts);

#ifdef SWIFT_DEBUG_CHECKS

  /* Verify that the links are correct */
  part_verify_links(s->parts, s->gparts, s->sparts, s->nr_parts, s->nr_gparts,
                    s->nr_sparts, e->verbose);
#endif

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Writes a snapshot with the current state of the engine
 *
 * @param e The #engine.
 */
void engine_dump_snapshot(struct engine *e) {

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);

  /* Be verbose about this */
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Dumping snapshot at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Dumping snapshot at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#else
  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Dumping snapshot at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Dumping snapshot at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

/* Dump... */
#if defined(HAVE_HDF5)
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
  write_output_parallel(e, e->snapshot_base_name, e->internal_units,
                        e->snapshot_units, e->nodeID, e->nr_nodes,
                        MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  write_output_serial(e, e->snapshot_base_name, e->internal_units,
                      e->snapshot_units, e->nodeID, e->nr_nodes, MPI_COMM_WORLD,
                      MPI_INFO_NULL);
#endif
#else
  write_output_single(e, e->snapshot_base_name, e->internal_units,
                      e->snapshot_units);
#endif
#endif

  /* Flag that we dumped a snapshot */
  e->step_props |= engine_step_prop_snapshot;

  clocks_gettime(&time2);
  if (e->verbose)
    message("writing particle properties took %.3f %s.",
            (float)clocks_diff(&time1, &time2), clocks_getunit());
}

#ifdef HAVE_SETAFFINITY
/**
 * @brief Returns the initial affinity the main thread is using.
 */
cpu_set_t *engine_entry_affinity(void) {

  static int use_entry_affinity = 0;
  static cpu_set_t entry_affinity;

  if (!use_entry_affinity) {
    pthread_t engine = pthread_self();
    pthread_getaffinity_np(engine, sizeof(entry_affinity), &entry_affinity);
    use_entry_affinity = 1;
  }

  return &entry_affinity;
}
#endif

/**
 * @brief  Ensure the NUMA node on which we initialise (first touch) everything
 * doesn't change before engine_init allocates NUMA-local workers.
 */
void engine_pin(void) {

#ifdef HAVE_SETAFFINITY
  cpu_set_t *entry_affinity = engine_entry_affinity();
  int pin;
  for (pin = 0; pin < CPU_SETSIZE && !CPU_ISSET(pin, entry_affinity); ++pin)
    ;

  cpu_set_t affinity;
  CPU_ZERO(&affinity);
  CPU_SET(pin, &affinity);
  if (sched_setaffinity(0, sizeof(affinity), &affinity) != 0) {
    error("failed to set engine's affinity");
  }
#else
  error("SWIFT was not compiled with support for pinning.");
#endif
}

/**
 * @brief Unpins the main thread.
 */
void engine_unpin(void) {
#ifdef HAVE_SETAFFINITY
  pthread_t main_thread = pthread_self();
  cpu_set_t *entry_affinity = engine_entry_affinity();
  pthread_setaffinity_np(main_thread, sizeof(*entry_affinity), entry_affinity);
#else
  error("SWIFT was not compiled with support for pinning.");
#endif
}

/**
 * @brief init an engine struct with the necessary properties for the
 *        simulation.
 *
 * Note do not use when restarting. Engine initialisation
 * is completed by a call to engine_config().
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param params The parsed parameter file.
 * @param Ngas total number of gas particles in the simulation.
 * @param Ngparts total number of gravity particles in the simulation.
 * @param Nstars total number of star particles in the simulation.
 * @param policy The queuing policy to use.
 * @param verbose Is this #engine talkative ?
 * @param reparttype What type of repartition algorithm are we using ?
 * @param internal_units The system of units used internally.
 * @param physical_constants The #phys_const used for this run.
 * @param cosmo The #cosmology used for this run.
 * @param hydro The #hydro_props used for this run.
 * @param gravity The #gravity_props used for this run.
 * @param mesh The #pm_mesh used for the long-range periodic forces.
 * @param potential The properties of the external potential.
 * @param cooling_func The properties of the cooling function.
 * @param chemistry The chemistry information.
 * @param sourceterms The properties of the source terms function.
 */
void engine_init(struct engine *e, struct space *s, struct swift_params *params,
                 long long Ngas, long long Ngparts, long long Nstars,
                 int policy, int verbose, struct repartition *reparttype,
                 const struct unit_system *internal_units,
                 const struct phys_const *physical_constants,
                 struct cosmology *cosmo, const struct hydro_props *hydro,
                 struct gravity_props *gravity, struct pm_mesh *mesh,
                 const struct external_potential *potential,
                 const struct cooling_function_data *cooling_func,
                 const struct chemistry_global_data *chemistry,
                 struct sourceterms *sourceterms) {

  /* Clean-up everything */
  bzero(e, sizeof(struct engine));

  /* Store the all values in the fields of the engine. */
  e->s = s;
  e->policy = policy;
  e->step = 0;
  e->total_nr_parts = Ngas;
  e->total_nr_gparts = Ngparts;
  e->total_nr_sparts = Nstars;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->reparttype = reparttype;
  e->ti_old = 0;
  e->ti_current = 0;
  e->time_step = 0.;
  e->time_base = 0.;
  e->time_base_inv = 0.;
  e->time_begin = 0.;
  e->time_end = 0.;
  e->max_active_bin = num_time_bins;
  e->min_active_bin = 1;
  e->internal_units = internal_units;
  e->a_first_snapshot =
      parser_get_opt_param_double(params, "Snapshots:scale_factor_first", 0.1);
  e->time_first_snapshot =
      parser_get_opt_param_double(params, "Snapshots:time_first", 0.);
  e->delta_time_snapshot =
      parser_get_param_double(params, "Snapshots:delta_time");
  e->ti_next_snapshot = 0;
  parser_get_param_string(params, "Snapshots:basename", e->snapshot_base_name);
  e->snapshot_compression =
      parser_get_opt_param_int(params, "Snapshots:compression", 0);
  e->snapshot_label_delta =
      parser_get_opt_param_int(params, "Snapshots:label_delta", 1);
  e->snapshot_units = (struct unit_system *)malloc(sizeof(struct unit_system));
  units_init_default(e->snapshot_units, params, "Snapshots", internal_units);
  e->snapshot_output_count = 0;
  e->dt_min = parser_get_param_double(params, "TimeIntegration:dt_min");
  e->dt_max = parser_get_param_double(params, "TimeIntegration:dt_max");
  e->dt_max_RMS_displacement = FLT_MAX;
  e->max_RMS_displacement_factor = parser_get_opt_param_double(
      params, "TimeIntegration:max_dt_RMS_factor", 0.25);
  e->a_first_statistics =
      parser_get_opt_param_double(params, "Statistics:scale_factor_first", 0.1);
  e->time_first_statistics =
      parser_get_opt_param_double(params, "Statistics:time_first", 0.);
  e->delta_time_statistics =
      parser_get_param_double(params, "Statistics:delta_time");
  e->ti_next_stats = 0;
  e->verbose = verbose;
  e->wallclock_time = 0.f;
  e->physical_constants = physical_constants;
  e->cosmology = cosmo;
  e->hydro_properties = hydro;
  e->gravity_properties = gravity;
  e->mesh = mesh;
  e->external_potential = potential;
  e->cooling_func = cooling_func;
  e->chemistry = chemistry;
  e->sourceterms = sourceterms;
  e->parameter_file = params;
  e->cell_loc = NULL;
#ifdef WITH_MPI
  e->cputime_last_step = 0;
  e->last_repartition = 0;
#endif

  /* Make the space link back to the engine. */
  s->e = e;

  /* Setup the timestep if non-cosmological */
  if (!(e->policy & engine_policy_cosmology)) {
    e->time_begin =
        parser_get_param_double(params, "TimeIntegration:time_begin");
    e->time_end = parser_get_param_double(params, "TimeIntegration:time_end");
    e->time_old = e->time_begin;
    e->time = e->time_begin;

    e->time_base = (e->time_end - e->time_begin) / max_nr_timesteps;
    e->time_base_inv = 1.0 / e->time_base;
    e->ti_current = 0;
  } else {

    e->time_begin = e->cosmology->time_begin;
    e->time_end = e->cosmology->time_end;
    e->time_old = e->time_begin;
    e->time = e->time_begin;

    /* Copy the relevent information from the cosmology model */
    e->time_base = e->cosmology->time_base;
    e->time_base_inv = e->cosmology->time_base_inv;
    e->ti_current = 0;
  }

  engine_init_output_lists(e, params);
}

/**
 * @brief configure an engine with the given number of threads, queues
 *        and core affinity. Also initialises the scheduler and opens various
 *        output files, computes the next timestep and initialises the
 *        threadpool.
 *
 * Assumes the engine is correctly initialised i.e. is restored from a restart
 * file or has been setup by engine_init(). When restarting any output log
 * files are positioned so that further output is appended. Note that
 * parameters are not read from the engine, just the parameter file, this
 * allows values derived in this function to be changed between runs.
 * When not restarting params should be the same as given to engine_init().
 *
 * @param restart true when restarting the application.
 * @param e The #engine.
 * @param params The parsed parameter file.
 * @param nr_nodes The number of MPI ranks.
 * @param nodeID The MPI rank of this node.
 * @param nr_threads The number of threads per MPI rank.
 * @param with_aff use processor affinity, if supported.
 * @param verbose Is this #engine talkative ?
 * @param restart_file The name of our restart file.
 */
void engine_config(int restart, struct engine *e, struct swift_params *params,
                   int nr_nodes, int nodeID, int nr_threads, int with_aff,
                   int verbose, const char *restart_file) {

  /* Store the values and initialise global fields. */
  e->nodeID = nodeID;
  e->nr_threads = nr_threads;
  e->nr_nodes = nr_nodes;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->forcerebuild = 1;
  e->forcerepart = 0;
  e->restarting = restart;
  e->step_props = engine_step_prop_none;
  e->links = NULL;
  e->nr_links = 0;
  e->file_stats = NULL;
  e->file_timesteps = NULL;
  e->verbose = verbose;
  e->wallclock_time = 0.f;
  e->restart_dump = 0;
  e->restart_file = restart_file;
  e->restart_next = 0;
  e->restart_dt = 0;
  e->timeFirstSTFOutput = 0;
  engine_rank = nodeID;

  /* Initialise VELOCIraptor. */
  if (e->policy & engine_policy_structure_finding) {
    parser_get_param_string(params, "StructureFinding:basename",
                            e->stfBaseName);
    e->timeFirstSTFOutput =
        parser_get_param_double(params, "StructureFinding:time_first");
    e->a_first_stf = parser_get_opt_param_double(
        params, "StructureFinding:scale_factor_first", 0.1);
    e->stf_output_freq_format =
        parser_get_param_int(params, "StructureFinding:output_time_format");
    if (e->stf_output_freq_format == STEPS) {
      e->deltaStepSTF =
          parser_get_param_int(params, "StructureFinding:delta_step");
    } else if (e->stf_output_freq_format == TIME) {
      e->deltaTimeSTF =
          parser_get_param_double(params, "StructureFinding:delta_time");
    } else
      error(
          "Invalid flag (%d) set for output time format of structure finding.",
          e->stf_output_freq_format);

    /* overwrite input if outputlist */
    if (e->output_list_stf) e->stf_output_freq_format = TIME;
  }

  /* Get the number of queues */
  int nr_queues =
      parser_get_opt_param_int(params, "Scheduler:nr_queues", nr_threads);
  if (nr_queues <= 0) nr_queues = e->nr_threads;
  if (nr_queues != nr_threads)
    message("Number of task queues set to %d", nr_queues);
  e->s->nr_queues = nr_queues;

/* Deal with affinity. For now, just figure out the number of cores. */
#if defined(HAVE_SETAFFINITY)
  const int nr_cores = sysconf(_SC_NPROCESSORS_ONLN);
  cpu_set_t *entry_affinity = engine_entry_affinity();
  const int nr_affinity_cores = CPU_COUNT(entry_affinity);

  if (nr_cores > CPU_SETSIZE) /* Unlikely, except on e.g. SGI UV. */
    error("must allocate dynamic cpu_set_t (too many cores per node)");

  char *buf = (char *)malloc((nr_cores + 1) * sizeof(char));
  buf[nr_cores] = '\0';
  for (int j = 0; j < nr_cores; ++j) {
    /* Reversed bit order from convention, but same as e.g. Intel MPI's
     * I_MPI_PIN_DOMAIN explicit mask: left-to-right, LSB-to-MSB. */
    buf[j] = CPU_ISSET(j, entry_affinity) ? '1' : '0';
  }

  if (verbose && with_aff) message("Affinity at entry: %s", buf);

  int *cpuid = NULL;
  cpu_set_t cpuset;

  if (with_aff) {

    cpuid = (int *)malloc(nr_affinity_cores * sizeof(int));

    int skip = 0;
    for (int k = 0; k < nr_affinity_cores; k++) {
      int c;
      for (c = skip; c < CPU_SETSIZE && !CPU_ISSET(c, entry_affinity); ++c)
        ;
      cpuid[k] = c;
      skip = c + 1;
    }

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
    if ((e->policy & engine_policy_cputight) != engine_policy_cputight) {

      if (numa_available() >= 0) {
        if (nodeID == 0) message("prefer NUMA-distant CPUs");

        /* Get list of numa nodes of all available cores. */
        int *nodes = (int *)malloc(nr_affinity_cores * sizeof(int));
        int nnodes = 0;
        for (int i = 0; i < nr_affinity_cores; i++) {
          nodes[i] = numa_node_of_cpu(cpuid[i]);
          if (nodes[i] > nnodes) nnodes = nodes[i];
        }
        nnodes += 1;

        /* Count cores per node. */
        int *core_counts = (int *)malloc(nnodes * sizeof(int));
        for (int i = 0; i < nr_affinity_cores; i++) {
          core_counts[nodes[i]] = 0;
        }
        for (int i = 0; i < nr_affinity_cores; i++) {
          core_counts[nodes[i]] += 1;
        }

        /* Index cores within each node. */
        int *core_indices = (int *)malloc(nr_affinity_cores * sizeof(int));
        for (int i = nr_affinity_cores - 1; i >= 0; i--) {
          core_indices[i] = core_counts[nodes[i]];
          core_counts[nodes[i]] -= 1;
        }

        /* Now sort so that we pick adjacent cpuids from different nodes
         * by sorting internal node core indices. */
        int done = 0;
        while (!done) {
          done = 1;
          for (int i = 1; i < nr_affinity_cores; i++) {
            if (core_indices[i] < core_indices[i - 1]) {
              int t = cpuid[i - 1];
              cpuid[i - 1] = cpuid[i];
              cpuid[i] = t;

              t = core_indices[i - 1];
              core_indices[i - 1] = core_indices[i];
              core_indices[i] = t;
              done = 0;
            }
          }
        }

        free(nodes);
        free(core_counts);
        free(core_indices);
      }
    }
#endif
  } else {
    if (nodeID == 0) message("no processor affinity used");

  } /* with_aff */

  /* Avoid (unexpected) interference between engine and runner threads. We can
   * do this once we've made at least one call to engine_entry_affinity and
   * maybe numa_node_of_cpu(sched_getcpu()), even if the engine isn't already
   * pinned. */
  if (with_aff) engine_unpin();
#endif

  if (with_aff && nodeID == 0) {
#ifdef HAVE_SETAFFINITY
#ifdef WITH_MPI
    printf("[%04i] %s engine_init: cpu map is [ ", nodeID,
           clocks_get_timesincestart());
#else
    printf("%s engine_init: cpu map is [ ", clocks_get_timesincestart());
#endif
    for (int i = 0; i < nr_affinity_cores; i++) printf("%i ", cpuid[i]);
    printf("].\n");
#endif
  }

  /* Are we doing stuff in parallel? */
  if (nr_nodes > 1) {
#ifndef WITH_MPI
    error("SWIFT was not compiled with MPI support.");
#else
    e->policy |= engine_policy_mpi;
    if ((e->proxies = (struct proxy *)calloc(sizeof(struct proxy),
                                             engine_maxproxies)) == NULL)
      error("Failed to allocate memory for proxies.");
    e->nr_proxies = 0;
#endif
  }

  /* Open some files */
  if (e->nodeID == 0) {

    /* When restarting append to these files. */
    const char *mode;
    if (restart)
      mode = "a";
    else
      mode = "w";

    char energyfileName[200] = "";
    parser_get_opt_param_string(params, "Statistics:energy_file_name",
                                energyfileName,
                                engine_default_energy_file_name);
    sprintf(energyfileName + strlen(energyfileName), ".txt");
    e->file_stats = fopen(energyfileName, mode);

    if (!restart) {
      fprintf(
          e->file_stats,
          "#%14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s %14s "
          "%14s %14s %14s %14s %14s %14s\n",
          "Time", "Mass", "E_tot", "E_kin", "E_int", "E_pot", "E_pot_self",
          "E_pot_ext", "E_radcool", "Entropy", "p_x", "p_y", "p_z", "ang_x",
          "ang_y", "ang_z", "com_x", "com_y", "com_z");
      fflush(e->file_stats);
    }

    char timestepsfileName[200] = "";
    parser_get_opt_param_string(params, "Statistics:timestep_file_name",
                                timestepsfileName,
                                engine_default_timesteps_file_name);

    sprintf(timestepsfileName + strlen(timestepsfileName), "_%d.txt",
            nr_nodes * nr_threads);
    e->file_timesteps = fopen(timestepsfileName, mode);

    if (!restart) {
      fprintf(
          e->file_timesteps,
          "# Host: %s\n# Branch: %s\n# Revision: %s\n# Compiler: %s, "
          "Version: %s \n# "
          "Number of threads: %d\n# Number of MPI ranks: %d\n# Hydrodynamic "
          "scheme: %s\n# Hydrodynamic kernel: %s\n# No. of neighbours: %.2f "
          "+/- %.4f\n# Eta: %f\n# Config: %s\n# CFLAGS: %s\n",
          hostname(), git_branch(), git_revision(), compiler_name(),
          compiler_version(), e->nr_threads, e->nr_nodes, SPH_IMPLEMENTATION,
          kernel_name, e->hydro_properties->target_neighbours,
          e->hydro_properties->delta_neighbours,
          e->hydro_properties->eta_neighbours, configuration_options(),
          compilation_cflags());

      fprintf(e->file_timesteps,
              "# Step Properties: Rebuild=%d, Redistribute=%d, Repartition=%d, "
              "Statistics=%d, Snapshot=%d, Restarts=%d\n",
              engine_step_prop_rebuild, engine_step_prop_redistribute,
              engine_step_prop_repartition, engine_step_prop_statistics,
              engine_step_prop_snapshot, engine_step_prop_restarts);

      fprintf(e->file_timesteps,
              "# %6s %14s %14s %10s %14s %9s %12s %12s %12s %16s [%s] %6s\n",
              "Step", "Time", "Scale-factor", "Redshift", "Time-step",
              "Time-bins", "Updates", "g-Updates", "s-Updates",
              "Wall-clock time", clocks_getunit(), "Props");
      fflush(e->file_timesteps);
    }
  }

  /* Print policy */
  engine_print_policy(e);

  /* Print information about the hydro scheme */
  if (e->policy & engine_policy_hydro)
    if (e->nodeID == 0) hydro_props_print(e->hydro_properties);

  /* Print information about the gravity scheme */
  if (e->policy & engine_policy_self_gravity)
    if (e->nodeID == 0) gravity_props_print(e->gravity_properties);

  /* Check we have sensible time bounds */
  if (e->time_begin >= e->time_end)
    error(
        "Final simulation time (t_end = %e) must be larger than the start time "
        "(t_beg = %e)",
        e->time_end, e->time_begin);

  /* Check we have sensible time-step values */
  if (e->dt_min > e->dt_max)
    error(
        "Minimal time-step size (%e) must be smaller than maximal time-step "
        "size (%e)",
        e->dt_min, e->dt_max);

  /* Info about time-steps */
  if (e->nodeID == 0) {
    message("Absolute minimal timestep size: %e", e->time_base);

    float dt_min = e->time_end - e->time_begin;
    while (dt_min > e->dt_min) dt_min /= 2.f;

    message("Minimal timestep size (on time-line): %e", dt_min);

    float dt_max = e->time_end - e->time_begin;
    while (dt_max > e->dt_max) dt_max /= 2.f;

    message("Maximal timestep size (on time-line): %e", dt_max);
  }

  if (e->dt_min < e->time_base && e->nodeID == 0)
    error(
        "Minimal time-step size smaller than the absolute possible minimum "
        "dt=%e",
        e->time_base);

  if (!(e->policy & engine_policy_cosmology))
    if (e->dt_max > (e->time_end - e->time_begin) && e->nodeID == 0)
      error("Maximal time-step size larger than the simulation run time t=%e",
            e->time_end - e->time_begin);

  /* Deal with outputs */
  if (e->policy & engine_policy_cosmology) {

    if (e->delta_time_snapshot <= 1.)
      error("Time between snapshots (%e) must be > 1.", e->delta_time_snapshot);

    if (e->delta_time_statistics <= 1.)
      error("Time between statistics (%e) must be > 1.",
            e->delta_time_statistics);

    if (e->a_first_snapshot < e->cosmology->a_begin)
      error(
          "Scale-factor of first snapshot (%e) must be after the simulation "
          "start a=%e.",
          e->a_first_snapshot, e->cosmology->a_begin);

    if (e->a_first_statistics < e->cosmology->a_begin)
      error(
          "Scale-factor of first stats output (%e) must be after the "
          "simulation start a=%e.",
          e->a_first_statistics, e->cosmology->a_begin);

    if ((e->policy & engine_policy_structure_finding) &&
        (e->stf_output_freq_format == TIME)) {

      if (e->deltaTimeSTF <= 1.)
        error("Time between STF (%e) must be > 1.", e->deltaTimeSTF);

      if (e->a_first_stf < e->cosmology->a_begin)
        error(
            "Scale-factor of first stf output (%e) must be after the "
            "simulation "
            "start a=%e.",
            e->a_first_stf, e->cosmology->a_begin);
    }
  } else {

    if (e->delta_time_snapshot <= 0.)
      error("Time between snapshots (%e) must be positive.",
            e->delta_time_snapshot);

    if (e->delta_time_statistics <= 0.)
      error("Time between statistics (%e) must be positive.",
            e->delta_time_statistics);

    /* Find the time of the first output */
    if (e->time_first_snapshot < e->time_begin)
      error(
          "Time of first snapshot (%e) must be after the simulation start "
          "t=%e.",
          e->time_first_snapshot, e->time_begin);

    if (e->time_first_statistics < e->time_begin)
      error(
          "Time of first stats output (%e) must be after the simulation start "
          "t=%e.",
          e->time_first_statistics, e->time_begin);

    if ((e->policy & engine_policy_structure_finding) &&
        (e->stf_output_freq_format == TIME)) {

      if (e->deltaTimeSTF <= 0.)
        error("Time between STF (%e) must be positive.", e->deltaTimeSTF);

      if (e->timeFirstSTFOutput < e->time_begin)
        error("Time of first STF (%e) must be after the simulation start t=%e.",
              e->timeFirstSTFOutput, e->time_begin);
    }
  }

  if (e->policy & engine_policy_structure_finding) {
    /* Find the time of the first stf output */
    if (e->stf_output_freq_format == TIME) {
      engine_compute_next_stf_time(e);
      message("Next STF step will be: %lld", e->ti_nextSTF);
    }
  }

  /* Get the total mass */
  e->total_mass = 0.;
  for (size_t i = 0; i < e->s->nr_gparts; ++i)
    e->total_mass += e->s->gparts[i].mass;

/* Reduce the total mass */
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->total_mass, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  /* Find the time of the first snapshot  output */
  engine_compute_next_snapshot_time(e);

  /* Find the time of the first statistics output */
  engine_compute_next_statistics_time(e);

  /* Whether restarts are enabled. Yes by default. Can be changed on restart. */
  e->restart_dump = parser_get_opt_param_int(params, "Restarts:enable", 1);

  /* Whether to save backup copies of the previous restart files. */
  e->restart_save = parser_get_opt_param_int(params, "Restarts:save", 1);

  /* Whether restarts should be dumped on exit. Not by default. Can be changed
   * on restart. */
  e->restart_onexit = parser_get_opt_param_int(params, "Restarts:onexit", 0);

  /* Hours between restart dumps. Can be changed on restart. */
  float dhours =
      parser_get_opt_param_float(params, "Restarts:delta_hours", 6.0);
  if (e->nodeID == 0) {
    if (e->restart_dump)
      message("Restarts will be dumped every %f hours", dhours);
    else
      message("WARNING: restarts will not be dumped");

    if (e->verbose && e->restart_onexit)
      message("Restarts will be dumped after the final step");
  }

  /* Internally we use ticks, so convert into a delta ticks. Assumes we can
   * convert from ticks into milliseconds. */
  e->restart_dt = clocks_to_ticks(dhours * 60.0 * 60.0 * 1000.0);

  /* The first dump will happen no sooner than restart_dt ticks in the
   * future. */
  e->restart_next = getticks() + e->restart_dt;

/* Construct types for MPI communications */
#ifdef WITH_MPI
  part_create_mpi_types();
  stats_create_MPI_type();
#endif

  /* Initialise the collection group. */
  collectgroup_init();

  /* Initialize the threadpool. */
  threadpool_init(&e->threadpool, e->nr_threads);

  /* First of all, init the barrier and lock it. */
  if (swift_barrier_init(&e->wait_barrier, NULL, e->nr_threads + 1) != 0 ||
      swift_barrier_init(&e->run_barrier, NULL, e->nr_threads + 1) != 0)
    error("Failed to initialize barrier.");

  /* Expected average for tasks per cell. If set to zero we use a heuristic
   * guess based on the numbers of cells and how many tasks per cell we expect.
   * On restart this number cannot be estimated (no cells yet), so we recover
   * from the end of the dumped run. Can be changed on restart.
   */
  e->tasks_per_cell =
      parser_get_opt_param_int(params, "Scheduler:tasks_per_cell", 0);
  int maxtasks = 0;
  if (restart)
    maxtasks = e->restart_max_tasks;
  else
    maxtasks = engine_estimate_nr_tasks(e);

  /* Init the scheduler. */
  scheduler_init(&e->sched, e->s, maxtasks, nr_queues,
                 (e->policy & scheduler_flag_steal), e->nodeID, &e->threadpool);

  /* Maximum size of MPI task messages, in KB, that should not be buffered,
   * that is sent using MPI_Issend, not MPI_Isend. 4Mb by default. Can be
   * changed on restart.
   */
  e->sched.mpi_message_limit =
      parser_get_opt_param_int(params, "Scheduler:mpi_message_limit", 4) * 1024;

  /* Allocate and init the threads. */
  if (posix_memalign((void **)&e->runners, SWIFT_CACHE_ALIGNMENT,
                     e->nr_threads * sizeof(struct runner)) != 0)
    error("Failed to allocate threads array.");
  for (int k = 0; k < e->nr_threads; k++) {
    e->runners[k].id = k;
    e->runners[k].e = e;
    if (pthread_create(&e->runners[k].thread, NULL, &runner_main,
                       &e->runners[k]) != 0)
      error("Failed to create runner thread.");

    /* Try to pin the runner to a given core */
    if (with_aff &&
        (e->policy & engine_policy_setaffinity) == engine_policy_setaffinity) {
#if defined(HAVE_SETAFFINITY)

      /* Set a reasonable queue ID. */
      int coreid = k % nr_affinity_cores;
      e->runners[k].cpuid = cpuid[coreid];

      if (nr_queues < e->nr_threads)
        e->runners[k].qid = cpuid[coreid] * nr_queues / nr_affinity_cores;
      else
        e->runners[k].qid = k;

      /* Set the cpu mask to zero | e->id. */
      CPU_ZERO(&cpuset);
      CPU_SET(cpuid[coreid], &cpuset);

      /* Apply this mask to the runner's pthread. */
      if (pthread_setaffinity_np(e->runners[k].thread, sizeof(cpu_set_t),
                                 &cpuset) != 0)
        error("Failed to set thread affinity.");

#else
      error("SWIFT was not compiled with affinity enabled.");
#endif
    } else {
      e->runners[k].cpuid = k;
      e->runners[k].qid = k * nr_queues / e->nr_threads;
    }

    /* Allocate particle caches. */
    e->runners[k].ci_gravity_cache.count = 0;
    e->runners[k].cj_gravity_cache.count = 0;
    gravity_cache_init(&e->runners[k].ci_gravity_cache, space_splitsize);
    gravity_cache_init(&e->runners[k].cj_gravity_cache, space_splitsize);
#ifdef WITH_VECTORIZATION
    e->runners[k].ci_cache.count = 0;
    e->runners[k].cj_cache.count = 0;
    cache_init(&e->runners[k].ci_cache, CACHE_SIZE);
    cache_init(&e->runners[k].cj_cache, CACHE_SIZE);
#endif

    if (verbose) {
      if (with_aff)
        message("runner %i on cpuid=%i with qid=%i.", e->runners[k].id,
                e->runners[k].cpuid, e->runners[k].qid);
      else
        message("runner %i using qid=%i no cpuid.", e->runners[k].id,
                e->runners[k].qid);
    }
  }

/* Free the affinity stuff */
#if defined(HAVE_SETAFFINITY)
  if (with_aff) {
    free(cpuid);
  }
  free(buf);
#endif

  /* Wait for the runner threads to be in place. */
  swift_barrier_wait(&e->wait_barrier);
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
    for (int k = 0; k <= engine_maxpolicy; k++)
      if (e->policy & (1 << k)) printf(" '%s' ", engine_policy_names[k + 1]);
    printf(" ]\n");
    fflush(stdout);
  }
#else
  printf("%s engine_policy: engine policies are [ ",
         clocks_get_timesincestart());
  for (int k = 0; k <= engine_maxpolicy; k++)
    if (e->policy & (1 << k)) printf(" '%s' ", engine_policy_names[k + 1]);
  printf(" ]\n");
  fflush(stdout);
#endif
}

/**
 * @brief Computes the next time (on the time line) for a dump
 *
 * @param e The #engine.
 */
void engine_compute_next_snapshot_time(struct engine *e) {
  /* Do outputlist file case */
  if (e->output_list_snapshots) {
    output_list_read_next_time(e->output_list_snapshots, e, "snapshots",
                               &e->ti_next_snapshot);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_snapshot;
  else
    time_end = e->time_end + e->delta_time_snapshot;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_snapshot;
  else
    time = e->time_first_snapshot;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_snapshot = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_snapshot = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_snapshot > e->ti_current) break;

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_snapshot;
    else
      time += e->delta_time_snapshot;
  }

  /* Deal with last snapshot */
  if (e->ti_next_snapshot >= max_nr_timesteps) {
    e->ti_next_snapshot = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_snapshot_time =
          exp(e->ti_next_snapshot * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next snapshot time set to a=%e.", next_snapshot_time);
    } else {
      const double next_snapshot_time =
          e->ti_next_snapshot * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next snapshot time set to t=%e.", next_snapshot_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for a statistics dump
 *
 * @param e The #engine.
 */
void engine_compute_next_statistics_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_stats) {
    output_list_read_next_time(e->output_list_stats, e, "stats",
                               &e->ti_next_stats);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_statistics;
  else
    time_end = e->time_end + e->delta_time_statistics;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_statistics;
  else
    time = e->time_first_statistics;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_stats = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_stats = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_stats > e->ti_current) break;

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_statistics;
    else
      time += e->delta_time_statistics;
  }

  /* Deal with last statistics */
  if (e->ti_next_stats >= max_nr_timesteps) {
    e->ti_next_stats = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_statistics_time =
          exp(e->ti_next_stats * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time for stats set to a=%e.",
                next_statistics_time);
    } else {
      const double next_statistics_time =
          e->ti_next_stats * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time for stats set to t=%e.",
                next_statistics_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for structure finding
 *
 * @param e The #engine.
 */
void engine_compute_next_stf_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_stf) {
    output_list_read_next_time(e->output_list_stf, e, "stf", &e->ti_nextSTF);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->deltaTimeSTF;
  else
    time_end = e->time_end + e->deltaTimeSTF;

  /* Find next snasphot above current time */
  double time = e->timeFirstSTFOutput;

  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_nextSTF = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_nextSTF = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_nextSTF > e->ti_current) break;

    if (e->policy & engine_policy_cosmology)
      time *= e->deltaTimeSTF;
    else
      time += e->deltaTimeSTF;
  }

  /* Deal with last snapshot */
  if (e->ti_nextSTF >= max_nr_timesteps) {
    e->ti_nextSTF = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const float next_snapshot_time =
          exp(e->ti_nextSTF * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time set to a=%e.", next_snapshot_time);
    } else {
      const float next_snapshot_time =
          e->ti_nextSTF * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time set to t=%e.", next_snapshot_time);
    }
  }
}

/**
 * @brief Initialize all the output_list required by the engine
 *
 * @param e The #engine.
 * @param params The #swift_params.
 */
void engine_init_output_lists(struct engine *e, struct swift_params *params) {
  /* Deal with snapshots */
  double snaps_time_first;
  e->output_list_snapshots = NULL;
  output_list_init(&e->output_list_snapshots, e, "Snapshots",
                   &e->delta_time_snapshot, &snaps_time_first);

  if (e->output_list_snapshots) {
    if (e->policy & engine_policy_cosmology)
      e->a_first_snapshot = snaps_time_first;
    else
      e->time_first_snapshot = snaps_time_first;
  }

  /* Deal with stats */
  double stats_time_first;
  e->output_list_stats = NULL;
  output_list_init(&e->output_list_stats, e, "Statistics",
                   &e->delta_time_statistics, &stats_time_first);

  if (e->output_list_stats) {
    if (e->policy & engine_policy_cosmology)
      e->a_first_statistics = stats_time_first;
    else
      e->time_first_statistics = stats_time_first;
  }

  /* Deal with stf */
  double stf_time_first;
  e->output_list_stf = NULL;
  output_list_init(&e->output_list_stf, e, "StructureFinding", &e->deltaTimeSTF,
                   &stf_time_first);

  if (e->output_list_stf) {
    if (e->policy & engine_policy_cosmology)
      e->a_first_stf = stf_time_first;
    else
      e->timeFirstSTFOutput = stf_time_first;
  }
}

/**
 * @brief Computes the maximal time-step allowed by the max RMS displacement
 * condition.
 *
 * @param e The #engine.
 */
void engine_recompute_displacement_constraint(struct engine *e) {

  /* Get the cosmological information */
  const struct cosmology *cosmo = e->cosmology;
  const float Om = cosmo->Omega_m;
  const float Ob = cosmo->Omega_b;
  const float rho_crit = cosmo->critical_density;
  const float a = cosmo->a;

  /* Start by reducing the minimal mass of each particle type */
  float min_mass[swift_type_count] = {e->s->min_part_mass,
                                      e->s->min_gpart_mass,
                                      FLT_MAX,
                                      FLT_MAX,
                                      e->s->min_spart_mass,
                                      FLT_MAX};
#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the minimal mass collection worked */
  float min_part_mass_check = FLT_MAX;
  for (size_t i = 0; i < e->s->nr_parts; ++i)
    min_part_mass_check =
        min(min_part_mass_check, hydro_get_mass(&e->s->parts[i]));
  if (min_part_mass_check != min_mass[swift_type_gas])
    error("Error collecting minimal mass of gas particles.");
#endif

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, min_mass, swift_type_count, MPI_FLOAT, MPI_MIN,
                MPI_COMM_WORLD);
#endif

  /* Do the same for the velocity norm sum */
  float vel_norm[swift_type_count] = {e->s->sum_part_vel_norm,
                                      e->s->sum_gpart_vel_norm,
                                      0.f,
                                      0.f,
                                      e->s->sum_spart_vel_norm,
                                      0.f};
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, vel_norm, swift_type_count, MPI_FLOAT, MPI_SUM,
                MPI_COMM_WORLD);
#endif

  /* Get the counts of each particle types */
  const long long total_nr_dm_gparts =
      e->total_nr_gparts - e->total_nr_parts - e->total_nr_sparts;
  float count_parts[swift_type_count] = {(float)e->total_nr_parts,
                                         (float)total_nr_dm_gparts,
                                         0.f,
                                         0.f,
                                         (float)e->total_nr_sparts,
                                         0.f};

  /* Count of particles for the two species */
  const float N_dm = count_parts[1];
  const float N_b = count_parts[0] + count_parts[4];

  /* Peculiar motion norm for the two species */
  const float vel_norm_dm = vel_norm[1];
  const float vel_norm_b = vel_norm[0] + vel_norm[4];

  /* Mesh forces smoothing scale */
  float r_s;
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic == 1)
    r_s = e->mesh->r_s;
  else
    r_s = FLT_MAX;

  float dt_dm = FLT_MAX, dt_b = FLT_MAX;

  /* DM case */
  if (N_dm > 0.f) {

    /* Minimal mass for the DM */
    const float min_mass_dm = min_mass[1];

    /* Inter-particle sepration for the DM */
    const float d_dm = cbrtf(min_mass_dm / ((Om - Ob) * rho_crit));

    /* RMS peculiar motion for the DM */
    const float rms_vel_dm = vel_norm_dm / N_dm;

    /* Time-step based on maximum displacement */
    dt_dm = a * a * min(r_s, d_dm) / sqrtf(rms_vel_dm);
  }

  /* Baryon case */
  if (N_b > 0.f) {

    /* Minimal mass for the baryons */
    const float min_mass_b = min(min_mass[0], min_mass[4]);

    /* Inter-particle sepration for the baryons */
    const float d_b = cbrtf(min_mass_b / (Ob * rho_crit));

    /* RMS peculiar motion for the baryons */
    const float rms_vel_b = vel_norm_b / N_b;

    /* Time-step based on maximum displacement */
    dt_b = a * a * min(r_s, d_b) / sqrtf(rms_vel_b);
  }

  /* Use the minimum */
  const float dt = min(dt_dm, dt_b);

  /* Apply the dimensionless factor */
  e->dt_max_RMS_displacement = dt * e->max_RMS_displacement_factor;

  if (e->verbose)
    message("max_dt_RMS_displacement = %e", e->dt_max_RMS_displacement);
}

/**
 * @brief Frees up the memory allocated for this #engine
 */
void engine_clean(struct engine *e) {

  for (int i = 0; i < e->nr_threads; ++i) {
#ifdef WITH_VECTORIZATION
    cache_clean(&e->runners[i].ci_cache);
    cache_clean(&e->runners[i].cj_cache);
#endif
    gravity_cache_clean(&e->runners[i].ci_gravity_cache);
    gravity_cache_clean(&e->runners[i].cj_gravity_cache);
  }
  free(e->runners);
  free(e->snapshot_units);
  if (e->output_list_snapshots) {
    output_list_clean(e->output_list_snapshots);
    free(e->output_list_snapshots);
  }
  if (e->output_list_stats) {
    output_list_clean(e->output_list_stats);
    free(e->output_list_stats);
  }
  if (e->output_list_stf) {
    output_list_clean(e->output_list_stf);
    free(e->output_list_stf);
  }
  free(e->links);
  free(e->cell_loc);
  scheduler_clean(&e->sched);
  space_clean(e->s);
  threadpool_clean(&e->threadpool);
}

/**
 * @brief Write the engine struct and its contents to the given FILE as a
 * stream of bytes.
 *
 * @param e the engine
 * @param stream the file stream
 */
void engine_struct_dump(struct engine *e, FILE *stream) {

  /* Dump the engine. Save the current tasks_per_cell estimate. */
  e->restart_max_tasks = engine_estimate_nr_tasks(e);
  restart_write_blocks(e, sizeof(struct engine), 1, stream, "engine",
                       "engine struct");

  /* And all the engine pointed data, these use their own dump functions. */
  space_struct_dump(e->s, stream);
  units_struct_dump(e->internal_units, stream);
  units_struct_dump(e->snapshot_units, stream);
  cosmology_struct_dump(e->cosmology, stream);

#ifdef WITH_MPI
  /* Save the partition for restoration. */
  partition_store_celllist(e->s, e->reparttype);
  partition_struct_dump(e->reparttype, stream);
#endif

  phys_const_struct_dump(e->physical_constants, stream);
  hydro_props_struct_dump(e->hydro_properties, stream);
  gravity_props_struct_dump(e->gravity_properties, stream);
  pm_mesh_struct_dump(e->mesh, stream);
  potential_struct_dump(e->external_potential, stream);
  cooling_struct_dump(e->cooling_func, stream);
  chemistry_struct_dump(e->chemistry, stream);
  sourceterms_struct_dump(e->sourceterms, stream);
  parser_struct_dump(e->parameter_file, stream);
  if (e->output_list_snapshots)
    output_list_struct_dump(e->output_list_snapshots, stream);
  if (e->output_list_stats)
    output_list_struct_dump(e->output_list_stats, stream);
  if (e->output_list_stf) output_list_struct_dump(e->output_list_stf, stream);
}

/**
 * @brief Re-create an engine struct and its contents from the given FILE
 *        stream.
 *
 * @param e the engine
 * @param stream the file stream
 */
void engine_struct_restore(struct engine *e, FILE *stream) {

  /* Read the engine. */
  restart_read_blocks(e, sizeof(struct engine), 1, stream, NULL,
                      "engine struct");

  /* Re-initializations as necessary for our struct and its members. */
  e->sched.tasks = NULL;
  e->sched.tasks_ind = NULL;
  e->sched.tid_active = NULL;
  e->sched.size = 0;

  /* Now for the other pointers, these use their own restore functions. */
  /* Note all this memory leaks, but is used once. */
  struct space *s = (struct space *)malloc(sizeof(struct space));
  space_struct_restore(s, stream);
  e->s = s;
  s->e = e;

  struct unit_system *us =
      (struct unit_system *)malloc(sizeof(struct unit_system));
  units_struct_restore(us, stream);
  e->internal_units = us;

  us = (struct unit_system *)malloc(sizeof(struct unit_system));
  units_struct_restore(us, stream);
  e->snapshot_units = us;

  struct cosmology *cosmo =
      (struct cosmology *)malloc(sizeof(struct cosmology));
  cosmology_struct_restore(e->policy & engine_policy_cosmology, cosmo, stream);
  e->cosmology = cosmo;

#ifdef WITH_MPI
  struct repartition *reparttype =
      (struct repartition *)malloc(sizeof(struct repartition));
  partition_struct_restore(reparttype, stream);
  e->reparttype = reparttype;
#endif

  struct phys_const *physical_constants =
      (struct phys_const *)malloc(sizeof(struct phys_const));
  phys_const_struct_restore(physical_constants, stream);
  e->physical_constants = physical_constants;

  struct hydro_props *hydro_properties =
      (struct hydro_props *)malloc(sizeof(struct hydro_props));
  hydro_props_struct_restore(hydro_properties, stream);
  e->hydro_properties = hydro_properties;

  struct gravity_props *gravity_properties =
      (struct gravity_props *)malloc(sizeof(struct gravity_props));
  gravity_props_struct_restore(gravity_properties, stream);
  e->gravity_properties = gravity_properties;

  struct pm_mesh *mesh = (struct pm_mesh *)malloc(sizeof(struct pm_mesh));
  pm_mesh_struct_restore(mesh, stream);
  e->mesh = mesh;

  struct external_potential *external_potential =
      (struct external_potential *)malloc(sizeof(struct external_potential));
  potential_struct_restore(external_potential, stream);
  e->external_potential = external_potential;

  struct cooling_function_data *cooling_func =
      (struct cooling_function_data *)malloc(
          sizeof(struct cooling_function_data));
  cooling_struct_restore(cooling_func, stream);
  e->cooling_func = cooling_func;

  struct chemistry_global_data *chemistry =
      (struct chemistry_global_data *)malloc(
          sizeof(struct chemistry_global_data));
  chemistry_struct_restore(chemistry, stream);
  e->chemistry = chemistry;

  struct sourceterms *sourceterms =
      (struct sourceterms *)malloc(sizeof(struct sourceterms));
  sourceterms_struct_restore(sourceterms, stream);
  e->sourceterms = sourceterms;

  struct swift_params *parameter_file =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_struct_restore(parameter_file, stream);
  e->parameter_file = parameter_file;

  if (e->output_list_snapshots) {
    struct output_list *output_list_snapshots =
        (struct output_list *)malloc(sizeof(struct output_list));
    output_list_struct_restore(output_list_snapshots, stream);
    e->output_list_snapshots = output_list_snapshots;
  }

  if (e->output_list_stats) {
    struct output_list *output_list_stats =
        (struct output_list *)malloc(sizeof(struct output_list));
    output_list_struct_restore(output_list_stats, stream);
    e->output_list_stats = output_list_stats;
  }

  if (e->output_list_stf) {
    struct output_list *output_list_stf =
        (struct output_list *)malloc(sizeof(struct output_list));
    output_list_struct_restore(output_list_stf, stream);
    e->output_list_stf = output_list_stf;
  }

  /* Want to force a rebuild before using this engine. Wait to repartition.*/
  e->forcerebuild = 1;
  e->forcerepart = 0;
}
