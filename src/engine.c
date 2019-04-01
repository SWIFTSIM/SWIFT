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
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "error.h"
#include "gravity.h"
#include "gravity_cache.h"
#include "hydro.h"
#include "logger.h"
#include "logger_io.h"
#include "map.h"
#include "memswap.h"
#include "memuse.h"
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
#include "star_formation.h"
#include "stars_io.h"
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
                                     "temperature",
                                     "cooling",
                                     "stars",
                                     "fof search",
                                     "structure finding",
                                     "star formation",
                                     "feedback",
                                     "time-step limiter"};

/** The rank of the engine as a global variable (for messages). */
int engine_rank;

/** The current step of the engine as a global variable (for messages). */
int engine_current_step;

/**
 * @brief Data collected from the cells at the end of a time-step
 */
struct end_of_step_data {

  size_t updated, g_updated, s_updated;
  size_t inhibited, g_inhibited, s_inhibited;
  integertime_t ti_hydro_end_min, ti_hydro_end_max, ti_hydro_beg_max;
  integertime_t ti_gravity_end_min, ti_gravity_end_max, ti_gravity_beg_max;
  integertime_t ti_stars_end_min;
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
    error(
        "Link table overflow. Increase the value of "
        "`Scheduler:links_per_tasks`.");
  }
  struct link *res = &e->links[ind];

  /* Set it atomically. */
  res->t = t;
  res->next = atomic_swap(l, res);
}

#ifdef WITH_MPI
/**
 * Do the exchange of one type of particles with all the other nodes.
 *
 * @param label a label for the memory allocations of this particle type.
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
static void *engine_do_redistribute(const char *label, int *counts, char *parts,
                                    size_t new_nr_parts, size_t sizeofparts,
                                    size_t alignsize, MPI_Datatype mpi_type,
                                    int nr_nodes, int nodeID) {

  /* Allocate a new particle array with some extra margin */
  char *parts_new = NULL;
  if (swift_memalign(
          label, (void **)&parts_new, alignsize,
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
            if (parts[k + offset].gpart->id_or_neg_offset > 0)                 \
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
      else if (s->gparts[k].type == swift_type_stars) {

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
  struct xpart *xparts = s->xparts;
  struct part *parts = s->parts;
  struct gpart *gparts = s->gparts;
  struct spart *sparts = s->sparts;
  ticks tic = getticks();

  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;

  /* Start by moving inhibited particles to the end of the arrays */
  for (size_t k = 0; k < nr_parts; /* void */) {
    if (parts[k].time_bin == time_bin_inhibited) {
      nr_parts -= 1;

      /* Swap the particle */
      memswap(&parts[k], &parts[nr_parts], sizeof(struct part));

      /* Swap the xpart */
      memswap(&xparts[k], &xparts[nr_parts], sizeof(struct xpart));

      /* Swap the link with the gpart */
      if (parts[k].gpart != NULL) {
        parts[k].gpart->id_or_neg_offset = -k;
      }
      if (parts[nr_parts].gpart != NULL) {
        parts[nr_parts].gpart->id_or_neg_offset = -nr_parts;
      }
    } else {
      k++;
    }
  }

  /* Now move inhibited star particles to the end of the arrays */
  for (size_t k = 0; k < nr_sparts; /* void */) {
    if (sparts[k].time_bin == time_bin_inhibited) {
      nr_sparts -= 1;

      /* Swap the particle */
      memswap(&s->sparts[k], &s->sparts[nr_sparts], sizeof(struct spart));

      /* Swap the link with the gpart */
      if (s->sparts[k].gpart != NULL) {
        s->sparts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->sparts[nr_sparts].gpart != NULL) {
        s->sparts[nr_sparts].gpart->id_or_neg_offset = -nr_sparts;
      }
    } else {
      k++;
    }
  }

  /* Finally do the same with the gravity particles */
  for (size_t k = 0; k < nr_gparts; /* void */) {
    if (gparts[k].time_bin == time_bin_inhibited) {
      nr_gparts -= 1;

      /* Swap the particle */
      memswap(&s->gparts[k], &s->gparts[nr_gparts], sizeof(struct gpart));

      /* Swap the link with part/spart */
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_stars) {
        s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }
      if (s->gparts[nr_gparts].type == swift_type_gas) {
        s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      } else if (s->gparts[nr_gparts].type == swift_type_stars) {
        s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      }
    } else {
      k++;
    }
  }

  /* Now we are ready to deal with real particles and can start the exchange. */

  /* Allocate temporary arrays to store the counts of particles to be sent
   * and the destination of each particle */
  int *counts;
  if ((counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate counts temporary buffer.");

  int *dest;
  if ((dest = (int *)swift_malloc("dest", sizeof(int) * nr_parts)) == NULL)
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
                 nr_parts, sizeof(struct part), 0, &redist_data);

  /* Sort the particles according to their cell index. */
  if (nr_parts > 0)
    space_parts_sort(s->parts, s->xparts, dest, &counts[nodeID * nr_nodes],
                     nr_nodes, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < nr_parts; k++) {
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
  if (nr_parts > 0 && nr_gparts > 0) {

    struct savelink_mapper_data savelink_data;
    savelink_data.nr_nodes = nr_nodes;
    savelink_data.counts = counts;
    savelink_data.parts = (void *)parts;
    savelink_data.nodeID = nodeID;
    threadpool_map(&e->threadpool, engine_redistribute_savelink_mapper_part,
                   nodes, nr_nodes, sizeof(int), 0, &savelink_data);
  }
  swift_free("dest", dest);

  /* Get destination of each s-particle */
  int *s_counts;
  if ((s_counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate s_counts temporary buffer.");

  int *s_dest;
  if ((s_dest = (int *)swift_malloc("s_dest", sizeof(int) * nr_sparts)) == NULL)
    error("Failed to allocate s_dest temporary buffer.");

  redist_data.counts = s_counts;
  redist_data.dest = s_dest;
  redist_data.base = (void *)sparts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_spart, sparts,
                 nr_sparts, sizeof(struct spart), 0, &redist_data);

  /* Sort the particles according to their cell index. */
  if (nr_sparts > 0)
    space_sparts_sort(s->sparts, s_dest, &s_counts[nodeID * nr_nodes], nr_nodes,
                      0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < nr_sparts; k++) {
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
  if (nr_sparts > 0) {

    struct savelink_mapper_data savelink_data;
    savelink_data.nr_nodes = nr_nodes;
    savelink_data.counts = s_counts;
    savelink_data.parts = (void *)sparts;
    savelink_data.nodeID = nodeID;
    threadpool_map(&e->threadpool, engine_redistribute_savelink_mapper_spart,
                   nodes, nr_nodes, sizeof(int), 0, &savelink_data);
  }
  swift_free("s_dest", s_dest);

  /* Get destination of each g-particle */
  int *g_counts;
  if ((g_counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate g_gcount temporary buffer.");

  int *g_dest;
  if ((g_dest = (int *)swift_malloc("g_dest", sizeof(int) * nr_gparts)) == NULL)
    error("Failed to allocate g_dest temporary buffer.");

  redist_data.counts = g_counts;
  redist_data.dest = g_dest;
  redist_data.base = (void *)gparts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_gpart, gparts,
                 nr_gparts, sizeof(struct gpart), 0, &redist_data);

  /* Sort the gparticles according to their cell index. */
  if (nr_gparts > 0)
    space_gparts_sort(s->gparts, s->parts, s->sparts, g_dest,
                      &g_counts[nodeID * nr_nodes], nr_nodes);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < nr_gparts; k++) {
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

  swift_free("g_dest", g_dest);

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
  size_t nr_parts_new = 0, nr_gparts_new = 0, nr_sparts_new = 0;
  for (int k = 0; k < nr_nodes; k++)
    nr_parts_new += counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_gparts_new += g_counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_sparts_new += s_counts[k * nr_nodes + nodeID];

  /* Now exchange the particles, type by type to keep the memory required
   * under control. */

  /* SPH particles. */
  void *new_parts = engine_do_redistribute(
      "parts", counts, (char *)s->parts, nr_parts_new, sizeof(struct part),
      part_align, part_mpi_type, nr_nodes, nodeID);
  swift_free("parts", s->parts);
  s->parts = (struct part *)new_parts;
  s->nr_parts = nr_parts_new;
  s->size_parts = engine_redistribute_alloc_margin * nr_parts_new;

  /* Extra SPH particle properties. */
  new_parts = engine_do_redistribute(
      "xparts", counts, (char *)s->xparts, nr_parts_new, sizeof(struct xpart),
      xpart_align, xpart_mpi_type, nr_nodes, nodeID);
  swift_free("xparts", s->xparts);
  s->xparts = (struct xpart *)new_parts;

  /* Gravity particles. */
  new_parts = engine_do_redistribute(
      "gparts", g_counts, (char *)s->gparts, nr_gparts_new,
      sizeof(struct gpart), gpart_align, gpart_mpi_type, nr_nodes, nodeID);
  swift_free("gparts", s->gparts);
  s->gparts = (struct gpart *)new_parts;
  s->nr_gparts = nr_gparts_new;
  s->size_gparts = engine_redistribute_alloc_margin * nr_gparts_new;

  /* Star particles. */
  new_parts = engine_do_redistribute(
      "sparts", s_counts, (char *)s->sparts, nr_sparts_new,
      sizeof(struct spart), spart_align, spart_mpi_type, nr_nodes, nodeID);
  swift_free("sparts", s->sparts);
  s->sparts = (struct spart *)new_parts;
  s->nr_sparts = nr_sparts_new;
  s->size_sparts = engine_redistribute_alloc_margin * nr_sparts_new;

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
  for (size_t k = 0; k < nr_parts_new; k++) {
    const int cid = cell_getid(s->cdim, s->parts[k].x[0] * s->iwidth[0],
                               s->parts[k].x[1] * s->iwidth[1],
                               s->parts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received particle (%zu) that does not belong here (nodeID=%i).", k,
            cells[cid].nodeID);
  }
  for (size_t k = 0; k < nr_gparts_new; k++) {
    const int cid = cell_getid(s->cdim, s->gparts[k].x[0] * s->iwidth[0],
                               s->gparts[k].x[1] * s->iwidth[1],
                               s->gparts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received g-particle (%zu) that does not belong here (nodeID=%i).",
            k, cells[cid].nodeID);
  }
  for (size_t k = 0; k < nr_sparts_new; k++) {
    const int cid = cell_getid(s->cdim, s->sparts[k].x[0] * s->iwidth[0],
                               s->sparts[k].x[1] * s->iwidth[1],
                               s->sparts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received s-particle (%zu) that does not belong here (nodeID=%i).",
            k, cells[cid].nodeID);
  }

  /* Verify that the links are correct */
  part_verify_links(s->parts, s->gparts, s->sparts, nr_parts_new, nr_gparts_new,
                    nr_sparts_new, e->verbose);
#endif

  /* Be verbose about what just happened. */
  if (e->verbose) {
    int my_cells = 0;
    for (int k = 0; k < nr_cells; k++)
      if (cells[k].nodeID == nodeID) my_cells += 1;
    message("node %i now has %zu parts, %zu sparts and %zu gparts in %i cells.",
            nodeID, nr_parts_new, nr_sparts_new, nr_gparts_new, my_cells);
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

#if defined(WITH_MPI) && (defined(HAVE_PARMETIS) || defined(HAVE_METIS))

  ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Be verbose about this. */
  if (e->nodeID == 0 || e->verbose) message("repartitioning space");
  fflush(stdout);

  /* Check that all cells have been drifted to the current time */
  space_check_drift_point(e->s, e->ti_current, /*check_multipoles=*/0);
#endif

  /* Clear the repartition flag. */
  e->forcerepart = 0;

  /* Nothing to do if only using a single node. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (e->nr_nodes == 1) return;

  /* Generate the fixed costs include file. */
  if (e->step > 3 && e->reparttype->trigger <= 1.f) {
    task_dump_stats("partition_fixed_costs.h", e, /* header = */ 1,
                    /* allranks = */ 1);
  }

  /* Do the repartitioning. */
  partition_repartition(e->reparttype, e->nodeID, e->nr_nodes, e->s,
                        e->sched.tasks, e->sched.nr_tasks);

  /* Partitioning requires copies of the particles, so we need to reduce the
   * memory in use to the minimum, we can free the sorting indices and the
   * tasks as these will be regenerated at the next rebuild. Also the foreign
   * particle arrays can go as these will be regenerated in proxy exchange. */

  /* Sorting indices. */
  if (e->s->cells_top != NULL) space_free_cells(e->s);

  /* Task arrays. */
  scheduler_free_tasks(&e->sched);

  /* Foreign parts. */
  space_free_foreign_parts(e->s);

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
    error("SWIFT was not compiled with MPI and METIS or ParMETIS support.");

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

  const ticks tic = getticks();

  /* Do nothing if there have not been enough steps since the last repartition
   * as we don't want to repeat this too often or immediately after a
   * repartition step. Also nothing to do when requested. */
  if (e->step - e->last_repartition >= 2 &&
      e->reparttype->type != REPART_NONE) {

    /* If we have fixed costs available and this is step 2 or we are forcing
     * repartitioning then we do a fixed costs one now. */
    if (e->reparttype->trigger > 1 ||
        (e->step == 2 && e->reparttype->use_fixed_costs)) {

      if (e->reparttype->trigger > 1) {
        if ((e->step % (int)e->reparttype->trigger) == 0) e->forcerepart = 1;
      } else {
        e->forcerepart = 1;
      }
      e->reparttype->use_ticks = 0;

    } else {

      /* It is only worth checking the CPU loads when we have processed a
       * significant number of all particles as we require all tasks to have
       * timings. */
      if ((e->updates > 1 &&
           e->updates >= e->total_nr_parts * e->reparttype->minfrac) ||
          (e->g_updates > 1 &&
           e->g_updates >= e->total_nr_gparts * e->reparttype->minfrac)) {

        /* Should we are use the task timings or fixed costs. */
        if (e->reparttype->use_fixed_costs > 1) {
          e->reparttype->use_ticks = 0;
        } else {
          e->reparttype->use_ticks = 1;
        }

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
          double abs_trigger = fabs(e->reparttype->trigger);
          if (((maxtime - mintime) / mean) > abs_trigger) {
            if (e->verbose)
              message("trigger fraction %.3f > %.3f will repartition",
                      (maxtime - mintime) / mean, abs_trigger);
            e->forcerepart = 1;
          } else {
            if (e->verbose)
              message("trigger fraction %.3f =< %.3f will not repartition",
                      (maxtime - mintime) / mean, abs_trigger);
          }
        }
      }

      /* All nodes do this together. */
      MPI_Bcast(&e->forcerepart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Remember we did this. */
    if (e->forcerepart) e->last_repartition = e->step;
  }

  /* We always reset CPU time for next check, unless it will not be used. */
  if (e->reparttype->type != REPART_NONE)
    e->cputime_last_step = clocks_get_cputime_used();

  if (e->verbose)
    message("took %.3f %s", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#endif
}

/**
 * @brief Exchange cell structures with other nodes.
 *
 * @param e The #engine.
 */
void engine_exchange_cells(struct engine *e) {

#ifdef WITH_MPI

  const int with_gravity = e->policy & engine_policy_self_gravity;
  const ticks tic = getticks();

  /* Exchange the cell structure with neighbouring ranks. */
  proxy_cells_exchange(e->proxies, e->nr_proxies, e->s, with_gravity);

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
 *        parts reside (i.e. the current number of local #part).
 * @param ind_part The foreign #cell ID of each part.
 * @param Npart The number of stray parts, contains the number of parts received
 *        on return.
 * @param offset_gparts The index in the gparts array as of which the foreign
 *        parts reside (i.e. the current number of local #gpart).
 * @param ind_gpart The foreign #cell ID of each gpart.
 * @param Ngpart The number of stray gparts, contains the number of gparts
 *        received on return.
 * @param offset_sparts The index in the sparts array as of which the foreign
 *        parts reside (i.e. the current number of local #spart).
 * @param ind_spart The foreign #cell ID of each spart.
 * @param Nspart The number of stray sparts, contains the number of sparts
 *        received on return.
 *
 * Note that this function does not mess-up the linkage between parts and
 * gparts, i.e. the received particles have correct linkeage.
 */
void engine_exchange_strays(struct engine *e, const size_t offset_parts,
                            const int *ind_part, size_t *Npart,
                            const size_t offset_gparts, const int *ind_gpart,
                            size_t *Ngpart, const size_t offset_sparts,
                            const int *ind_spart, size_t *Nspart) {

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

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_part[k] == -1) continue;

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

#ifdef SWIFT_DEBUG_CHECKS
    if (s->parts[offset_parts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the part and xpart into the proxy. */
    proxy_parts_load(&e->proxies[pid], &s->parts[offset_parts + k],
                     &s->xparts[offset_parts + k], 1);
  }

  /* Put the sparts into the corresponding proxies. */
  for (size_t k = 0; k < *Nspart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_spart[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_spart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lld, x=[%e,%e,%e].",
          node_id, s->sparts[offset_sparts + k].id,
          s->sparts[offset_sparts + k].x[0], s->sparts[offset_sparts + k].x[1],
          s->sparts[offset_sparts + k].x[2]);
    }

    /* Re-link the associated gpart with the buffer offset of the spart. */
    if (s->sparts[offset_sparts + k].gpart != NULL) {
      s->sparts[offset_sparts + k].gpart->id_or_neg_offset =
          -e->proxies[pid].nr_sparts_out;
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->sparts[offset_sparts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the spart into the proxy */
    proxy_sparts_load(&e->proxies[pid], &s->sparts[offset_sparts + k], 1);
  }

  /* Put the gparts into the corresponding proxies. */
  for (size_t k = 0; k < *Ngpart; k++) {

    /* Ignore the particles we want to get rid of (inhibited, ...). */
    if (ind_gpart[k] == -1) continue;

    /* Get the target node and proxy ID. */
    const int node_id = e->s->cells_top[ind_gpart[k]].nodeID;
    if (node_id < 0 || node_id >= e->nr_nodes)
      error("Bad node ID %i.", node_id);
    const int pid = e->proxy_ind[node_id];
    if (pid < 0) {
      error(
          "Do not have a proxy for the requested nodeID %i for part with "
          "id=%lli, x=[%e,%e,%e].",
          node_id, s->gparts[offset_gparts + k].id_or_neg_offset,
          s->gparts[offset_gparts + k].x[0], s->gparts[offset_gparts + k].x[1],
          s->gparts[offset_gparts + k].x[2]);
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (s->gparts[offset_gparts + k].time_bin == time_bin_inhibited)
      error("Attempting to exchange an inhibited particle");
#endif

    /* Load the gpart into the proxy */
    proxy_gparts_load(&e->proxies[pid], &s->gparts[offset_gparts + k], 1);
  }

  /* Launch the proxies. */
  MPI_Request reqs_in[4 * engine_maxproxies];
  MPI_Request reqs_out[4 * engine_maxproxies];
  for (int k = 0; k < e->nr_proxies; k++) {
    proxy_parts_exchange_first(&e->proxies[k]);
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
    proxy_parts_exchange_second(&e->proxies[pid]);
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
    if (swift_memalign("parts", (void **)&parts_new, part_align,
                       sizeof(struct part) * s->size_parts) != 0 ||
        swift_memalign("xparts", (void **)&xparts_new, xpart_align,
                       sizeof(struct xpart) * s->size_parts) != 0)
      error("Failed to allocate new part data.");
    memcpy(parts_new, s->parts, sizeof(struct part) * offset_parts);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * offset_parts);
    swift_free("parts", s->parts);
    swift_free("xparts", s->xparts);
    s->parts = parts_new;
    s->xparts = xparts_new;

    /* Reset the links */
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
    if (swift_memalign("sparts", (void **)&sparts_new, spart_align,
                       sizeof(struct spart) * s->size_sparts) != 0)
      error("Failed to allocate new spart data.");
    memcpy(sparts_new, s->sparts, sizeof(struct spart) * offset_sparts);
    swift_free("sparts", s->sparts);
    s->sparts = sparts_new;

    /* Reset the links */
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
    if (swift_memalign("gparts", (void **)&gparts_new, gpart_align,
                       sizeof(struct gpart) * s->size_gparts) != 0)
      error("Failed to allocate new gpart data.");
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * offset_gparts);
    swift_free("gparts", s->gparts);
    s->gparts = gparts_new;

    /* Reset the links */
    for (size_t k = 0; k < offset_gparts; k++) {
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_stars) {
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
        } else if (gp->type == swift_type_stars) {
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

  ticks tic = getticks();

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
   * We use our home-made reduction operation that simply performs a XOR
   * operation on the multipoles. Since only local multipoles are non-zero and
   * each multipole is only present once, the bit-by-bit XOR will
   * create the desired result.
   */
  int err = MPI_Allreduce(MPI_IN_PLACE, e->s->multipoles_top, e->s->nr_cells,
                          multipole_mpi_type, multipole_mpi_reduce_op,
                          MPI_COMM_WORLD);
  if (err != MPI_SUCCESS)
    mpi_error(err, "Failed to all-reduce the top-level multipoles.");

#ifdef SWIFT_DEBUG_CHECKS
  long long counter = 0;

  /* Let's check that what we received makes sense */
  for (int i = 0; i < e->s->nr_cells; ++i) {
    const struct gravity_tensors *m = &e->s->multipoles_top[i];
    counter += m->m_pole.num_gpart;
    if (m->m_pole.num_gpart < 0) {
      error("m->m_pole.num_gpart is negative: %lld", m->m_pole.num_gpart);
    }
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
    error(
        "Total particles in multipoles inconsistent with engine.\n "
        "  counter = %lld, nr_gparts = %lld",
        counter, e->total_nr_gparts);
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

void engine_exchange_proxy_multipoles(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();

  /* Start by counting the number of cells to send and receive */
  int count_send_cells = 0;
  int count_recv_cells = 0;
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
      count_recv_cells += p->cells_in[k]->mpi.pcell_size;

    for (int k = 0; k < p->nr_cells_out; k++)
      count_send_cells += p->cells_out[k]->mpi.pcell_size;
  }

  /* Allocate the buffers for the packed data */
  struct gravity_tensors *buffer_send = NULL;
  if (swift_memalign("send_gravity_tensors", (void **)&buffer_send,
                     SWIFT_CACHE_ALIGNMENT,
                     count_send_cells * sizeof(struct gravity_tensors)) != 0)
    error("Unable to allocate memory for multipole transactions");

  struct gravity_tensors *buffer_recv = NULL;
  if (swift_memalign("recv_gravity_tensors", (void **)&buffer_recv,
                     SWIFT_CACHE_ALIGNMENT,
                     count_recv_cells * sizeof(struct gravity_tensors)) != 0)
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

      const int num_elements = p->cells_in[k]->mpi.pcell_size;

      /* Receive everything */
      MPI_Irecv(&buffer_recv[this_recv], num_elements, multipole_mpi_type,
                p->cells_in[k]->nodeID, p->cells_in[k]->mpi.tag, MPI_COMM_WORLD,
                &requests[this_request]);

      /* Move to the next slot in the buffers */
      this_recv += num_elements;
      this_request++;
    }

    /* Loop over the proxies to issue the sends. */
    for (int k = 0; k < p->nr_cells_out; k++) {

      /* Number of multipoles in this cell hierarchy */
      const int num_elements = p->cells_out[k]->mpi.pcell_size;

      /* Let's pack everything recursively */
      cell_pack_multipoles(p->cells_out[k], &buffer_send[this_send]);

      /* Send everything (note the use of cells_in[0] to get the correct node
       * ID. */
      MPI_Isend(&buffer_send[this_send], num_elements, multipole_mpi_type,
                p->cells_in[0]->nodeID, p->cells_out[k]->mpi.tag,
                MPI_COMM_WORLD, &requests[this_request]);

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

      const int num_elements = p->cells_in[k]->mpi.pcell_size;

#ifdef SWIFT_DEBUG_CHECKS

      /* Check that the first element (top-level cell's multipole) matches what
       * we received */
      if (p->cells_in[k]->grav.multipole->m_pole.num_gpart !=
          buffer_recv[this_recv].m_pole.num_gpart)
        error("Current: M_000=%e num_gpart=%lld\n New: M_000=%e num_gpart=%lld",
              p->cells_in[k]->grav.multipole->m_pole.M_000,
              p->cells_in[k]->grav.multipole->m_pole.num_gpart,
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
 * @brief Allocate memory for the foreign particles.
 *
 * We look into the proxies for cells that have tasks and count
 * the number of particles in these cells. We then allocate
 * memory and link all the cells that have tasks and all cells
 * deeper in the tree.
 *
 * @param e The #engine.
 */
void engine_allocate_foreign_particles(struct engine *e) {

#ifdef WITH_MPI

  const int nr_proxies = e->nr_proxies;
  struct space *s = e->s;
  ticks tic = getticks();

  /* Count the number of particles we need to import and re-allocate
     the buffer if needed. */
  size_t count_parts_in = 0, count_gparts_in = 0, count_sparts_in = 0;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro) {
        count_parts_in += cell_count_parts_for_tasks(e->proxies[k].cells_in[j]);
      }

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity) {
        count_gparts_in +=
            cell_count_gparts_for_tasks(e->proxies[k].cells_in[j]);
      }

      /* For stars, we just use the numbers in the top-level cells */
      count_sparts_in += e->proxies[k].cells_in[j]->stars.count;
    }
  }

  if (e->verbose)
    message("Counting number of foreign particles took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Allocate space for the foreign particles we will receive */
  if (count_parts_in > s->size_parts_foreign) {
    if (s->parts_foreign != NULL)
      swift_free("sparts_foreign", s->parts_foreign);
    s->size_parts_foreign = engine_foreign_alloc_margin * count_parts_in;
    if (swift_memalign("parts_foreign", (void **)&s->parts_foreign, part_align,
                       sizeof(struct part) * s->size_parts_foreign) != 0)
      error("Failed to allocate foreign part data.");
  }

  /* Allocate space for the foreign particles we will receive */
  if (count_gparts_in > s->size_gparts_foreign) {
    if (s->gparts_foreign != NULL)
      swift_free("gparts_foreign", s->gparts_foreign);
    s->size_gparts_foreign = engine_foreign_alloc_margin * count_gparts_in;
    if (swift_memalign("gparts_foreign", (void **)&s->gparts_foreign,
                       gpart_align,
                       sizeof(struct gpart) * s->size_gparts_foreign) != 0)
      error("Failed to allocate foreign gpart data.");
  }

  /* Allocate space for the foreign particles we will receive */
  if (count_sparts_in > s->size_sparts_foreign) {
    if (s->sparts_foreign != NULL)
      swift_free("sparts_foreign", s->sparts_foreign);
    s->size_sparts_foreign = engine_foreign_alloc_margin * count_sparts_in;
    if (swift_memalign("sparts_foreign", (void **)&s->sparts_foreign,
                       spart_align,
                       sizeof(struct spart) * s->size_sparts_foreign) != 0)
      error("Failed to allocate foreign spart data.");
  }

  if (e->verbose)
    message("Allocating %zd/%zd/%zd foreign part/gpart/spart (%zd/%zd/%zd MB)",
            s->size_parts_foreign, s->size_gparts_foreign,
            s->size_sparts_foreign,
            s->size_parts_foreign * sizeof(struct part) / (1024 * 1024),
            s->size_gparts_foreign * sizeof(struct gpart) / (1024 * 1024),
            s->size_sparts_foreign * sizeof(struct spart) / (1024 * 1024));

  /* Unpack the cells and link to the particle data. */
  struct part *parts = s->parts_foreign;
  struct gpart *gparts = s->gparts_foreign;
  struct spart *sparts = s->sparts_foreign;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro) {

        const size_t count_parts =
            cell_link_foreign_parts(e->proxies[k].cells_in[j], parts);
        parts = &parts[count_parts];
      }

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity) {

        const size_t count_gparts =
            cell_link_foreign_gparts(e->proxies[k].cells_in[j], gparts);
        gparts = &gparts[count_gparts];
      }

      /* For stars, we just use the numbers in the top-level cells */
      cell_link_sparts(e->proxies[k].cells_in[j], sparts);
      sparts = &sparts[e->proxies[k].cells_in[j]->stars.count];
    }
  }

  /* Update the counters */
  s->nr_parts_foreign = parts - s->parts_foreign;
  s->nr_gparts_foreign = gparts - s->gparts_foreign;
  s->nr_sparts_foreign = sparts - s->sparts_foreign;

  if (e->verbose)
    message("Recursively linking foreign arrays took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Prints the number of tasks in the engine
 *
 * @param e The #engine.
 */
void engine_print_task_counts(const struct engine *e) {

  const ticks tic = getticks();
  const struct scheduler *sched = &e->sched;
  const int nr_tasks = sched->nr_tasks;
  const struct task *const tasks = sched->tasks;

  /* Global tasks and cells when using MPI. */
#ifdef WITH_MPI
  if (e->nodeID == 0 && e->total_nr_tasks > 0)
    printf(
        "[%04i] %s engine_print_task_counts: System total: %lld,"
        " no. cells: %lld\n",
        e->nodeID, clocks_get_timesincestart(), e->total_nr_tasks,
        e->total_nr_cells);
  fflush(stdout);
#endif

  /* Report value that can be used to estimate the task_per_cells parameter. */
  float tasks_per_cell = (float)nr_tasks / (float)e->s->tot_cells;

#ifdef WITH_MPI
  message("Total = %d (per cell = %.2f)", nr_tasks, tasks_per_cell);

  /* And the system maximum on rank 0, only after first step, increase by our
   * margin to allow for some variation in repartitioning. */
  if (e->nodeID == 0 && e->total_nr_tasks > 0) {
    message("Total = %d (maximum per cell = %.2f)", nr_tasks,
            e->tasks_per_cell_max * engine_tasks_per_cell_margin);
  }

#else
  message("Total = %d (per cell = %.2f)", nr_tasks, tasks_per_cell);
#endif
  fflush(stdout);

  /* Count and print the number of each task type. */
  int counts[task_type_count + 1];
  for (int k = 0; k <= task_type_count; k++) counts[k] = 0;
  for (int k = 0; k < nr_tasks; k++) {
    if (tasks[k].skip)
      counts[task_type_count] += 1;
    else
      counts[(int)tasks[k].type] += 1;
  }

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
 * If e->tasks_per_cell is set greater than 0.0 then that value is used
 * as the estimate of the average number of tasks per cell,
 * otherwise we attempt an estimate.
 *
 * @param e the #engine
 *
 * @return the estimated total number of tasks
 */
int engine_estimate_nr_tasks(const struct engine *e) {

  float tasks_per_cell = e->tasks_per_cell;
  if (tasks_per_cell > 0.0f) {
    if (e->verbose)
      message("tasks per cell given as: %.2f, so maximum tasks: %d",
              e->tasks_per_cell, (int)(e->s->tot_cells * tasks_per_cell));
    return (int)(e->s->tot_cells * tasks_per_cell);
  }

  /* Our guess differs depending on the types of tasks we are using, but we
   * basically use a formula <n1>*ntopcells + <n2>*(totcells - ntopcells).
   * Where <n1> is the expected maximum tasks per top-level/super cell, and
   * <n2> the expected maximum tasks for all other cells. These should give
   * a safe upper limit. */
  int n1 = 0;
  int n2 = 0;
  if (e->policy & engine_policy_hydro) {
    /* 2 self (density, force), 1 sort, 26/2 density pairs
       26/2 force pairs, 1 drift, 3 ghosts, 2 kicks, 1 time-step,
       1 end_force, 2 extra space
     */
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
  if (e->policy & engine_policy_limiter) {
    n1 += 18;
    n2 += 1;
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
    /* Cooling task + extra space */
    n1 += 2;
  }
  if (e->policy & engine_policy_star_formation) {
    n1 += 1;
  }
  if (e->policy & engine_policy_stars) {
    /* 2 self (density, feedback), 1 sort, 26/2 density pairs
       26/2 feedback pairs, 1 drift, 3 ghosts, 2 kicks, 1 time-step,
       1 end_force, 2 extra space
     */
    n1 += 37;
    n2 += 2;
#ifdef WITH_MPI
    n1 += 6;
#endif
  }
  if (e->policy & engine_policy_fof) {
    n1 += 2;
  }
#if defined(WITH_LOGGER)
  /* each cell logs its particles */
  n1 += 1;
#endif

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
    int nparts = c->hydro.count + c->grav.count + c->stars.count;
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

  float ntasks = n1 * ntop + n2 * (ncells - ntop);
  if (ncells > 0) tasks_per_cell = ceil(ntasks / ncells);

  if (tasks_per_cell < 1.0f) tasks_per_cell = 1.0f;
  if (e->verbose)
    message("tasks per cell estimated as: %.2f, maximum tasks: %d",
            tasks_per_cell, (int)(ncells * tasks_per_cell));

  return (int)(ncells * tasks_per_cell);
}

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 * @param repartitioned Did we just redistribute?
 * @param clean_smoothing_length_values Are we cleaning up the values of
 * the smoothing lengths before building the tasks ?
 */
void engine_rebuild(struct engine *e, int repartitioned,
                    int clean_smoothing_length_values) {

  const ticks tic = getticks();

  /* Clear the forcerebuild flag, whatever it was. */
  e->forcerebuild = 0;
  e->restarting = 0;

  /* Re-build the space. */
  space_rebuild(e->s, repartitioned, e->verbose);

  /* Report the number of cells and memory */
  if (e->verbose)
    message(
        "Nr. of top-level cells: %d Nr. of local cells: %d memory use: %zd MB.",
        e->s->nr_cells, e->s->tot_cells,
        (e->s->nr_cells + e->s->tot_cells) * sizeof(struct cell) /
            (1024 * 1024));

  /* Report the number of multipoles and memory */
  if (e->verbose && (e->policy & engine_policy_self_gravity))
    message(
        "Nr. of top-level mpoles: %d Nr. of local mpoles: %d memory use: %zd "
        "MB.",
        e->s->nr_cells, e->s->tot_cells,
        (e->s->nr_cells + e->s->tot_cells) * sizeof(struct gravity_tensors) /
            (1024 * 1024));

  const ticks tic2 = getticks();

  /* Update the global counters of particles */
  long long num_particles[3] = {
      (long long)(e->s->nr_parts - e->s->nr_extra_parts),
      (long long)(e->s->nr_gparts - e->s->nr_extra_gparts),
      (long long)(e->s->nr_sparts - e->s->nr_extra_sparts)};
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, num_particles, 3, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  e->total_nr_parts = num_particles[0];
  e->total_nr_gparts = num_particles[1];
  e->total_nr_sparts = num_particles[2];

  /* Flag that there are no inhibited particles */
  e->nr_inhibited_parts = 0;
  e->nr_inhibited_gparts = 0;
  e->nr_inhibited_sparts = 0;

  if (e->verbose)
    message("updating particle counts took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Re-compute the mesh forces */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
    pm_mesh_compute_potential(e->mesh, e->s, &e->threadpool, e->verbose);

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
  if (e->policy & engine_policy_self_gravity) engine_exchange_top_multipoles(e);

  engine_exchange_cells(e);
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Let's check that what we received makes sense */
  if (e->policy & engine_policy_self_gravity) {
    long long counter = 0;

    for (int i = 0; i < e->s->nr_cells; ++i) {
      const struct gravity_tensors *m = &e->s->multipoles_top[i];
      counter += m->m_pole.num_gpart;
    }
    if (counter != e->total_nr_gparts)
      error("Total particles in multipoles inconsistent with engine");
  }
#endif

  /* Re-build the tasks. */
  engine_maketasks(e);

  /* Make the list of top-level cells that have tasks */
  space_list_useful_top_level_cells(e->s);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);

  if (e->policy & engine_policy_self_gravity) {
    for (int k = 0; k < e->s->nr_local_cells; k++)
      cell_check_foreign_multipole(&e->s->cells_top[e->s->local_cells_top[k]]);
  }
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
  int repartitioned = 0;

  /* Unskip active tasks and check for rebuild */
  if (!e->forcerebuild && !e->forcerepart && !e->restarting) engine_unskip(e);

  const ticks tic3 = getticks();

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->forcerebuild, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#endif

  if (e->verbose)
    message("Communicating rebuild flag took %.3f %s.",
            clocks_from_ticks(getticks() - tic3), clocks_getunit());

  /* Do we need repartitioning ? */
  if (e->forcerepart) {

    /* Let's start by drifting everybody to the current time */
    engine_drift_all(e, /*drift_mpole=*/0);
    drifted_all = 1;

    /* And repartition */
    engine_repartition(e);
    repartitioned = 1;
  }

  /* Do we need rebuilding ? */
  if (e->forcerebuild) {

    /* Let's start by drifting everybody to the current time */
    if (!e->restarting && !drifted_all) engine_drift_all(e, /*drift_mpole=*/0);

    /* And rebuild */
    engine_rebuild(e, repartitioned, 0);
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
 * @brief Recursive function gathering end-of-step data.
 *
 * We recurse until we encounter a timestep or time-step MPI recv task
 * as the values will have been set at that level. We then bring these
 * values upwards.
 *
 * @param c The #cell to recurse into.
 * @param e The #engine.
 */
void engine_collect_end_of_step_recurse(struct cell *c,
                                        const struct engine *e) {

/* Skip super-cells (Their values are already set) */
#ifdef WITH_MPI
  if (c->timestep != NULL || c->mpi.recv_ti != NULL) return;
#else
  if (c->timestep != NULL) return;
#endif /* WITH_MPI */

  /* Counters for the different quantities. */
  size_t updated = 0, g_updated = 0, s_updated = 0;
  size_t inhibited = 0, g_inhibited = 0, s_inhibited = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps;

  /* Collect the values from the progeny. */
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp != NULL &&
        (cp->hydro.count > 0 || cp->grav.count > 0 || cp->stars.count > 0)) {

      /* Recurse */
      engine_collect_end_of_step_recurse(cp, e);

      /* And update */
      ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, cp->hydro.ti_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);

      ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
      ti_gravity_end_max = max(ti_gravity_end_max, cp->grav.ti_end_max);
      ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);

      ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);

      updated += cp->hydro.updated;
      g_updated += cp->grav.updated;
      s_updated += cp->stars.updated;

      inhibited += cp->hydro.inhibited;
      g_inhibited += cp->grav.inhibited;
      s_inhibited += cp->stars.inhibited;

      /* Collected, so clear for next time. */
      cp->hydro.updated = 0;
      cp->grav.updated = 0;
      cp->stars.updated = 0;
    }
  }

  /* Store the collected values in the cell. */
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_end_max = ti_gravity_end_max;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->hydro.updated = updated;
  c->grav.updated = g_updated;
  c->stars.updated = s_updated;
  c->hydro.inhibited = inhibited;
  c->grav.inhibited = g_inhibited;
  c->stars.inhibited = s_inhibited;
}

/**
 * @brief Mapping function to collect the data from the end of the step
 *
 * This function will call a recursive function on all the top-level cells
 * to collect the information we are after.
 *
 * @param map_data The list of cells with tasks on this node.
 * @param num_elements The number of elements in the list this thread will work
 * on.
 * @param extra_data The #engine.
 */
void engine_collect_end_of_step_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct end_of_step_data *data = (struct end_of_step_data *)extra_data;
  const struct engine *e = data->e;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;

  /* Local collectible */
  size_t updated = 0, g_updated = 0, s_updated = 0;
  size_t inhibited = 0, g_inhibited = 0, s_inhibited = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];

    if (c->hydro.count > 0 || c->grav.count > 0 || c->stars.count > 0) {

      /* Make the top-cells recurse */
      engine_collect_end_of_step_recurse(c, e);

      /* And aggregate */
      if (c->hydro.ti_end_min > e->ti_current)
        ti_hydro_end_min = min(ti_hydro_end_min, c->hydro.ti_end_min);
      ti_hydro_end_max = max(ti_hydro_end_max, c->hydro.ti_end_max);
      ti_hydro_beg_max = max(ti_hydro_beg_max, c->hydro.ti_beg_max);

      if (c->grav.ti_end_min > e->ti_current)
        ti_gravity_end_min = min(ti_gravity_end_min, c->grav.ti_end_min);
      ti_gravity_end_max = max(ti_gravity_end_max, c->grav.ti_end_max);
      ti_gravity_beg_max = max(ti_gravity_beg_max, c->grav.ti_beg_max);

      if (c->stars.ti_end_min > e->ti_current)
        ti_stars_end_min = min(ti_stars_end_min, c->stars.ti_end_min);

      updated += c->hydro.updated;
      g_updated += c->grav.updated;
      s_updated += c->stars.updated;

      inhibited += c->hydro.inhibited;
      g_inhibited += c->grav.inhibited;
      s_inhibited += c->stars.inhibited;

      /* Collected, so clear for next time. */
      c->hydro.updated = 0;
      c->grav.updated = 0;
      c->stars.updated = 0;
    }
  }

  /* Let's write back to the global data.
   * We use the space lock to garanty single access*/
  if (lock_lock(&s->lock) == 0) {
    data->updated += updated;
    data->g_updated += g_updated;
    data->s_updated += s_updated;

    data->inhibited += inhibited;
    data->g_inhibited += g_inhibited;
    data->s_inhibited += s_inhibited;

    if (ti_hydro_end_min > e->ti_current)
      data->ti_hydro_end_min = min(ti_hydro_end_min, data->ti_hydro_end_min);
    data->ti_hydro_end_max = max(ti_hydro_end_max, data->ti_hydro_end_max);
    data->ti_hydro_beg_max = max(ti_hydro_beg_max, data->ti_hydro_beg_max);

    if (ti_gravity_end_min > e->ti_current)
      data->ti_gravity_end_min =
          min(ti_gravity_end_min, data->ti_gravity_end_min);
    data->ti_gravity_end_max =
        max(ti_gravity_end_max, data->ti_gravity_end_max);
    data->ti_gravity_beg_max =
        max(ti_gravity_beg_max, data->ti_gravity_beg_max);

    if (ti_stars_end_min > e->ti_current)
      data->ti_stars_end_min = min(ti_stars_end_min, data->ti_stars_end_min);
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
  struct space *s = e->s;
  struct end_of_step_data data;
  data.updated = 0, data.g_updated = 0, data.s_updated = 0;
  data.inhibited = 0, data.g_inhibited = 0, data.s_inhibited = 0;
  data.ti_hydro_end_min = max_nr_timesteps, data.ti_hydro_end_max = 0,
  data.ti_hydro_beg_max = 0;
  data.ti_gravity_end_min = max_nr_timesteps, data.ti_gravity_end_max = 0,
  data.ti_gravity_beg_max = 0;
  data.e = e;

  /* Collect information from the local top-level cells */
  threadpool_map(&e->threadpool, engine_collect_end_of_step_mapper,
                 s->local_cells_with_tasks_top, s->nr_local_cells_with_tasks,
                 sizeof(int), 0, &data);

  /* Store the local number of inhibited particles */
  s->nr_inhibited_parts = data.inhibited;
  s->nr_inhibited_gparts = data.g_inhibited;
  s->nr_inhibited_sparts = data.s_inhibited;

  /* Store these in the temporary collection group. */
  collectgroup1_init(
      &e->collect_group1, data.updated, data.g_updated, data.s_updated,
      data.inhibited, data.g_inhibited, data.s_inhibited, data.ti_hydro_end_min,
      data.ti_hydro_end_max, data.ti_hydro_beg_max, data.ti_gravity_end_min,
      data.ti_gravity_end_max, data.ti_gravity_beg_max, e->forcerebuild,
      e->s->tot_cells, e->sched.nr_tasks,
      (float)e->sched.nr_tasks / (float)e->s->tot_cells);

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
    out_ll[0] = data.updated;
    out_ll[1] = data.g_updated;
    out_ll[2] = data.s_updated;
    if (MPI_Allreduce(out_ll, in_ll, 3, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.updated)
      error("Failed to get same updated, is %lld, should be %lld", in_ll[0],
            e->collect_group1.updated);
    if (in_ll[1] != (long long)e->collect_group1.g_updated)
      error("Failed to get same g_updated, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_updated);
    if (in_ll[2] != (long long)e->collect_group1.s_updated)
      error("Failed to get same s_updated, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_updated);

    out_ll[0] = data.inhibited;
    out_ll[1] = data.g_inhibited;
    out_ll[2] = data.s_inhibited;
    if (MPI_Allreduce(out_ll, in_ll, 3, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.inhibited)
      error("Failed to get same inhibited, is %lld, should be %lld", in_ll[0],
            e->collect_group1.inhibited);
    if (in_ll[1] != (long long)e->collect_group1.g_inhibited)
      error("Failed to get same g_inhibited, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_inhibited);
    if (in_ll[2] != (long long)e->collect_group1.s_inhibited)
      error("Failed to get same s_inhibited, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_inhibited);

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
  space_check_drift_point(e->s, e->ti_current, /*chek_mpoles=*/0);

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
        t->type == task_type_drift_spart || t->type == task_type_kick1 ||
        t->type == task_type_kick2 || t->type == task_type_timestep ||
        t->type == task_type_timestep_limiter ||
        t->subtype == task_subtype_force ||
        t->subtype == task_subtype_limiter || t->subtype == task_subtype_grav ||
        t->type == task_type_end_hydro_force ||
        t->type == task_type_end_grav_force ||
        t->type == task_type_grav_long_range || t->type == task_type_grav_mm ||
        t->type == task_type_grav_down || t->type == task_type_grav_down_in ||
        t->type == task_type_drift_gpart_out || t->type == task_type_cooling ||
        t->type == task_type_stars_in || t->type == task_type_stars_out ||
        t->type == task_type_star_formation ||
        t->type == task_type_extra_ghost ||
        t->subtype == task_subtype_gradient ||
        t->subtype == task_subtype_stars_feedback ||
        t->subtype == task_subtype_tend || t->subtype == task_subtype_rho ||
        t->subtype == task_subtype_gpart)
      t->skip = 1;
  }

  /* Run through the cells and clear some flags. */
  space_map_cells_pre(e->s, 1, cell_clear_drift_flags, NULL);
  space_map_cells_pre(e->s, 1, cell_clear_limiter_flags, NULL);
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
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_drift_spart)
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
 * @param compute_init_accel Are we computing the initial acceleration of
 *particles?
 */
void engine_init_particles(struct engine *e, int flag_entropy_ICs,
                           int clean_h_values, int compute_init_accel) {

  struct space *s = e->s;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_props_update(e->gravity_properties, e->cosmology);

  /* Udpate the hydro properties */
  if (e->policy & engine_policy_hydro)
    hydro_props_update(e->hydro_properties, e->gravity_properties,
                       e->cosmology);

  /* Start by setting the particles in a good state */
  if (e->nodeID == 0) message("Setting particles to a valid state...");
  engine_first_init_particles(e);

  if (e->nodeID == 0) message("Computing initial gas densities.");

  /* Construct all cells and tasks to start everything */
  engine_rebuild(e, 0, clean_h_values);

  /* No time integration. We just want the density and ghosts */
  engine_skip_force_and_kick(e);

  /* Activate the self and pair FOF tasks and skip all other tasks. */
  for (int i = 0; i < e->sched.nr_tasks; i++) {

    struct task *t = &e->sched.tasks[i];

    t->skip = 1;

    if (t->type == task_type_fof_self || t->type == task_type_fof_pair) {
      scheduler_activate(&e->sched, t);
    }
  }

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

  /* Init the particle data (by hand). */
  space_init_parts(s, e->verbose);
  space_init_gparts(s, e->verbose);
  space_init_sparts(s, e->verbose);

  /* Update the cooling function */
  if ((e->policy & engine_policy_cooling) ||
      (e->policy & engine_policy_temperature))
    cooling_update(e->cosmology, e->cooling_func, e->s);

#ifdef WITH_LOGGER
  /* Mark the first time step in the particle logger file. */
  logger_log_timestamp(e->logger, e->ti_current, e->time,
                       &e->logger->timestamp_offset);
  /* Make sure that we have enough space in the particle logger file
   * to store the particles in current time step. */
  logger_ensure_size(e->logger, e->total_nr_parts, e->total_nr_gparts, 0);
#endif

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

  /* Exit if only computing the FOF. */
  if (!compute_init_accel) return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].grav.multipole->m_pole.num_gpart;
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

  /* Construct all cells again for a new round (need to update h_max) */
  engine_rebuild(e, 0, 0);

  /* No drift this time */
  engine_skip_drift(e);

  /* Init the particle data (by hand). */
  space_init_parts(e->s, e->verbose);
  space_init_gparts(e->s, e->verbose);
  space_init_sparts(e->s, e->verbose);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Run the brute-force gravity calculation for some gparts */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_compute(e->s, e);
#endif

  scheduler_write_dependencies(&e->sched, e->verbose);
  space_write_cell_hierarchy(e->s);
  if (e->nodeID == 0) scheduler_write_task_level(&e->sched);

  /* Run the 0th time-step */
  TIMER_TIC2;
  engine_launch(e);
  TIMER_TOC2(timer_runners);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Check the accuracy of the gravity calculation */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_check(e->s, e, 1e-1);
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure all woken-up particles have been processed */
  space_check_limiter(e->s);
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

      /* Ignore fake buffer particles for on-the-fly creation */
      if (s->parts[k].time_bin == time_bin_not_created) continue;

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

      /* Ignore fake buffer particles for on-the-fly creation */
      if (s->gparts[k].time_bin == time_bin_not_created) continue;

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
      if (c->nodeID == engine_rank && c->hydro.count > 0) {
        float part_h_max = c->hydro.parts[0].h;
        for (int k = 1; k < c->hydro.count; k++) {
          if (c->hydro.parts[k].h > part_h_max)
            part_h_max = c->hydro.parts[k].h;
        }
        c->hydro.h_max = max(part_h_max, c->hydro.h_max);
      }
    }
  }

  if (s->cells_top != NULL && s->nr_sparts > 0) {
    for (int i = 0; i < s->nr_cells; i++) {
      struct cell *c = &s->cells_top[i];
      if (c->nodeID == engine_rank && c->stars.count > 0) {
        float spart_h_max = c->stars.parts[0].h;
        for (int k = 1; k < c->stars.count; k++) {
          if (c->stars.parts[k].h > spart_h_max)
            spart_h_max = c->stars.parts[k].h;
        }
        c->stars.h_max = max(spart_h_max, c->stars.h_max);
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

  e->tic_step = getticks();

  if (e->nodeID == 0) {

    /* Print some information to the screen */
    printf(
        "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %21.3f "
        "%6d\n",
        e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
        e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
        e->s_updates, e->wallclock_time, e->step_props);
#ifdef SWIFT_DEBUG_CHECKS
    fflush(stdout);
#endif

    if (!e->restarting)
      fprintf(
          e->file_timesteps,
          "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %21.3f "
          "%6d\n",
          e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
          e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
          e->s_updates, e->wallclock_time, e->step_props);
#ifdef SWIFT_DEBUG_CHECKS
    fflush(e->file_timesteps);
#endif
  }

  /* We need some cells to exist but not the whole task stuff. */
  if (e->restarting) space_rebuild(e->s, 0, e->verbose);

  /* Move forward in time */
  e->ti_old = e->ti_current;
  e->ti_current = e->ti_end_min;
  e->max_active_bin = get_max_active_bin(e->ti_end_min);
  e->min_active_bin = get_min_active_bin(e->ti_current, e->ti_old);
  e->step += 1;
  engine_current_step = e->step;
  e->step_props = engine_step_prop_none;

  /* When restarting, move everyone to the current time. */
  if (e->restarting) engine_drift_all(e, /*drift_mpole=*/1);

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

  /* Update the cooling function */
  if ((e->policy & engine_policy_cooling) ||
      (e->policy & engine_policy_temperature))
    cooling_update(e->cosmology, e->cooling_func, e->s);

  /*****************************************************/
  /* OK, we now know what the next end of time-step is */
  /*****************************************************/

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_props_update(e->gravity_properties, e->cosmology);

  /* Udpate the hydro properties */
  if (e->policy & engine_policy_hydro)
    hydro_props_update(e->hydro_properties, e->gravity_properties,
                       e->cosmology);

  /* Trigger a tree-rebuild if we passed the frequency threshold */
  if ((e->policy & engine_policy_self_gravity) &&
      ((double)e->g_updates_since_rebuild >
       ((double)e->total_nr_gparts) * e->gravity_properties->rebuild_frequency))
    e->forcerebuild = 1;

#ifdef WITH_LOGGER
  /* Mark the current time step in the particle logger file. */
  logger_log_timestamp(e->logger, e->ti_current, e->time,
                       &e->logger->timestamp_offset);
  /* Make sure that we have enough space in the particle logger file
   * to store the particles in current time step. */
  logger_ensure_size(e->logger, e->total_nr_parts, e->total_nr_gparts, 0);
#endif

  /* Are we drifting everything (a la Gadget/GIZMO) ? */
  if (e->policy & engine_policy_drift_all && !e->forcerebuild)
    engine_drift_all(e, /*drift_mpole=*/1);

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
    // dumpCells("cells", 1, 0, 0, 0, e->s, e->nodeID, e->step);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].grav.multipole->m_pole.num_gpart;
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

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure all woken-up particles have been processed */
  space_check_limiter(e->s);
#endif

  /* Collect information about the next time-step */
  engine_collect_end_of_step(e, 1);
  e->forcerebuild = e->collect_group1.forcerebuild;
  e->updates_since_rebuild += e->collect_group1.updated;
  e->g_updates_since_rebuild += e->collect_group1.g_updated;
  e->s_updates_since_rebuild += e->collect_group1.s_updated;

#ifdef SWIFT_DEBUG_CHECKS
  if (e->ti_end_min == e->ti_current && e->ti_end_min < max_nr_timesteps)
    error("Obtained a time-step of size 0");
#endif

  /********************************************************/
  /* OK, we are done with the regular stuff. Time for i/o */
  /********************************************************/

  /* Create a restart file if needed. */
  engine_dump_restarts(e, 0, e->restart_onexit && engine_is_done(e));

  engine_check_for_dumps(e);

  TIMER_TOC2(timer_step);

  clocks_gettime(&time2);
  e->wallclock_time = (float)clocks_diff(&time1, &time2);

  /* Time in ticks at the end of this step. */
  e->toc_step = getticks();
}

/**
 * @brief Check whether any kind of i/o has to be performed during this
 * step.
 *
 * This includes snapshots, stats and halo finder. We also handle the case
 * of multiple outputs between two steps.
 *
 * @param e The #engine.
 */
void engine_check_for_dumps(struct engine *e) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_stf = (e->policy & engine_policy_structure_finding);

  /* What kind of output are we getting? */
  enum output_type {
    output_none,
    output_snapshot,
    output_statistics,
    output_stf
  };

  /* What kind of output do we want? And at which time ?
   * Find the earliest output (amongst all kinds) that takes place
   * before the next time-step */
  enum output_type type = output_none;
  integertime_t ti_output = max_nr_timesteps;

  /* Save some statistics ? */
  if (e->ti_end_min > e->ti_next_stats && e->ti_next_stats > 0) {
    if (e->ti_next_stats < ti_output) {
      ti_output = e->ti_next_stats;
      type = output_statistics;
    }
  }

  /* Do we want a snapshot? */
  if (e->ti_end_min > e->ti_next_snapshot && e->ti_next_snapshot > 0) {
    if (e->ti_next_snapshot < ti_output) {
      ti_output = e->ti_next_snapshot;
      type = output_snapshot;
    }
  }

  /* Do we want to perform structure finding? */
  if (with_stf) {
    if (e->ti_end_min > e->ti_next_stf && e->ti_next_stf > 0) {
      if (e->ti_next_stf < ti_output) {
        ti_output = e->ti_next_stf;
        type = output_stf;
      }
    }
  }

  /* Store information before attempting extra dump-related drifts */
  const integertime_t ti_current = e->ti_current;
  const timebin_t max_active_bin = e->max_active_bin;
  const double time = e->time;

  while (type != output_none) {

    /* Let's fake that we are at the dump time */
    e->ti_current = ti_output;
    e->max_active_bin = 0;
    if (with_cosmology) {
      cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
      e->time = e->cosmology->time;
    } else {
      e->time = ti_output * e->time_base + e->time_begin;
    }

    /* Drift everyone */
    engine_drift_all(e, /*drift_mpole=*/0);

    /* Write some form of output */
    switch (type) {
      case output_snapshot:

        /* Do we want a corresponding VELOCIraptor output? */
        if (with_stf && e->snapshot_invoke_stf) {

#ifdef HAVE_VELOCIRAPTOR
          velociraptor_invoke(e, /*linked_with_snap=*/1);
          e->step_props |= engine_step_prop_stf;
#else
          error(
              "Asking for a VELOCIraptor output but SWIFT was compiled without "
              "the interface!");
#endif
        }

        /* Perform a FOF search. */
        if (e->run_fof) {

          fof_search_tree(e->s);
          e->run_fof = 0;
        }

          /* Dump... */
#ifdef WITH_LOGGER
        /* Write a file containing the offsets in the particle logger. */
        engine_dump_index(e);
#else
        engine_dump_snapshot(e);
#endif

        /* Free the memory allocated for VELOCIraptor i/o. */
        if (with_stf && e->snapshot_invoke_stf) {
#ifdef HAVE_VELOCIRAPTOR
          swift_free("gpart_group_data", e->s->gpart_group_data);
          e->s->gpart_group_data = NULL;
#endif
        }

        /* ... and find the next output time */
        engine_compute_next_snapshot_time(e);
        break;

      case output_statistics:

        /* Dump */
        engine_print_stats(e);

        /* and move on */
        engine_compute_next_statistics_time(e);

        break;

      case output_stf:

#ifdef HAVE_VELOCIRAPTOR
        /* Unleash the raptor! */
        velociraptor_invoke(e, /*linked_with_snap=*/0);
        e->step_props |= engine_step_prop_stf;

        /* ... and find the next output time */
        engine_compute_next_stf_time(e);
#else
        error(
            "Asking for a VELOCIraptor output but SWIFT was compiled without "
            "the interface!");
#endif
        break;

      default:
        error("Invalid dump type");
    }

    /* Perform a FOF search. */
    // if(run_fof) {

    //  // MATTHIEU: Add a drift_all here. And check the order with the order
    //  i/o
    //  // options.
    //
    //  fof_search_tree(e->s);

    //  run_fof = 0;
    //}

    /* We need to see whether whether we are in the pathological case
     * where there can be another dump before the next step. */

    type = output_none;
    ti_output = max_nr_timesteps;

    /* Save some statistics ? */
    if (e->ti_end_min > e->ti_next_stats && e->ti_next_stats > 0) {
      if (e->ti_next_stats < ti_output) {
        ti_output = e->ti_next_stats;
        type = output_statistics;
      }
    }

    /* Do we want a snapshot? */
    if (e->ti_end_min > e->ti_next_snapshot && e->ti_next_snapshot > 0) {
      if (e->ti_next_snapshot < ti_output) {
        ti_output = e->ti_next_snapshot;
        type = output_snapshot;
      }
    }

    /* Do we want to perform structure finding? */
    if (with_stf) {
      if (e->ti_end_min > e->ti_next_stf && e->ti_next_stf > 0) {
        if (e->ti_next_stf < ti_output) {
          ti_output = e->ti_next_stf;
          type = output_stf;
        }
      }
    }

  } /* While loop over output types */

  /* Restore the information we stored */
  e->ti_current = ti_current;
  if (e->policy & engine_policy_cosmology)
    cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
  e->max_active_bin = max_active_bin;
  e->time = time;
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
      if (!drifted_all) engine_drift_all(e, /*drift_mpole=*/1);
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
  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if ((e->policy & engine_policy_hydro && cell_is_active_hydro(c, e)) ||
        (e->policy & engine_policy_self_gravity &&
         cell_is_active_gravity(c, e)) ||
        (e->policy & engine_policy_external_gravity &&
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
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit_inv = e->gravity_properties->theta_crit_inv;
  const double theta_crit2 = e->gravity_properties->theta_crit2;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from a shifted CoM to centre of cell */
  const double delta_CoM = engine_max_proxy_centre_frac * r_diag;

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = r_diag + 2. * delta_CoM;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Compute how many cells away we need to walk */
  int delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max * theta_crit_inv;
    delta_cells = (int)(distance / cells[0].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
    }
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d top-level cells away (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Loop over all its neighbours neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Early abort  */
              if (cid >= cjd) continue;

              /* Early abort (both same node) */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID == nodeID)
                continue;

              /* Early abort (both foreign node) */
              if (cells[cid].nodeID != nodeID && cells[cjd].nodeID != nodeID)
                continue;

              int proxy_type = 0;

              /* In the hydro case, only care about direct neighbours */
              if (with_hydro) {

                // MATTHIEU: to do: Write a better expression for the
                // non-periodic case.

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1)))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* First just add the direct neighbours. Then look for
                   some further out if the opening angle demands it */

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1))) {

                  proxy_type |= (int)proxy_cell_type_gravity;
                } else {

                  /* We don't have multipoles yet (or there CoMs) so we will
                     have to cook up something based on cell locations only. We
                     hence need an upper limit on the distance that the CoMs in
                     those cells could have. We then can decide whether we are
                     too close for an M2L interaction and hence require a proxy
                     as this pair of cells cannot rely on just an M2L
                     calculation. */

                  /* Minimal distance between any two points in the cells */
                  const double min_dist_centres2 = cell_min_dist2_same_size(
                      &cells[cid], &cells[cjd], periodic, dim);

                  /* Let's now assume the CoMs will shift a bit */
                  const double min_dist_CoM =
                      sqrt(min_dist_centres2) - 2. * delta_CoM;
                  const double min_dist_CoM2 = min_dist_CoM * min_dist_CoM;

                  /* Are we beyond the distance where the truncated forces are 0
                   * but not too far such that M2L can be used? */
                  if (periodic) {

                    if ((min_dist_CoM2 < max_mesh_dist2) &&
                        (!gravity_M2L_accept(r_max, r_max, theta_crit2,
                                             min_dist_CoM2)))
                      proxy_type |= (int)proxy_cell_type_gravity;

                  } else {

                    if (!gravity_M2L_accept(r_max, r_max, theta_crit2,
                                            min_dist_CoM2))
                      proxy_type |= (int)proxy_cell_type_gravity;
                  }
                }
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cells[cjd].nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cjd].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cjd].nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                cells[cid].mpi.sendto |= (1ULL << proxy_id);
              }

              /* Same for the symmetric case? */
              if (cells[cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cells[cid].nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cid].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
                  proxy_id = e->nr_proxies;
                  e->nr_proxies += 1;

                  /* Check the maximal proxy limit */
                  if ((size_t)proxy_id > 8 * sizeof(long long))
                    error(
                        "Created more than %zd proxies. cell.mpi.sendto will "
                        "overflow.",
                        8 * sizeof(long long));
                }

                /* Add the cell to the proxy */
                proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cells[cjd].mpi.sendto |= (1ULL << proxy_id);
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
  const ticks tic = getticks();

  struct space *s = e->s;

  /* Do the initial partition of the cells. */
  partition_initial_partition(initial_partition, e->nodeID, e->nr_nodes, s);

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Re-allocate the local parts. */
  if (e->verbose)
    message("Re-allocating parts array from %zu to %zu.", s->size_parts,
            (size_t)(s->nr_parts * engine_redistribute_alloc_margin));
  s->size_parts = s->nr_parts * engine_redistribute_alloc_margin;
  struct part *parts_new = NULL;
  struct xpart *xparts_new = NULL;
  if (swift_memalign("parts", (void **)&parts_new, part_align,
                     sizeof(struct part) * s->size_parts) != 0 ||
      swift_memalign("xparts", (void **)&xparts_new, xpart_align,
                     sizeof(struct xpart) * s->size_parts) != 0)
    error("Failed to allocate new part data.");

  if (s->nr_parts > 0) {
    memcpy(parts_new, s->parts, sizeof(struct part) * s->nr_parts);
    memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->nr_parts);
  }
  swift_free("parts", s->parts);
  swift_free("xparts", s->xparts);
  s->parts = parts_new;
  s->xparts = xparts_new;

  /* Re-link the gparts to their parts. */
  if (s->nr_parts > 0 && s->nr_gparts > 0)
    part_relink_gparts_to_parts(s->parts, s->nr_parts, 0);

  /* Re-allocate the local sparts. */
  if (e->verbose)
    message("Re-allocating sparts array from %zu to %zu.", s->size_sparts,
            (size_t)(s->nr_sparts * engine_redistribute_alloc_margin));
  s->size_sparts = s->nr_sparts * engine_redistribute_alloc_margin;
  struct spart *sparts_new = NULL;
  if (swift_memalign("sparts", (void **)&sparts_new, spart_align,
                     sizeof(struct spart) * s->size_sparts) != 0)
    error("Failed to allocate new spart data.");

  if (s->nr_sparts > 0)
    memcpy(sparts_new, s->sparts, sizeof(struct spart) * s->nr_sparts);
  swift_free("sparts", s->sparts);
  s->sparts = sparts_new;

  /* Re-link the gparts to their sparts. */
  if (s->nr_sparts > 0 && s->nr_gparts > 0)
    part_relink_gparts_to_sparts(s->sparts, s->nr_sparts, 0);

  /* Re-allocate the local gparts. */
  if (e->verbose)
    message("Re-allocating gparts array from %zu to %zu.", s->size_gparts,
            (size_t)(s->nr_gparts * engine_redistribute_alloc_margin));
  s->size_gparts = s->nr_gparts * engine_redistribute_alloc_margin;
  struct gpart *gparts_new = NULL;
  if (swift_memalign("gparts", (void **)&gparts_new, gpart_align,
                     sizeof(struct gpart) * s->size_gparts) != 0)
    error("Failed to allocate new gpart data.");

  if (s->nr_gparts > 0)
    memcpy(gparts_new, s->gparts, sizeof(struct gpart) * s->nr_gparts);
  swift_free("gparts", s->gparts);
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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

#ifdef DEBUG_INTERACTIONS_STARS
/**
 * @brief Exchange the feedback counters between stars
 * @param e The #engine.
 */
void engine_collect_stars_counter(struct engine *e) {

#ifdef WITH_MPI
  if (e->total_nr_sparts > 1e5) {
    message("WARNING: too many sparts, skipping exchange of counters");
    return;
  }

  /* Get number of sparticles for each rank */
  size_t *n_sparts = (size_t *)malloc(e->nr_nodes * sizeof(size_t));

  int err = MPI_Allgather(&e->s->nr_sparts_foreign, 1, MPI_UNSIGNED_LONG,
                          n_sparts, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) error("Communication failed");

  /* Compute derivated quantities */
  int total = 0;
  int *n_sparts_int = (int *)malloc(e->nr_nodes * sizeof(int));
  int *displs = (int *)malloc(e->nr_nodes * sizeof(int));
  for (int i = 0; i < e->nr_nodes; i++) {
    displs[i] = total;
    total += n_sparts[i];
    n_sparts_int[i] = n_sparts[i];
  }

  /* Get all sparticles */
  struct spart *sparts =
      (struct spart *)swift_malloc("sparts", total * sizeof(struct spart));
  err = MPI_Allgatherv(e->s->sparts_foreign, e->s->nr_sparts_foreign,
                       spart_mpi_type, sparts, n_sparts_int, displs,
                       spart_mpi_type, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) error("Communication failed");

  /* Reset counters */
  for (size_t i = 0; i < e->s->nr_sparts_foreign; i++) {
    e->s->sparts_foreign[i].num_ngb_force = 0;
  }

  /* Update counters */
  struct spart *local_sparts = e->s->sparts;
  for (size_t i = 0; i < e->s->nr_sparts; i++) {
    const long long id_i = local_sparts[i].id;

    for (int j = 0; j < total; j++) {
      const long long id_j = sparts[j].id;

      if (id_j == id_i) {
        if (j >= displs[engine_rank] &&
            j < displs[engine_rank] + n_sparts_int[engine_rank]) {
          error(
              "Found a local spart in foreign cell ID=%lli: j=%i, displs=%i, "
              "n_sparts=%i",
              id_j, j, displs[engine_rank], n_sparts_int[engine_rank]);
        }

        local_sparts[i].num_ngb_force += sparts[j].num_ngb_force;
      }
    }
  }

  free(n_sparts);
  free(n_sparts_int);
  swift_free("sparts", sparts);
#endif
}

#endif

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
  space_check_drift_point(e->s, e->ti_current, /* check_mpole=*/0);

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

#ifdef DEBUG_INTERACTIONS_STARS
  engine_collect_stars_counter(e);
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

/**
 * @brief Writes an index file with the current state of the engine
 *
 * @param e The #engine.
 */
void engine_dump_index(struct engine *e) {

#if defined(WITH_LOGGER)
  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Writing index at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Writing index at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }

  /* Dump... */
  write_index_single(e, e->logger->base_name, e->internal_units,
                     e->snapshot_units);

  /* Flag that we dumped a snapshot */
  e->step_props |= engine_step_prop_logger_index;

  clocks_gettime(&time2);
  if (e->verbose)
    message("writing particle indices took %.3f %s.",
            (float)clocks_diff(&time1, &time2), clocks_getunit());
#else
  error("SWIFT was not compiled with the logger");
#endif
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
 * @param entropy_floor The #entropy_floor_properties for this run.
 * @param gravity The #gravity_props used for this run.
 * @param stars The #stars_props used for this run.
 * @param mesh The #pm_mesh used for the long-range periodic forces.
 * @param potential The properties of the external potential.
 * @param cooling_func The properties of the cooling function.
 * @param starform The #star_formation model of this run.
 * @param chemistry The chemistry information.
 */
void engine_init(struct engine *e, struct space *s, struct swift_params *params,
                 long long Ngas, long long Ngparts, long long Nstars,
                 int policy, int verbose, struct repartition *reparttype,
                 const struct unit_system *internal_units,
                 const struct phys_const *physical_constants,
                 struct cosmology *cosmo, struct hydro_props *hydro,
                 const struct entropy_floor_properties *entropy_floor,
                 struct gravity_props *gravity, const struct stars_props *stars,
                 struct pm_mesh *mesh,
                 const struct external_potential *potential,
                 struct cooling_function_data *cooling_func,
                 const struct star_formation *starform,
                 const struct chemistry_global_data *chemistry) {

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
  e->snapshot_int_time_label_on =
      parser_get_opt_param_int(params, "Snapshots:int_time_label_on", 0);
  e->snapshot_invoke_stf =
      parser_get_opt_param_int(params, "Snapshots:invoke_stf", 0);
  e->snapshot_units = (struct unit_system *)malloc(sizeof(struct unit_system));
  units_init_default(e->snapshot_units, params, "Snapshots", internal_units);
  e->snapshot_output_count = 0;
  e->stf_output_count = 0;
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
  e->entropy_floor = entropy_floor;
  e->gravity_properties = gravity;
  e->stars_properties = stars;
  e->mesh = mesh;
  e->external_potential = potential;
  e->cooling_func = cooling_func;
  e->star_formation = starform;
  e->chemistry = chemistry;
  e->parameter_file = params;
#ifdef WITH_MPI
  e->cputime_last_step = 0;
  e->last_repartition = 0;
#endif
  e->total_nr_cells = 0;
  e->total_nr_tasks = 0;

#if defined(WITH_LOGGER)
  e->logger = (struct logger *)malloc(sizeof(struct logger));
  logger_init(e->logger, params);
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

  /* Initialise VELOCIraptor output. */
  if (e->policy & engine_policy_structure_finding) {
    parser_get_param_string(params, "StructureFinding:basename",
                            e->stf_base_name);
    parser_get_param_string(params, "StructureFinding:config_file_name",
                            e->stf_config_file_name);

    e->time_first_stf_output =
        parser_get_opt_param_double(params, "StructureFinding:time_first", 0.);
    e->a_first_stf_output = parser_get_opt_param_double(
        params, "StructureFinding:scale_factor_first", 0.1);
    e->delta_time_stf =
        parser_get_opt_param_double(params, "StructureFinding:delta_time", -1.);
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
  e->run_fof = 0;
  engine_rank = nodeID;

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
              "Statistics=%d, Snapshot=%d, Restarts=%d STF=%d, logger=%d\n",
              engine_step_prop_rebuild, engine_step_prop_redistribute,
              engine_step_prop_repartition, engine_step_prop_statistics,
              engine_step_prop_snapshot, engine_step_prop_restarts,
              engine_step_prop_stf, engine_step_prop_logger_index);

      fprintf(e->file_timesteps,
              "# %6s %14s %12s %12s %14s %9s %12s %12s %12s %16s [%s] %6s\n",
              "Step", "Time", "Scale-factor", "Redshift", "Time-step",
              "Time-bins", "Updates", "g-Updates", "s-Updates",
              "Wall-clock time", clocks_getunit(), "Props");
      fflush(e->file_timesteps);
    }
  }

  /* Print policy */
  engine_print_policy(e);

  /* Print information about the hydro scheme */
  if (e->policy & engine_policy_hydro) {
    if (e->nodeID == 0) hydro_props_print(e->hydro_properties);
    if (e->nodeID == 0) entropy_floor_print(e->entropy_floor);
  }

  /* Print information about the gravity scheme */
  if (e->policy & engine_policy_self_gravity)
    if (e->nodeID == 0) gravity_props_print(e->gravity_properties);

  if (e->policy & engine_policy_stars)
    if (e->nodeID == 0) stars_props_print(e->stars_properties);

  /* Check we have sensible time bounds */
  if (e->time_begin >= e->time_end)
    error(
        "Final simulation time (t_end = %e) must be larger than the start time "
        "(t_beg = %e)",
        e->time_end, e->time_begin);

  /* Check we don't have inappropriate time labels */
  if ((e->snapshot_int_time_label_on == 1) && (e->time_end <= 1.f))
    error("Snapshot integer time labels enabled but end time <= 1");

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

    if (e->policy & engine_policy_structure_finding) {

      if (e->delta_time_stf == -1. && !e->snapshot_invoke_stf)
        error("A value for `StructureFinding:delta_time` must be specified");

      if (e->delta_time_stf <= 1. && e->delta_time_stf != -1.)
        error("Time between STF (%e) must be > 1.", e->delta_time_stf);

      if (e->a_first_stf_output < e->cosmology->a_begin)
        error(
            "Scale-factor of first stf output (%e) must be after the "
            "simulation start a=%e.",
            e->a_first_stf_output, e->cosmology->a_begin);
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

    if (e->policy & engine_policy_structure_finding) {

      if (e->delta_time_stf == -1. && !e->snapshot_invoke_stf)
        error("A value for `StructureFinding:delta_time` must be specified");

      if (e->delta_time_stf <= 0. && e->delta_time_stf != -1.)
        error("Time between STF (%e) must be positive.", e->delta_time_stf);

      if (e->time_first_stf_output < e->time_begin)
        error("Time of first STF (%e) must be after the simulation start t=%e.",
              e->time_first_stf_output, e->time_begin);
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

#if defined(WITH_LOGGER)
  if (e->nodeID == 0)
    message(
        "WARNING: There is currently no way of predicting the output "
        "size, please use it carefully");
#endif

  /* Find the time of the first snapshot output */
  engine_compute_next_snapshot_time(e);

  /* Find the time of the first statistics output */
  engine_compute_next_statistics_time(e);

  /* Find the time of the first stf output */
  if (e->policy & engine_policy_structure_finding) {
    engine_compute_next_stf_time(e);
  }

  /* Check that we are invoking VELOCIraptor only if we have it */
  if (e->snapshot_invoke_stf &&
      !(e->policy & engine_policy_structure_finding)) {
    error(
        "Invoking VELOCIraptor after snapshots but structure finding wasn't "
        "activated at runtime (Use --velociraptor).");
  }

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
  multipole_create_mpi_types();
  stats_create_mpi_type();
  proxy_create_mpi_type();
  task_create_mpi_comms();
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
   * from the end of the dumped run. Can be changed on restart. */
  e->tasks_per_cell =
      parser_get_opt_param_float(params, "Scheduler:tasks_per_cell", 0.0);
  e->tasks_per_cell_max = 0.0f;

  float maxtasks = 0;
  if (restart)
    maxtasks = e->restart_max_tasks;
  else
    maxtasks = engine_estimate_nr_tasks(e);

  /* Estimated number of links per tasks */
  e->links_per_tasks =
      parser_get_opt_param_int(params, "Scheduler:links_per_tasks", 10);

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
  if (swift_memalign("runners", (void **)&e->runners, SWIFT_CACHE_ALIGNMENT,
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

#ifdef WITH_LOGGER
  /* Write the particle logger header */
  logger_write_file_header(e->logger, e);
#endif

  /* Initialise the structure finder */
#ifdef HAVE_VELOCIRAPTOR
  if (e->policy & engine_policy_structure_finding) velociraptor_init(e);
#endif

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

  int found_snapshot_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_snapshot = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_snapshot = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_snapshot > e->ti_current) {
      found_snapshot_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_snapshot;
    else
      time += e->delta_time_snapshot;
  }

  /* Deal with last snapshot */
  if (!found_snapshot_time) {
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

  int found_stats_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_stats = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_stats = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_stats > e->ti_current) {
      found_stats_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_statistics;
    else
      time += e->delta_time_statistics;
  }

  /* Deal with last statistics */
  if (!found_stats_time) {
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
    output_list_read_next_time(e->output_list_stf, e, "stf", &e->ti_next_stf);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_stf;
  else
    time_end = e->time_end + e->delta_time_stf;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_stf_output;
  else
    time = e->time_first_stf_output;

  int found_stf_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_stf = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_stf = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_stf > e->ti_current) {
      found_stf_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_stf;
    else
      time += e->delta_time_stf;
  }

  /* Deal with last snapshot */
  if (!found_stf_time) {
    e->ti_next_stf = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const float next_stf_time =
          exp(e->ti_next_stf * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next VELOCIraptor time set to a=%e.", next_stf_time);
    } else {
      const float next_stf_time = e->ti_next_stf * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next VELOCIraptor time set to t=%e.", next_stf_time);
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
  output_list_init(&e->output_list_stf, e, "StructureFinding",
                   &e->delta_time_stf, &stf_time_first);

  if (e->output_list_stf) {
    if (e->policy & engine_policy_cosmology)
      e->a_first_stf_output = stf_time_first;
    else
      e->time_first_stf_output = stf_time_first;
  }
}

/**
 * @brief Computes the maximal time-step allowed by the max RMS displacement
 * condition.
 *
 * @param e The #engine.
 */
void engine_recompute_displacement_constraint(struct engine *e) {

  const ticks tic = getticks();

  /* Get the cosmological information */
  const struct cosmology *cosmo = e->cosmology;
  const float Om = cosmo->Omega_m;
  const float Ob = cosmo->Omega_b;
  const float H0 = cosmo->H0;
  const float a = cosmo->a;
  const float G_newton = e->physical_constants->const_newton_G;
  const float rho_crit0 = 3.f * H0 * H0 / (8.f * M_PI * G_newton);

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
  for (size_t i = 0; i < e->s->nr_parts; ++i) {
    if (e->s->parts[i].time_bin >= num_time_bins) continue;
    min_part_mass_check =
        min(min_part_mass_check, hydro_get_mass(&e->s->parts[i]));
  }
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
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
    r_s = e->mesh->r_s;
  else
    r_s = FLT_MAX;

  float dt_dm = FLT_MAX, dt_b = FLT_MAX;

  /* DM case */
  if (N_dm > 0.f) {

    /* Minimal mass for the DM */
    const float min_mass_dm = min_mass[1];

    /* Inter-particle sepration for the DM */
    const float d_dm = cbrtf(min_mass_dm / ((Om - Ob) * rho_crit0));

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
    const float d_b = cbrtf(min_mass_b / (Ob * rho_crit0));

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

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
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
  swift_free("runners", e->runners);
  free(e->snapshot_units);

  output_list_clean(&e->output_list_snapshots);
  output_list_clean(&e->output_list_stats);
  output_list_clean(&e->output_list_stf);

  swift_free("links", e->links);
#if defined(WITH_LOGGER)
  logger_clean(e->logger);
  free(e->logger);
#endif
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
  entropy_floor_struct_dump(e->entropy_floor, stream);
  gravity_props_struct_dump(e->gravity_properties, stream);
  stars_props_struct_dump(e->stars_properties, stream);
  pm_mesh_struct_dump(e->mesh, stream);
  potential_struct_dump(e->external_potential, stream);
  cooling_struct_dump(e->cooling_func, stream);
  starformation_struct_dump(e->star_formation, stream);
  chemistry_struct_dump(e->chemistry, stream);
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

  struct entropy_floor_properties *entropy_floor =
      (struct entropy_floor_properties *)malloc(
          sizeof(struct entropy_floor_properties));
  entropy_floor_struct_restore(entropy_floor, stream);
  e->entropy_floor = entropy_floor;

  struct gravity_props *gravity_properties =
      (struct gravity_props *)malloc(sizeof(struct gravity_props));
  gravity_props_struct_restore(gravity_properties, stream);
  e->gravity_properties = gravity_properties;

  struct stars_props *stars_properties =
      (struct stars_props *)malloc(sizeof(struct stars_props));
  stars_props_struct_restore(stars_properties, stream);
  e->stars_properties = stars_properties;

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
  cooling_struct_restore(cooling_func, stream, e->cosmology);
  e->cooling_func = cooling_func;

  struct star_formation *star_formation =
      (struct star_formation *)malloc(sizeof(struct star_formation));
  starformation_struct_restore(star_formation, stream);
  e->star_formation = star_formation;

  struct chemistry_global_data *chemistry =
      (struct chemistry_global_data *)malloc(
          sizeof(struct chemistry_global_data));
  chemistry_struct_restore(chemistry, stream);
  e->chemistry = chemistry;

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

#ifdef EOS_PLANETARY
  eos_init(&eos, e->physical_constants, e->snapshot_units, e->parameter_file);
#endif

  /* Want to force a rebuild before using this engine. Wait to repartition.*/
  e->forcerebuild = 1;
  e->forcerepart = 0;
}
