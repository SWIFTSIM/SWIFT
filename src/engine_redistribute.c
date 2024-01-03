/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "memswap.h"

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
 * @param syncredist whether to use slower more memory friendly synchronous
 *                   exchanges.
 *
 * @result new particle data constructed from all the exchanges with the
 *         given alignment.
 */
static void *engine_do_redistribute(const char *label, int *counts, char *parts,
                                    size_t new_nr_parts, size_t sizeofparts,
                                    size_t alignsize, MPI_Datatype mpi_type,
                                    int nr_nodes, int nodeID, int syncredist) {

  /* Allocate a new particle array with some extra margin */
  char *parts_new = NULL;
  if (swift_memalign(
          label, (void **)&parts_new, alignsize,
          sizeofparts * new_nr_parts * engine_redistribute_alloc_margin) != 0)
    error("Failed to allocate new particle data.");

  if (syncredist) {

    /* Slow synchronous redistribute,. */
    size_t offset_send = 0, offset_recv = 0;

    /* Only send and receive only "chunk" particles per request.
     * Fixing the message size to 2GB. */
    const int chunk = INT_MAX / sizeofparts;
    int res = 0;
    for (int k = 0; k < nr_nodes; k++) {
      int kk = k;

      /* Rank 0 decides the index of sending node */
      MPI_Bcast(&kk, 1, MPI_INT, 0, MPI_COMM_WORLD);

      int ind_recv = kk * nr_nodes + nodeID;

      if (nodeID == kk) {

        /*  Send out our particles. */
        offset_send = 0;
        for (int j = 0; j < nr_nodes; j++) {

          int ind_send = kk * nr_nodes + j;

          /*  Just copy our own parts */
          if (counts[ind_send] > 0) {
            if (j == nodeID) {
              memcpy(&parts_new[offset_recv * sizeofparts],
                     &parts[offset_send * sizeofparts],
                     sizeofparts * counts[ind_recv]);
              offset_send += counts[ind_send];
              offset_recv += counts[ind_recv];
            } else {
              for (int i = 0, n = 0; i < counts[ind_send]; n++) {

                /* Count and index, with chunk parts at most. */
                size_t sendc = min(chunk, counts[ind_send] - i);
                size_t sendo = offset_send + i;

                res = MPI_Send(&parts[sendo * sizeofparts], sendc, mpi_type, j,
                               n, MPI_COMM_WORLD);
                if (res != MPI_SUCCESS) {
                  mpi_error(res, "Failed to send parts to node %i from %i.", j,
                            nodeID);
                }
                i += sendc;
              }
              offset_send += counts[ind_send];
            }
          }
        }
      } else {
        /*  Listen for sends from kk. */
        if (counts[ind_recv] > 0) {
          for (int i = 0, n = 0; i < counts[ind_recv]; n++) {
            /* Count and index, with +chunk parts at most. */
            size_t recvc = min(chunk, counts[ind_recv] - i);
            size_t recvo = offset_recv + i;

            MPI_Status status;
            res = MPI_Recv(&parts_new[recvo * sizeofparts], recvc, mpi_type, kk,
                           n, MPI_COMM_WORLD, &status);
            if (res != MPI_SUCCESS) {
              mpi_error(res, "Failed to recv of parts from node %i to %i.", kk,
                        nodeID);
            }
            i += recvc;
          }
          offset_recv += counts[ind_recv];
        }
      }
    }

  } else {
    /* Asynchronous redistribute, can take a lot of memory. */

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
                MPI_Isend(&parts[offset_send * sizeofparts], sending, mpi_type,
                          k, ind_send, MPI_COMM_WORLD, &reqs[2 * k + 0]);
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
            int res = MPI_Irecv(&parts_new[offset_recv * sizeofparts],
                                receiving, mpi_type, k, ind_recv,
                                MPI_COMM_WORLD, &reqs[2 * k + 1]);
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
  }

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
        if (parts[k].x[j] == s->dim[j]) parts[k].x[j] = 0.0;               \
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
void ENGINE_REDISTRIBUTE_DEST_MAPPER(part);

/**
 * @brief Accumulate the counts of star particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * spart version.
 */
void ENGINE_REDISTRIBUTE_DEST_MAPPER(spart);

/**
 * @brief Accumulate the counts of gravity particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * gpart version.
 */
void ENGINE_REDISTRIBUTE_DEST_MAPPER(gpart);

/**
 * @brief Accumulate the counts of black holes particles per cell.
 * Threadpool helper for accumulating the counts of particles per cell.
 *
 * bpart version.
 */
void ENGINE_REDISTRIBUTE_DEST_MAPPER(bpart);

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
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(part, 1);
#else
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(part, 0);
#endif

/**
 * @brief Save position of spart-gpart links.
 * Threadpool helper for accumulating the counts of particles per cell.
 */
#ifdef SWIFT_DEBUG_CHECKS
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(spart, 1);
#else
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(spart, 0);
#endif

/**
 * @brief Save position of bpart-gpart links.
 * Threadpool helper for accumulating the counts of particles per cell.
 */
#ifdef SWIFT_DEBUG_CHECKS
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(bpart, 1);
#else
void ENGINE_REDISTRIBUTE_SAVELINK_MAPPER(bpart, 0);
#endif

#endif /* savelink_mapper_data */

#ifdef WITH_MPI /* relink_mapper_data */

/* Support for relinking parts, gparts, sparts and bparts after moving between
 * nodes. */
struct relink_mapper_data {
  int nodeID;
  int nr_nodes;
  int *counts;
  int *s_counts;
  int *g_counts;
  int *b_counts;
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
void engine_redistribute_relink_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  int *nodes = (int *)map_data;
  struct relink_mapper_data *mydata = (struct relink_mapper_data *)extra_data;

  int nodeID = mydata->nodeID;
  int nr_nodes = mydata->nr_nodes;
  int *counts = mydata->counts;
  int *g_counts = mydata->g_counts;
  int *s_counts = mydata->s_counts;
  int *b_counts = mydata->b_counts;
  struct space *s = mydata->s;

  for (int i = 0; i < num_elements; i++) {

    int node = nodes[i];

    /* Get offsets to correct parts of the counts arrays for this node. */
    size_t offset_parts = 0;
    size_t offset_gparts = 0;
    size_t offset_sparts = 0;
    size_t offset_bparts = 0;
    for (int n = 0; n < node; n++) {
      int ind_recv = n * nr_nodes + nodeID;
      offset_parts += counts[ind_recv];
      offset_gparts += g_counts[ind_recv];
      offset_sparts += s_counts[ind_recv];
      offset_bparts += b_counts[ind_recv];
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

      /* Does this gpart have a black hole partner ? */
      else if (s->gparts[k].type == swift_type_black_hole) {

        const ptrdiff_t partner_index =
            offset_bparts - s->gparts[k].id_or_neg_offset;

        /* Re-link */
        s->gparts[k].id_or_neg_offset = -partner_index;
        s->bparts[partner_index].gpart = &s->gparts[k];
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
 * 5) Asynchronous or synchronous communications are issued to transfer the
 * data.
 *
 *
 * @param e The #engine.
 */
void engine_redistribute(struct engine *e) {

#ifdef WITH_MPI
#ifdef SWIFT_DEBUG_CHECKS
  const int nr_sinks_new = 0;
#endif
  if (e->policy & engine_policy_sinks) {
    error("Not implemented yet");
  }

  const int nr_nodes = e->nr_nodes;
  const int nodeID = e->nodeID;
  struct space *s = e->s;
  struct cell *cells = s->cells_top;
  const int nr_cells = s->nr_cells;
  struct xpart *xparts = s->xparts;
  struct part *parts = s->parts;
  struct gpart *gparts = s->gparts;
  struct spart *sparts = s->sparts;
  struct bpart *bparts = s->bparts;
  ticks tic = getticks();

  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;
  size_t nr_bparts = s->nr_bparts;

  /* Start by moving inhibited particles to the end of the arrays */
  for (size_t k = 0; k < nr_parts; /* void */) {
    if (parts[k].time_bin == time_bin_inhibited ||
        parts[k].time_bin == time_bin_not_created) {
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
    if (sparts[k].time_bin == time_bin_inhibited ||
        sparts[k].time_bin == time_bin_not_created) {
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

  /* Now move inhibited black hole particles to the end of the arrays */
  for (size_t k = 0; k < nr_bparts; /* void */) {
    if (bparts[k].time_bin == time_bin_inhibited ||
        bparts[k].time_bin == time_bin_not_created) {
      nr_bparts -= 1;

      /* Swap the particle */
      memswap(&s->bparts[k], &s->bparts[nr_bparts], sizeof(struct bpart));

      /* Swap the link with the gpart */
      if (s->bparts[k].gpart != NULL) {
        s->bparts[k].gpart->id_or_neg_offset = -k;
      }
      if (s->bparts[nr_bparts].gpart != NULL) {
        s->bparts[nr_bparts].gpart->id_or_neg_offset = -nr_bparts;
      }
    } else {
      k++;
    }
  }

  /* Finally do the same with the gravity particles */
  for (size_t k = 0; k < nr_gparts; /* void */) {
    if (gparts[k].time_bin == time_bin_inhibited ||
        gparts[k].time_bin == time_bin_not_created) {
      nr_gparts -= 1;

      /* Swap the particle */
      memswap_unaligned(&s->gparts[k], &s->gparts[nr_gparts],
                        sizeof(struct gpart));

      /* Swap the link with part/spart */
      if (s->gparts[k].type == swift_type_gas) {
        s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_stars) {
        s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      } else if (s->gparts[k].type == swift_type_black_hole) {
        s->bparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
      }

      if (s->gparts[nr_gparts].type == swift_type_gas) {
        s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      } else if (s->gparts[nr_gparts].type == swift_type_stars) {
        s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
            &s->gparts[nr_gparts];
      } else if (s->gparts[nr_gparts].type == swift_type_black_hole) {
        s->bparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
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
                 nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                 &redist_data);

  /* Sort the particles according to their cell index. */
  if (nr_parts > 0)
    space_parts_sort(s->parts, s->xparts, dest, &counts[nodeID * nr_nodes],
                     nr_nodes, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < nr_parts; k++) {
    const struct part *p = &s->parts[k];

    if (p->time_bin == time_bin_inhibited)
      error("Inhibited particle found after sorting!");

    if (p->time_bin == time_bin_not_created)
      error("Inhibited particle found after sorting!");

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
                   nodes, nr_nodes, sizeof(int), threadpool_auto_chunk_size,
                   &savelink_data);
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
                 nr_sparts, sizeof(struct spart), threadpool_auto_chunk_size,
                 &redist_data);

  /* Sort the particles according to their cell index. */
  if (nr_sparts > 0)
    space_sparts_sort(s->sparts, s_dest, &s_counts[nodeID * nr_nodes], nr_nodes,
                      0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < nr_sparts; k++) {
    const struct spart *sp = &s->sparts[k];

    if (sp->time_bin == time_bin_inhibited)
      error("Inhibited particle found after sorting!");

    if (sp->time_bin == time_bin_not_created)
      error("Inhibited particle found after sorting!");

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
                   nodes, nr_nodes, sizeof(int), threadpool_auto_chunk_size,
                   &savelink_data);
  }
  swift_free("s_dest", s_dest);

  /* Get destination of each b-particle */
  int *b_counts;
  if ((b_counts = (int *)calloc(sizeof(int), nr_nodes * nr_nodes)) == NULL)
    error("Failed to allocate b_counts temporary buffer.");

  int *b_dest;
  if ((b_dest = (int *)swift_malloc("b_dest", sizeof(int) * nr_bparts)) == NULL)
    error("Failed to allocate b_dest temporary buffer.");

  redist_data.counts = b_counts;
  redist_data.dest = b_dest;
  redist_data.base = (void *)bparts;

  threadpool_map(&e->threadpool, engine_redistribute_dest_mapper_bpart, bparts,
                 nr_bparts, sizeof(struct bpart), threadpool_auto_chunk_size,
                 &redist_data);

  /* Sort the particles according to their cell index. */
  if (nr_bparts > 0)
    space_bparts_sort(s->bparts, b_dest, &b_counts[nodeID * nr_nodes], nr_nodes,
                      0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the bpart have been sorted correctly. */
  for (size_t k = 0; k < nr_bparts; k++) {
    const struct bpart *bp = &s->bparts[k];

    if (bp->time_bin == time_bin_inhibited)
      error("Inhibited particle found after sorting!");

    if (bp->time_bin == time_bin_not_created)
      error("Inhibited particle found after sorting!");

    /* New cell index */
    const int new_cid =
        cell_getid(s->cdim, bp->x[0] * s->iwidth[0], bp->x[1] * s->iwidth[1],
                   bp->x[2] * s->iwidth[2]);

    /* New cell of this bpart */
    const struct cell *c = &s->cells_top[new_cid];
    const int new_node = c->nodeID;

    if (b_dest[k] != new_node)
      error("bpart's new node index not matching sorted index.");

    if (bp->x[0] < c->loc[0] || bp->x[0] > c->loc[0] + c->width[0] ||
        bp->x[1] < c->loc[1] || bp->x[1] > c->loc[1] + c->width[1] ||
        bp->x[2] < c->loc[2] || bp->x[2] > c->loc[2] + c->width[2])
      error("bpart not sorted into the right top-level cell!");
  }
#endif

  /* We need to re-link the gpart partners of bparts. */
  if (nr_bparts > 0) {

    struct savelink_mapper_data savelink_data;
    savelink_data.nr_nodes = nr_nodes;
    savelink_data.counts = b_counts;
    savelink_data.parts = (void *)bparts;
    savelink_data.nodeID = nodeID;
    threadpool_map(&e->threadpool, engine_redistribute_savelink_mapper_bpart,
                   nodes, nr_nodes, sizeof(int), threadpool_auto_chunk_size,
                   &savelink_data);
  }
  swift_free("b_dest", b_dest);

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
                 nr_gparts, sizeof(struct gpart), threadpool_auto_chunk_size,
                 &redist_data);

  /* Sort the gparticles according to their cell index. */
  if (nr_gparts > 0)
    space_gparts_sort(s->gparts, s->parts, s->sinks, s->sparts, s->bparts,
                      g_dest, &g_counts[nodeID * nr_nodes], nr_nodes);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < nr_gparts; k++) {
    const struct gpart *gp = &s->gparts[k];

    if (gp->time_bin == time_bin_inhibited)
      error("Inhibited particle found after sorting!");

    if (gp->time_bin == time_bin_not_created)
      error("Inhibited particle found after sorting!");

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

  /* Get all the g_counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, g_counts, nr_nodes * nr_nodes, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce gparticle transfer counts.");

  /* Get all the s_counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, s_counts, nr_nodes * nr_nodes, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce sparticle transfer counts.");

  /* Get all the b_counts from all the nodes. */
  if (MPI_Allreduce(MPI_IN_PLACE, b_counts, nr_nodes * nr_nodes, MPI_INT,
                    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to allreduce bparticle transfer counts.");

  /* Report how many particles will be moved. */
  if (e->verbose) {
    if (e->nodeID == 0) {
      size_t total = 0, g_total = 0, s_total = 0, b_total = 0;
      size_t unmoved = 0, g_unmoved = 0, s_unmoved = 0, b_unmoved = 0;
      for (int p = 0, r = 0; p < nr_nodes; p++) {
        for (int n = 0; n < nr_nodes; n++) {
          total += counts[r];
          g_total += g_counts[r];
          s_total += s_counts[r];
          b_total += b_counts[r];
          if (p == n) {
            unmoved += counts[r];
            g_unmoved += g_counts[r];
            s_unmoved += s_counts[r];
            b_unmoved += b_counts[r];
          }
          r++;
        }
      }
      if (total > 0)
        message("%zu of %zu (%.2f%%) of particles moved", total - unmoved,
                total, 100.0 * (double)(total - unmoved) / (double)total);
      if (g_total > 0)
        message("%zu of %zu (%.2f%%) of g-particles moved", g_total - g_unmoved,
                g_total,
                100.0 * (double)(g_total - g_unmoved) / (double)g_total);
      if (s_total > 0)
        message("%zu of %zu (%.2f%%) of s-particles moved", s_total - s_unmoved,
                s_total,
                100.0 * (double)(s_total - s_unmoved) / (double)s_total);
      if (b_total > 0)
        message("%ld of %ld (%.2f%%) of b-particles moved", b_total - b_unmoved,
                b_total,
                100.0 * (double)(b_total - b_unmoved) / (double)b_total);
    }
  }

  /* Now each node knows how many parts, sparts, bparts, and gparts will be
   * transferred to every other node. Get the new numbers of particles for this
   * node. */
  size_t nr_parts_new = 0, nr_gparts_new = 0, nr_sparts_new = 0,
         nr_bparts_new = 0;
  for (int k = 0; k < nr_nodes; k++)
    nr_parts_new += counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_gparts_new += g_counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_sparts_new += s_counts[k * nr_nodes + nodeID];
  for (int k = 0; k < nr_nodes; k++)
    nr_bparts_new += b_counts[k * nr_nodes + nodeID];

#ifdef WITH_CSDS
  const int initial_redistribute = e->ti_current == 0;
  if (!initial_redistribute && e->policy & engine_policy_csds) {
    /* Log the particles before sending them out */
    size_t part_offset = 0;
    size_t spart_offset = 0;
    size_t gpart_offset = 0;
    size_t bpart_offset = 0;

    for (int i = 0; i < nr_nodes; i++) {
      const size_t c_ind = engine_rank * nr_nodes + i;

      /* No need to log the local particles. */
      if (i == engine_rank) {
        part_offset += counts[c_ind];
        spart_offset += s_counts[c_ind];
        gpart_offset += g_counts[c_ind];
        bpart_offset += b_counts[c_ind];
        continue;
      }

      /* Log the hydro parts. */
      csds_log_parts(e->csds, &parts[part_offset], &xparts[part_offset],
                     counts[c_ind], e, /* log_all_fields */ 1,
                     csds_flag_mpi_exit, i);

      /* Log the stellar parts. */
      csds_log_sparts(e->csds, &sparts[spart_offset], s_counts[c_ind], e,
                      /* log_all_fields */ 1, csds_flag_mpi_exit, i);

      /* Log the gparts */
      csds_log_gparts(e->csds, &gparts[gpart_offset], g_counts[c_ind], e,
                      /* log_all_fields */ 1, csds_flag_mpi_exit, i);

      /* Log the bparts */
      if (b_counts[c_ind] > 0) {
        error("TODO");
      }

      /* Update the counters */
      part_offset += counts[c_ind];
      spart_offset += s_counts[c_ind];
      gpart_offset += g_counts[c_ind];
      bpart_offset += b_counts[c_ind];
    }
  }
#endif

  /* Now exchange the particles, type by type to keep the memory required
   * under control. */

  /* SPH particles. */
  void *new_parts = engine_do_redistribute(
      "parts", counts, (char *)s->parts, nr_parts_new, sizeof(struct part),
      part_align, part_mpi_type, nr_nodes, nodeID, e->syncredist);
  swift_free("parts", s->parts);
  s->parts = (struct part *)new_parts;
  s->nr_parts = nr_parts_new;
  s->size_parts = engine_redistribute_alloc_margin * nr_parts_new;

  /* Extra SPH particle properties. */
  new_parts = engine_do_redistribute(
      "xparts", counts, (char *)s->xparts, nr_parts_new, sizeof(struct xpart),
      xpart_align, xpart_mpi_type, nr_nodes, nodeID, e->syncredist);
  swift_free("xparts", s->xparts);
  s->xparts = (struct xpart *)new_parts;

  /* Gravity particles. */
  new_parts =
      engine_do_redistribute("gparts", g_counts, (char *)s->gparts,
                             nr_gparts_new, sizeof(struct gpart), gpart_align,
                             gpart_mpi_type, nr_nodes, nodeID, e->syncredist);
  swift_free("gparts", s->gparts);
  s->gparts = (struct gpart *)new_parts;
  s->nr_gparts = nr_gparts_new;
  s->size_gparts = engine_redistribute_alloc_margin * nr_gparts_new;

  /* Star particles. */
  new_parts =
      engine_do_redistribute("sparts", s_counts, (char *)s->sparts,
                             nr_sparts_new, sizeof(struct spart), spart_align,
                             spart_mpi_type, nr_nodes, nodeID, e->syncredist);
  swift_free("sparts", s->sparts);
  s->sparts = (struct spart *)new_parts;
  s->nr_sparts = nr_sparts_new;
  s->size_sparts = engine_redistribute_alloc_margin * nr_sparts_new;

  /* Black holes particles. */
  new_parts =
      engine_do_redistribute("bparts", b_counts, (char *)s->bparts,
                             nr_bparts_new, sizeof(struct bpart), bpart_align,
                             bpart_mpi_type, nr_nodes, nodeID, e->syncredist);
  swift_free("bparts", s->bparts);
  s->bparts = (struct bpart *)new_parts;
  s->nr_bparts = nr_bparts_new;
  s->size_bparts = engine_redistribute_alloc_margin * nr_bparts_new;

  /* All particles have now arrived. Time for some final operations on the
     stuff we just received */

#ifdef WITH_CSDS
  if (!initial_redistribute && e->policy & engine_policy_csds) {
    size_t part_offset = 0;
    size_t spart_offset = 0;
    size_t gpart_offset = 0;
    size_t bpart_offset = 0;

    for (int i = 0; i < nr_nodes; i++) {
      const size_t c_ind = i * nr_nodes + engine_rank;

      /* No need to log the local particles. */
      if (i == engine_rank) {
        part_offset += counts[c_ind];
        spart_offset += s_counts[c_ind];
        gpart_offset += g_counts[c_ind];
        bpart_offset += b_counts[c_ind];
        continue;
      }

      /* Log the hydro parts. */
      csds_log_parts(e->csds, &s->parts[part_offset], &s->xparts[part_offset],
                     counts[c_ind], e,
                     /* log_all_fields */ 1, csds_flag_mpi_enter, i);

      /* Log the stellar parts. */
      csds_log_sparts(e->csds, &s->sparts[spart_offset], s_counts[c_ind], e,
                      /* log_all_fields */ 1, csds_flag_mpi_enter, i);

      /* Log the gparts */
      csds_log_gparts(e->csds, &s->gparts[gpart_offset], g_counts[c_ind], e,
                      /* log_all_fields */ 1, csds_flag_mpi_enter, i);

      /* Log the bparts */
      if (b_counts[c_ind] > 0) {
        error("TODO");
      }

      /* Update the counters */
      part_offset += counts[c_ind];
      spart_offset += s_counts[c_ind];
      gpart_offset += g_counts[c_ind];
      bpart_offset += b_counts[c_ind];
    }
  }
#endif

  /* Restore the part<->gpart and spart<->gpart links.
   * Generate indices and counts for threadpool tasks. Note we process a node
   * at a time. */
  struct relink_mapper_data relink_data;
  relink_data.s = s;
  relink_data.counts = counts;
  relink_data.g_counts = g_counts;
  relink_data.s_counts = s_counts;
  relink_data.b_counts = b_counts;
  relink_data.nodeID = nodeID;
  relink_data.nr_nodes = nr_nodes;

  threadpool_map(&e->threadpool, engine_redistribute_relink_mapper, nodes,
                 nr_nodes, sizeof(int), 1, &relink_data);
  free(nodes);

  /* Clean up the counts now we are done. */
  free(counts);
  free(g_counts);
  free(s_counts);
  free(b_counts);

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
  for (size_t k = 0; k < nr_bparts_new; k++) {
    const int cid = cell_getid(s->cdim, s->bparts[k].x[0] * s->iwidth[0],
                               s->bparts[k].x[1] * s->iwidth[1],
                               s->bparts[k].x[2] * s->iwidth[2]);
    if (cells[cid].nodeID != nodeID)
      error("Received b-particle (%zu) that does not belong here (nodeID=%i).",
            k, cells[cid].nodeID);
  }

  /* Verify that the links are correct */
  part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                    nr_parts_new, nr_gparts_new, nr_sinks_new, nr_sparts_new,
                    nr_bparts_new, e->verbose);

#endif

  /* Be verbose about what just happened. */
  if (e->verbose) {
    int my_cells = 0;
    for (int k = 0; k < nr_cells; k++)
      if (cells[k].nodeID == nodeID) my_cells += 1;
    message(
        "node %i now has %zu parts, %zu sparts, %zu bparts and %zu gparts in "
        "%i cells.",
        nodeID, nr_parts_new, nr_sparts_new, nr_bparts_new, nr_gparts_new,
        my_cells);
  }

  /* Flag that we do not have any extra particles any more */
  s->nr_extra_parts = 0;
  s->nr_extra_gparts = 0;
  s->nr_extra_sparts = 0;
  s->nr_extra_bparts = 0;

  /* Flag that a redistribute has taken place */
  e->step_props |= engine_step_prop_redistribute;

  ticks tock = getticks();
  e->local_treebuild_time += clocks_diff_ticks(tock, tic);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(tock - tic), clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
