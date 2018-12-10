/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#endif

/* This object's header. */
#include "proxy.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "space.h"

#ifdef WITH_MPI
/* MPI data type for the communications */
MPI_Datatype pcell_mpi_type;
#endif

/**
 * @brief Exchange tags between nodes.
 *
 * Note that this function assumes that the cell structures have already
 * been exchanged, e.g. via #proxy_cells_exchange.
 *
 * @param proxies The list of #proxy that will send/recv tags
 * @param num_proxies The number of proxies.
 * @param s The space into which the tags will be unpacked.
 */
void proxy_tags_exchange(struct proxy *proxies, int num_proxies,
                         struct space *s) {

#ifdef WITH_MPI

  ticks tic2 = getticks();

  /* Run through the cells and get the size of the tags that will be sent off.
   */
  int count_out = 0;
  int offset_out[s->nr_cells];
  for (int k = 0; k < s->nr_cells; k++) {
    offset_out[k] = count_out;
    if (s->cells_top[k].mpi.sendto) {
      count_out += s->cells_top[k].mpi.pcell_size;
    }
  }

  /* Run through the proxies and get the count of incoming tags. */
  int count_in = 0;
  int offset_in[s->nr_cells];
  for (int k = 0; k < num_proxies; k++) {
    for (int j = 0; j < proxies[k].nr_cells_in; j++) {
      offset_in[proxies[k].cells_in[j] - s->cells_top] = count_in;
      count_in += proxies[k].cells_in[j]->mpi.pcell_size;
    }
  }

  /* Allocate the tags. */
  int *tags_in = NULL;
  int *tags_out = NULL;
  if (posix_memalign((void **)&tags_in, SWIFT_CACHE_ALIGNMENT,
                     sizeof(int) * count_in) != 0 ||
      posix_memalign((void **)&tags_out, SWIFT_CACHE_ALIGNMENT,
                     sizeof(int) * count_out) != 0)
    error("Failed to allocate tags buffers.");

  /* Pack the local tags. */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].mpi.sendto) {
      cell_pack_tags(&s->cells_top[k], &tags_out[offset_out[k]]);
    }
  }

  if (s->e->verbose)
    message("Cell pack tags took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Allocate the incoming and outgoing request handles. */
  int num_reqs_out = 0;
  int num_reqs_in = 0;
  for (int k = 0; k < num_proxies; k++) {
    num_reqs_in += proxies[k].nr_cells_in;
    num_reqs_out += proxies[k].nr_cells_out;
  }
  MPI_Request *reqs_in = NULL;
  int *cids_in = NULL;
  if ((reqs_in = (MPI_Request *)malloc(sizeof(MPI_Request) *
                                       (num_reqs_in + num_reqs_out))) == NULL ||
      (cids_in = (int *)malloc(sizeof(int) * (num_reqs_in + num_reqs_out))) ==
          NULL)
    error("Failed to allocate MPI_Request arrays.");
  MPI_Request *reqs_out = &reqs_in[num_reqs_in];
  int *cids_out = &cids_in[num_reqs_in];

  /* Emit the sends and recvs. */
  for (int send_rid = 0, recv_rid = 0, k = 0; k < num_proxies; k++) {
    for (int j = 0; j < proxies[k].nr_cells_in; j++) {
      const int cid = proxies[k].cells_in[j] - s->cells_top;
      cids_in[recv_rid] = cid;
      int err = MPI_Irecv(
          &tags_in[offset_in[cid]], proxies[k].cells_in[j]->mpi.pcell_size,
          MPI_INT, proxies[k].nodeID, cid, MPI_COMM_WORLD, &reqs_in[recv_rid]);
      if (err != MPI_SUCCESS) mpi_error(err, "Failed to irecv tags.");
      recv_rid += 1;
    }
    for (int j = 0; j < proxies[k].nr_cells_out; j++) {
      const int cid = proxies[k].cells_out[j] - s->cells_top;
      cids_out[send_rid] = cid;
      int err = MPI_Isend(
          &tags_out[offset_out[cid]], proxies[k].cells_out[j]->mpi.pcell_size,
          MPI_INT, proxies[k].nodeID, cid, MPI_COMM_WORLD, &reqs_out[send_rid]);
      if (err != MPI_SUCCESS) mpi_error(err, "Failed to isend tags.");
      send_rid += 1;
    }
  }

  tic2 = getticks();

  /* Wait for each recv and unpack the tags into the local cells. */
  for (int k = 0; k < num_reqs_in; k++) {
    int pid = MPI_UNDEFINED;
    MPI_Status status;
    if (MPI_Waitany(num_reqs_in, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    const int cid = cids_in[pid];
    cell_unpack_tags(&tags_in[offset_in[cid]], &s->cells_top[cid]);
  }

  if (s->e->verbose)
    message("Cell unpack tags took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Wait for all the sends to have completed. */
  if (MPI_Waitall(num_reqs_out, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Clean up. */
  free(tags_in);
  free(tags_out);
  free(reqs_in);
  free(cids_in);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchange cells with a remote node, first part.
 *
 * The first part of the transaction sends the local cell count and the packed
 * #pcell array to the destination node, and enqueues an @c MPI_Irecv for
 * the foreign cell counts.
 *
 * @param p The #proxy.
 */
void proxy_cells_exchange_first(struct proxy *p) {

#ifdef WITH_MPI

  /* Get the number of pcells we will need to send. */
  p->size_pcells_out = 0;
  for (int k = 0; k < p->nr_cells_out; k++)
    p->size_pcells_out += p->cells_out[k]->mpi.pcell_size;

  /* Send the number of pcells. */
  int err = MPI_Isend(&p->size_pcells_out, 1, MPI_INT, p->nodeID,
                      p->mynodeID * proxy_tag_shift + proxy_tag_count,
                      MPI_COMM_WORLD, &p->req_cells_count_out);
  if (err != MPI_SUCCESS) mpi_error(err, "Failed to isend nr of pcells.");
  // message( "isent pcell count (%i) from node %i to node %i." ,
  // p->size_pcells_out , p->mynodeID , p->nodeID ); fflush(stdout);

  /* Allocate and fill the pcell buffer. */
  if (p->pcells_out != NULL) free(p->pcells_out);
  if (posix_memalign((void **)&p->pcells_out, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct pcell) * p->size_pcells_out) != 0)
    error("Failed to allocate pcell_out buffer.");
  for (int ind = 0, k = 0; k < p->nr_cells_out; k++) {
    memcpy(&p->pcells_out[ind], p->cells_out[k]->mpi.pcell,
           sizeof(struct pcell) * p->cells_out[k]->mpi.pcell_size);
    ind += p->cells_out[k]->mpi.pcell_size;
  }

  /* Send the pcell buffer. */
  err = MPI_Isend(p->pcells_out, p->size_pcells_out, pcell_mpi_type, p->nodeID,
                  p->mynodeID * proxy_tag_shift + proxy_tag_cells,
                  MPI_COMM_WORLD, &p->req_cells_out);

  if (err != MPI_SUCCESS) mpi_error(err, "Failed to pcell_out buffer.");
  // message( "isent pcells (%i) from node %i to node %i." , p->size_pcells_out
  // , p->mynodeID , p->nodeID ); fflush(stdout);

  /* Receive the number of pcells. */
  err = MPI_Irecv(&p->size_pcells_in, 1, MPI_INT, p->nodeID,
                  p->nodeID * proxy_tag_shift + proxy_tag_count, MPI_COMM_WORLD,
                  &p->req_cells_count_in);
  if (err != MPI_SUCCESS) mpi_error(err, "Failed to irecv nr of pcells.");
    // message( "irecv pcells count on node %i from node %i." , p->mynodeID ,
    // p->nodeID ); fflush(stdout);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchange cells with a remote node, second part.
 *
 * Once the incomming cell count has been received, allocate a buffer
 * for the foreign packed #pcell array and emit the @c MPI_Irecv for
 * it.
 *
 * @param p The #proxy.
 */
void proxy_cells_exchange_second(struct proxy *p) {

#ifdef WITH_MPI

  /* Re-allocate the pcell_in buffer. */
  if (p->pcells_in != NULL) free(p->pcells_in);
  if (posix_memalign((void **)&p->pcells_in, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct pcell) * p->size_pcells_in) != 0)
    error("Failed to allocate pcell_in buffer.");

  /* Receive the particle buffers. */
  int err = MPI_Irecv(p->pcells_in, p->size_pcells_in, pcell_mpi_type,
                      p->nodeID, p->nodeID * proxy_tag_shift + proxy_tag_cells,
                      MPI_COMM_WORLD, &p->req_cells_in);

  if (err != MPI_SUCCESS) mpi_error(err, "Failed to irecv part data.");
    // message( "irecv pcells (%i) on node %i from node %i." , p->size_pcells_in
    // , p->mynodeID , p->nodeID ); fflush(stdout);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

#ifdef WITH_MPI

void proxy_cells_count_mapper(void *map_data, int num_elements,
                              void *extra_data) {
  struct cell *cells = (struct cell *)map_data;

  for (int k = 0; k < num_elements; k++) {
    if (cells[k].mpi.sendto) cells[k].mpi.pcell_size = cell_getsize(&cells[k]);
  }
}

struct pack_mapper_data {
  struct space *s;
  int *offset;
  struct pcell *pcells;
  int with_gravity;
};

void proxy_cells_pack_mapper(void *map_data, int num_elements,
                             void *extra_data) {
  struct cell *cells = (struct cell *)map_data;
  struct pack_mapper_data *data = (struct pack_mapper_data *)extra_data;

  for (int k = 0; k < num_elements; k++) {
    if (cells[k].mpi.sendto) {
      ptrdiff_t ind = &cells[k] - data->s->cells_top;
      cells[k].mpi.pcell = &data->pcells[data->offset[ind]];
      cell_pack(&cells[k], cells[k].mpi.pcell, data->with_gravity);
    }
  }
}

void proxy_cells_exchange_first_mapper(void *map_data, int num_elements,
                                       void *extra_data) {
  struct proxy *proxies = (struct proxy *)map_data;

  for (int k = 0; k < num_elements; k++) {
    proxy_cells_exchange_first(&proxies[k]);
  }
}

struct wait_and_unpack_mapper_data {
  struct space *s;
  int num_proxies;
  MPI_Request *reqs_in;
  struct proxy *proxies;
  int with_gravity;
  swift_lock_type lock;
};

void proxy_cells_wait_and_unpack_mapper(void *unused_map_data, int num_elements,
                                        void *extra_data) {

  // MATTHIEU: This is currently unused. Scalar (non-threadpool) version is
  // faster but we still need to explore why this happens.

  struct wait_and_unpack_mapper_data *data =
      (struct wait_and_unpack_mapper_data *)extra_data;

  for (int k = 0; k < num_elements; k++) {
    int pid = MPI_UNDEFINED;
    MPI_Status status;
    int res;

    /* We need a lock to prevent concurrent calls to MPI_Waitany on
       the same array of requests since this is not supported in the MPI
       standard (v3.1). This is not really a problem since the threads
       would block inside MPI_Waitany anyway. */
    lock_lock(&data->lock);
    if ((res = MPI_Waitany(data->num_proxies, data->reqs_in, &pid, &status)) !=
            MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      mpi_error(res, "MPI_Waitany failed.");
    if (lock_unlock(&data->lock) != 0) {
      error("Failed to release lock.");
    }

    // message( "cell data from proxy %i has arrived." , pid );
    for (int count = 0, j = 0; j < data->proxies[pid].nr_cells_in; j++)
      count += cell_unpack(&data->proxies[pid].pcells_in[count],
                           data->proxies[pid].cells_in[j], data->s,
                           data->with_gravity);
  }
}

#endif  // WITH_MPI

/**
 * @brief Exchange the cell structures with all proxies.
 *
 * @param proxies The list of #proxy that will send/recv cells.
 * @param num_proxies The number of proxies.
 * @param s The space into which the particles will be unpacked.
 * @param with_gravity Are we running with gravity and hence need
 *      to exchange multipoles?
 */
void proxy_cells_exchange(struct proxy *proxies, int num_proxies,
                          struct space *s, const int with_gravity) {

#ifdef WITH_MPI

  MPI_Request *reqs;
  if ((reqs = (MPI_Request *)malloc(sizeof(MPI_Request) * 2 * num_proxies)) ==
      NULL)
    error("Failed to allocate request buffers.");
  MPI_Request *reqs_in = reqs;
  MPI_Request *reqs_out = &reqs[num_proxies];

  ticks tic2 = getticks();

  /* Run through the cells and get the size of the ones that will be sent off.
   */
  threadpool_map(&s->e->threadpool, proxy_cells_count_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), /*chunk=*/0,
                 /*extra_data=*/NULL);
  int count_out = 0;
  int offset[s->nr_cells];
  for (int k = 0; k < s->nr_cells; k++) {
    offset[k] = count_out;
    if (s->cells_top[k].mpi.sendto) count_out += s->cells_top[k].mpi.pcell_size;
  }

  if (s->e->verbose)
    message("Counting cells to send took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Allocate the pcells. */
  struct pcell *pcells = NULL;
  if (posix_memalign((void **)&pcells, SWIFT_CACHE_ALIGNMENT,
                     sizeof(struct pcell) * count_out) != 0)
    error("Failed to allocate pcell buffer.");

  tic2 = getticks();

  /* Pack the cells. */
  struct pack_mapper_data data = {s, offset, pcells, with_gravity};
  threadpool_map(&s->e->threadpool, proxy_cells_pack_mapper, s->cells_top,
                 s->nr_cells, sizeof(struct cell), /*chunk=*/0, &data);

  if (s->e->verbose)
    message("Packing cells took %.3f %s.", clocks_from_ticks(getticks() - tic2),
            clocks_getunit());

  /* Launch the first part of the exchange. */
  threadpool_map(&s->e->threadpool, proxy_cells_exchange_first_mapper, proxies,
                 num_proxies, sizeof(struct proxy), /*chunk=*/0,
                 /*extra_data=*/NULL);
  for (int k = 0; k < num_proxies; k++) {
    reqs_in[k] = proxies[k].req_cells_count_in;
    reqs_out[k] = proxies[k].req_cells_count_out;
  }

  /* Wait for each count to come in and start the recv. */
  for (int k = 0; k < num_proxies; k++) {
    int pid = MPI_UNDEFINED;
    MPI_Status status;
    if (MPI_Waitany(num_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "request from proxy %i has arrived." , pid );
    proxy_cells_exchange_second(&proxies[pid]);
  }

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(num_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Set the requests for the cells. */
  for (int k = 0; k < num_proxies; k++) {
    reqs_in[k] = proxies[k].req_cells_in;
    reqs_out[k] = proxies[k].req_cells_out;
  }

  tic2 = getticks();

  /* Wait for each pcell array to come in from the proxies. */
  for (int k = 0; k < num_proxies; k++) {
    int pid = MPI_UNDEFINED;
    MPI_Status status;
    if (MPI_Waitany(num_proxies, reqs_in, &pid, &status) != MPI_SUCCESS ||
        pid == MPI_UNDEFINED)
      error("MPI_Waitany failed.");
    // message( "cell data from proxy %i has arrived." , pid );
    for (int count = 0, j = 0; j < proxies[pid].nr_cells_in; j++)
      count += cell_unpack(&proxies[pid].pcells_in[count],
                           proxies[pid].cells_in[j], s, with_gravity);
  }

  if (s->e->verbose)
    message("Un-packing cells took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

  /* Wait for all the sends to have finished too. */
  if (MPI_Waitall(num_proxies, reqs_out, MPI_STATUSES_IGNORE) != MPI_SUCCESS)
    error("MPI_Waitall on sends failed.");

  /* Clean up. */
  free(reqs);
  free(pcells);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Add a cell to the given proxy's input list.
 *
 * @param p The #proxy.
 * @param c The #cell.
 * @param type Why is this cell in the proxy (hdro, gravity, ...) ?
 */
void proxy_addcell_in(struct proxy *p, struct cell *c, int type) {

  if (type == proxy_cell_type_none) error("Invalid type for proxy");

  /* Check if the cell is already registered with the proxy. */
  for (int k = 0; k < p->nr_cells_in; k++)
    if (p->cells_in[k] == c) {

      /* Update the type */
      p->cells_in_type[k] |= type;
      return;
    }

  /* Do we need to grow the number of in cells? */
  if (p->nr_cells_in == p->size_cells_in) {

    p->size_cells_in *= proxy_buffgrow;

    struct cell **temp_cell;
    if ((temp_cell = (struct cell **)malloc(sizeof(struct cell *) *
                                            p->size_cells_in)) == NULL)
      error("Failed to allocate incoming cell list.");
    memcpy(temp_cell, p->cells_in, sizeof(struct cell *) * p->nr_cells_in);
    free(p->cells_in);
    p->cells_in = temp_cell;

    int *temp_type;
    if ((temp_type = (int *)malloc(sizeof(int) * p->size_cells_in)) == NULL)
      error("Failed to allocate incoming cell type list.");
    memcpy(temp_type, p->cells_in_type, sizeof(int) * p->nr_cells_in);
    free(p->cells_in_type);
    p->cells_in_type = temp_type;
  }

  /* Add the cell. */
  p->cells_in[p->nr_cells_in] = c;
  p->cells_in_type[p->nr_cells_in] = type;
  p->nr_cells_in += 1;
}

/**
 * @brief Add a cell to the given proxy's output list.
 *
 * @param p The #proxy.
 * @param c The #cell.
 * @param type Why is this cell in the proxy (hdro, gravity, ...) ?
 */
void proxy_addcell_out(struct proxy *p, struct cell *c, int type) {

  if (type == proxy_cell_type_none) error("Invalid type for proxy");

  /* Check if the cell is already registered with the proxy. */
  for (int k = 0; k < p->nr_cells_out; k++)
    if (p->cells_out[k] == c) {

      /* Update the type */
      p->cells_out_type[k] |= type;
      return;
    }

  /* Do we need to grow the number of out cells? */
  if (p->nr_cells_out == p->size_cells_out) {
    p->size_cells_out *= proxy_buffgrow;

    struct cell **temp_cell;
    if ((temp_cell = (struct cell **)malloc(sizeof(struct cell *) *
                                            p->size_cells_out)) == NULL)
      error("Failed to allocate outgoing cell list.");
    memcpy(temp_cell, p->cells_out, sizeof(struct cell *) * p->nr_cells_out);
    free(p->cells_out);
    p->cells_out = temp_cell;

    int *temp_type;
    if ((temp_type = (int *)malloc(sizeof(int) * p->size_cells_out)) == NULL)
      error("Failed to allocate outgoing cell type list.");
    memcpy(temp_type, p->cells_out_type, sizeof(int) * p->nr_cells_out);
    free(p->cells_out_type);
    p->cells_out_type = temp_type;
  }

  /* Add the cell. */
  p->cells_out[p->nr_cells_out] = c;
  p->cells_out_type[p->nr_cells_out] = type;
  p->nr_cells_out += 1;
}

/**
 * @brief Exchange particles with a remote node.
 *
 * @param p The #proxy.
 */
void proxy_parts_exchange_first(struct proxy *p) {

#ifdef WITH_MPI

  /* Send the number of particles. */
  p->buff_out[0] = p->nr_parts_out;
  p->buff_out[1] = p->nr_gparts_out;
  p->buff_out[2] = p->nr_sparts_out;
  if (MPI_Isend(p->buff_out, 3, MPI_INT, p->nodeID,
                p->mynodeID * proxy_tag_shift + proxy_tag_count, MPI_COMM_WORLD,
                &p->req_parts_count_out) != MPI_SUCCESS)
    error("Failed to isend nr of parts.");
  /* message( "isent particle counts [%i, %i] from node %i to node %i." ,
  p->buff_out[0], p->buff_out[1], p->mynodeID , p->nodeID ); fflush(stdout); */

  /* Send the particle buffers. */
  if (p->nr_parts_out > 0) {
    if (MPI_Isend(p->parts_out, p->nr_parts_out, part_mpi_type, p->nodeID,
                  p->mynodeID * proxy_tag_shift + proxy_tag_parts,
                  MPI_COMM_WORLD, &p->req_parts_out) != MPI_SUCCESS ||
        MPI_Isend(p->xparts_out, p->nr_parts_out, xpart_mpi_type, p->nodeID,
                  p->mynodeID * proxy_tag_shift + proxy_tag_xparts,
                  MPI_COMM_WORLD, &p->req_xparts_out) != MPI_SUCCESS)
      error("Failed to isend part data.");
    // message( "isent particle data (%i) to node %i." , p->nr_parts_out ,
    // p->nodeID ); fflush(stdout);
    /*for (int k = 0; k < p->nr_parts_out; k++)
      message("sending particle %lli, x=[%.3e %.3e %.3e], h=%.3e, to node %i.",
              p->parts_out[k].id, p->parts_out[k].x[0], p->parts_out[k].x[1],
              p->parts_out[k].x[2], p->parts_out[k].h, p->nodeID);*/
  }
  if (p->nr_gparts_out > 0) {
    if (MPI_Isend(p->gparts_out, p->nr_gparts_out, gpart_mpi_type, p->nodeID,
                  p->mynodeID * proxy_tag_shift + proxy_tag_gparts,
                  MPI_COMM_WORLD, &p->req_gparts_out) != MPI_SUCCESS)
      error("Failed to isend gpart data.");
    // message( "isent gpart data (%i) to node %i." , p->nr_parts_out ,
    // p->nodeID ); fflush(stdout);
  }

  if (p->nr_sparts_out > 0) {
    if (MPI_Isend(p->sparts_out, p->nr_sparts_out, spart_mpi_type, p->nodeID,
                  p->mynodeID * proxy_tag_shift + proxy_tag_sparts,
                  MPI_COMM_WORLD, &p->req_sparts_out) != MPI_SUCCESS)
      error("Failed to isend spart data.");
    // message( "isent gpart data (%i) to node %i." , p->nr_parts_out ,
    // p->nodeID ); fflush(stdout);
  }

  /* Receive the number of particles. */
  if (MPI_Irecv(p->buff_in, 3, MPI_INT, p->nodeID,
                p->nodeID * proxy_tag_shift + proxy_tag_count, MPI_COMM_WORLD,
                &p->req_parts_count_in) != MPI_SUCCESS)
    error("Failed to irecv nr of parts.");

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

void proxy_parts_exchange_second(struct proxy *p) {

#ifdef WITH_MPI

  /* Unpack the incomming parts counts. */
  p->nr_parts_in = p->buff_in[0];
  p->nr_gparts_in = p->buff_in[1];
  p->nr_sparts_in = p->buff_in[2];

  /* Is there enough space in the buffers? */
  if (p->nr_parts_in > p->size_parts_in) {
    do {
      p->size_parts_in *= proxy_buffgrow;
    } while (p->nr_parts_in > p->size_parts_in);
    free(p->parts_in);
    free(p->xparts_in);
    if ((p->parts_in = (struct part *)malloc(sizeof(struct part) *
                                             p->size_parts_in)) == NULL ||
        (p->xparts_in = (struct xpart *)malloc(sizeof(struct xpart) *
                                               p->size_parts_in)) == NULL)
      error("Failed to re-allocate parts_in buffers.");
  }
  if (p->nr_gparts_in > p->size_gparts_in) {
    do {
      p->size_gparts_in *= proxy_buffgrow;
    } while (p->nr_gparts_in > p->size_gparts_in);
    free(p->gparts_in);
    if ((p->gparts_in = (struct gpart *)malloc(sizeof(struct gpart) *
                                               p->size_gparts_in)) == NULL)
      error("Failed to re-allocate gparts_in buffers.");
  }
  if (p->nr_sparts_in > p->size_sparts_in) {
    do {
      p->size_sparts_in *= proxy_buffgrow;
    } while (p->nr_sparts_in > p->size_sparts_in);
    free(p->sparts_in);
    if ((p->sparts_in = (struct spart *)malloc(sizeof(struct spart) *
                                               p->size_sparts_in)) == NULL)
      error("Failed to re-allocate sparts_in buffers.");
  }

  /* Receive the particle buffers. */
  if (p->nr_parts_in > 0) {
    if (MPI_Irecv(p->parts_in, p->nr_parts_in, part_mpi_type, p->nodeID,
                  p->nodeID * proxy_tag_shift + proxy_tag_parts, MPI_COMM_WORLD,
                  &p->req_parts_in) != MPI_SUCCESS ||
        MPI_Irecv(p->xparts_in, p->nr_parts_in, xpart_mpi_type, p->nodeID,
                  p->nodeID * proxy_tag_shift + proxy_tag_xparts,
                  MPI_COMM_WORLD, &p->req_xparts_in) != MPI_SUCCESS)
      error("Failed to irecv part data.");
    // message( "irecv particle data (%i) from node %i." , p->nr_parts_in ,
    // p->nodeID ); fflush(stdout);
  }
  if (p->nr_gparts_in > 0) {
    if (MPI_Irecv(p->gparts_in, p->nr_gparts_in, gpart_mpi_type, p->nodeID,
                  p->nodeID * proxy_tag_shift + proxy_tag_gparts,
                  MPI_COMM_WORLD, &p->req_gparts_in) != MPI_SUCCESS)
      error("Failed to irecv gpart data.");
    // message( "irecv gpart data (%i) from node %i." , p->nr_gparts_in ,
    // p->nodeID ); fflush(stdout);
  }
  if (p->nr_sparts_in > 0) {
    if (MPI_Irecv(p->sparts_in, p->nr_sparts_in, spart_mpi_type, p->nodeID,
                  p->nodeID * proxy_tag_shift + proxy_tag_sparts,
                  MPI_COMM_WORLD, &p->req_sparts_in) != MPI_SUCCESS)
      error("Failed to irecv spart data.");
    // message( "irecv gpart data (%i) from node %i." , p->nr_gparts_in ,
    // p->nodeID ); fflush(stdout);
  }

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Load parts onto a proxy for exchange.
 *
 * @param p The #proxy.
 * @param parts Pointer to an array of #part to send.
 * @param xparts Pointer to an array of #xpart to send.
 * @param N The number of parts.
 */
void proxy_parts_load(struct proxy *p, const struct part *parts,
                      const struct xpart *xparts, int N) {

  /* Is there enough space in the buffer? */
  if (p->nr_parts_out + N > p->size_parts_out) {
    do {
      p->size_parts_out *= proxy_buffgrow;
    } while (p->nr_parts_out + N > p->size_parts_out);
    struct part *tp = NULL;
    struct xpart *txp = NULL;
    if ((tp = (struct part *)malloc(sizeof(struct part) * p->size_parts_out)) ==
            NULL ||
        (txp = (struct xpart *)malloc(sizeof(struct xpart) *
                                      p->size_parts_out)) == NULL)
      error("Failed to re-allocate parts_out buffers.");
    memcpy(tp, p->parts_out, sizeof(struct part) * p->nr_parts_out);
    memcpy(txp, p->xparts_out, sizeof(struct xpart) * p->nr_parts_out);
    free(p->parts_out);
    free(p->xparts_out);
    p->parts_out = tp;
    p->xparts_out = txp;
  }

  /* Copy the parts and xparts data to the buffer. */
  memcpy(&p->parts_out[p->nr_parts_out], parts, sizeof(struct part) * N);
  memcpy(&p->xparts_out[p->nr_parts_out], xparts, sizeof(struct xpart) * N);

  /* Increase the counters. */
  p->nr_parts_out += N;
}

/**
 * @brief Load gparts onto a proxy for exchange.
 *
 * @param p The #proxy.
 * @param gparts Pointer to an array of #gpart to send.
 * @param N The number of gparts.
 */
void proxy_gparts_load(struct proxy *p, const struct gpart *gparts, int N) {

  /* Is there enough space in the buffer? */
  if (p->nr_gparts_out + N > p->size_gparts_out) {
    do {
      p->size_gparts_out *= proxy_buffgrow;
    } while (p->nr_gparts_out + N > p->size_gparts_out);
    struct gpart *tp;
    if ((tp = (struct gpart *)malloc(sizeof(struct gpart) *
                                     p->size_gparts_out)) == NULL)
      error("Failed to re-allocate gparts_out buffers.");
    memcpy(tp, p->gparts_out, sizeof(struct gpart) * p->nr_gparts_out);
    free(p->gparts_out);
    p->gparts_out = tp;
  }

  /* Copy the parts and xparts data to the buffer. */
  memcpy(&p->gparts_out[p->nr_gparts_out], gparts, sizeof(struct gpart) * N);

  /* Increase the counters. */
  p->nr_gparts_out += N;
}

/**
 * @brief Load sparts onto a proxy for exchange.
 *
 * @param p The #proxy.
 * @param sparts Pointer to an array of #spart to send.
 * @param N The number of sparts.
 */
void proxy_sparts_load(struct proxy *p, const struct spart *sparts, int N) {

  /* Is there enough space in the buffer? */
  if (p->nr_sparts_out + N > p->size_sparts_out) {
    do {
      p->size_sparts_out *= proxy_buffgrow;
    } while (p->nr_sparts_out + N > p->size_sparts_out);
    struct spart *tp;
    if ((tp = (struct spart *)malloc(sizeof(struct spart) *
                                     p->size_sparts_out)) == NULL)
      error("Failed to re-allocate sparts_out buffers.");
    memcpy(tp, p->sparts_out, sizeof(struct spart) * p->nr_sparts_out);
    free(p->sparts_out);
    p->sparts_out = tp;
  }

  /* Copy the parts and xparts data to the buffer. */
  memcpy(&p->sparts_out[p->nr_sparts_out], sparts, sizeof(struct spart) * N);

  /* Increase the counters. */
  p->nr_sparts_out += N;
}

/**
 * @brief Initialize the given proxy.
 *
 * @param p The #proxy.
 * @param mynodeID The node this proxy is running on.
 * @param nodeID The node with which this proxy will communicate.
 */
void proxy_init(struct proxy *p, int mynodeID, int nodeID) {

  /* Set the nodeID. */
  p->mynodeID = mynodeID;
  p->nodeID = nodeID;

  /* Allocate the cell send and receive buffers, if needed. */
  if (p->cells_in == NULL) {
    p->size_cells_in = proxy_buffinit;
    if ((p->cells_in =
             (struct cell **)malloc(sizeof(void *) * p->size_cells_in)) == NULL)
      error("Failed to allocate cells_in buffer.");
    if ((p->cells_in_type = (int *)malloc(sizeof(int) * p->size_cells_in)) ==
        NULL)
      error("Failed to allocate cells_in_type buffer.");
  }
  p->nr_cells_in = 0;
  if (p->cells_out == NULL) {
    p->size_cells_out = proxy_buffinit;
    if ((p->cells_out = (struct cell **)malloc(sizeof(void *) *
                                               p->size_cells_out)) == NULL)
      error("Failed to allocate cells_out buffer.");
    if ((p->cells_out_type = (int *)malloc(sizeof(int) * p->size_cells_out)) ==
        NULL)
      error("Failed to allocate cells_out_type buffer.");
  }
  p->nr_cells_out = 0;

  /* Allocate the part send and receive buffers, if needed. */
  if (p->parts_in == NULL) {
    p->size_parts_in = proxy_buffinit;
    if ((p->parts_in = (struct part *)malloc(sizeof(struct part) *
                                             p->size_parts_in)) == NULL ||
        (p->xparts_in = (struct xpart *)malloc(sizeof(struct xpart) *
                                               p->size_parts_in)) == NULL)
      error("Failed to allocate parts_in buffers.");
  }
  p->nr_parts_in = 0;
  if (p->parts_out == NULL) {
    p->size_parts_out = proxy_buffinit;
    if ((p->parts_out = (struct part *)malloc(sizeof(struct part) *
                                              p->size_parts_out)) == NULL ||
        (p->xparts_out = (struct xpart *)malloc(sizeof(struct xpart) *
                                                p->size_parts_out)) == NULL)
      error("Failed to allocate parts_out buffers.");
  }
  p->nr_parts_out = 0;

  /* Allocate the gpart send and receive buffers, if needed. */
  if (p->gparts_in == NULL) {
    p->size_gparts_in = proxy_buffinit;
    if ((p->gparts_in = (struct gpart *)malloc(sizeof(struct gpart) *
                                               p->size_gparts_in)) == NULL)
      error("Failed to allocate gparts_in buffers.");
  }
  p->nr_gparts_in = 0;
  if (p->gparts_out == NULL) {
    p->size_gparts_out = proxy_buffinit;
    if ((p->gparts_out = (struct gpart *)malloc(sizeof(struct gpart) *
                                                p->size_gparts_out)) == NULL)
      error("Failed to allocate gparts_out buffers.");
  }
  p->nr_gparts_out = 0;

  /* Allocate the spart send and receive buffers, if needed. */
  if (p->sparts_in == NULL) {
    p->size_sparts_in = proxy_buffinit;
    if ((p->sparts_in = (struct spart *)malloc(sizeof(struct spart) *
                                               p->size_sparts_in)) == NULL)
      error("Failed to allocate sparts_in buffers.");
  }
  p->nr_sparts_in = 0;
  if (p->sparts_out == NULL) {
    p->size_sparts_out = proxy_buffinit;
    if ((p->sparts_out = (struct spart *)malloc(sizeof(struct spart) *
                                                p->size_sparts_out)) == NULL)
      error("Failed to allocate sparts_out buffers.");
  }
  p->nr_sparts_out = 0;
}

/**
 * @brief Registers the MPI types for the proxy cells.
 */
void proxy_create_mpi_type(void) {

#ifdef WITH_MPI
  if (MPI_Type_contiguous(sizeof(struct pcell) / sizeof(unsigned char),
                          MPI_BYTE, &pcell_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&pcell_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for parts.");
  }
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
