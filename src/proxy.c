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
#include "error.h"

/**
 * @brief Exchange cells with a remote node.
 *
 * @param p The #proxy.
 */
void proxy_cells_exch1(struct proxy *p) {

#ifdef WITH_MPI

  /* Get the number of pcells we will need to send. */
  p->size_pcells_out = 0;
  for (int k = 0; k < p->nr_cells_out; k++)
    p->size_pcells_out += p->cells_out[k]->pcell_size;

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
    memcpy(&p->pcells_out[ind], p->cells_out[k]->pcell,
           sizeof(struct pcell) * p->cells_out[k]->pcell_size);
    ind += p->cells_out[k]->pcell_size;
  }

  /* Send the pcell buffer. */
  err = MPI_Isend(p->pcells_out, sizeof(struct pcell) * p->size_pcells_out,
                  MPI_BYTE, p->nodeID,
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

void proxy_cells_exch2(struct proxy *p) {

#ifdef WITH_MPI

  /* Re-allocate the pcell_in buffer. */
  if (p->pcells_in != NULL) free(p->pcells_in);
  if (posix_memalign((void **)&p->pcells_in, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct pcell) * p->size_pcells_in) != 0)
    error("Failed to allocate pcell_in buffer.");

  /* Receive the particle buffers. */
  int err = MPI_Irecv(p->pcells_in, sizeof(struct pcell) * p->size_pcells_in,
                      MPI_BYTE, p->nodeID,
                      p->nodeID * proxy_tag_shift + proxy_tag_cells,
                      MPI_COMM_WORLD, &p->req_cells_in);

  if (err != MPI_SUCCESS) mpi_error(err, "Failed to irecv part data.");
    // message( "irecv pcells (%i) on node %i from node %i." , p->size_pcells_in
    // , p->mynodeID , p->nodeID ); fflush(stdout);

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
void proxy_parts_exch1(struct proxy *p) {

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

void proxy_parts_exch2(struct proxy *p) {

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
