/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

/* Standard headers */
#include <stdlib.h>
#include <limits.h>

/* Local headers */
#include "error.h"

/**
 * @brief Given an array of structs of size element_size, send 
 * nr_send[i] elements to each node i. Allocates the receive
 * buffer recvbuf to the appropriate size and returns its size
 * in nr_recv_tot.
 *
 * @param nr_send Number of elements to send to each node
 * @param nr_recv Number of elements to receive from each node
 * @param sendbuf The elements to send
 * @param recvbuf The output buffer
 *
 */
void exchange_structs(size_t *nr_send, void *sendbuf,
                      size_t *nr_recv, void *recvbuf,
                      size_t element_size) {

#ifdef WITH_MPI

  /* Determine rank, number of ranks */
  int nr_nodes, nodeID;
  MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &nodeID);

  /* Compute send offsets */
  size_t *send_offset = malloc(nr_nodes * sizeof(size_t));
  send_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    send_offset[i] = send_offset[i - 1] + nr_send[i - 1];
  }

  /* Compute receive offsets */
  size_t *recv_offset = malloc(nr_nodes * sizeof(size_t));
  recv_offset[0] = 0;
  for (int i = 1; i < nr_nodes; i += 1) {
    recv_offset[i] = recv_offset[i - 1] + nr_recv[i - 1];
  }

  /* Allocate request objects (one send and receive per node) */
  MPI_Request *request = malloc(2 * sizeof(MPI_Request) * nr_nodes);

  /* Make type to communicate the struct */
  MPI_Datatype value_mpi_type;
  if (MPI_Type_contiguous(element_size, MPI_BYTE, &value_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&value_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for struct to exchange.");
  }

  /*
   * Post the send operations. This is an alltoallv really but
   * we want to avoid the limits imposed by int counts and offsets
   * in MPI_Alltoallv.
   */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_send[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_send[i] > INT_MAX) error("exchange_structs() fails if nr_send > INT_MAX!");

      char *buf = (char *) sendbuf;
      MPI_Isend(&(buf[send_offset[i]*element_size]), (int)nr_send[i],
                value_mpi_type, i, 0, MPI_COMM_WORLD, &(request[i]));
    } else {
      request[i] = MPI_REQUEST_NULL;
    }
  }

  /* Post the receives */
  for (int i = 0; i < nr_nodes; i += 1) {
    if (nr_recv[i] > 0) {

      /* TODO: handle very large messages */
      if(nr_recv[i] > INT_MAX) error("exchange_structs() fails if nr_recv > INT_MAX!");

      char *buf = (char *) recvbuf;
      MPI_Irecv(&(buf[recv_offset[i]*element_size]), (int)nr_recv[i],
                value_mpi_type, i, 0, MPI_COMM_WORLD,
                &(request[i + nr_nodes]));
    } else {
      request[i + nr_nodes] = MPI_REQUEST_NULL;
    }
  }

  /* Wait for everything to complete */
  MPI_Waitall(2 * nr_nodes, request, MPI_STATUSES_IGNORE);

  /* Done with the MPI type */
  MPI_Type_free(&value_mpi_type);

  /* Tidy up */
  free(recv_offset);
  free(send_offset);
  free(request);
#else
  error("should only be called in MPI mode");
#endif
}
