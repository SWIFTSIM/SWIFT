/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CUSTOM_MPI_H
#define SWIFT_CUSTOM_MPI_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <limits.h>
#include <stdlib.h>

/* Local headers. */
#include "inline.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>

/**
 * @brief Creates a datatype describing a contiguous run of more elements than
 * fit in an int.
 *
 * The run is split into full chunks of INT_MAX elements followed by the
 * remainder, if any, glued together in a struct. Counts that fit in an int
 * simply yield a contiguous type.
 *
 * @param total_count Number of elements.
 * @param oldtype Datatype of the elements.
 * @param newtype (return) The new, uncommitted, datatype. MPI_DATATYPE_NULL if
 * total_count is 0.
 */
INLINE static int create_large_count_type(const size_t total_count,
                                          MPI_Datatype oldtype,
                                          MPI_Datatype *newtype) {

  if (total_count == 0) {
    *newtype = MPI_DATATYPE_NULL;
    return MPI_SUCCESS;
  }

  /* If it fits into a standard int, use a normal contiguous block */
  if (total_count <= INT_MAX) {
    return MPI_Type_contiguous((int)total_count, oldtype, newtype);
  }

  /* Split the count into chunks of INT_MAX plus a remainder */
  const size_t chunks = total_count / INT_MAX;
  const int remainder = (int)(total_count % INT_MAX);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(oldtype, &lb, &extent);

  /* At most two blocks: the full chunks and, if there is one, the remainder
   * placed right after them. */
  const int count = (remainder > 0) ? 2 : 1;
  int blocklengths[2] = {1, 1};
  MPI_Aint displacements[2] = {0, (MPI_Aint)chunks * INT_MAX * extent};
  MPI_Datatype types[2] = {MPI_DATATYPE_NULL, MPI_DATATYPE_NULL};

  MPI_Type_vector((int)chunks, INT_MAX, INT_MAX, oldtype, &types[0]);
  if (remainder > 0) MPI_Type_contiguous(remainder, oldtype, &types[1]);

  int status = MPI_Type_create_struct(count, blocklengths, displacements, types,
                                      newtype);

  /* The struct holds its own references to the sub-types */
  MPI_Type_free(&types[0]);
  if (remainder > 0) MPI_Type_free(&types[1]);

  return status;
}

/**
 * @brief Implementation of MPI_Allgatherv which allows for size_t arguments
 *
 * Every rank sends its data to every rank with point-to-point messages. The
 * displacements are applied as plain pointer arithmetic, so derived datatypes
 * are only ever built for the (rare) counts that do not fit in an int.
 *
 * @param sendbuf Starting address of send buffer.
 * @param sendcount Number of elements in send buffer.
 * @param sendtype Data type of send buffer elements.
 * @param recvbuf Address of receive buffer.
 * @param recvcounts size_t array (of length group size) containing the number
 * of elements that are to be received from each process.
 * @param displs size_t array (of length group size). Entry i specifies the
 * displacement (relative to recvbuf ) at which to place the incoming data from
 * process i.
 * @param recvtype Data type of receive buffer elements.
 * @param comm MPI communicator.
 */
INLINE static int swift_mpi_allgatherv_sizet(
    const void *sendbuf, size_t sendcount, MPI_Datatype sendtype, void *recvbuf,
    const size_t recvcounts[], const size_t displs[], MPI_Datatype recvtype,
    MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(recvtype, &lb, &extent);

  /* Per-rank receive types (only used for counts that do not fit in an int)
   * and one request per receive and per send. */
  MPI_Datatype *recv_types =
      (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));
  MPI_Request *requests = (MPI_Request *)malloc(2 * size * sizeof(MPI_Request));

  if (!recv_types || !requests) {
    free(recv_types);
    free(requests);
    return MPI_ERR_INTERN;
  }

  /* Nothing committed yet: the clean-up below relies on these being set. */
  MPI_Datatype send_type = MPI_DATATYPE_NULL;
  for (int i = 0; i < size; ++i) recv_types[i] = MPI_DATATYPE_NULL;

  int req_count = 0;
  int status = MPI_SUCCESS;

  /* Post the non-blocking receives, each straight into its slot */
  for (int i = 0; i < size; ++i) {
    if (recvcounts[i] == 0) continue;

    char *const dest = (char *)recvbuf + (MPI_Aint)displs[i] * extent;

    if (recvcounts[i] <= INT_MAX) {
      MPI_Irecv(dest, (int)recvcounts[i], recvtype, i, 0, comm,
                &requests[req_count++]);
    } else {
      MPI_Datatype layout;
      status = create_large_count_type(recvcounts[i], recvtype, &layout);
      if (status != MPI_SUCCESS) break;
      recv_types[i] = layout;
      MPI_Type_commit(&recv_types[i]);
      MPI_Irecv(dest, 1, recv_types[i], i, 0, comm, &requests[req_count++]);
    }
  }

  /* Post the non-blocking sends of the same buffer to every rank */
  if (status == MPI_SUCCESS && sendcount > 0) {
    if (sendcount <= INT_MAX) {
      for (int i = 0; i < size; ++i)
        MPI_Isend(sendbuf, (int)sendcount, sendtype, i, 0, comm,
                  &requests[req_count++]);
    } else {
      MPI_Datatype layout;
      status = create_large_count_type(sendcount, sendtype, &layout);
      if (status == MPI_SUCCESS) {
        send_type = layout;
        MPI_Type_commit(&send_type);
        for (int i = 0; i < size; ++i)
          MPI_Isend(sendbuf, 1, send_type, i, 0, comm, &requests[req_count++]);
      }
    }
  }

  /* Complete network operations and perform structural deallocations */
  if (req_count > 0) {
    if (status == MPI_SUCCESS) {
      status = MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
    } else {

      /* We bailed out whilst setting up the transfers. Only receives have been
       * posted at this point and they will never be matched, so cancel them
       * and reap them all before the request buffer disappears. */
      for (int i = 0; i < req_count; ++i) MPI_Cancel(&requests[i]);
      MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
    }
  }

  /* Garbage collection of custom committed structures */
  if (send_type != MPI_DATATYPE_NULL) MPI_Type_free(&send_type);
  for (int i = 0; i < size; ++i) {
    if (recv_types[i] != MPI_DATATYPE_NULL) MPI_Type_free(&recv_types[i]);
  }

  free(recv_types);
  free(requests);
  return status;
}

#endif /* WITH_MPI */

#endif /* SWIFT_CUSTOM_MPI_H */
