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
#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>

/**
 * @brief Create MPI data types for > 2^31 bytes transfers.
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

  /* Build the structural definitions for the large layout */
  const int count = (remainder > 0) ? 2 : 1;
  int *blocklengths = (int *)malloc(count * sizeof(int));
  MPI_Aint *displacements = (MPI_Aint *)malloc(count * sizeof(MPI_Aint));
  MPI_Datatype *types = (MPI_Datatype *)malloc(count * sizeof(MPI_Datatype));

  if (!blocklengths || !displacements || !types) {
    free(blocklengths);
    free(displacements);
    free(types);
    return MPI_ERR_INTERN;
  }

  /* Create a vector representing the massive blocks of full INT_MAX chunks */
  MPI_Datatype chunks_type;
  MPI_Type_vector((int)chunks, INT_MAX, INT_MAX, oldtype, &chunks_type);

  blocklengths[0] = 1;
  displacements[0] = 0;
  types[0] = chunks_type;

  /* If a remainder exists, tie it to the end of the full chunks */
  if (remainder > 0) {
    MPI_Datatype remainder_type;
    MPI_Type_contiguous(remainder, oldtype, &remainder_type);

    blocklengths[1] = 1;
    displacements[1] = (MPI_Aint)chunks * INT_MAX * extent;
    types[1] = remainder_type;
  }

  /* Merge everything seamlessly using an absolute struct layout */
  int status = MPI_Type_create_struct(count, blocklengths, displacements, types,
                                      newtype);

  /* Clean up temporary sub-types
   * Note: they are retained inside the struct wrapper */
  MPI_Type_free(&chunks_type);
  if (remainder > 0) {
    MPI_Type_free(&types[1]);
  }

  free(blocklengths);
  free(displacements);
  free(types);
  return status;
}

/**
 * @brief Implementation of MPI_Allgatherv which allows for size_t arguments
 *
 * Code in part inspired by Gemini AI.
 *
 * @param sendbuf Starting address of send buffer.
 * @param sendcount Number of elements in send buffer.
 * @param sendtype Data type of send buffer elements.
 * @param recvbuf Address of receive buffer.
 * @param recvcounts Integer array (of length group size) containing the number
 * of elements that are to be received from each process.
 * @param displs Integer array (of length group size). Entry i specifies the
 * displacement (relative to recvbuf ) at which to place the incoming data from
 * process i.
 * @param recvtype Data type of receive buffer elements.
 * @param comm MPI communicator.
 */
int MPI_Allgatherv_sizet(const void *sendbuf, size_t sendcount,
                         MPI_Datatype sendtype, void *recvbuf,
                         const size_t recvcounts[], const size_t displs[],
                         MPI_Datatype recvtype, MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(recvtype, &lb, &extent);

  /* Create custom types */
  MPI_Datatype *send_type_mpi = (MPI_Datatype *)malloc(sizeof(MPI_Datatype));
  MPI_Datatype *recv_types_mpi =
      (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));
  MPI_Request *requests = (MPI_Request *)malloc(2 * size * sizeof(MPI_Request));

  if (!send_type_mpi || !recv_types_mpi || !requests) {
    free(send_type_mpi);
    free(recv_types_mpi);
    free(requests);
    return MPI_ERR_INTERN;
  }

  int req_count = 0;
  int status = MPI_SUCCESS;

  /*vPost Non-Blocking Receives with custom structural byte offsets */
  for (int i = 0; i < size; ++i) {
    recv_types_mpi[i] = MPI_DATATYPE_NULL;
    if (recvcounts[i] == 0) continue;

    MPI_Datatype data_layout;
    status = create_large_count_type(recvcounts[i], recvtype, &data_layout);
    if (status != MPI_SUCCESS) break;

    /* Apply 64-bit byte displacements directly to the base buffer pointer */
    MPI_Aint byte_disp = (MPI_Aint)displs[i] * extent;
    MPI_Type_create_hindexed_block(1, 1, &byte_disp, data_layout,
                                   &recv_types_mpi[i]);
    MPI_Type_commit(&recv_types_mpi[i]);
    MPI_Type_free(&data_layout);

    MPI_Irecv(recvbuf, 1, recv_types_mpi[i], i, 0, comm,
              &requests[req_count++]);
  }

  /* Post Non-Blocking Sends using a localized large-count structure */
  if (status == MPI_SUCCESS && sendcount > 0) {
    status = create_large_count_type(sendcount, sendtype, send_type_mpi);
    if (status == MPI_SUCCESS) {
      MPI_Type_commit(send_type_mpi);
      for (int i = 0; i < size; ++i) {
        MPI_Isend(sendbuf, 1, *send_type_mpi, i, 0, comm,
                  &requests[req_count++]);
      }
    }
  }

  /* Complete network operations and perform structural deallocations */
  if (status == MPI_SUCCESS && req_count > 0) {
    status = MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
  }

  /* Garbage collection of custom committed structures */
  if (sendcount > 0 && *send_type_mpi != MPI_DATATYPE_NULL) {
    MPI_Type_free(send_type_mpi);
  }
  for (int i = 0; i < size; ++i) {
    if (recv_types_mpi[i] != MPI_DATATYPE_NULL) {
      MPI_Type_free(&recv_types_mpi[i]);
    }
  }

  free(send_type_mpi);
  free(recv_types_mpi);
  free(requests);
  return status;
}

#endif /* WITH_MPI */

#endif /* SWIFT_CUSTOM_MPI_H */
