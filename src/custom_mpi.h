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

int MPI_Allgatherv_sizet(const void *sendbuf, size_t sendcount,
                         MPI_Datatype sendtype, void *recvbuf,
                         const size_t recvcounts[], const size_t displs[],
                         MPI_Datatype recvtype, MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  /* Get the byte size of the receive datatype to calculate 64-bit byte
   * displacements */
  MPI_Aint lb, extent;
  MPI_Type_get_extent(recvtype, &lb, &extent);

  /* Standard MPI-3 point-to-point element counts are limited to 'int'. */
  if (sendcount > INT_MAX) {
    return MPI_ERR_COUNT;
  }

  /* Allocate tracking arrays for active transactions
   * Each rank can have at most 1 send and 1 receive request */
  MPI_Request *requests = (MPI_Request *)malloc(2 * size * sizeof(MPI_Request));
  MPI_Datatype *recvtypes_mpi =
      (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));

  if (!requests || !recvtypes_mpi) {
    free(requests);
    free(recvtypes_mpi);
    return MPI_ERR_INTERN;
  }

  /* Initialize all custom datatype slots to null */
  for (int i = 0; i < size; ++i) {
    recvtypes_mpi[i] = MPI_DATATYPE_NULL;
  }

  int req_count = 0;
  int status = MPI_SUCCESS;

  /* Post non-blocking receives ONLY for ranks sending > 0 elements */
  for (int i = 0; i < size; ++i) {
    if (recvcounts[i] == 0) {
      continue; /* Skip completely to avoid passing 0 blocklength to MPI_Type */
    }

    if (recvcounts[i] > INT_MAX) {
      status = MPI_ERR_COUNT;
      break;
    }

    /* Calculate 64-bit byte displacement */
    MPI_Aint byte_disp = (MPI_Aint)displs[i] * extent;
    MPI_Type_create_hindexed_block(1, (int)recvcounts[i], &byte_disp, recvtype,
                                   &recvtypes_mpi[i]);
    MPI_Type_commit(&recvtypes_mpi[i]);

    /* Receive directly into the base of 'recvbuf' */
    MPI_Irecv(recvbuf, 1, recvtypes_mpi[i], i, 0, comm, &requests[req_count++]);
  }

  /* If an error occurred during the receive setup, clean up and exit */
  if (status != MPI_SUCCESS) {
    for (int i = 0; i < size; ++i) {
      if (recvtypes_mpi[i] != MPI_DATATYPE_NULL) {
        MPI_Type_free(&recvtypes_mpi[i]);
      }
    }
    free(requests);
    free(recvtypes_mpi);
    return status;
  }

  /* Post non-blocking sends ONLY if this process actually has data to share */
  if (sendcount > 0) {
    for (int i = 0; i < size; ++i) {
      MPI_Isend(sendbuf, (int)sendcount, sendtype, i, 0, comm,
                &requests[req_count++]);
    }
  }

  /* Wait for all active communications to finish */
  if (req_count > 0) {
    status = MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
  }

  /* Free up memory allocations and clean up custom types */
  for (int i = 0; i < size; ++i) {
    if (recvtypes_mpi[i] != MPI_DATATYPE_NULL) {
      MPI_Type_free(&recvtypes_mpi[i]);
    }
  }

  free(requests);
  free(recvtypes_mpi);

  return status;
}

#endif /* WITH_MPI */

#endif /* SWIFT_CUSTOM_MPI_H */
