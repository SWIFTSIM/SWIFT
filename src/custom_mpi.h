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

  // Get the byte size of the receive datatype to calculate 64-bit byte
  // displacements
  MPI_Aint lb, extent;
  MPI_Type_get_extent(recvtype, &lb, &extent);

  // Standard MPI-3 limits local element counts in type creation to 'int'.
  if (sendcount > INT_MAX) {
    return MPI_ERR_COUNT;
  }

  // Allocate arrays for MPI_Alltoallw.
  // Under MPI-3, displacements in MPI_Alltoallw are 'int' arrays containing
  // byte offsets.
  int *sendcounts_mpi = (int *)malloc(size * sizeof(int));
  int *senddispls_mpi = (int *)malloc(size * sizeof(int));
  MPI_Datatype *sendtypes_mpi =
      (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));

  int *recvcounts_mpi = (int *)malloc(size * sizeof(int));
  int *recvdispls_mpi = (int *)malloc(size * sizeof(int));
  MPI_Datatype *recvtypes_mpi =
      (MPI_Datatype *)malloc(size * sizeof(MPI_Datatype));

  // Handle allocation failures gracefully
  if (!sendcounts_mpi || !senddispls_mpi || !sendtypes_mpi || !recvcounts_mpi ||
      !recvdispls_mpi || !recvtypes_mpi) {
    free(sendcounts_mpi);
    free(senddispls_mpi);
    free(sendtypes_mpi);
    free(recvcounts_mpi);
    free(recvdispls_mpi);
    free(recvtypes_mpi);
    return MPI_ERR_INTERN;
  }

  // Initialize tracking arrays
  for (int i = 0; i < size; ++i) {
    sendcounts_mpi[i] = 0;
    senddispls_mpi[i] = 0;
    sendtypes_mpi[i] = sendtype;

    recvcounts_mpi[i] =
        1;  // We receive exactly 1 custom composite datatype from each rank
    recvdispls_mpi[i] = 0;  // The actual memory offset is already hardcoded
                            // inside recvtypes_mpi
  }

  // Every process pushes its localized data slot to all other nodes
  sendcounts_mpi[rank] = (int)sendcount;

  // Create customized datatypes for every incoming slot to safely use 64-bit
  // displacements
  int status = MPI_SUCCESS;
  for (int i = 0; i < size; ++i) {
    if (recvcounts[i] > INT_MAX) {
      status = MPI_ERR_COUNT;
      // Clean up previously successfully created types before escaping
      for (int j = 0; j < i; ++j) {
        MPI_Type_free(&recvtypes_mpi[j]);
      }
      break;
    }

    // Safely project the size_t element displacement onto a byte offset
    // (MPI_Aint)
    MPI_Aint byte_disp = (MPI_Aint)displs[i] * extent;

    // MPI_Type_create_hindexed_block takes an integer blocklength,
    // but crucially uses an MPI_Aint (64-bit) for byte displacements.
    MPI_Type_create_hindexed_block(1, (int)recvcounts[i], &byte_disp, recvtype,
                                   &recvtypes_mpi[i]);
    MPI_Type_commit(&recvtypes_mpi[i]);
  }

  if (status == MPI_SUCCESS) {
    // MPI_Alltoallw accepts 'int[]' for displacements in MPI-3.
    // It reads 'recvdispls_mpi' as a zero-offset baseline because our
    // custom datatypes already have the true 64-bit offset embedded inside
    // them.
    status = MPI_Alltoallw(sendbuf, sendcounts_mpi, senddispls_mpi,
                           sendtypes_mpi, recvbuf, recvcounts_mpi,
                           recvdispls_mpi, recvtypes_mpi, comm);

    // Standard cleanup of the temporary MPI committed objects
    for (int i = 0; i < size; ++i) {
      MPI_Type_free(&recvtypes_mpi[i]);
    }
  }

  // Release allocated heap buffers
  free(sendcounts_mpi);
  free(senddispls_mpi);
  free(sendtypes_mpi);
  free(recvcounts_mpi);
  free(recvdispls_mpi);
  free(recvtypes_mpi);

  return status;
}

#endif /* WITH_MPI */

#endif /* SWIFT_CUSTOM_MPI_H */
