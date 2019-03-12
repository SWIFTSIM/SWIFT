/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* This object's header. */
#include "multipole.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef WITH_MPI

/* MPI data type for the multipole transfer and reduction */
MPI_Datatype multipole_mpi_type;
MPI_Op multipole_mpi_reduce_op;

/**
 * @brief Apply a bit-by-bit XOR operattion on #gravity_tensors (i.e. does
 * a^=b).
 *
 * @param a The #gravity_tensors to add to.
 * @param b The #gravity_tensors to add.
 */
void gravity_binary_xor(struct gravity_tensors *a,
                        const struct gravity_tensors *b) {

  char *aa = (char *)a;
  const char *bb = (const char *)b;

  for (size_t i = 0; i < sizeof(struct gravity_tensors); ++i) {
    aa[i] ^= bb[i];
  }
}

/**
 * @brief MPI reduction function for the #gravity_tensors.
 *
 * @param invec Array of #gravity_tensors to read.
 * @param inoutvec Array of #gravity_tensors to read and do the reduction into.
 * @param len The length of the array.
 * @param datatype The MPI type this function acts upon (unused).
 */
void gravity_tensors_mpi_reduce(void *invec, void *inoutvec, int *len,
                                MPI_Datatype *datatype) {

  for (int i = 0; i < *len; ++i) {
    gravity_binary_xor(&((struct gravity_tensors *)inoutvec)[i],
                       &((const struct gravity_tensors *)invec)[i]);
  }
}

void multipole_create_mpi_types(void) {

  /* Create the datatype for multipoles */
  /* We just consider each structure to be a byte field disregarding their */
  /* detailed content */
  if (MPI_Type_contiguous(
          sizeof(struct gravity_tensors) / sizeof(unsigned char), MPI_BYTE,
          &multipole_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&multipole_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for multipole.");
  }

  /* And the reduction operator */
  MPI_Op_create(gravity_tensors_mpi_reduce, 1, &multipole_mpi_reduce_op);
}

#endif
