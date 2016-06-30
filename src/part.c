/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "error.h"
#include "part.h"

/**
 * @brief Re-link the #gpart%s associated with the list of #part%s.
 *
 * @param parts The list of #part.
 * @param N The number of particles to re-link;
 * @param offset The offset of #part%s relative to the global parts list.
 */
void part_relink_gparts(struct part *parts, size_t N, ptrdiff_t offset) {
  for (size_t k = 0; k < N; k++) {
    if (parts[k].gpart) {
      parts[k].gpart->id_or_neg_offset = -(k + offset);
    }
  }
}

/**
 * @brief Re-link the #gpart%s associated with the list of #part%s.
 *
 * @param gparts The list of #gpart.
 * @param N The number of particles to re-link;
 * @param parts The global part array in which to find the #gpart offsets.
 */
void part_relink_parts(struct gpart *gparts, size_t N, struct part *parts) {
  for (size_t k = 0; k < N; k++) {
    if (gparts[k].id_or_neg_offset <= 0) {
      parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
    }
  }
}

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
MPI_Datatype part_mpi_type;
MPI_Datatype xpart_mpi_type;
MPI_Datatype gpart_mpi_type;
#endif

#ifdef WITH_MPI
/**
 * @brief Registers MPI particle types.
 */
void part_create_mpi_types() {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the part structure. */
  if (MPI_Type_contiguous(sizeof(struct part) / sizeof(unsigned char), MPI_BYTE,
                          &part_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&part_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for parts.");
  }
  if (MPI_Type_contiguous(sizeof(struct xpart) / sizeof(unsigned char),
                          MPI_BYTE, &xpart_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&xpart_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for xparts.");
  }
  if (MPI_Type_contiguous(sizeof(struct gpart) / sizeof(unsigned char),
                          MPI_BYTE, &gpart_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&gpart_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for gparts.");
  }
}
#endif
