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
#include "part.h"

#ifdef WITH_MPI
/**
 * @brief Registers and returns an MPI type for the particles
 *
 * @param part_type The type container
 */
void part_create_mpi_type(MPI_Datatype* part_type) {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the part structure. */
  MPI_Type_contiguous(sizeof(struct part) / sizeof(unsigned char), MPI_BYTE,
                      part_type);
  MPI_Type_commit(part_type);
}

/**
 * @brief Registers and returns an MPI type for the xparticles
 *
 * @param xpart_type The type container
 */
void xpart_create_mpi_type(MPI_Datatype* xpart_type) {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the part structure. */
  MPI_Type_contiguous(sizeof(struct xpart) / sizeof(unsigned char), MPI_BYTE,
                      xpart_type);
  MPI_Type_commit(xpart_type);
}

/**
 * @brief Registers and returns an MPI type for the gparticles
 *
 * @param gpart_type The type container
 */
void gpart_create_mpi_type(MPI_Datatype* gpart_type) {

  /* This is not the recommended way of doing this.
     One should define the structure field by field
     But as long as we don't do serialization via MPI-IO
     we don't really care.
     Also we would have to modify this function everytime something
     is added to the part structure. */
  MPI_Type_contiguous(sizeof(struct gpart) / sizeof(unsigned char), MPI_BYTE,
                      gpart_type);
  MPI_Type_commit(gpart_type);
}

#endif
