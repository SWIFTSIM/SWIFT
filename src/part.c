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
#include "multipole.h"

/* Local headers */
#include "error.h"
#include "part.h"

/**
 * @brief Re-link the #gpart%s associated with the list of #part%s.
 *
 * @param parts The list of #part.
 * @param N The number of particles to re-link;
 * @param offset The offset of #part%s relative to the global parts list.
 */
void part_relink_gparts_to_parts(struct part *parts, size_t N,
                                 ptrdiff_t offset) {
  for (size_t k = 0; k < N; k++) {
    if (parts[k].gpart) {
      parts[k].gpart->id_or_neg_offset = -(k + offset);
    }
  }
}

/**
 * @brief Re-link the #gpart%s associated with the list of #spart%s.
 *
 * @param sparts The list of #spart.
 * @param N The number of s-particles to re-link;
 * @param offset The offset of #spart%s relative to the global sparts list.
 */
void part_relink_gparts_to_sparts(struct spart *sparts, size_t N,
                                  ptrdiff_t offset) {
  for (size_t k = 0; k < N; k++) {
    if (sparts[k].gpart) {
      sparts[k].gpart->id_or_neg_offset = -(k + offset);
    }
  }
}

/**
 * @brief Re-link the #part%s associated with the list of #gpart%s.
 *
 * @param gparts The list of #gpart.
 * @param N The number of particles to re-link;
 * @param parts The global #part array in which to find the #gpart offsets.
 */
void part_relink_parts_to_gparts(struct gpart *gparts, size_t N,
                                 struct part *parts) {
  for (size_t k = 0; k < N; k++) {
    if (gparts[k].type == swift_type_gas) {
      parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
    }
  }
}

/**
 * @brief Re-link the #spart%s associated with the list of #gpart%s.
 *
 * @param gparts The list of #gpart.
 * @param N The number of particles to re-link;
 * @param sparts The global #spart array in which to find the #gpart offsets.
 */
void part_relink_sparts_to_gparts(struct gpart *gparts, size_t N,
                                  struct spart *sparts) {
  for (size_t k = 0; k < N; k++) {
    if (gparts[k].type == swift_type_stars) {
      sparts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
    }
  }
}

/**
 * @brief Re-link both the #part%s and #spart%s associated with the list of
 * #gpart%s.
 *
 * @param gparts The list of #gpart.
 * @param N The number of particles to re-link;
 * @param parts The global #part array in which to find the #gpart offsets.
 * @param sparts The global #spart array in which to find the #gpart offsets.
 */
void part_relink_all_parts_to_gparts(struct gpart *gparts, size_t N,
                                     struct part *parts, struct spart *sparts) {
  for (size_t k = 0; k < N; k++) {
    if (gparts[k].type == swift_type_gas) {
      parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
    } else if (gparts[k].type == swift_type_stars) {
      sparts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
    }
  }
}

/**
 * @brief Verifies that the #gpart, #part and #spart are correctly linked
 * together
 * and that the particle poisitions match.
 *
 * This is a debugging function.
 *
 * @param parts The #part array.
 * @param gparts The #gpart array.
 * @param sparts The #spart array.
 * @param nr_parts The number of #part in the array.
 * @param nr_gparts The number of #gpart in the array.
 * @param nr_sparts The number of #spart in the array.
 * @param verbose Do we report verbosely in case of success ?
 */
void part_verify_links(struct part *parts, struct gpart *gparts,
                       struct spart *sparts, size_t nr_parts, size_t nr_gparts,
                       size_t nr_sparts, int verbose) {

  ticks tic = getticks();

  for (size_t k = 0; k < nr_gparts; ++k) {

    /* We have a real DM particle */
    if (gparts[k].type == swift_type_dark_matter &&
        gparts[k].time_bin != time_bin_not_created) {

      /* Check that it's not linked */
      if (gparts[k].id_or_neg_offset <= 0)
        error("DM gpart particle linked to something !");
    }

    /* We have a gas particle */
    else if (gparts[k].type == swift_type_gas) {

      /* Check that it is linked */
      if (gparts[k].id_or_neg_offset > 0)
        error("Gas gpart not linked to anything!");

      /* Find its link */
      const struct part *part = &parts[-gparts[k].id_or_neg_offset];

      /* Check the reverse link */
      if (part->gpart != &gparts[k]) error("Linking problem!");

      /* Check that the particles are at the same place */
      if (gparts[k].x[0] != part->x[0] || gparts[k].x[1] != part->x[1] ||
          gparts[k].x[2] != part->x[2])
        error(
            "Linked particles are not at the same position!\n"
            "gp->x=[%e %e %e] p->x=[%e %e %e] diff=[%e %e %e]",
            gparts[k].x[0], gparts[k].x[1], gparts[k].x[2], part->x[0],
            part->x[1], part->x[2], gparts[k].x[0] - part->x[0],
            gparts[k].x[1] - part->x[1], gparts[k].x[2] - part->x[2]);

      /* Check that the particles are at the same time */
      if (gparts[k].time_bin != part->time_bin)
        error("Linked particles are not at the same time !");
    }

    else if (gparts[k].type == swift_type_stars) {

      /* Check that it is linked */
      if (gparts[k].id_or_neg_offset > 0)
        error("Stars gpart not linked to anything !");

      /* Find its link */
      const struct spart *spart = &sparts[-gparts[k].id_or_neg_offset];

      /* Check the reverse link */
      if (spart->gpart != &gparts[k]) error("Linking problem !");

      /* Check that the particles are at the same place */
      if (gparts[k].x[0] != spart->x[0] || gparts[k].x[1] != spart->x[1] ||
          gparts[k].x[2] != spart->x[2])
        error(
            "Linked particles are not at the same position !\n"
            "gp->x=[%e %e %e] sp->x=[%e %e %e] diff=[%e %e %e]",
            gparts[k].x[0], gparts[k].x[1], gparts[k].x[2], spart->x[0],
            spart->x[1], spart->x[2], gparts[k].x[0] - spart->x[0],
            gparts[k].x[1] - spart->x[1], gparts[k].x[2] - spart->x[2]);

      /* Check that the particles are at the same time */
      if (gparts[k].time_bin != spart->time_bin)
        error("Linked particles are not at the same time !");
    }
  }

  /* Now check that all parts are linked */
  for (size_t k = 0; k < nr_parts; ++k) {

    /* Ok, there is a link */
    if (parts[k].gpart != NULL) {

      /* Check the link */
      if (parts[k].gpart->id_or_neg_offset != -(ptrdiff_t)k) {
        error("Linking problem !");
      }

      /* Check that the particles are at the same place */
      if (parts[k].x[0] != parts[k].gpart->x[0] ||
          parts[k].x[1] != parts[k].gpart->x[1] ||
          parts[k].x[2] != parts[k].gpart->x[2])
        error("Linked particles are not at the same position !");

      /* Check that the particles are at the same time */
      if (parts[k].time_bin != parts[k].gpart->time_bin)
        error("Linked particles are not at the same time !");
    }
  }

  /* Now check that all sparts are linked */
  for (size_t k = 0; k < nr_sparts; ++k) {

    /* Ok, there is a link */
    if (sparts[k].gpart != NULL) {

      /* Check the link */
      if (sparts[k].gpart->id_or_neg_offset != -(ptrdiff_t)k) {
        error("Linking problem !");

        /* Check that the particles are at the same place */
        if (sparts[k].x[0] != sparts[k].gpart->x[0] ||
            sparts[k].x[1] != sparts[k].gpart->x[1] ||
            sparts[k].x[2] != sparts[k].gpart->x[2])
          error("Linked particles are not at the same position !");

        /* Check that the particles are at the same time */
        if (sparts[k].time_bin != sparts[k].gpart->time_bin)
          error("Linked particles are not at the same time !");
      }
    }
  }

  if (verbose) message("All links OK");
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

#ifdef WITH_MPI
/* MPI data type for the particle transfers */
MPI_Datatype part_mpi_type;
MPI_Datatype xpart_mpi_type;
MPI_Datatype gpart_mpi_type;
MPI_Datatype spart_mpi_type;

/**
 * @brief Registers MPI particle types.
 */
void part_create_mpi_types(void) {

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
  if (MPI_Type_contiguous(sizeof(struct spart) / sizeof(unsigned char),
                          MPI_BYTE, &spart_mpi_type) != MPI_SUCCESS ||
      MPI_Type_commit(&spart_mpi_type) != MPI_SUCCESS) {
    error("Failed to create MPI type for sparts.");
  }
}
#endif
