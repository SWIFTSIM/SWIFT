/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

#if defined(WITH_LOGGER)

/* Some standard headers. */
#include "common_io.h"

#include <errno.h>
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

/* This object's header. */
#include "logger_io.h"

/* Local includes. */
#include "chemistry_io.h"
#include "common_io.h"
#include "cooling_io.h"
#include "dimension.h"
#include "engine.h"
#include "error.h"
#include "gravity_io.h"
#include "gravity_properties.h"
#include "hydro_io.h"
#include "hydro_properties.h"
#include "io_properties.h"
#include "kernel_hydro.h"
#include "parallel_io.h"
#include "part.h"
#include "serial_io.h"
#include "single_io.h"
#include "stars_io.h"
#include "threadpool.h"
#include "tracers_io.h"
#include "units.h"
#include "version.h"
#include "xmf.h"

/**
 * @brief Writes the data array in the index file.
 *
 * @param e The #engine we are writing from.
 * @param f The file to use.
 * @param props The #io_props array.
 * @param n_props The number of element in @props.
 * @param N The number of particles to write.
 */
void write_index_array(const struct engine* e, FILE* f, struct io_props* props,
                       size_t n_props, size_t N) {

  /* Check that the assumptions are corrects */
  if (n_props != 2)
    error("Not implemented: The index file can only write two props.");

  /* Get a few variables */
  const size_t type_size0 = io_sizeof_type(props[0].type) * props[0].dimension;
  const size_t type_size1 = io_sizeof_type(props[1].type) * props[1].dimension;
  const size_t type_size = type_size0 + type_size1;

  /* Convert FILE to int */
  int fd = fileno(f);
  if (fd == -1) {
    error("Failed to get the integer descriptor");
  }

  /* Get a few variables for the mapping */
  char* data;
  const size_t offset = ftell(f);
  const size_t count = N * type_size + offset;

  /* Truncate the file to the correct length. */
  if (ftruncate(fd, count) != 0) {
    error("Failed to truncate dump file (%s).", strerror(errno));
  }

  /* Map the file */
  if ((data = mmap(NULL, count, PROT_WRITE, MAP_SHARED, fd, 0)) == MAP_FAILED) {
    error("Failed to allocate map of size %zi bytes (%s).", count,
          strerror(errno));
  }

  /* Copy the data into the file */
  char* first = data + offset;
  char* second = data + type_size0 + offset;
  for (size_t i = 0; i < N; i++) {
    memcpy(first + i * type_size, props[0].field + i * props[0].partSize,
           type_size0);
    memcpy(second + i * type_size, props[1].field + i * props[1].partSize,
           type_size1);
  }

  /* Unmap the data in memory. */
  if (munmap(data, count) != 0) {
    error("Failed to unmap dump data (%s).", strerror(errno));
  }

  /* Move the file position */
  fseek(f, count, SEEK_SET);
}

/**
 * @brief Write the history (created or deleted) for all the particles type.
 *
 * @param history The list of history to write.
 * @param e The #engine.
 * @param f The opened file to use.
 */
void logger_write_history(struct logger_history* history, struct engine* e,
                          FILE* f) {

  /* Write the number of particles. */
  uint64_t size[swift_type_count];
  for (int i = 0; i < swift_type_count; i++) {
    size[i] = history[i].size;
  }
  fwrite(size, sizeof(uint64_t), swift_type_count, f);

  /* Write the data */
  for (int i = 0; i < swift_type_count; i++) {
    logger_history_write(&history[i], e, f);
  }
}

/**
 * @brief Writes a logger index file
 *
 * @param log The #logger_writer.
 * @param e The engine containing all the system.
 *
 * Creates an output file and writes the offset and id of particles
 * contained in the engine. If such a file already exists, it is erased and
 * replaced by the new one.
 *
 * An index file is constructed by writing first a few variables (e.g. time,
 * number of particles, if the file is sorted, ...) and then an array of index
 * and offset for each particle type.
 *
 * Calls #error() if an error occurs.
 *
 */
void logger_write_index_file(struct logger_writer* log, struct engine* e) {

  struct part* parts = e->s->parts;
  struct xpart* xparts = e->s->xparts;
  struct gpart* gparts = e->s->gparts;
  struct spart* sparts = e->s->sparts;

  /* Number of particles currently in the arrays */
  const size_t Ntot = e->s->nr_gparts;
  const size_t Ngas = e->s->nr_parts;
  const size_t Nstars = e->s->nr_sparts;
  /* const size_t Nblackholes = e->s->nr_bparts; */

  /* Number of particles that we will write */
  const size_t Ntot_written =
      e->s->nr_gparts - e->s->nr_inhibited_gparts - e->s->nr_extra_gparts;
  const size_t Ngas_written =
      e->s->nr_parts - e->s->nr_inhibited_parts - e->s->nr_extra_parts;
  const size_t Nstars_written =
      e->s->nr_sparts - e->s->nr_inhibited_sparts - e->s->nr_extra_sparts;
  const size_t Nblackholes_written =
      e->s->nr_bparts - e->s->nr_inhibited_bparts - e->s->nr_extra_bparts;
  const size_t Nbaryons_written =
      Ngas_written + Nstars_written + Nblackholes_written;
  const size_t Ndm_written =
      Ntot_written > 0 ? Ntot_written - Nbaryons_written : 0;

  /* Format things in a Gadget-friendly array */
  uint64_t N_total[swift_type_count] = {
      (uint64_t)Ngas_written,   (uint64_t)Ndm_written,        0, 0,
      (uint64_t)Nstars_written, (uint64_t)Nblackholes_written};

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%.100s_%04i_%04i.index",
           e->logger->base_name, engine_rank, log->index_file_number);
  log->index_file_number++;

  /* Open file (include reading for mmap) */
  FILE* f = NULL;
  f = fopen(fileName, "w+b");

  if (f == NULL) {
    error("Failed to open file %s", fileName);
  }

  /* Write double time */
  fwrite(&e->time, sizeof(double), 1, f);

  /* Write integer time */
  fwrite(&e->ti_current, sizeof(integertime_t), 1, f);

  /* Write number of particles */
  fwrite(N_total, sizeof(uint64_t), swift_type_count, f);

  /* Write if the file is sorted */
  const char sorted = 0;
  fwrite(&sorted, sizeof(char), 1, f);

  /* Ensure the data to be aligned */
  size_t cur_pos = ftell(f);
  size_t d_align = ((cur_pos + 7) & ~7) - cur_pos;
  if (d_align > 0) {
    long int tmp = 0;
    /* Fill the memory with garbage */
    fwrite(&tmp, d_align, 1, f);
  }

  /* Loop over all particle types */
  for (int ptype = 0; ptype < swift_type_count; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N_total[ptype] == 0) continue;

    /* Number of properties (the code cannot deal with more than two props
       per particle type) */
    size_t N = 0;
    int num_fields = 0;
    struct io_props list[2];

    struct part* parts_written = NULL;
    struct xpart* xparts_written = NULL;
    struct gpart* gparts_written = NULL;
    struct velociraptor_gpart_data* gpart_group_data_written = NULL;
    struct spart* sparts_written = NULL;
    struct bpart* bparts_written = NULL;

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case swift_type_gas:
        if (Ngas == Ngas_written) {

          /* No inhibted particles: easy case */
          N = Ngas;
          num_fields += hydro_write_index(parts, xparts, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ngas_written;

          /* Allocate temporary arrays */
          if (swift_memalign("parts_written", (void**)&parts_written,
                             part_align,
                             Ngas_written * sizeof(struct part)) != 0)
            error("Error while allocating temporary memory for parts");
          if (swift_memalign("xparts_written", (void**)&xparts_written,
                             xpart_align,
                             Ngas_written * sizeof(struct xpart)) != 0)
            error("Error while allocating temporary memory for xparts");

          /* Collect the particles we want to write */
          io_collect_parts_to_write(parts, xparts, parts_written,
                                    xparts_written, Ngas, Ngas_written);

          /* Select the fields to write */
          num_fields += hydro_write_index(parts_written, xparts_written, list);
        }
        break;

      case swift_type_dark_matter:
        if (Ntot == Ndm_written) {

          /* This is a DM-only run without inhibited particles */
          N = Ntot;
          num_fields += darkmatter_write_index(gparts, list);
        } else {

          /* Ok, we need to fish out the particles we want */
          N = Ndm_written;

          /* Allocate temporary array */
          if (swift_memalign("gparts_written", (void**)&gparts_written,
                             gpart_align,
                             Ndm_written * sizeof(struct gpart)) != 0)
            error("Error while allocating temporary memory for gparts");

          /* Collect the non-inhibited DM particles from gpart */
          const int with_stf = 0;
          io_collect_gparts_to_write(gparts, e->s->gpart_group_data,
                                     gparts_written, gpart_group_data_written,
                                     Ntot, Ndm_written, with_stf);

          /* Select the fields to write */
          num_fields += darkmatter_write_index(gparts_written, list);
        }
        break;

      case swift_type_stars:
        if (Nstars == Nstars_written) {

          /* No inhibted particles: easy case */
          N = Nstars;
          num_fields += stars_write_index(sparts, list);

        } else {

          /* Ok, we need to fish out the particles we want */
          N = Nstars_written;

          /* Allocate temporary arrays */
          if (swift_memalign("sparts_written", (void**)&sparts_written,
                             spart_align,
                             Nstars_written * sizeof(struct spart)) != 0)
            error("Error while allocating temporary memory for sparts");

          /* Collect the particles we want to write */
          io_collect_sparts_to_write(sparts, sparts_written, Nstars,
                                     Nstars_written);

          /* Select the fields to write */
          num_fields += stars_write_index(sparts_written, list);
        }
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    if (num_fields != 2) {
      error(
          "The code expects only two fields per particle type for the logger");
    }

    /* Write ids */
    write_index_array(e, f, list, num_fields, N);

    /* Free temporary arrays */
    if (parts_written) swift_free("parts_written", parts_written);
    if (xparts_written) swift_free("xparts_written", xparts_written);
    if (gparts_written) swift_free("gparts_written", gparts_written);
    if (gpart_group_data_written)
      swift_free("gpart_group_written", gpart_group_data_written);
    if (sparts_written) swift_free("sparts_written", sparts_written);
    if (bparts_written) swift_free("bparts_written", bparts_written);
  }

  /* Write the particles created */
  logger_write_history(log->history_new, e, f);

  /* Write the particles removed */
  logger_write_history(log->history_removed, e, f);

  /* Close file */
  fclose(f);
}

/**
 * @brief Write the parameters into a yaml file.
 *
 * @params log The #logger.
 * @params e The #engine.
 */
void logger_write_description(struct logger_writer* log, struct engine* e) {
  /* Only the master writes the description */
  if (engine_rank != 0) {
    return;
  }
  /* const struct unit_system *internal_units = e->internal_units; */
  /* const struct unit_system *snapshot_units = e->snapshot_units; */

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "%.100s.yml", e->logger->base_name);

  /* Open file */
  FILE* f = NULL;
  f = fopen(fileName, "wb");

  if (f == NULL) {
    error("Failed to open file %s", fileName);
  }

  /* TODO Write stuff */

  /* Close file */
  fclose(f);
}

#endif /* WITH_LOGGER */
