/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

#ifdef HAVE_POSIX_FALLOCATE /* Are we on a sensible platform? */
#ifdef WITH_LOGGER

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Define the particles first */
#include "part.h"

/* This object's header. */
#include "logger.h"

/* Local headers. */
#include "atomic.h"
#include "dump.h"
#include "engine.h"
#include "error.h"
#include "gravity_logger.h"
#include "hydro_io.h"
#include "stars_io.h"
#include "units.h"

/*
 * Thoses are definitions from the format and therefore should not be changed!
 */
/* Number of bytes for a mask. */
// TODO change this to number of bits
#define logger_mask_size 2

/* Number of bits for record header. */
#define logger_header_bytes 8

/* Number bytes for an offset. */
#define logger_offset_size logger_header_bytes - logger_mask_size

/* Number of bytes for the file format information. */
#define logger_format_size 20

/* Number of bytes for the labels in the header. */
#define logger_label_size 20

char logger_file_format[logger_format_size] = "SWIFT_LOGGER";

/*
 * The two following defines need to correspond to the list's order
 * in logger_init_masks.
 */
/* Index of the special flags in the list of masks */
#define logger_index_special_flags 0
/* Index of the timestamp in the list of masks */
#define logger_index_timestamp 1

/**
 * @brief Write the header of a record (offset + mask).
 *
 * This is maybe broken for big(?) endian.
 *
 * @param buff The buffer where to write the mask and offset.
 * @param mask The mask to write inside the buffer.
 * @param offset The offset of the previous record.
 * @param offset_new The offset of the current record.
 *
 * @return updated buff
 */
char *logger_write_record_header(char *buff, const unsigned int *mask,
                                 const size_t *offset,
                                 const size_t offset_new) {
  /* write mask. */
  memcpy(buff, mask, logger_mask_size);
  buff += logger_mask_size;

  /* write offset. */
  uint64_t diff_offset = offset_new - *offset;
  memcpy(buff, &diff_offset, logger_offset_size);
  buff += logger_offset_size;

  return buff;
}

/**
 * @brief Write to the dump.
 *
 * @param d #dump file
 * @param offset (return) offset of the data
 * @param size number of bytes to write
 * @param p pointer to the data
 */
void logger_write_data(struct dump *d, size_t *offset, size_t size,
                       const void *p) {
  /* get buffer. */
  char *buff = dump_get(d, size, offset);

  /* write data to the buffer. */
  memcpy(buff, p, size);

  /* Update offset to end of record. */
  *offset += size;
}

/**
 * @brief log all particles in the engine.
 *
 * @param log The #logger_writer
 * @param e The #engine
 */
void logger_log_all_particles(struct logger_writer *log,
                              const struct engine *e) {

  /* Ensure that enough space is available. */
  logger_ensure_size(log, e->s->nr_parts, e->s->nr_gparts, e->s->nr_sparts);

  /* some constants. */
  const struct space *s = e->s;

  /* log the parts. */
  logger_log_parts(log, s->parts, s->xparts, s->nr_parts, e,
                   /* log_all_fields= */ 1, /* flag= */ 0, /* flag_data= */ 0);

  /* log the gparts */
  logger_log_gparts(log, s->gparts, s->nr_gparts, e,
                    /* log_all_fields= */ 1, /* flag= */ 0,
                    /* flag_data= */ 0);

  /* log the parts */
  logger_log_sparts(log, s->sparts, s->nr_sparts, e,
                    /* log_all_fields= */ 1, /* flag= */ 0,
                    /* flag_data= */ 0);

  if (e->total_nr_bparts > 0) error("Not implemented");
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #logger_writer
 * @param p The #part to copy.
 * @param xp The #xpart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data for the special flags.
 */
void logger_copy_part_fields(const struct logger_writer *log,
                             const struct part *p, const struct xpart *xp,
                             const struct engine *e, unsigned int mask,
                             size_t *offset, size_t offset_new, char *buff,
                             const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = logger_write_record_header(buff, &mask, offset, offset_new);

  /* Write the hydro fields */
  buff = hydro_logger_write_particle(log->mask_data_pointers.hydro, p, xp,
                                     &mask, buff);

  /* Special flags */
  if (mask & log->logger_mask_data[logger_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->logger_mask_data[logger_index_special_flags].size);
    buff += log->logger_mask_data[logger_index_special_flags].size;
    mask &= ~log->logger_mask_data[logger_index_special_flags].mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in parts. %u", mask);
  }
#endif
}

/**
 * @brief Dump a #part to the log.
 *
 * @param log The #logger_writer
 * @param p The #part to dump.
 * @param xp The #xpart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_part(struct logger_writer *log, const struct part *p,
                     struct xpart *xp, const struct engine *e,
                     const int log_all_fields,
                     const enum logger_special_flags flag,
                     const int flag_data) {

  logger_log_parts(log, p, xp, /* count= */ 1, e, log_all_fields, flag,
                   flag_data);
}

/**
 * @brief Dump a group of #part to the log.
 *
 * @param log The #logger_writer.
 * @param p The #part to dump.
 * @param xp The #xpart to dump.
 * @param count The number of particle to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_parts(struct logger_writer *log, const struct part *p,
                      struct xpart *xp, int count, const struct engine *e,
                      const int log_all_fields,
                      const enum logger_special_flags flag,
                      const int flag_data) {

  /* Build the special flag */
  const uint32_t special_flags = logger_pack_flags_and_data(flag, flag_data);

  /* Compute the size of the buffer. */
  size_t size_total = 0;
  if (log_all_fields) {
    size_total = count * (log->max_size_record_part + logger_header_bytes);
  } else {
    for (int i = 0; i < count; i++) {
      unsigned int mask = 0;
      size_t size = 0;
      hydro_logger_compute_size_and_mask(log->mask_data_pointers.hydro, &p[i],
                                         &xp[i], log_all_fields, &size, &mask);
      size_total += size + logger_header_bytes;
    }
  }

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);

  /* Write the particles */
  for (int i = 0; i < count; i++) {
    /* Get the masks */
    size_t size = 0;
    unsigned int mask = 0;
    hydro_logger_compute_size_and_mask(log->mask_data_pointers.hydro, &p[i],
                                       &xp[i], log_all_fields, &size, &mask);
    size += logger_header_bytes;

    if (special_flags != 0) {
      mask |= log->logger_mask_data[logger_index_special_flags].mask;
    }

    /* Copy everything into the buffer */
    logger_copy_part_fields(log, &p[i], &xp[i], e, mask,
                            &xp[i].logger_data.last_offset, offset_new, buff,
                            special_flags);

    /* Write the particle into the history if needed. */
    if (flag & logger_flag_create || flag & logger_flag_mpi_enter) {
      logger_history_log(&log->history_new[swift_type_gas], p->id,
                         xp->logger_data.last_offset);

    } else if (flag & logger_flag_change_type || flag & logger_flag_delete ||
               flag & logger_flag_mpi_exit) {
      logger_history_log(&log->history_removed[swift_type_gas], p->id,
                         xp->logger_data.last_offset);
    }

    /* Update the pointers */
    xp[i].logger_data.last_offset = offset_new;
    xp[i].logger_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #logger_writer.
 * @param sp The #spart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data for the special flags.
 */
void logger_copy_spart_fields(const struct logger_writer *log,
                              const struct spart *sp, const struct engine *e,
                              unsigned int mask, size_t *offset,
                              size_t offset_new, char *buff,
                              const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = logger_write_record_header(buff, &mask, offset, offset_new);

  /* Write the stellar fields */
  buff = stars_logger_write_particle(log->mask_data_pointers.stars, sp, &mask,
                                     buff);

  /* Special flags */
  if (mask & log->logger_mask_data[logger_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->logger_mask_data[logger_index_special_flags].size);
    buff += log->logger_mask_data[logger_index_special_flags].size;
    mask &= ~log->logger_mask_data[logger_index_special_flags].mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in sparts. %u", mask);
  }
#endif
}

/**
 * @brief Dump a #spart to the log.
 *
 * @param log The #logger_writer
 * @param sp The #spart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_spart(struct logger_writer *log, struct spart *sp,
                      const struct engine *e, const int log_all_fields,
                      const enum logger_special_flags flag,
                      const int flag_data) {

  logger_log_sparts(log, sp, /* count */ 1, e, log_all_fields, flag, flag_data);
}

/**
 * @brief Dump a group of #spart to the log.
 *
 * @param log The #logger_writer
 * @param sp The #spart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param count The number of particle to dump.
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_sparts(struct logger_writer *log, struct spart *sp, int count,
                       const struct engine *e, const int log_all_fields,
                       const enum logger_special_flags flag,
                       const int flag_data) {
  /* Build the special flag */
  const uint32_t special_flags = logger_pack_flags_and_data(flag, flag_data);

  /* Compute the size of the buffer. */
  size_t size_total = 0;
  if (log_all_fields) {
    size_total = count * (log->max_size_record_spart + logger_header_bytes);
  } else {
    for (int i = 0; i < count; i++) {
      unsigned int mask = 0;
      size_t size = 0;
      stars_logger_compute_size_and_mask(log->mask_data_pointers.stars, &sp[i],
                                         log_all_fields, &size, &mask);
      size_total += size + logger_header_bytes;
    }
  }

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);

  for (int i = 0; i < count; i++) {
    /* Get the masks */
    size_t size = 0;
    unsigned int mask = 0;
    stars_logger_compute_size_and_mask(log->mask_data_pointers.stars, &sp[i],
                                       log_all_fields, &size, &mask);
    size += logger_header_bytes;

    if (special_flags != 0) {
      mask |= log->logger_mask_data[logger_index_special_flags].mask;
    }

    /* Copy everything into the buffer */
    logger_copy_spart_fields(log, &sp[i], e, mask,
                             &sp[i].logger_data.last_offset, offset_new, buff,
                             special_flags);

    /* Write the particle into the history if needed. */
    if (flag & logger_flag_create || flag & logger_flag_mpi_enter) {
      logger_history_log(&log->history_new[swift_type_stars], sp->id,
                         sp->logger_data.last_offset);

    } else if (flag & logger_flag_change_type || flag & logger_flag_delete ||
               flag & logger_flag_mpi_exit) {
      logger_history_log(&log->history_removed[swift_type_stars], sp->id,
                         sp->logger_data.last_offset);
    }

    /* Update the pointers */
    sp[i].logger_data.last_offset = offset_new;
    sp[i].logger_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }
}

/**
 * @brief Copy the particle fields into a given buffer.
 *
 * @param log The #logger_writer.
 * @param gp The #gpart to copy.
 * @param e The #engine.
 * @param mask The mask for the fields to write.
 * @param offset The offset to the previous log.
 * @param offset_new The offset of the current record.
 * @param buff The buffer to use when writing.
 * @param special_flags The data of the special flag.
 */
void logger_copy_gpart_fields(const struct logger_writer *log,
                              const struct gpart *gp, const struct engine *e,
                              unsigned int mask, size_t *offset,
                              size_t offset_new, char *buff,
                              const uint32_t special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (mask == 0) {
    error("You should always log at least one field.");
  }
#endif

  /* Write the header. */
  buff = logger_write_record_header(buff, &mask, offset, offset_new);

  /* Write the hydro fields */
  buff = gravity_logger_write_particle(log->mask_data_pointers.gravity, gp,
                                       &mask, buff);

  /* Special flags */
  if (mask & log->logger_mask_data[logger_index_special_flags].mask) {
    memcpy(buff, &special_flags,
           log->logger_mask_data[logger_index_special_flags].size);
    buff += log->logger_mask_data[logger_index_special_flags].size;
    mask &= ~log->logger_mask_data[logger_index_special_flags].mask;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (mask) {
    error("Requested logging of values not present in gparts. %u", mask);
  }
#endif
}

/**
 * @brief Dump a #gpart to the log.
 *
 * @param log The #logger_writer
 * @param p The #gpart to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_gpart(struct logger_writer *log, struct gpart *p,
                      const struct engine *e, const int log_all_fields,
                      const enum logger_special_flags flag,
                      const int flag_data) {
  logger_log_gparts(log, p, /* count */ 1, e, log_all_fields, flag, flag_data);
}

/**
 * @brief Dump a group of #gpart to the log.
 *
 * @param log The #logger_writer
 * @param p The #gpart to dump.
 * @param count The number of particle to dump.
 * @param e The #engine.
 * @param log_all_fields Should we log all the fields?
 * @param flag The value of the special flags.
 * @param flag_data The data to write for the flag.
 */
void logger_log_gparts(struct logger_writer *log, struct gpart *p, int count,
                       const struct engine *e, const int log_all_fields,
                       const enum logger_special_flags flag,
                       const int flag_data) {
  /* Build the special flag */
  const uint32_t special_flags = logger_pack_flags_and_data(flag, flag_data);

  /* Compute the size of the buffer. */
  size_t size_total = 0;
  if (log_all_fields) {
    size_total = count * (log->max_size_record_gpart + logger_header_bytes);
  } else {
    for (int i = 0; i < count; i++) {
      /* Log only the dark matter */
      if (p[i].type != swift_type_dark_matter) continue;

      unsigned int mask = 0;
      size_t size = 0;
      gravity_logger_compute_size_and_mask(log->mask_data_pointers.gravity,
                                           &p[i], log_all_fields, &size, &mask);
      size_total += size + logger_header_bytes;
    }
  }

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size_total, &offset_new);

  for (int i = 0; i < count; i++) {
    /* Log only the dark matter */
    if (p[i].type != swift_type_dark_matter) continue;

    /* Get the masks */
    size_t size = 0;
    unsigned int mask = 0;
    gravity_logger_compute_size_and_mask(log->mask_data_pointers.gravity, &p[i],
                                         log_all_fields, &size, &mask);
    size += logger_header_bytes;

    if (special_flags != 0) {
      mask |= log->logger_mask_data[logger_index_special_flags].mask;
    }

    /* Copy everything into the buffer */
    logger_copy_gpart_fields(log, &p[i], e, mask, &p[i].logger_data.last_offset,
                             offset_new, buff, special_flags);

    /* Write the particle into the history if needed. */
    if (flag & logger_flag_create || flag & logger_flag_mpi_enter) {
      logger_history_log(&log->history_new[swift_type_dark_matter],
                         p->id_or_neg_offset, p->logger_data.last_offset);

    } else if (flag & logger_flag_change_type || flag & logger_flag_delete ||
               flag & logger_flag_mpi_exit) {
      logger_history_log(&log->history_removed[swift_type_dark_matter],
                         p->id_or_neg_offset, p->logger_data.last_offset);
    }

    /* Update the pointers */
    p[i].logger_data.last_offset = offset_new;
    p[i].logger_data.steps_since_last_output = 0;
    buff += size;
    offset_new += size;
  }
}

/**
 * @brief write a timestamp
 *
 * @param log The #logger_writer
 * @param timestamp time to write
 * @param time time or scale factor
 * @param offset Pointer to the offset of the previous log of this particle;
 * (return) offset of this log.
 */
void logger_log_timestamp(struct logger_writer *log, integertime_t timestamp,
                          double time, size_t *offset) {
  struct dump *dump = &log->dump;
  /* Start by computing the size of the message. */
  const int size =
      log->logger_mask_data[logger_index_timestamp].size + logger_header_bytes;

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  unsigned int mask = log->logger_mask_data[logger_index_timestamp].mask;
  buff = logger_write_record_header(buff, &mask, offset, offset_new);

  /* Store the timestamp. */
  memcpy(buff, &timestamp, sizeof(integertime_t));
  buff += sizeof(integertime_t);

  /* Store the time. */
  memcpy(buff, &time, sizeof(double));

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief Ensure that the buffer is large enough for a step.
 *
 * Check if logger parameters are large enough to write all particles
 * and ensure that enough space is available in the buffer.
 *
 * @param log The #logger_writer
 * @param total_nr_parts total number of part
 * @param total_nr_gparts total number of gpart
 * @param total_nr_sparts total number of spart
 */
void logger_ensure_size(struct logger_writer *log, size_t total_nr_parts,
                        size_t total_nr_gparts, size_t total_nr_sparts) {

  /* count part memory */
  size_t limit = 0;

  /* count part memory */
  limit += total_nr_parts;

  /* count gpart memory */
  limit += total_nr_gparts;

  /* count spart memory. */
  limit += total_nr_sparts;

  // TODO improve estimate with the size of each particle
  limit *= log->max_record_size;

  /* ensure enough space in dump */
  dump_ensure(&log->dump, limit, log->buffer_scale * limit);
}

/** @brief Generate the name of the dump files
 *
 * @param log The #logger_writer.
 * @param filename The filename of the dump file.
 */
void logger_get_dump_name(struct logger_writer *log, char *filename) {
  sprintf(filename, "%s_%04i.dump", log->base_name, engine_rank);
}

/**
 * @brief Initialize the variable logger_mask_data.
 *
 * @param log The #logger_writer.
 * @param e The #engine.
 */
void logger_init_masks(struct logger_writer *log, const struct engine *e) {
  /* Set the pointers to 0 */
  log->mask_data_pointers.hydro = NULL;
  log->mask_data_pointers.stars = NULL;
  log->mask_data_pointers.gravity = NULL;

  struct mask_data list[100];
  int num_fields = 0;

  /* The next fields must be the two first ones. */
  /* Add the special flags (written manually => no need of offset) */
  if (logger_index_special_flags != 0) {
    error("Expecting the special flags to be the first element.");
  }
  list[logger_index_special_flags] =
      logger_create_mask_entry("SpecialFlags", sizeof(int));
  num_fields += 1;

  /* Add the timestamp */
  if (logger_index_timestamp != 1) {
    error("Expecting the timestamp to be the first element.");
  }
  list[logger_index_timestamp] = logger_create_mask_entry(
      "Timestamp", sizeof(integertime_t) + sizeof(double));
  list[num_fields].type = mask_type_timestep;  // flag it as timestamp
  num_fields += 1;

  // TODO add chemistry, cooling, ... + xpart + spart

  /* Get all the fields that need to be written for the hydro. */
  struct mask_data *tmp = &list[num_fields];

  /* Set the mask_data_pointers */
  log->mask_data_pointers.hydro = tmp;

  /* Set the masks */
  int tmp_num_fields = hydro_logger_writer_populate_mask_data(tmp);
  /* Set the particle type */
  for (int i = 0; i < tmp_num_fields; i++) {
    tmp[i].type = mask_type_gas;
  }
  num_fields += tmp_num_fields;

  /* Get all the fields that need to be written for the stars. */
  tmp = &list[num_fields];

  /* Set the mask_data_pointers */
  log->mask_data_pointers.stars = tmp;

  /* Set the masks */
  tmp_num_fields = stars_logger_writer_populate_mask_data(tmp);
  /* Set the particle type */
  for (int i = 0; i < tmp_num_fields; i++) {
    tmp[i].type = mask_type_stars;
  }
  num_fields += tmp_num_fields;

  /* Get all the fields that need to be written for the gravity. */
  tmp = &list[num_fields];

  /* Set the mask_data_pointers */
  log->mask_data_pointers.gravity = tmp;

  /* Set the masks */
  tmp_num_fields = gravity_logger_writer_populate_mask_data(tmp);
  /* Set the particle type */
  for (int i = 0; i < tmp_num_fields; i++) {
    tmp[i].type = mask_type_dark_matter;
  }
  num_fields += tmp_num_fields;

  /* Set the masks and ensure to have only one for the common fields
     (e.g. Coordinates).
     Initially we have (Name, mask, part_type):
      - Coordinates, 0, 0
      - Velocity, 0, 0
      - Coordinates, 0, 1

      And get:
      - Coordinates, 1, 0
      - Velocity, 2, 0
      - Coordinates, 1, 1
  */
  int mask = 0;
  for (int i = 0; i < num_fields; i++) {
    /* Skip the elements already processed. */
    if (list[i].mask != 0) {
      continue;
    }
    const char *name = list[i].name;
    list[i].mask = 1 << mask;

    /* Check if the field exists in the other particle type. */
    for (int j = i + 1; j < num_fields; j++) {
      /* Check if the name is the same */
      if (strcmp(name, list[j].name) == 0) {
        /* Check if the data are the same */
        if (list[i].size != list[j].size) {
          error("Found two same fields but with different data size (%s).",
                name);
        }

        list[j].mask = 1 << mask;
      }
    }
    mask += 1;
  }

  /* Check that we have enough available flags. */
  if (mask >= 8 * logger_mask_size) {
    error(
        "Not enough available flags for all the fields. "
        "Please reduce the number of output fields.");
  }

  /* Save the data */
  size_t size_list = sizeof(struct mask_data) * num_fields;
  log->logger_mask_data = (struct mask_data *)malloc(size_list);
  memcpy(log->logger_mask_data, list, size_list);

  /* Update the pointers */
  if (log->mask_data_pointers.hydro != NULL) {
    log->mask_data_pointers.hydro =
        log->logger_mask_data + (log->mask_data_pointers.hydro - list);
  }
  if (log->mask_data_pointers.stars != NULL) {
    log->mask_data_pointers.stars =
        log->logger_mask_data + (log->mask_data_pointers.stars - list);
  }
  if (log->mask_data_pointers.gravity != NULL) {
    log->mask_data_pointers.gravity =
        log->logger_mask_data + (log->mask_data_pointers.gravity - list);
  }

  /* Compute the maximal size of the records. */
  log->max_size_record_part = 0;
  for (int i = 0; i < hydro_logger_field_count; i++) {
    log->max_size_record_part += log->mask_data_pointers.hydro[i].size;
  }

  log->max_size_record_gpart = 0;
  for (int i = 0; i < gravity_logger_field_count; i++) {
    log->max_size_record_gpart += log->mask_data_pointers.gravity[i].size;
  }

  log->max_size_record_spart = 0;
  for (int i = 0; i < stars_logger_field_count; i++) {
    log->max_size_record_spart += log->mask_data_pointers.stars[i].size;
  }

  /* Set the counter */
  log->logger_count_mask = num_fields;

#ifdef SWIFT_DEBUG_CHECKS
  message("The logger contains the following masks:");
  for (int i = 0; i < log->logger_count_mask; i++) {
    message("%20s:\t mask=%03u\t size=%i", log->logger_mask_data[i].name,
            log->logger_mask_data[i].mask, log->logger_mask_data[i].size);
  }
#endif
}

/**
 * @brief intialize the logger structure
 *
 * @param log The #logger_writer
 * @param e The #engine.
 * @param params The #swift_params
 */
void logger_init(struct logger_writer *log, const struct engine *e,
                 struct swift_params *params) {
  /* read parameters. */
  log->delta_step = parser_get_param_int(params, "Logger:delta_step");
  size_t buffer_size =
      parser_get_opt_param_float(params, "Logger:initial_buffer_size", 0.5) *
      1e9;
  log->buffer_scale =
      parser_get_opt_param_float(params, "Logger:buffer_scale", 10);
  parser_get_param_string(params, "Logger:basename", log->base_name);

  log->index.mem_frac =
      parser_get_opt_param_float(params, "Logger:index_mem_frac", 0.05);

  /* Initialize the logger_mask_data */
  logger_init_masks(log, e);

  /* set initial value of parameters. */
  log->timestamp_offset = 0;
  log->index.dump_size_last_output = 0;
  log->index_file_number = 0;

  /* generate dump filename. */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  logger_get_dump_name(log, logger_name_file);

  /* Compute max size for a particle record. */
  int max_size = logger_offset_size + logger_mask_size;

  /* Loop over all fields except timestamp. */
  for (int i = 0; i < log->logger_count_mask; i++) {
    /* Skip the timestamp */
    if (i == logger_index_timestamp) continue;

    max_size += log->logger_mask_data[i].size;
  }
  log->max_record_size = max_size;

  /* init dump. */
  dump_init(&log->dump, logger_name_file, buffer_size);

  /* Read the maximal size of the history. */
  const float max_memory_size =
      parser_get_opt_param_float(params, "Logger:maximal_memory_size", 1.);
  log->maximal_size_history =
      max_memory_size / sizeof(struct logger_index_data);

  if (e->nodeID == 0) {
    message("Maximal memory size for the logger history: %g GB",
            max_memory_size);
  }

  /* initialize the history */
  for (int i = 0; i < swift_type_count; i++) {
    logger_history_init(&log->history_removed[i]);
    logger_history_init(&log->history_new[i]);
  }
}

/**
 * @brief Close dump file and desallocate memory
 *
 * @param log The #logger_writer
 */
void logger_free(struct logger_writer *log) {
  dump_close(&log->dump);

  free(log->logger_mask_data);
  log->logger_mask_data = NULL;
  log->logger_count_mask = 0;

  for (int i = 0; i < swift_type_count; i++) {
    logger_history_free(&log->history_new[i]);
    logger_history_free(&log->history_removed[i]);
  }
}

/**
 * @brief Write a file header to a logger file
 *
 * @param log The #logger_writer
 *
 */
void logger_write_file_header(struct logger_writer *log) {

  /* get required variables. */
  struct dump *dump = &log->dump;

  uint64_t file_offset = dump->file_offset;

  if (file_offset != 0)
    error(
        "The logger is not empty."
        "This function should be called before writing anything in the logger");

  /* Write format information. */
  logger_write_data(dump, &file_offset, logger_format_size,
                    &logger_file_format);

  /* Write the major version number. */
  int major = logger_major_version;
  logger_write_data(dump, &file_offset, sizeof(int), &major);

  /* Write the minor version number. */
  int minor = logger_minor_version;
  logger_write_data(dump, &file_offset, sizeof(int), &minor);

  /* write offset direction. */
  const int reversed = 0;
  logger_write_data(dump, &file_offset, sizeof(int), &reversed);

  /* placeholder to write the offset of the first log here. */
  char *skip_header = dump_get(dump, logger_offset_size, &file_offset);

  /* write number of bytes used for names. */
  const unsigned int label_size = logger_label_size;
  logger_write_data(dump, &file_offset, sizeof(unsigned int), &label_size);

  /* placeholder to write the number of unique masks. */
  char *skip_unique_masks = dump_get(dump, sizeof(unsigned int), &file_offset);

  /* write masks. */
  // loop over all mask type.
  unsigned int unique_mask = 0;
  for (int i = 0; i < log->logger_count_mask; i++) {
    /* Check if the mask was not already written */
    int is_written = 0;
    for (int j = 0; j < i; j++) {
      if (log->logger_mask_data[i].mask == log->logger_mask_data[j].mask) {
        is_written = 1;
        break;
      }
    }

    if (is_written) {
      continue;
    }

    unique_mask += 1;

    // mask name.
    logger_write_data(dump, &file_offset, logger_label_size,
                      &log->logger_mask_data[i].name);

    // mask size.
    logger_write_data(dump, &file_offset, sizeof(unsigned int),
                      &log->logger_mask_data[i].size);
  }
  memcpy(skip_unique_masks, &unique_mask, sizeof(unsigned int));

  /* last step: write first offset. */
  memcpy(skip_header, &file_offset, logger_offset_size);
}

/**
 * @brief read record header
 *
 * @param buff The reading buffer
 * @param mask The mask to read
 * @param offset (return) the offset pointed by this record (absolute)
 * @param cur_offset The current record offset
 *
 * @return Number of bytes read
 */
__attribute__((always_inline)) INLINE static int logger_read_record_header(
    const char *buff, unsigned int *mask, size_t *offset, size_t cur_offset) {
  memcpy(mask, buff, logger_mask_size);
  buff += logger_mask_size;

  *offset = 0;
  memcpy(offset, buff, logger_offset_size);
  *offset = cur_offset - *offset;

  return logger_mask_size + logger_offset_size;
}

/**
 * @brief Read a logger message and store the data in a #part.
 *
 * @param p The #part in which to store the values.
 * @param offset Pointer to the offset of the logger message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded logger message.
 *
 * @return The mask containing the values read.
 */
int logger_read_part(const struct logger_writer *log, struct part *p,
                     size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_record_header(buff, &mask, offset, cur_offset);

  for (int i = 0; i < log->logger_count_mask; i++) {
    if ((mask & log->logger_mask_data[i].mask) &&
        (log->logger_mask_data[i].type == mask_type_gas)) {

      const char *name = log->logger_mask_data[i].name;
      if (strcmp("Coordinates", name) == 0) {
        memcpy(p->x, buff, 3 * sizeof(double));
        buff += 3 * sizeof(double);
      } else if (strcmp("Velocities", name) == 0) {
        memcpy(p->v, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("Accelerations", name) == 0) {
        memcpy(p->a_hydro, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("SmoothingLengths", name) == 0) {
        memcpy(&p->h, buff, sizeof(float));
        buff += sizeof(float);
      }
#if defined(GADGET2_SPH)
      else if (strcmp("Entropies", name) == 0) {
        memcpy(&p->entropy, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("Masses", name) == 0) {
        memcpy(&p->mass, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("Densities", name) == 0) {
        memcpy(&p->rho, buff, sizeof(float));
        buff += sizeof(float);
      } else if (strcmp("ParticleIDs", name) == 0) {
        memcpy(&p->id, buff, sizeof(long long));
        buff += sizeof(long long);
      }
#endif
      else {
        error("Field '%s' not found", name);
      }
    }
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a logger message and store the data in a #gpart.
 *
 * @param p The #gpart in which to store the values.
 * @param offset Pointer to the offset of the logger message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded logger message.
 *
 * @return The mask containing the values read.
 */
int logger_read_gpart(const struct logger_writer *log, struct gpart *p,
                      size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_record_header(buff, &mask, offset, cur_offset);

  for (int i = 0; i < log->logger_count_mask; i++) {
    if ((mask & log->logger_mask_data[i].mask) &&
        (log->logger_mask_data[i].type == mask_type_dark_matter)) {

      const char *name = log->logger_mask_data[i].name;
      if (strcmp("Coordinates", name) == 0) {
        memcpy(p->x, buff, 3 * sizeof(double));
        buff += 3 * sizeof(double);
      } else if (strcmp("Velocities", name) == 0) {
        memcpy(p->v_full, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("Accelerations", name) == 0) {
        memcpy(p->a_grav, buff, 3 * sizeof(float));
        buff += 3 * sizeof(float);
      } else if (strcmp("ParticleIDs", name) == 0) {
        memcpy(&p->id_or_neg_offset, buff, sizeof(long long));
        buff += sizeof(long long);
      } else if (strcmp("Masses", name) == 0) {
        memcpy(&p->mass, buff, sizeof(float));
        buff += sizeof(float);
      } else {
        error("Field '%s' not found", name);
      }
    }
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a logger message for a timestamp.
 *
 * @param log The #logger_writer.
 * @param t The timestamp in which to store the value.
 * @param time The time in which to store the value.
 * @param offset Pointer to the offset of the logger message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded logger message.
 *
 * @return The mask containing the values read.
 */
int logger_read_timestamp(const struct logger_writer *log, integertime_t *t,
                          double *time, size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_record_header(buff, &mask, offset, cur_offset);

  /* We are only interested in timestamps. */
  if (!(mask & log->logger_mask_data[logger_index_timestamp].mask))
    error("Trying to read timestamp from a particle.");

  /* Make sure we don't have extra fields. */
  if (mask != log->logger_mask_data[logger_index_timestamp].mask)
    error("Timestamp message contains extra fields.");

  /* Copy the timestamp value from the buffer. */
  memcpy(t, buff, sizeof(integertime_t));
  buff += sizeof(integertime_t);

  /* Copy the timestamp value from the buffer. */
  memcpy(time, buff, sizeof(double));

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Write a swift_params struct to the given FILE as a stream of bytes.
 *
 * @param log the struct
 * @param stream the file stream
 */
void logger_struct_dump(const struct logger_writer *log, FILE *stream) {
  restart_write_blocks((void *)log, sizeof(struct logger_writer), 1, stream,
                       "logger", "logger");

  /* Write the masks */
  restart_write_blocks((void *)log->logger_mask_data, sizeof(struct mask_data),
                       log->logger_count_mask, stream, "logger_masks",
                       "logger_masks");

  /* Dump the logger mpi history */
  for (int i = 0; i < swift_type_count; i++) {
    logger_history_dump(&log->history_new[i], stream);
    logger_history_dump(&log->history_removed[i], stream);
  }
}

/**
 * @brief Restore a logger struct from the given FILE as a stream of
 * bytes.
 *
 * @param logger the struct
 * @param stream the file stream
 */
void logger_struct_restore(struct logger_writer *log, FILE *stream) {
  /* Read the block */
  restart_read_blocks((void *)log, sizeof(struct logger_writer), 1, stream,
                      NULL, "logger");

  /* Read the masks */
  const struct mask_data *old_logger_mask_data = log->logger_mask_data;
  log->logger_mask_data = (struct mask_data *)malloc(sizeof(struct mask_data) *
                                                     log->logger_count_mask);

  restart_read_blocks((void *)log->logger_mask_data, sizeof(struct mask_data),
                      log->logger_count_mask, stream, NULL, "logger_masks");

  /* Restore the pointers */
  log->mask_data_pointers.hydro =
      log->logger_mask_data +
      (log->mask_data_pointers.hydro - old_logger_mask_data);
  log->mask_data_pointers.gravity =
      log->logger_mask_data +
      (log->mask_data_pointers.gravity - old_logger_mask_data);
  log->mask_data_pointers.stars =
      log->logger_mask_data +
      (log->mask_data_pointers.stars - old_logger_mask_data);

  /* Restart the dump file. */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  logger_get_dump_name(log, logger_name_file);

  dump_restart(&log->dump, logger_name_file);

  /* Restore the logger mpi history */
  for (int i = 0; i < swift_type_count; i++) {
    logger_history_restore(&log->history_new[i], stream);
    logger_history_restore(&log->history_removed[i], stream);
  }
}

#endif /* WITH_LOGGER */

#endif /* HAVE_POSIX_FALLOCATE */
