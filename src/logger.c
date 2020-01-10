/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* This object's header. */
#include "logger.h"

/* Local headers. */
#include "atomic.h"
#include "dump.h"
#include "engine.h"
#include "error.h"
#include "part.h"
#include "units.h"

/*
 * Thoses are definitions from the format and therefore should not be changed!
 */
/* Number of bytes for a mask. */
// TODO change this to number of bits
#define logger_mask_size 2

/* Number of bits for chunk header. */
#define logger_header_bytes 8

/* Number bytes for an offset. */
#define logger_offset_size logger_header_bytes - logger_mask_size

/* Number of bytes for the file format information. */
#define logger_format_size 20

/* Number of bytes for the labels in the header. */
#define logger_label_size 20

char logger_file_format[logger_format_size] = "SWIFT_LOGGER";

const struct mask_data logger_mask_data[logger_count_mask] = {
    /* Particle's position. */
    {3 * sizeof(double), 1 << logger_x, "positions"},
    /* Particle's velocity. */
    {3 * sizeof(float), 1 << logger_v, "velocities"},
    /* Particle's acceleration. */
    {3 * sizeof(float), 1 << logger_a, "accelerations"},
    /* Particle's entropy. */
    {sizeof(float), 1 << logger_u, "entropy"},
    /* Particle's smoothing length. */
    {sizeof(float), 1 << logger_h, "smoothing length"},
    /* Particle's density. */
    {sizeof(float), 1 << logger_rho, "density"},
    /* Particle's constants: mass (float) and ID (long long). */
    {sizeof(float) + sizeof(long long), 1 << logger_consts, "consts"},
    /* Flag for special cases (e.g. change of MPI rank, star formation, ...) */
    {sizeof(int), 1 << logger_special_flags, "special flags"},
    /* Simulation time stamp: integertime and double time (e.g. scale
       factor or time). */
    {sizeof(integertime_t) + sizeof(double), 1 << logger_timestamp,
     "timestamp"},
};

/**
 * @brief Write the header of a chunk (offset + mask).
 *
 * This is maybe broken for big(?) endian.
 *
 * @param buff The writing buffer
 * @param mask The mask to write
 * @param offset The old offset
 * @param offset_new The new offset
 *
 * @return updated buff
 */
char *logger_write_chunk_header(char *buff, const unsigned int *mask,
                                const size_t *offset, const size_t offset_new) {
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

  /* Update offset to end of chunk. */
  *offset += size;
}

/**
 * @brief Compute the size of a message given its mask.
 *
 * @param mask The mask that will be used to dump a #part or #gpart.
 *
 * @return The size of the logger message in bytes.
 */
int logger_compute_chunk_size(unsigned int mask) {

  /* Start with 8 bytes for the header. */
  int size = logger_mask_size + logger_offset_size;

  /* Is this a particle or a timestep? */
  if (mask & logger_mask_data[logger_timestamp].mask) {

    /* The timestamp should not contain any other bits. */
    if (mask != logger_mask_data[logger_timestamp].mask)
      error("Timestamps should not include any other data.");

    /* A timestamp consists of an unsigned long long int. */
    size += logger_mask_data[logger_timestamp].size;

  } else {

    for (int i = 0; i < logger_count_mask; i++) {
      if (mask & logger_mask_data[i].mask) {
        size += logger_mask_data[i].size;
      }
    }
  }

  return size;
}

/**
 * @brief log all particles in the engine.
 *
 * @param log The #logger_writer
 * @param e The #engine
 */
void logger_log_all(struct logger_writer *log, const struct engine *e) {

  /* Ensure that enough space is available. */
  logger_ensure_size(log, e->s->nr_parts, e->s->nr_gparts, e->s->nr_sparts);

  /* some constants. */
  const struct space *s = e->s;
  const unsigned int mask_hydro =
      logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
      logger_mask_data[logger_a].mask | logger_mask_data[logger_u].mask |
      logger_mask_data[logger_h].mask | logger_mask_data[logger_rho].mask |
      logger_mask_data[logger_consts].mask;

  /* loop over all parts. */
  for (size_t i = 0; i < s->nr_parts; i++) {
    logger_log_part(log, &s->parts[i], mask_hydro,
                    &s->xparts[i].logger_data.last_offset,
                    /* Special flags */ 0);
    s->xparts[i].logger_data.steps_since_last_output = 0;
  }

  const unsigned int mask_grav =
      logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
      logger_mask_data[logger_a].mask | logger_mask_data[logger_consts].mask;

  /* loop over all gparts */
  for (size_t i = 0; i < s->nr_gparts; i++) {
    /* Log only the dark matter */
    if (s->gparts[i].type != swift_type_dark_matter) continue;

    logger_log_gpart(log, &s->gparts[i], mask_grav,
                     &s->gparts[i].logger_data.last_offset,
                     /* Special flags */ 0);
    s->gparts[i].logger_data.steps_since_last_output = 0;
  }

  const unsigned int mask_stars = logger_mask_data[logger_x].mask |
                                  logger_mask_data[logger_v].mask |
                                  logger_mask_data[logger_consts].mask;

  /* loop over all sparts */
  for (size_t i = 0; i < s->nr_sparts; i++) {
    logger_log_spart(log, &s->sparts[i], mask_stars,
                     &s->sparts[i].logger_data.last_offset,
                     /* Special flags */ 0);
    s->sparts[i].logger_data.steps_since_last_output = 0;
  }

  if (e->total_nr_bparts > 0) error("Not implemented");
}

/**
 * @brief Dump a #part to the log.
 *
 * @param log The #logger_writer
 * @param p The #part to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle;
 * @param special_flags The value of the special flag.
 * (return) offset of this log.
 */
void logger_log_part(struct logger_writer *log, const struct part *p,
                     unsigned int mask, size_t *offset,
                     const int special_flags) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_data[logger_timestamp].mask)
    error("You should not log particles as timestamps.");

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size, &offset_new);

  /* Write the header. */
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Particle position as three doubles. */
  if (mask & logger_mask_data[logger_x].mask) {
    memcpy(buff, p->x, logger_mask_data[logger_x].size);
    buff += logger_mask_data[logger_x].size;
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_data[logger_v].mask) {
    memcpy(buff, p->v, logger_mask_data[logger_v].size);
    buff += logger_mask_data[logger_v].size;
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_data[logger_a].mask) {
    memcpy(buff, p->a_hydro, logger_mask_data[logger_a].size);
    buff += logger_mask_data[logger_a].size;
  }

#if defined(GADGET2_SPH)

  /* Particle internal energy as a single float. */
  if (mask & logger_mask_data[logger_u].mask) {
    memcpy(buff, &p->entropy, logger_mask_data[logger_u].size);
    buff += logger_mask_data[logger_u].size;
  }

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_data[logger_h].mask) {
    memcpy(buff, &p->h, logger_mask_data[logger_h].size);
    buff += logger_mask_data[logger_h].size;
  }

  /* Particle density as a single float. */
  if (mask & logger_mask_data[logger_rho].mask) {
    memcpy(buff, &p->rho, logger_mask_data[logger_rho].size);
    buff += logger_mask_data[logger_rho].size;
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_data[logger_consts].mask) {
    // TODO make it dependent of logger_mask_data
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
    const int64_t id = p->id;
    memcpy(buff, &id, sizeof(int64_t));
    buff += sizeof(int64_t);
  }

#endif

  /* Special flags */
  if (mask & logger_mask_data[logger_special_flags].mask) {
    memcpy(buff, &special_flags, logger_mask_data[logger_special_flags].size);
    buff += logger_mask_data[logger_special_flags].size;
  }

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief Dump a #spart to the log.
 *
 * @param log The #logger_writer
 * @param sp The #spart to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle;
 * @param special_flags The value of the special flag.
 * (return) offset of this log.
 */
void logger_log_spart(struct logger_writer *log, const struct spart *sp,
                      unsigned int mask, size_t *offset,
                      const int special_flags) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_data[logger_timestamp].mask)
    error("You should not log particles as timestamps.");

  /* Make sure we're not looging fields not supported by gparts. */
  if (mask &
      (logger_mask_data[logger_u].mask | logger_mask_data[logger_rho].mask |
       logger_mask_data[logger_a].mask))
    error("Can't log SPH quantities for sparts.");

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size, &offset_new);

  /* Write the header. */
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Particle position as three doubles. */
  if (mask & logger_mask_data[logger_x].mask) {
    memcpy(buff, sp->x, logger_mask_data[logger_x].size);
    buff += logger_mask_data[logger_x].size;
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_data[logger_v].mask) {
    memcpy(buff, sp->v, logger_mask_data[logger_v].size);
    buff += logger_mask_data[logger_v].size;
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_data[logger_consts].mask) {
    // TODO make it dependent of logger_mask_data
    memcpy(buff, &sp->mass, sizeof(float));
    buff += sizeof(float);
    const int64_t id = sp->id;
    memcpy(buff, &id, sizeof(int64_t));
    buff += sizeof(int64_t);
  }

  /* Special flags */
  if (mask & logger_mask_data[logger_special_flags].mask) {
    memcpy(buff, &special_flags, logger_mask_data[logger_special_flags].size);
    buff += logger_mask_data[logger_special_flags].size;
  }

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief Dump a #gpart to the log.
 *
 * @param log The #logger_writer
 * @param p The #gpart to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle;
 * @param special_flags The value of the special flags.
 * (return) offset of this log.
 */
void logger_log_gpart(struct logger_writer *log, const struct gpart *p,
                      unsigned int mask, size_t *offset,
                      const int special_flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->id_or_neg_offset < 0) {
    error("Cannot log a gpart attached to another particle");
  }
#endif

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_data[logger_timestamp].mask)
    error("You should not log particles as timestamps.");

  /* Make sure we're not looging fields not supported by gparts. */
  if (mask &
      (logger_mask_data[logger_u].mask | logger_mask_data[logger_rho].mask))
    error("Can't log SPH quantities for gparts.");

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(&log->dump, size, &offset_new);

  /* Write the header. */
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Particle position as three doubles. */
  if (mask & logger_mask_data[logger_x].mask) {
    memcpy(buff, p->x, logger_mask_data[logger_x].size);
    buff += logger_mask_data[logger_x].size;
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_data[logger_v].mask) {
    memcpy(buff, p->v_full, logger_mask_data[logger_v].size);
    buff += logger_mask_data[logger_v].size;
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_data[logger_a].mask) {
    memcpy(buff, p->a_grav, logger_mask_data[logger_a].size);
    buff += logger_mask_data[logger_a].size;
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_data[logger_consts].mask) {
    // TODO make it dependent of logger_mask_data.
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
    const int64_t id = p->id_or_neg_offset;
    memcpy(buff, &id, sizeof(int64_t));
    buff += sizeof(int64_t);
  }

  /* Special flags */
  if (mask & logger_mask_data[logger_special_flags].mask) {
    memcpy(buff, &special_flags, logger_mask_data[logger_special_flags].size);
    buff += logger_mask_data[logger_special_flags].size;
  }

  /* Update the log message offset. */
  *offset = offset_new;
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
      logger_compute_chunk_size(logger_mask_data[logger_timestamp].mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  unsigned int mask = logger_mask_data[logger_timestamp].mask;
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Store the timestamp. */
  // TODO make it dependent of logger_mask_data.
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
  limit *= log->max_chunk_size;

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
 * @brief intialize the logger structure
 *
 * @param log The #logger_writer
 * @param params The #swift_params
 */
void logger_init(struct logger_writer *log, struct swift_params *params) {
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

  /* set initial value of parameters. */
  log->timestamp_offset = 0;
  log->index.dump_size_last_output = 0;

  /* generate dump filename. */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  logger_get_dump_name(log, logger_name_file);

  /* Compute max size for a particle chunk. */
  int max_size = logger_offset_size + logger_mask_size;

  /* Loop over all fields except timestamp. */
  for (int i = 0; i < logger_count_mask - 1; i++) {
    max_size += logger_mask_data[i].size;
  }
  log->max_chunk_size = max_size;

  /* init dump. */
  dump_init(&log->dump, logger_name_file, buffer_size);
}

/**
 * @brief Close dump file and desallocate memory
 *
 * @param log The #logger_writer
 */
void logger_free(struct logger_writer *log) { dump_close(&log->dump); }

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

  /* write number of masks. */
  const unsigned int count_mask = logger_count_mask;
  logger_write_data(dump, &file_offset, sizeof(unsigned int), &count_mask);

  /* write masks. */
  // loop over all mask type.
  for (int i = 0; i < logger_count_mask; i++) {
    // mask name.
    logger_write_data(dump, &file_offset, logger_label_size,
                      &logger_mask_data[i].name);

    // mask size.
    logger_write_data(dump, &file_offset, sizeof(unsigned int),
                      &logger_mask_data[i].size);
  }

  /* last step: write first offset. */
  memcpy(skip_header, &file_offset, logger_offset_size);
}

/**
 * @brief read chunk header
 *
 * @param buff The reading buffer
 * @param mask The mask to read
 * @param offset (return) the offset pointed by this chunk (absolute)
 * @param cur_offset The current chunk offset
 *
 * @return Number of bytes read
 */
__attribute__((always_inline)) INLINE static int logger_read_chunk_header(
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
int logger_read_part(struct part *p, size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_chunk_header(buff, &mask, offset, cur_offset);

  /* We are only interested in particle data. */
  if (mask & logger_mask_data[logger_timestamp].mask)
    error("Trying to read timestamp as particle.");

  /* Particle position as three doubles. */
  if (mask & logger_mask_data[logger_x].mask) {
    memcpy(p->x, buff, logger_mask_data[logger_x].size);
    buff += logger_mask_data[logger_x].size;
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_data[logger_v].mask) {
    memcpy(p->v, buff, logger_mask_data[logger_v].size);
    buff += logger_mask_data[logger_v].size;
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_data[logger_a].mask) {
    memcpy(p->a_hydro, buff, logger_mask_data[logger_a].size);
    buff += logger_mask_data[logger_a].size;
  }

#if defined(GADGET2_SPH)

  /* Particle internal energy as a single float. */
  if (mask & logger_mask_data[logger_u].mask) {
    memcpy(&p->entropy, buff, logger_mask_data[logger_u].size);
    buff += logger_mask_data[logger_u].size;
  }

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_data[logger_h].mask) {
    memcpy(&p->h, buff, logger_mask_data[logger_h].size);
    buff += logger_mask_data[logger_h].size;
  }

  /* Particle density as a single float. */
  if (mask & logger_mask_data[logger_rho].mask) {
    memcpy(&p->rho, buff, logger_mask_data[logger_rho].size);
    buff += logger_mask_data[logger_rho].size;
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_data[logger_rho].mask) {
    // TODO make it dependent of logger_mask_data.
    memcpy(&p->mass, buff, sizeof(float));
    buff += sizeof(float);
    int64_t id = 0;
    memcpy(&id, buff, sizeof(int64_t));
    p->id = id;
    buff += sizeof(int64_t);
  }

#endif

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
int logger_read_gpart(struct gpart *p, size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_chunk_header(buff, &mask, offset, cur_offset);

  /* We are only interested in particle data. */
  if (mask & logger_mask_data[logger_timestamp].mask)
    error("Trying to read timestamp as particle.");

  /* We can't store all part fields in a gpart. */
  if (mask &
      (logger_mask_data[logger_u].mask | logger_mask_data[logger_rho].mask))
    error("Trying to read SPH quantities into a gpart.");

  /* Particle position as three doubles. */
  if (mask & logger_mask_data[logger_x].mask) {
    memcpy(p->x, buff, logger_mask_data[logger_x].size);
    buff += logger_mask_data[logger_x].size;
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_data[logger_v].mask) {
    memcpy(p->v_full, buff, logger_mask_data[logger_v].size);
    buff += logger_mask_data[logger_v].size;
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_data[logger_a].mask) {
    memcpy(p->a_grav, buff, logger_mask_data[logger_a].size);
    buff += logger_mask_data[logger_a].size;
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_data[logger_rho].mask) {
    // TODO make it dependent of logger_mask_data
    memcpy(&p->mass, buff, sizeof(float));
    buff += sizeof(float);
    int64_t id = p->id_or_neg_offset;
    memcpy(&id, buff, sizeof(int64_t));
    buff += sizeof(int64_t);
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a logger message for a timestamp.
 *
 * @param t The timestamp in which to store the value.
 * @param time The time in which to store the value.
 * @param offset Pointer to the offset of the logger message in the buffer,
 *        will be overwritten with the offset of the previous message.
 * @param buff Pointer to the start of an encoded logger message.
 *
 * @return The mask containing the values read.
 */
int logger_read_timestamp(unsigned long long int *t, double *time,
                          size_t *offset, const char *buff) {

  /* Jump to the offset. */
  buff = &buff[*offset];

  /* Start by reading the logger mask for this entry. */
  const size_t cur_offset = *offset;
  unsigned int mask = 0;
  buff += logger_read_chunk_header(buff, &mask, offset, cur_offset);

  /* We are only interested in timestamps. */
  if (!(mask & logger_mask_data[logger_timestamp].mask))
    error("Trying to read timestamp from a particle.");

  /* Make sure we don't have extra fields. */
  if (mask != logger_mask_data[logger_timestamp].mask)
    error("Timestamp message contains extra fields.");

  /* Copy the timestamp value from the buffer. */
  // TODO make it dependent of logger_mask_data.
  memcpy(t, buff, sizeof(unsigned long long int));
  buff += sizeof(unsigned long long int);

  /* Copy the timestamp value from the buffer. */
  memcpy(time, buff, sizeof(double));

  /* Finally, return the mask of the values we just read. */
  return mask;
}


#ifdef WITH_MPI
/**
 * @brief Log all the particles leaving the current rank.
 *
 * @param log The #logger_writer.
 * @param nr_nodes Number of nodes used in the simulation.
 * @param sneding Are we sending the particles (or receiving)?
 * @param parts The list of #part.
 * @param nr_parts The number of parts.
 * @param count The number of parts in each ranks.
 * @param gparts The list of #gpart.
 * @param nr_gparts The number of gparts.
 * @param gcount The number of gparts in each ranks.
 * @param sparts The list of #spart.
 * @param nr_sparts The number of sparts.
 * @param s_counts The number of sparts in each ranks.
 * @param bparts The list of #bpart.
 * @param nr_bparts The number of bparts.
 * @param b_counts The number of bparts in each ranks.
 *
 */
void logger_log_repartition(
    struct logger_writer *log, int nr_nodes, int sending, struct part *parts,
    struct xpart *xparts, size_t nr_parts, int *counts,
    struct gpart *gparts, size_t nr_gparts, int *g_counts,
    struct spart *sparts, size_t nr_sparts, int *s_counts,
    struct bpart *bparts, size_t nr_bparts, int *b_counts) {

  size_t part_offset = 0;
  size_t spart_offset = 0;
  size_t gpart_offset = 0;
  size_t bpart_offset = 0;

  for(int i = 0; i < nr_nodes; i++) {
    const size_t c_ind = sending ? engine_rank * nr_nodes + i:
      i * nr_nodes + engine_rank;

    /* No need to log the local particles. */
    if (i == engine_rank) {
      part_offset += counts[c_ind];
      spart_offset += s_counts[c_ind];
      gpart_offset += g_counts[c_ind];
      bpart_offset += b_counts[c_ind];
      continue;
    }

    const int flag = logger_generate_flag(
      logger_flag_mpi | logger_flag_delete, i);

    const unsigned int mask_hydro =
      logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
      logger_mask_data[logger_a].mask | logger_mask_data[logger_u].mask |
      logger_mask_data[logger_h].mask | logger_mask_data[logger_rho].mask |
      logger_mask_data[logger_consts].mask |
      logger_mask_data[logger_special_flags].mask;

    /* Log the hydro parts. */
    for(int j = 0; j < counts[c_ind]; j++) {
      size_t ind = part_offset + j;
      logger_log_part(log, &parts[ind], mask_hydro,
                      &xparts[ind].logger_data.last_offset,
                      flag);
      xparts[ind].logger_data.steps_since_last_output = 0;
    }

    const unsigned int mask_stars = logger_mask_data[logger_x].mask |
      logger_mask_data[logger_v].mask |
      logger_mask_data[logger_consts].mask |
      logger_mask_data[logger_special_flags].mask;

    /* Log the stellar parts. */
    for(int j = 0; j < s_counts[c_ind]; j++) {
      size_t ind = spart_offset + j;
      logger_log_spart(log, &sparts[ind], mask_stars,
                       &sparts[ind].logger_data.last_offset,
                       flag);
      sparts[ind].logger_data.steps_since_last_output = 0;
    }

    const unsigned int mask_grav =
      logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
      logger_mask_data[logger_a].mask | logger_mask_data[logger_consts].mask |
      logger_mask_data[logger_special_flags].mask;

    /* Log the gparts */
    for(int j = 0; j < g_counts[c_ind]; j++) {
      size_t ind = gpart_offset + j;
      /* Log only the dark matter */
      if (gparts[ind].type != swift_type_dark_matter) continue;

      logger_log_gpart(log, &gparts[ind], mask_grav,
                       &gparts[ind].logger_data.last_offset,
                       flag);
      gparts[ind].logger_data.steps_since_last_output = 0;
    }

    /* Log the bparts */
    if (b_counts[c_ind] > 0) {
      error("TODO");
    }

    /* Update the counters */
    part_offset += counts[c_ind];
    spart_offset += s_counts[c_ind];
    gpart_offset += g_counts[c_ind];
    bpart_offset += b_counts[c_ind];
  }
}

/**
 * @brief Log all the particles arriving in the current rank.
 *
 * @param parts The list of #part.
 * @param nr_parts The number of parts.
 * @param count The number of parts in each ranks.
 * @param gparts The list of #gpart.
 * @param nr_gparts The number of gparts.
 * @param gcount The number of gparts in each ranks.
 * @param sparts The list of #spart.
 * @param nr_sparts The number of sparts.
 * @param s_counts The number of sparts in each ranks.
 * @param bparts The list of #bpart.
 * @param nr_bparts The number of bparts.
 * @param b_counts The number of bparts in each ranks.
 *
 */
void logger_log_recv_strays(
    struct logger_writer *log,
    struct part *parts, struct xpart *xparts, size_t nr_parts,
    struct gpart *gparts, size_t nr_gparts,
    struct spart *sparts, size_t nr_sparts,
    struct bpart *bparts, size_t nr_bparts,
    int node_id) {

  const int flag = logger_generate_flag(logger_flag_mpi | logger_flag_create,
                                        node_id);

  /* Log the gas particles */
  const unsigned int mask_hydro =
    logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
    logger_mask_data[logger_a].mask | logger_mask_data[logger_u].mask |
    logger_mask_data[logger_h].mask | logger_mask_data[logger_rho].mask |
    logger_mask_data[logger_consts].mask |
    logger_mask_data[logger_special_flags].mask;

  for(size_t i = 0; i < nr_parts; i++) {
    logger_log_part(log, &parts[i], mask_hydro,
                    &xparts[i].logger_data.last_offset,
                    flag);
    xparts[i].logger_data.steps_since_last_output = 0;
  }

  /* Log the stellar particles */
  const unsigned int mask_stars = logger_mask_data[logger_x].mask |
    logger_mask_data[logger_v].mask |
    logger_mask_data[logger_consts].mask |
    logger_mask_data[logger_special_flags].mask;
  for(size_t i = 0; i < nr_sparts; i++) {
    logger_log_spart(log, &sparts[i], mask_stars,
                     &sparts[i].logger_data.last_offset,
                     /* Special flags */ 0);
    sparts[i].logger_data.steps_since_last_output = 0;
  }


  /* Log the gparts */
  const unsigned int mask_grav =
    logger_mask_data[logger_x].mask | logger_mask_data[logger_v].mask |
    logger_mask_data[logger_a].mask | logger_mask_data[logger_consts].mask |
    logger_mask_data[logger_special_flags].mask;
  for(size_t i = 0; i < nr_gparts; i++) {
    /* Log only the dark matter */
    if (gparts[i].type != swift_type_dark_matter) continue;

    logger_log_gpart(log, &gparts[i], mask_grav,
                     &gparts[i].logger_data.last_offset,
                     /* Special flags */ 0);
    gparts[i].logger_data.steps_since_last_output = 0;
  }

  /* Log the bparts */
  if (nr_bparts > 0) {
    error("TODO");
  }
}

#endif

/**
 * @brief Write a swift_params struct to the given FILE as a stream of bytes.
 *
 * @param log the struct
 * @param stream the file stream
 */
void logger_struct_dump(const struct logger_writer *log, FILE *stream) {
  restart_write_blocks((void *)log, sizeof(struct logger_writer), 1, stream,
                       "logger", "logger");
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

  /* generate dump filename */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  logger_get_dump_name(log, logger_name_file);

  dump_restart(&log->dump, logger_name_file);
}

#endif /* WITH_LOGGER */

#endif /* HAVE_POSIX_FALLOCATE */
