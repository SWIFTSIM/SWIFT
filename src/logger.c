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
/* number of bytes for a mask */
#define logger_mask_size 1

/* number bytes for an offset */
#define logger_offset_size 7

/* number of bytes for the version information */
#define logger_version_size 20

/* number of bytes for the labels in the header */
#define logger_label_size 20

/* number of bytes for the number in the header */
#define logger_number_size 4

char logger_version[logger_version_size] = "0.1";

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
  /* write mask */
  memcpy(buff, mask, logger_mask_size);
  buff += logger_mask_size;

  /* write offset */
  size_t diff_offset = offset_new - *offset;
  memcpy(buff, &diff_offset, logger_offset_size);
  buff += logger_offset_size;

  return buff;
}

/**
 * @brief Write to the dump
 *
 * @param d #dump file
 * @param offset (return) offset of the data
 * @param size number of bytes to write
 * @param p pointer to the data
 */
void logger_write_data(struct dump *d, size_t *offset, size_t size,
                       const void *p) {
  /* get buffer */
  char *buff = dump_get(d, size, offset);

  /* write data to the buffer */
  memcpy(buff, p, size);

  /* Update offset to end of chunk */
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
  if (mask & logger_mask_timestamp) {

    /* The timestamp should not contain any other bits. */
    if (mask != logger_mask_timestamp)
      error("Timestamps should not include any other data.");

    /* A timestamp consists of an unsigned long long int. */
    size += sizeof(unsigned long long int);
    size += sizeof(double);

  } else {

    /* Particle position as three doubles. */
    if (mask & logger_mask_x) size += 3 * sizeof(double);

    /* Particle velocity as three floats. */
    if (mask & logger_mask_v) size += 3 * sizeof(float);

    /* Particle accelleration as three floats. */
    if (mask & logger_mask_a) size += 3 * sizeof(float);

    /* Particle internal energy as a single float. */
    if (mask & logger_mask_u) size += sizeof(float);

    /* Particle smoothing length as a single float. */
    if (mask & logger_mask_h) size += sizeof(float);

    /* Particle density as a single float. */
    if (mask & logger_mask_rho) size += sizeof(float);

    /* Particle constants, which is a bit more complicated. */
    if (mask & logger_mask_rho) {
      size += sizeof(float) +     // mass
              sizeof(long long);  // id
    }
  }

  return size;
}

/**
 * @brief log all particles in the engine.
 *
 * @param log The #logger
 * @param e The #engine
 */
void logger_log_all(struct logger *log, const struct engine *e) {

  /* Ensure that enough space is available */
  logger_ensure_size(log, e->total_nr_parts, e->total_nr_gparts, 0);
#ifdef SWIFT_DEBUG_CHECKS
  message("Need to implement stars");
#endif

  /* some constants */
  const struct space *s = e->s;
  const unsigned int mask = logger_mask_x | logger_mask_v | logger_mask_a |
                            logger_mask_u | logger_mask_h | logger_mask_rho |
                            logger_mask_consts;

  /* loop over all parts */
  for (long long i = 0; i < e->total_nr_parts; i++) {
    logger_log_part(log, &s->parts[i], mask,
                    &s->xparts[i].logger_data.last_offset);
    s->xparts[i].logger_data.steps_since_last_output = 0;
  }

  /* loop over all gparts */
  if (e->total_nr_gparts > 0) error("Not implemented");

  /* loop over all sparts */
  // TODO
}

/**
 * @brief Dump a #part to the log.
 *
 * @param log The #logger
 * @param p The #part to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle;
 * (return) offset of this log.
 */
void logger_log_part(struct logger *log, const struct part *p,
                     unsigned int mask, size_t *offset) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_timestamp)
    error("You should not log particles as timestamps.");

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(log->dump, size, &offset_new);

  /* Write the header. */
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Particle position as three doubles. */
  if (mask & logger_mask_x) {
    memcpy(buff, p->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_v) {
    memcpy(buff, p->v, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_a) {
    memcpy(buff, p->a_hydro, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

#if defined(GADGET2_SPH)

  /* Particle internal energy as a single float. */
  if (mask & logger_mask_u) {
    memcpy(buff, &p->entropy, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_h) {
    memcpy(buff, &p->h, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle density as a single float. */
  if (mask & logger_mask_rho) {
    memcpy(buff, &p->rho, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_consts) {
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
    memcpy(buff, &p->id, sizeof(long long));
    buff += sizeof(long long);
  }

#endif

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief Dump a #gpart to the log.
 *
 * @param log The #logger
 * @param p The #gpart to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle;
 * (return) offset of this log.
 */
void logger_log_gpart(struct logger *log, const struct gpart *p,
                      unsigned int mask, size_t *offset) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_timestamp)
    error("You should not log particles as timestamps.");

  /* Make sure we're not looging fields not supported by gparts. */
  if (mask & (logger_mask_u | logger_mask_rho))
    error("Can't log SPH quantities for gparts.");

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(log->dump, size, &offset_new);

  /* Write the header. */
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Particle position as three doubles. */
  if (mask & logger_mask_x) {
    memcpy(buff, p->x, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_v) {
    memcpy(buff, p->v_full, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_a) {
    memcpy(buff, p->a_grav, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_consts) {
    memcpy(buff, &p->mass, sizeof(float));
    buff += sizeof(float);
    memcpy(buff, &p->id_or_neg_offset, sizeof(long long));
    buff += sizeof(long long);
  }

  /* Update the log message offset. */
  *offset = offset_new;
}

/**
 * @brief write a timestamp
 *
 * @param log The #logger
 * @param timestamp time to write
 * @param time time or scale factor
 * @param offset Pointer to the offset of the previous log of this particle;
 * (return) offset of this log.
 */
void logger_log_timestamp(struct logger *log, integertime_t timestamp,
                          double time, size_t *offset) {
  struct dump *dump = log->dump;

  /* Start by computing the size of the message. */
  const int size = logger_compute_chunk_size(logger_mask_timestamp);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  unsigned int mask = logger_mask_timestamp;
  buff = logger_write_chunk_header(buff, &mask, offset, offset_new);

  /* Store the timestamp. */
  memcpy(buff, &timestamp, sizeof(integertime_t));
  buff += sizeof(integertime_t);

  /* Store the time */
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
 * @param log The #logger
 * @param total_nr_parts total number of part
 * @param total_nr_gparts total number of gpart
 * @param total_nr_sparts total number of spart
 */
void logger_ensure_size(struct logger *log, size_t total_nr_parts,
                        size_t total_nr_gparts, size_t total_nr_sparts) {

  /* count part memory */
  size_t limit = log->max_chunk_size;

  limit *= total_nr_parts;

  /* count gpart memory */
  if (total_nr_gparts > 0) error("Not implemented");

  /* count spart memory */
  if (total_nr_sparts > 0) error("Not implemented");

  /* ensure enough space in dump */
  dump_ensure(log->dump, limit, log->buffer_scale * limit);
}

/**
 * @brief Generate a list of the mask size.
 *
 * @param log The #logger.
 * @param output (output) the mask sizes (length = number of masks)
 */
void logger_get_list_mask_size(const struct logger *log, int *output) {
  /* Position */
  output[0] = 3 * sizeof(double);
  /* Velocity */
  output[1] = 3 * sizeof(float);
  /* Acceleration */
  output[2] = 3 * sizeof(float);
  /* Internal Energy */
  output[3] = sizeof(float);
  /* Smoothing Length */
  output[4] = sizeof(float);
  /* Density */
  output[5] = sizeof(float);
  /* Constants */
  output[6] = sizeof(float) + sizeof(long long);
  /* Time stamp */
  output[7] = sizeof(integertime_t) + sizeof(double);

}

/**
 * @brief Generate a list of the mask names.
 *
 * @param log The #logger.
 * @param output (output) the mask name (length = number of masks * logger_label_size)
 */
void logger_get_list_mask_name(const struct logger *log, char *output) {
  /* set the mask names */
  char tmp[logger_label_size];
  strcpy(tmp, "position");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "velocity");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "acceleration");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "entropy");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "cutoff radius");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "density");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "consts");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;

  strcpy(tmp, "timestamp");
  memcpy(output, &tmp, logger_label_size);
  output += logger_label_size;
}


/**
 * @brief Generate a list of the masks.
 *
 * @param log The #logger.
 * @param output (output) the masks (length = number of masks)
 */
void logger_get_list_mask(const struct logger *log, int *output) {
  /* Position */
  output[0] = logger_mask_x;
  /* Velocity */
  output[1] = logger_mask_v;
  /* Acceleration */
  output[2] = logger_mask_a;
  /* Internal Energy */
  output[3] = logger_mask_u;
  /* Smoothing Length */
  output[4] = logger_mask_h;
  /* Density */
  output[5] = logger_mask_rho;
  /* Constants */
  output[6] = logger_mask_consts;
  /* Time stamp */
  output[7] = logger_mask_timestamp;
}

/**
 * @brief Compute the maximal size of a chunk.
 *
 * @param log The #logger.
 */
int logger_compute_max_chunk_size(const struct logger *log) {

  int *output = malloc(sizeof(int) * log->number_masks);
  if (output == NULL)
    error("Unable to allocate memory");
  logger_get_list_mask_size(log, output);

  int max_size = logger_offset_size + logger_mask_size;
  /* Loop over all fields except timestamp */
  for (int i = 0; i < log->number_masks - 1; i++) {
    max_size += output[i];
  }

  free(output);
  return max_size;
}

/**
 * @brief intialize the logger structure
 *
 * @param log The #logger
 * @param params The #swift_params
 */
void logger_init(struct logger *log, struct swift_params *params) {
  /* read parameters */
  log->delta_step = parser_get_param_int(params, "Logger:delta_step");
  size_t buffer_size =
      parser_get_opt_param_float(params, "Logger:initial_buffer_size", 0.5) *
      1e9;
  log->buffer_scale =
      parser_get_opt_param_float(params, "Logger:buffer_scale", 10);
  parser_get_param_string(params, "Logger:basename", log->base_name);

  /* set initial value of parameters */
  log->timestamp_offset = 0;

  /* generate dump filename */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  strcpy(logger_name_file, log->base_name);
  strcat(logger_name_file, ".dump");

  /* Define some constants */
  log->number_masks = 8;
  log->max_chunk_size = logger_compute_max_chunk_size(log);
  
  /* init dump */
  log->dump = malloc(sizeof(struct dump));
  struct dump *dump_file = log->dump;

  dump_init(dump_file, logger_name_file, buffer_size);
}

/**
 * @brief Close dump file and desallocate memory
 *
 * @param log The #logger
 */
void logger_clean(struct logger *log) {
  dump_close(log->dump);
  free(log->dump);
  log->dump = NULL;
}

/**
 * @brief Write a file header to a logger file
 *
 * @param log The #logger
 * @param dump The #dump in which to log the particle data.
 *
 */
void logger_write_file_header(struct logger *log, const struct engine *e) {

  /* get required variables */
  struct dump *dump = log->dump;

  size_t file_offset = dump->file_offset;

  if (file_offset != 0)
    error(
        "The logger is not empty."
        "This function should be called before writing anything in the logger");

  /* Write version information */
  logger_write_data(dump, &file_offset, logger_version_size, &logger_version);

  /* write offset direction */
  const int reversed = 0;
  logger_write_data(dump, &file_offset, logger_number_size,
                    &reversed);

  /* placeholder to write the offset of the first log here */
  char *skip_header = dump_get(dump, logger_offset_size, &file_offset);

  /* write number of bytes used for names */
  const int label_size = logger_label_size;
  logger_write_data(dump, &file_offset, logger_number_size,
                    &label_size);

  /* write number of masks */
  logger_write_data(dump, &file_offset, logger_number_size,
                    &log->number_masks);

  /* Get masks sizes */
  int *mask_sizes = malloc(sizeof(int) * log->number_masks);
  if (mask_sizes == NULL)
    error("Unable to allocate memory for mask sizes");

  logger_get_list_mask_size(log, mask_sizes);

  /* Get masks */
  int *masks = malloc(sizeof(int) * log->number_masks);
  if (masks == NULL)
    error("Unable to allocate memory for masks");
  
  logger_get_list_mask(log, masks);

  /* Get masks names */
  char *mask_names = malloc(logger_label_size * log->number_masks);
  if (mask_names == NULL)
    error("Unable to allocate memory");

  logger_get_list_mask_name(log, mask_names);
  
  /* write masks */
  // loop over all mask type
  for (int i = 0; i < log->number_masks; i++) {
    // mask name
    size_t j = i * logger_label_size;
    logger_write_data(dump, &file_offset, logger_label_size,
                      &mask_names[j]);

    // mask size
    logger_write_data(dump, &file_offset, logger_number_size,
                      &mask_sizes[i]);
  }
  /* cleanup memory */
  free(mask_sizes);
  free(masks);
  free(mask_names);
  
  /* last step: write first offset */
  memcpy(skip_header, &file_offset, logger_offset_size);
}

/**
 * @brief read chunk header
 *
 * @param buff The reading buffer
 * @param mask The mask to read
 * @param offset (return) the offset pointed by this chunk (absolute)
 * @param offset_cur The current chunk offset
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
  if (mask & logger_mask_timestamp)
    error("Trying to read timestamp as particle.");

  /* Particle position as three doubles. */
  if (mask & logger_mask_x) {
    memcpy(p->x, buff, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_v) {
    memcpy(p->v, buff, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_a) {
    memcpy(p->a_hydro, buff, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

#if defined(GADGET2_SPH)

  /* Particle internal energy as a single float. */
  if (mask & logger_mask_u) {
    memcpy(&p->entropy, buff, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_h) {
    memcpy(&p->h, buff, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle density as a single float. */
  if (mask & logger_mask_rho) {
    memcpy(&p->rho, buff, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_rho) {
    memcpy(&p->mass, buff, sizeof(float));
    buff += sizeof(float);
    memcpy(&p->id, buff, sizeof(long long));
    buff += sizeof(long long);
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
  if (mask & logger_mask_timestamp)
    error("Trying to read timestamp as particle.");

  /* We can't store all part fields in a gpart. */
  if (mask & (logger_mask_u | logger_mask_rho))
    error("Trying to read SPH quantities into a gpart.");

  /* Particle position as three doubles. */
  if (mask & logger_mask_x) {
    memcpy(p->x, buff, 3 * sizeof(double));
    buff += 3 * sizeof(double);
  }

  /* Particle velocity as three floats. */
  if (mask & logger_mask_v) {
    memcpy(p->v_full, buff, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle accelleration as three floats. */
  if (mask & logger_mask_a) {
    memcpy(p->a_grav, buff, 3 * sizeof(float));
    buff += 3 * sizeof(float);
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_rho) {
    memcpy(&p->mass, buff, sizeof(float));
    buff += sizeof(float);
    memcpy(&p->id_or_neg_offset, buff, sizeof(long long));
    buff += sizeof(long long);
  }

  /* Finally, return the mask of the values we just read. */
  return mask;
}

/**
 * @brief Read a logger message for a timestamp.
 *
 * @param t The timestamp in which to store the value.
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
  if (!(mask & logger_mask_timestamp))
    error("Trying to read timestamp from a particle.");

  /* Make sure we don't have extra fields. */
  if (mask != logger_mask_timestamp)
    error("Timestamp message contains extra fields.");

  /* Copy the timestamp value from the buffer. */
  memcpy(t, buff, sizeof(unsigned long long int));
  buff += sizeof(unsigned long long int);

  /* Copy the timestamp value from the buffer. */
  memcpy(time, buff, sizeof(double));

  /* Finally, return the mask of the values we just read. */
  return mask;
}

#endif /* WITH_LOGGER */

#endif /* HAVE_POSIX_FALLOCATE */
