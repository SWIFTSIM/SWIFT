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
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

/* This object's header. */
#include "logger.h"

/* Local headers. */
#include "atomic.h"
#include "dump.h"
#include "error.h"
#include "part.h"
#include "units.h"
#include "engine.h"

char logger_version[logger_version_size] = "0.1";


const unsigned int logger_datatype_size[logger_data_count] = {
  sizeof(int),
  sizeof(float),
  sizeof(double),
  sizeof(char),
  sizeof(long long),
  1,
};

/**
 * @brief write chunk header
 *
 * @param buff The writing buffer
 * @param mask The mask to write
 * @param offset The old offset
 * @param offset_new The new offset
 *
 * @return updated buff
 */
char *logger_write_chunk_header(char *buff, const unsigned int *mask, const size_t *offset, const size_t offset_new) {
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
 * @brief write a data to the file
 *
 * @param d #dump file
 * @param offset offset at which to write
 * @param size number of bytes to write
 * @param p pointer to the data
 */
void logger_write_data(struct dump *d, size_t *offset, const size_t size, void *const p)
{
  char *buff = dump_get(d, size, offset);
  memcpy(buff, p, size);
}

/**
 * @brief write a general data to the file
 *
 * write data in the following order: name, data type, data
 *
 * @param d #dump file
 * @param log #logger_parameters file format informations
 * @param offset offset at which to write (moved by the data size)
 * @param p pointer to the data
 * @param name data name (should be smaller than log->name)
 * @param data_type #logger_datatype to write
 */
void logger_write_general_data(struct dump *d, struct logger_parameters *params, size_t *offset,
			       const void *p, char* name, size_t data_type)
{
  char *buff;
  /* write name */
  buff = dump_get(d, params->name_size, offset);
  memcpy(buff, name, params->name_size);

  /* write data type */
  buff = dump_get(d, params->data_type_size, offset);
  memcpy(buff, &data_type, params->data_type_size);

  /* write value */
  if (data_type >= logger_data_count)
    error("Not implemented");
  size_t size = logger_data_size[data_type];
  
  buff = dump_get(d, size, offset);
  memcpy(buff, p, size);

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
  int size = 8;

  /* Is this a particle or a timestep? */
  if (mask & logger_mask_timestamp) {

    /* The timestamp should not contain any other bits. */
    if (mask != logger_mask_timestamp)
      error("Timestamps should not include any other data.");

    /* A timestamp consists of an unsigned long long int. */
    size += sizeof(unsigned long long int);

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
 * @brief log all particles
 *
 * @param log The #logger
 * @param e The #engine
 */
void logger_log_all(struct logger *log, struct engine *e) {
  struct space *s = e->s;
  const unsigned int mask = logger_mask_x | logger_mask_v | logger_mask_a |
    logger_mask_u | logger_mask_h | logger_mask_rho |
    logger_mask_consts;
  
  for(long long i=0; i < e->total_nr_parts; i++) {
    logger_log_part(log, &s->parts[i], mask, &s->xparts[i].logger_data.last_offset);
    s->xparts[i].logger_data.last_output = 0;
  }

  if (e->total_nr_gparts > 0)
    error("Not implemented");

}

/**
 * @brief Dump a #part to the log.
 *
 * @param log The #logger
 * @param p The #part to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle.
 */
void logger_log_part(struct logger *log, const struct part *p, const unsigned int mask, size_t *offset) {

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
 * @param offset Pointer to the offset of the previous log of this particle.
 */
void logger_log_gpart(struct logger *log, const struct gpart *p, const unsigned int mask, size_t *offset) {

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
 * @param offset In: previous offset, out: offset of this chunk
 */
void logger_log_timestamp(struct logger *log, integertime_t timestamp, size_t *offset) {
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

  /* Update the log message offset. */
  *offset = offset_new;
}


/**
 * @brief ensure that the buffer is large enough
 *
 * Check if logger parameters are large enough to write all particles
 * and ensure that enough space is available in the buffer
 *
 * @param log The #logger
 * @param total_nr_nparts total number of particle
 */
void logger_ensure_size(
    struct logger *log, size_t total_nr_parts,
    size_t total_nr_gparts, size_t total_nr_sparts) {
  size_t limit, i;
  struct logger_parameters *log_params = log->params;

  limit = log_params->offset_size + log_params->mask_size;

  for(i=0; i < log_params->nber_mask; i++) {
    if (log_params->masks[i] != logger_mask_timestamp)
      limit += log_params->masks_size[i];
  }

  limit *= total_nr_parts;

  if (log->buffer_size < limit) error("Need a larger logger size");

  if (total_nr_gparts > 0)
    error("Not implemented");

  if (total_nr_sparts > 0)
    error("Not implemented");
  
  dump_ensure(log->dump, log->buffer_size);
}

/**
 * @brief intialize the logger structure
 *
 * @param log The #logger
 * @param params The #swift_params
 * @param e The #engine
 */
void logger_init(struct logger *log, const struct swift_params *params, const struct engine *e) {
  /* read parameters */
  log->delta_step = parser_get_param_int(params, "Logger:delta_step");
  log->buffer_size = parser_get_param_float(params, "Logger:mmaped_buffer_size") * 1e9;
  parser_get_param_string(params, "Logger:basename", log->base_name);

  /* generate dump filename */
  char logger_name_file[PARSER_MAX_LINE_SIZE];
  strcpy(logger_name_file, log->base_name);
  strcat(logger_name_file, ".dump");

  /* init dump */
  log->dump = malloc(sizeof(struct dump));
  struct dump *dump_file = log->dump;

  /* init parameters */
  log->params = (struct logger_parameters*) malloc(sizeof(struct logger_parameters));
  logger_parameters_init(log->params);
    
  dump_init(dump_file, logger_name_file, log->buffer_size);
  dump_ensure(dump_file, log->buffer_size);
  log->timestamp_offset = 0;
}

/**
 * @brief Close dump file and desallocate memory
 *
 * @param log The #logger
 */
void logger_clean(struct logger *log) {
  dump_close(log->dump);
  free(log->params);
}

/**
 * @brief Write a file header to a logger file
 *
 * @param log The #logger
 * @param dump The #dump in which to log the particle data.
 *
 */
void logger_write_file_header(struct logger *log, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  message("writing header");
#endif
  
  size_t i;
  char *skip_header;
  size_t file_offset;

  const struct logger_parameters log_params = *log->params;
  struct dump *dump = log->dump;
  
  file_offset = dump->file_offset;
    
  if (file_offset != 0) error("Something was already written in the dump file");

  /* Write version information */
  logger_write_data(dump, &file_offset, logger_string_size, &logger_version);
 
  /* write number of bytes used for the offsets */
  logger_write_data(dump, &file_offset, logger_number_size, &log_params.offset_size);

  /* write offset direction */
  int reversed = 0;
  logger_write_data(dump, &file_offset, logger_datatype_size[logger_data_bool], &reversed);

  /* will write the offset of the first particle here */
  skip_header = dump_get(dump, log_params.offset_size, &file_offset);

  /* write number of bytes used for names */
  logger_write_data(dump, &file_offset, logger_number_size, &log_params.name_size);

  /* write number of bytes used for numbers */
  logger_write_data(dump, &file_offset, logger_number_size, &log_params.number_size);

  /* write number of bytes used for masks */
  logger_write_data(dump, &file_offset, logger_number_size, &log_params.mask_size);

  /* write number of masks */
  logger_write_data(dump, &file_offset, log_params.number_size, &log_params.nber_mask);
  
  /* write masks */
  // loop over all mask type
  for(i=0; i<log_params.nber_mask; i++) {
    // mask name
    size_t j = i * log_params.name_size;
    logger_write_data(dump, &file_offset, log_params.name_size, &log_params.masks_name[j]);
    
    // mask
    logger_write_data(dump, &file_offset, log_params.mask_size, &log_params.masks[i]);

    // mask size
    logger_write_data(dump, &file_offset, log_params.number_size, &log_params.masks_size[i]);
  }

  /* write mask data */
  /* loop over each mask and each data in this mask */
  /* write number of bytes for each field */
  /* write data type (float, double, ...) */
  /* write data name (mass, id, ...) */
  
  /* Write data */
  char *name = malloc(sizeof(char) * log_params.name_size);
  strcpy(name, "time_base");
  logger_write_general_data(dump, &log_params, &file_offset, &e->time_base,
			    name, logger_data_double);

  /* last step: write first offset */
  memcpy(skip_header, &file_offset, log_params.offset_size);

  /* free memory */
  free(name);

  dump_ensure(log->dump, log->buffer_size);
  log->timestamp_offset = 0;
}

/**
 * @brief initialize the #logger_parameters with the format informations
 *
 * @param log_params #logger_parameters to initialize
 */
void logger_parameters_init(struct logger_parameters* log_params) {
  log_params->name_size = 20;
  log_params->offset_size = 7;
  log_params->mask_size = 1;
  log_params->number_size = 1;
  log_params->data_type_size = 1;

  log_params->nber_mask = 8;

  char *cur_name;
  char tmp[log_params->name];
  size_t block_size;
  
  // masks value
  log_params->masks = malloc(sizeof(size_t)*log_params->nber_mask);
  log_params->masks[0] = logger_mask_x;
  log_params->masks[1] = logger_mask_v;
  log_params->masks[2] = logger_mask_a;
  log_params->masks[3] = logger_mask_u;
  log_params->masks[4] = logger_mask_h;
  log_params->masks[5] = logger_mask_rho;
  log_params->masks[6] = logger_mask_consts;
  log_params->masks[7] = logger_mask_timestamp;

  // masks name
  block_size = log_params->name_size * log_params->nber_mask;
  log_params->masks_name = malloc(block_size);
  cur_name = log_params->masks_name;

  strcpy(tmp, "position");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "velocity");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "acceleration");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "entropy");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "cutoff radius");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "density");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "consts");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  strcpy(tmp, "timestamp");
  memcpy(cur_name, &tmp, log_params->name_size);
  cur_name += log_params->name_size;

  log_params->masks_size = malloc(sizeof(size_t) * log_params->nber_mask);
  log_params->masks_size[0] = 3 * sizeof(double);
  log_params->masks_size[1] = 3 * sizeof(float);
  log_params->masks_size[2] = 3 * sizeof(float);
  log_params->masks_size[3] = sizeof(float);
  log_params->masks_size[4] = sizeof(float);
  log_params->masks_size[5] = sizeof(float);
  log_params->masks_size[6] = sizeof(float) + sizeof(long long);
  log_params->masks_size[7] = sizeof(integertime_t);

  // todo masks_type

}

/**
 * @brief free the memory allocated when initializing the #logger_parameters
 *
 * @param log_params #logger_parameters to clean
 */
void logger_parameters_free(struct logger_parameters* log_params) {
  free(log_params->masks);
  free(log_params->masks_name);
  free(log_params->masks_data_size);
}


/**
 * @brief read chunk header
 *
 * @param buff The reading buffer
 * @param mask The mask to read
 * @param offset Out: the offset pointed by this chunk (absolute)
 * @param offset_cur The current chunk offset
 *
 * @return Number of bytes read
 */
__attribute__((always_inline)) INLINE static int logger_read_chunk_header(const char *buff, unsigned int *mask, size_t *offset, size_t cur_offset) {
  memcpy(mask, buff, logger_size_mask);
  buff += logger_size_mask;

  *offset = 0;
  memcpy(offset, buff, logger_size_offset);
  *offset = cur_offset - *offset;
  
  return logger_size_mask + logger_size_offset;
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
int logger_read_timestamp(unsigned long long int *t, size_t *offset,
                          const char *buff) {

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

  /* Finally, return the mask of the values we just read. */
  return mask;
}

#endif /* WITH_LOGGER */

#endif /* HAVE_POSIX_FALLOCATE */

