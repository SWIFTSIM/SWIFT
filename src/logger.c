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


void logger_write_data(struct dump *d, size_t *offset, size_t size, void *p)
{
  char *buff;
  buff = dump_get(d, size, offset);
  memcpy(buff, p, size);
}

/**
 * @brief Compute the size of a message given its mask.
 *
 * @param mask The mask that will be used to dump a #part or #gpart.
 *
 * @return The size of the logger message in bytes.
 */
int logger_size(unsigned int mask) {

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
 * @brief Dump a #part to the log.
 *
 * @param p The #part to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle.
 * @param dump The #dump in which to log the particle data.
 */
void logger_log_part(struct part *p, unsigned int mask, size_t *offset,
                     struct dump *dump) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_timestamp)
    error("You should not log particles as timestamps.");

  /* Start by computing the size of the message. */
  const int size = logger_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  uint64_t temp = (((uint64_t)(offset_new - *offset)) & 0xffffffffffffffULL) |
                  ((uint64_t)mask << 56);
  memcpy(buff, &temp, 8);
  buff += 8;

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
 * @param p The #gpart to dump.
 * @param mask The mask of the data to dump.
 * @param offset Pointer to the offset of the previous log of this particle.
 * @param dump The #dump in which to log the particle data.
 */
void logger_log_gpart(struct gpart *p, unsigned int mask, size_t *offset,
                      struct dump *dump) {

  /* Make sure we're not writing a timestamp. */
  if (mask & logger_mask_timestamp)
    error("You should not log particles as timestamps.");

  /* Make sure we're not looging fields not supported by gparts. */
  if (mask & (logger_mask_u | logger_mask_rho))
    error("Can't log SPH quantities for gparts.");

  /* Start by computing the size of the message. */
  const int size = logger_size(mask);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  uint64_t temp = (((uint64_t)(offset_new - *offset)) & 0xffffffffffffffULL) |
                  ((uint64_t)mask << 56);
  memcpy(buff, &temp, 8);
  buff += 8;

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

void logger_log_timestamp(unsigned long long int timestamp, size_t *offset,
                          struct dump *dump) {

  /* Start by computing the size of the message. */
  const int size = logger_size(logger_mask_timestamp);

  /* Allocate a chunk of memory in the dump of the right size. */
  size_t offset_new;
  char *buff = (char *)dump_get(dump, size, &offset_new);

  /* Write the header. */
  uint64_t temp = (((uint64_t)(offset_new - *offset)) & 0xffffffffffffffULL) |
                  ((uint64_t)logger_mask_timestamp << 56);
  memcpy(buff, &temp, 8);
  buff += 8;

  /* Store the timestamp. */
  memcpy(buff, &timestamp, sizeof(unsigned long long int));

  /* Update the log message offset. */
  *offset = offset_new;
}


void logger_ensure_size(size_t total_nr_parts, size_t logger_size) {
  size_t limit, i;
  struct logger_const log_const;
  logger_const_init(&log_const);


  limit = log_const.offset + log_const.mask;

  for(i=0; i < log_const.nber_mask; i++) {
    limit += log_const.masks_size[i];
  }

  limit *= total_nr_parts;

  if (logger_size < limit) error("Need a larger logger size");
  
  logger_const_free(&log_const);
}

/**
 * @brief Write a file header to a logger file
 *
 * @param offset Pointer to the offset of the previous log of this particle.
 * @param dump The #dump in which to log the particle data.
 *
 */
void logger_write_file_header(struct dump *dump, struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  message("writing header");
#endif
  
  size_t i;
  char *skip_header;
  size_t *file_offset;

  struct logger_const log_const;
  logger_const_init(&log_const);

  file_offset = &dump->file_offset;
    
  if (*file_offset != 0) error("Something was already written in the dump file");

  /* Write version information */
  logger_write_data(dump, file_offset, LOGGER_VERSION_SIZE, LOGGER_VERSION);
 
  /* write number of bytes used for the offsets */
  logger_write_data(dump, file_offset, LOGGER_OFFSET_SIZE, &log_const.offset);

  /* will write the offset of the first particle here */
  skip_header = dump_get(dump, log_const.offset, file_offset);

  /* write number of bytes used for names */
  logger_write_data(dump, file_offset, LOGGER_NAME_SIZE, &log_const.name);

  /* write number of bytes used for numbers */
  logger_write_data(dump, file_offset, LOGGER_NBER_SIZE, &log_const.number);

  /* write number of bytes used for masks */
  logger_write_data(dump, file_offset, LOGGER_MASK_SIZE, &log_const.mask);

  /* write number of masks */
  logger_write_data(dump, file_offset, log_const.number, &log_const.nber_mask);
  
  /* write masks */
  // loop over all mask type
  for(i=0; i<log_const.nber_mask; i++) {
    // mask name
    size_t j = i * log_const.name;
    logger_write_data(dump, file_offset, log_const.name, &log_const.masks_name[j]);
    
    // mask
    logger_write_data(dump, file_offset, log_const.mask, &log_const.masks[i]);

    // mask size
    logger_write_data(dump, file_offset, log_const.number, &log_const.masks_size[i]);
  }

  /* Write data */
  

  /* last step */
  memcpy(skip_header, file_offset, log_const.offset);
  logger_const_free(&log_const);
}

void logger_const_init(struct logger_const* log_const) {
  log_const->name = 20;
  log_const->offset = 7;
  log_const->mask = 1;
  log_const->number = 1;

  log_const->nber_mask = 8;

  char *cur_name;
  char tmp[log_const->name];
  size_t block_size;
  
  // masks value
  log_const->masks = malloc(sizeof(size_t)*log_const->nber_mask);
  log_const->masks[0] = logger_mask_x;
  log_const->masks[1] = logger_mask_v;
  log_const->masks[2] = logger_mask_a;
  log_const->masks[3] = logger_mask_u;
  log_const->masks[4] = logger_mask_h;
  log_const->masks[5] = logger_mask_rho;
  log_const->masks[6] = logger_mask_consts;
  log_const->masks[7] = logger_mask_timestamp;

  // masks name
  block_size = log_const->name * log_const->nber_mask;
  log_const->masks_name = malloc(block_size);
  cur_name = log_const->masks_name;

  strcpy(tmp, "position");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "velocity");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "acceleration");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "entropy");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "cutoff radius");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "density");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "consts");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  strcpy(tmp, "timestamp");
  memcpy(cur_name, &tmp, log_const->name);
  cur_name += log_const->name;

  log_const->masks_size = malloc(sizeof(size_t) * log_const->nber_mask);
  log_const->masks_size[0] = 3 * sizeof(double);
  log_const->masks_size[1] = 3 * sizeof(float);
  log_const->masks_size[2] = 3 * sizeof(float);
  log_const->masks_size[3] = sizeof(float);
  log_const->masks_size[4] = sizeof(float);
  log_const->masks_size[5] = sizeof(float);
  log_const->masks_size[6] = sizeof(float) + sizeof(long long);
  log_const->masks_size[7] = 8;

  // todo masks_type

}

void logger_const_free(struct logger_const* log_const) {
  free(log_const->masks);
  free(log_const->masks_name);
  free(log_const->masks_size);
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
  uint64_t temp;
  memcpy(&temp, buff, 8);
  const int mask = temp >> 56;
  *offset -= temp & 0xffffffffffffffULL;
  buff += 8;

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
  uint64_t temp;
  memcpy(&temp, buff, 8);
  const int mask = temp >> 56;
  *offset -= temp & 0xffffffffffffffULL;
  buff += 8;

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
  uint64_t temp;
  memcpy(&temp, buff, 8);
  const int mask = temp >> 56;
  *offset -= temp & 0xffffffffffffffULL;
  buff += 8;

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

