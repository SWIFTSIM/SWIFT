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

/* Some standard headers. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "logger.h"

/* Local headers. */
#include "atomic.h"
#include "dump.h"
#include "error.h"
#include "part.h"

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
  char *buff = dump_get(dump, size, &offset_new);

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
  if (mask & logger_mask_rho) {
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
  char *buff = dump_get(dump, size, &offset_new);

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

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_h) {
    memcpy(buff, &p->epsilon, sizeof(float));
    buff += sizeof(float);
  }

  /* Particle constants, which is a bit more complicated. */
  if (mask & logger_mask_rho) {
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
  char *buff = dump_get(dump, size, &offset_new);

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

  /* Particle smoothing length as a single float. */
  if (mask & logger_mask_h) {
    memcpy(&p->epsilon, buff, sizeof(float));
    buff += sizeof(float);
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
