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
#ifndef SWIFT_LOGGER_H
#define SWIFT_LOGGER_H

#ifdef WITH_LOGGER

/* Includes. */
#include "common_io.h"
#include "engine.h"
#include "logger_struct.h"
#include "part.h"
#include "units.h"

/* Forward declaration */
struct dump;

/**
 * Logger entries contain messages representing the particle data at a given
 * point in time during the simulation.
 *
 * The logger messages always start with an 8-byte header structured as
 * follows:
 *
 *   data: [ mask |                     offset                     ]
 *   byte: [  01  |  02  |  03  |  04  |  05  |  06  |  07  |  08  ]
 *
 * I.e. a first "mask" byte followed by 7 "offset" bytes. The mask contains
 * information on what kind of data is packed after the header. The mask
 * bits correspond to the following data:
 *
 *   bit | name   | size | comment
 *   -------------------------------------------------------------------------
 *   0   | x      | 24   | The particle position, in absolute coordinates,
 *       |        |      | stored as three doubles.
 *   1   | v      | 12   | Particle velocity, stored as three floats.
 *   2   | a      | 12   | Particle acceleration, stored as three floats.
 *   3   | u      | 4    | Particle internal energy (or entropy, if Gadget-SPH
 *       |        |      | is used), stored as a single float.
 *   4   | h      | 4    | Particle smoothing length (or epsilon, if a gpart),
 *       |        |      | stored as a single float.
 *   5   | rho    | 4    | Particle density, stored as a single float.
 *   6   | consts | 12   | Particle constants, i.e. mass and ID.
 *   7   | time   | 8    | Timestamp, not associated with a particle, just
 *       |        |      | marks the transitions from one timestep to another.
 *
 * There is no distinction between gravity and SPH particles.
 *
 * The offset refers to the relative location of the previous message for the
 * same particle or for the previous timestamp (if mask bit 7 is set). I.e.
 * the previous log entry will be at the address of the current mask byte minus
 * the unsigned value stored in the offset. An offset equal to the chunk offset
 * indicated that this is the first message for the given particle/timestamp.
 */

/* Some constants. */
#define logger_mask_x 1
#define logger_mask_v 2
#define logger_mask_a 4
#define logger_mask_u 8
#define logger_mask_h 16
#define logger_mask_rho 32
#define logger_mask_consts 64
#define logger_mask_timestamp 128

#define logger_size_mask 1 // size of the mask
#define logger_size_offset 7 // size of the offset

/* header constants
 * Thoses are definitions from the format and therefore should not be changed!
 * Size in bytes
 */
#define LOGGER_VERSION_SIZE 20 // size of the version message
#define LOGGER_NAME_SIZE 2 // size of the labels size information
#define LOGGER_MASK_SIZE 1 // size of the masks size information
#define LOGGER_NBER_SIZE 1 // size of the number of elements information
#define LOGGER_OFFSET_SIZE 1// size of the offset size information
#define LOGGER_DATATYPE_SIZE 1 // size of the data type information

extern char LOGGER_VERSION[LOGGER_VERSION_SIZE];


struct logger_const {
  size_t name; // labels size
  size_t offset; // offset size
  size_t mask; // mask size
  size_t number; // number size
  size_t nber_mask; // number of different masks
  size_t *masks; // value of each masks (e.g. logger_mask_...)
  size_t *masks_size; // size of each mask
  char *masks_name; // label of each mask
  char *masks_type; // type of data (e.g. 'CHAR', 'INT', 'FLOAT')
};

enum logger_datatype {
  logger_data_int,
  logger_data_float,
  logger_data_double,
  logger_data_char,
  logger_data_longlong,
  logger_data_bool,
  logger_data_count /* should be last */
};

extern const unsigned int logger_data_size[];

/* Function prototypes. */
int logger_compute_chunk_size(unsigned int mask);
void logger_log_all(struct logger *log, struct engine *e);
void logger_log_part(struct logger *log, const struct part *p, const unsigned int mask, size_t *offset);
void logger_log_gpart(struct logger *log, const struct gpart *p, const unsigned int mask, size_t *offset);
void logger_init(struct logger *log, const struct swift_params *params, const struct engine *e);
void logger_clean(struct logger *log);
void logger_log_timestamp(struct logger *log, integertime_t t, size_t *offset);
void logger_ensure_size(struct logger *log, size_t total_nr_parts,
    size_t total_nr_gparts, size_t total_nr_sparts);
void logger_write_file_header(struct logger *log, const struct engine* e);

int logger_read_part(struct part *p, size_t *offset, const char *buff);
int logger_read_gpart(struct gpart *p, size_t *offset, const char *buff);
int logger_read_timestamp(unsigned long long int *t, size_t *offset,
                          const char *buff);

void logger_const_init(struct logger_const* log_const);
void logger_const_free(struct logger_const* log_const);

#endif /* WITH_LOGGER */

#endif /* SWIFT_LOGGER_H */
