/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_CSDS_H
#define SWIFT_CSDS_H

#include "../config.h"

#ifdef WITH_CSDS

/* Includes. */
#include "align.h"
#include "common_io.h"
#include "dump.h"
#include "error.h"
#include "inline.h"
#include "timeline.h"
#include "units.h"

/* Forward declaration. */
struct dump;
struct gpart;
struct part;
struct engine;

#define csds_major_version 1
#define csds_minor_version 3
/* Size of the strings. */
#define csds_string_length 200

/*
 * The two following defines need to correspond to the list's order
 * in csds_init_masks.
 */
/* Index of the special flags in the list of masks */
#define csds_index_special_flags 0
/* Index of the timestamp in the list of masks */
#define csds_index_timestamp 1

/**
 * Csds entries contain messages representing the particle data at a given
 * point in time during the simulation.
 *
 * The csds messages always start with an 8-byte header structured as
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
 * the unsigned value stored in the offset. An offset equal to the record offset
 * indicated that this is the first message for the given particle/timestamp.
 */

enum csds_special_flags {
  csds_flag_none = 0,        /* No flag */
  csds_flag_change_type = 1, /* Flag for a change of particle type */
  csds_flag_mpi_enter, /* Flag for a particle received from another  MPI rank
                        */
  csds_flag_mpi_exit,  /* Flag for a particle sent to another MPI rank */
  csds_flag_delete,    /* Flag for a deleted particle */
  csds_flag_create,    /* Flag for a created particle */
} __attribute__((packed));

/**
 * @brief structure containing global data for the particle csds.
 */
struct csds_writer {
  /* Number of particle updates between log entries. */
  short int delta_step;

  /* Csds basename. */
  char base_name[csds_string_length];

  /*  Dump file (In the reader, the dump is cleaned, therefore it is renamed
   * logfile). */
  struct dump dump;

  /* timestamp offset for csds. */
  size_t timestamp_offset;

  /* scaling factor when buffer is too small. */
  float buffer_scale;

  /* Size of a record if every mask are activated. */
  int max_record_size;

  /* Description of all the fields that can be written. */
  struct mask_data *csds_mask_data;

  /* Pointer to the variable csds_mask_data for each module. */
  struct {
    /* pointer for the hydro */
    struct mask_data *hydro;

    /* pointer for the chemistry */
    struct mask_data *chemistry_part;

    /* pointer for the chemistry */
    struct mask_data *chemistry_spart;

    /* pointer for the gravity */
    struct mask_data *gravity;

    /* pointer for the stars */
    struct mask_data *stars;

    /* pointer for the star formation */
    struct mask_data *star_formation;
  } mask_data_pointers;

  /* Number of elements in csds_mask_data. */
  int csds_count_mask;

  /* Maximum size for a hydro record. */
  int max_size_record_part;

  /* Maximum size for a gravity record. */
  int max_size_record_gpart;

  /* Maximum size for a star record. */
  int max_size_record_spart;

} SWIFT_STRUCT_ALIGN;

/* required structure for each particle type. */
struct csds_part_data {
  /* Number of particle updates since last output. */
  int steps_since_last_output;

  /* offset of last particle log entry. */
  uint64_t last_offset;
};

/* Function prototypes. */
void csds_log_all_particles(struct csds_writer *log, const struct engine *e,
                            int first_log);
void csds_log_part(struct csds_writer *log, const struct part *p,
                   struct xpart *xp, const struct engine *e,
                   const int log_all_fields, const enum csds_special_flags flag,
                   const int flag_data);
void csds_log_parts(struct csds_writer *log, const struct part *p,
                    struct xpart *xp, int count, const struct engine *e,
                    const int log_all_fields,
                    const enum csds_special_flags flag, const int flag_data);
void csds_log_spart(struct csds_writer *log, struct spart *p,
                    const struct engine *e, const int log_all_fields,
                    const enum csds_special_flags flag, const int flag_data);
void csds_log_sparts(struct csds_writer *log, struct spart *sp, int count,
                     const struct engine *e, const int log_all_fields,
                     const enum csds_special_flags flag, const int flag_data);
void csds_log_gpart(struct csds_writer *log, struct gpart *p,
                    const struct engine *e, const int log_all_fields,
                    const enum csds_special_flags flag, const int flag_data);
void csds_log_gparts(struct csds_writer *log, struct gpart *gp, int count,
                     const struct engine *e, const int log_all_fields,
                     const enum csds_special_flags flag, const int flag_data);
void csds_init(struct csds_writer *log, const struct engine *e,
               struct swift_params *params);
void csds_free(struct csds_writer *log);
void csds_log_timestamp(struct csds_writer *log, integertime_t t, double time,
                        size_t *offset);
void csds_ensure_size(struct csds_writer *log, size_t total_nr_parts,
                      size_t total_nr_gparts, size_t total_nr_sparts);
void csds_write_file_header(struct csds_writer *log);

int csds_read_part(const struct csds_writer *log, struct part *p,
                   size_t *offset, const char *buff);
int csds_read_gpart(const struct csds_writer *log, struct gpart *p,
                    size_t *offset, const char *buff);
int csds_read_timestamp(const struct csds_writer *log, integertime_t *t,
                        double *time, size_t *offset, const char *buff);
void csds_struct_dump(const struct csds_writer *log, FILE *stream);
void csds_struct_restore(struct csds_writer *log, FILE *stream);

/**
 * @brief Generate the data for the special flags.
 *
 * @param flag The special flag to use.
 * @param flag_data The data to write in the record.
 * @param type The type of the particle.
 */
INLINE static uint32_t csds_pack_flags_and_data(enum csds_special_flags flag,
                                                int flag_data,
                                                enum part_type type) {
#ifdef SWIFT_DEBUG_CHECKS
  if (flag & 0xFFFFFF00) {
    error(
        "The special flag in the particle CSDS cannot be larger than 1 "
        "byte.");
  }
  if (flag_data & ~0xFFFF) {
    error(
        "The data for the special flag in the particle CSDS cannot be larger "
        "than 2 bytes.");
  }
#endif
  return ((uint32_t)flag << (3 * 8)) | ((flag_data & 0xFFFF) << 8) |
         (type & 0xFF);
}

/**
 * @brief Initialize the csds data for a particle.
 *
 * @param csds The #csds_part_data.
 */
INLINE static void csds_part_data_init(struct csds_part_data *csds) {
  csds->last_offset = 0;
  csds->steps_since_last_output = 0;
}

/**
 * @brief Should this particle write its data now ?
 *
 * @param csds_data The #csds_part_data of a particle.
 * @param log The #csds_writer.
 *
 * @return 1 if the particle should be writen, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int csds_should_write(
    const struct csds_part_data *csds_data, const struct csds_writer *log) {

  return (csds_data->steps_since_last_output > log->delta_step);
}

#endif /* WITH_CSDS */

#endif /* SWIFT_CSDS_H */
