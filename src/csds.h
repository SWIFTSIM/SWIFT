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

#include <config.h>

#ifdef WITH_CSDS

/* Includes. */
#include "align.h"
#include "common_io.h"
#include "error.h"
#include "inline.h"
#include "timeline.h"
#include "units.h"

/* Include the CSDS */
#include "csds/src/logfile_writer.h"

/* Forward declaration. */
struct gpart;
struct part;
struct engine;

/**
 * Csds entries contain messages representing the particle data at a given
 * point in time during the simulation.
 *
 * The csds messages always start with an 8-byte header structured as
 * follows:
 *
 *   data: [ mask        |              offset                     ]
 *   byte: [  01  |  02  |  03  |  04  |  05  |  06  |  07  |  08  ]
 *
 * I.e. a first "mask" byte followed by 6 "offset" bytes. The mask contains
 * information on what kind of data is packed after the header. The mask
 * bits correspond to the following data:
 *
 * There is no distinction between gravity and SPH particles.
 *
 * The offset refers to the relative location of the previous message for the
 * same particle or for the previous timestamp. I.e.
 * the previous log entry will be at the address of the current mask byte minus
 * the unsigned value stored in the offset. An offset equal to the record offset
 * indicated that this is the first message for the given particle/timestamp.
 */

/**
 * @brief structure containing global data for the particle csds.
 */
struct csds_writer {
  /* Number of particle updates between log entries. */
  short int delta_step;

  /* Csds basename. */
  char base_name[CSDS_STRING_SIZE];

  /*  The logfile writer. */
  struct csds_logfile_writer logfile;

  /* timestamp offset for csds. */
  size_t timestamp_offset;

  /* scaling factor when buffer is too small. */
  float buffer_scale;

  /* Size of a record if every mask are activated. */
  int max_record_size;

  /* Description of all the fields that can be written. */
  struct csds_field *list_fields;

  /* Pointer to the variable list_fields for each module. */
  struct csds_field *field_pointers[swift_type_count];

  /* Number of fields for each particle type. */
  int number_fields[swift_type_count];

  /* Number of elements in list_fields. */
  int total_number_fields;

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
                            const enum csds_special_flags flag);
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
void csds_ensure_size(struct csds_writer *log, const struct engine *e);
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
