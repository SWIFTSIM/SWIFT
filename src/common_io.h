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
#ifndef SWIFT_COMMON_IO_H
#define SWIFT_COMMON_IO_H

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "units.h"

#define FIELD_BUFFER_SIZE 200
#define PARTICLE_GROUP_BUFFER_SIZE 50
#define FILENAME_BUFFER_SIZE 150
#define IO_BUFFER_ALIGNMENT 1024

/* Avoid cyclic inclusion problems */
struct part;
struct gpart;
struct spart;
struct io_props;
struct engine;
struct threadpool;

/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum IO_DATA_TYPE {
  INT,
  LONG,
  LONGLONG,
  UINT,
  ULONG,
  ULONGLONG,
  FLOAT,
  DOUBLE,
  CHAR
};

/**
 * @brief The different formats for when to run structure finding.
 *
 */
enum IO_STF_OUTPUT_FORMAT { STEPS = 0, TIME };

#if defined(HAVE_HDF5)

hid_t io_hdf5_type(enum IO_DATA_TYPE type);

void io_read_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                       void* data);

void io_write_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                        void* data, int num);

void io_write_attribute_d(hid_t grp, const char* name, double data);
void io_write_attribute_f(hid_t grp, const char* name, float data);
void io_write_attribute_i(hid_t grp, const char* name, int data);
void io_write_attribute_l(hid_t grp, const char* name, long data);
void io_write_attribute_s(hid_t grp, const char* name, const char* str);

void io_write_code_description(hid_t h_file);
void io_write_engine_policy(hid_t h_file, const struct engine* e);

void io_read_unit_system(hid_t h_file, struct unit_system* ic_units,
                         const struct unit_system* internal_units,
                         int mpi_rank);
void io_write_unit_system(hid_t h_grp, const struct unit_system* us,
                          const char* groupName);

void io_copy_temp_buffer(void* temp, const struct engine* e,
                         const struct io_props props, size_t N,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units);

#endif /* defined HDF5 */

size_t io_sizeof_type(enum IO_DATA_TYPE type);
int io_is_double_precision(enum IO_DATA_TYPE type);

void io_collect_dm_gparts(const struct gpart* const gparts, size_t Ntot,
                          struct gpart* const dmparts, size_t Ndm);
void io_prepare_dm_gparts(struct threadpool* tp, struct gpart* const gparts,
                          size_t Ndm);
void io_duplicate_hydro_gparts(struct threadpool* tp, struct part* const parts,
                               struct gpart* const gparts, size_t Ngas,
                               size_t Ndm);
void io_duplicate_star_gparts(struct threadpool* tp, struct spart* const sparts,
                              struct gpart* const gparts, size_t Nstars,
                              size_t Ndm);

void io_check_output_fields(const struct swift_params* params,
                            const long long N_total[3]);

void io_write_output_field_parameter(const char* filename);

#endif /* SWIFT_COMMON_IO_H */
