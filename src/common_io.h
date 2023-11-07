/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl).
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
#include <config.h>

/* Local includes. */
#include "part_type.h"

#define FIELD_BUFFER_SIZE 64
#define DESCRIPTION_BUFFER_SIZE 600
#define PARTICLE_GROUP_BUFFER_SIZE 50
#define FILENAME_BUFFER_SIZE 150
#define IO_BUFFER_ALIGNMENT 1024

/* Avoid cyclic inclusion problems */
struct cell;
struct space;
struct part;
struct gpart;
struct velociraptor_gpart_data;
struct spart;
struct bpart;
struct xpart;
struct sink;
struct io_props;
struct engine;
struct threadpool;
struct output_list;
struct output_options;
struct unit_system;

/**
 * @brief The different types of data used in the GADGET IC files.
 *
 * (This is admittedly a poor substitute to C++ templates...)
 */
enum IO_DATA_TYPE {
  INT,
  LONG,
  LONGLONG,
  UINT8,
  UINT,
  UINT64,
  ULONG,
  ULONGLONG,
  FLOAT,
  DOUBLE,
  CHAR,
  SIZE_T,
};

#if defined(HAVE_HDF5)

/* Library header */
#include <hdf5.h>

hid_t io_hdf5_type(enum IO_DATA_TYPE type);

hsize_t io_get_number_element_in_attribute(hid_t attr);
hsize_t io_get_number_element_in_dataset(hid_t dataset);
void io_read_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                       void* data);
void io_read_attribute_graceful(hid_t grp, const char* name,
                                enum IO_DATA_TYPE type, void* data);
void io_assert_valid_header_cosmology(hid_t h_grp, double a);

void io_read_array_attribute(hid_t grp, const char* name,
                             enum IO_DATA_TYPE type, void* data,
                             hsize_t number_element);
void io_read_array_dataset(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                           void* data, hsize_t number_element);
void io_write_attribute(hid_t grp, const char* name, enum IO_DATA_TYPE type,
                        const void* data, int num);

void io_write_attribute_d(hid_t grp, const char* name, double data);
void io_write_attribute_f(hid_t grp, const char* name, float data);
void io_write_attribute_i(hid_t grp, const char* name, int data);
void io_write_attribute_l(hid_t grp, const char* name, long data);
void io_write_attribute_ll(hid_t grp, const char* name, long long data);
void io_write_attribute_s(hid_t grp, const char* name, const char* str);

void io_write_meta_data(hid_t h_file, const struct engine* e,
                        const struct unit_system* internal_units,
                        const struct unit_system* snapshot_units,
                        const int fof);

void io_write_code_description(hid_t h_file);
void io_write_engine_policy(hid_t h_file, const struct engine* e);
void io_write_part_type_names(hid_t h_grp);

void io_write_cell_offsets(hid_t h_grp, const int cdim[3], const double dim[3],
                           const struct cell* cells_top, const int nr_cells,
                           const double width[3], const int nodeID,
                           const int distributed,
                           const int subsample[swift_type_count],
                           const float subsample_fraction[swift_type_count],
                           const int snap_num,
                           const long long global_counts[swift_type_count],
                           const long long global_offsets[swift_type_count],
                           const int to_write[swift_type_count],
                           const int num_fields[swift_type_count],
                           const struct unit_system* internal_units,
                           const struct unit_system* snapshot_units);

void io_read_unit_system(hid_t h_file, struct unit_system* ic_units,
                         const struct unit_system* internal_units,
                         int mpi_rank);
void io_write_unit_system(hid_t h_grp, const struct unit_system* us,
                          const char* groupName);

void io_copy_temp_buffer(void* temp, const struct engine* e,
                         const struct io_props props, size_t N,
                         const struct unit_system* internal_units,
                         const struct unit_system* snapshot_units);

#endif /* HAVE_HDF5 */

size_t io_sizeof_type(enum IO_DATA_TYPE type);
int io_is_double_precision(enum IO_DATA_TYPE type);

long long io_count_gas_to_write(const struct space* s, const int subsample,
                                const float subsample_ratio,
                                const int snap_num);

long long io_count_dark_matter_to_write(const struct space* s,
                                        const int subsample,
                                        const float subsample_ratio,
                                        const int snap_num);

long long io_count_background_dark_matter_to_write(const struct space* s,
                                                   const int subsample,
                                                   const float subsample_ratio,
                                                   const int snap_num);

long long io_count_stars_to_write(const struct space* s, const int subsample,
                                  const float subsample_ratio,
                                  const int snap_num);

long long io_count_sinks_to_write(const struct space* s, const int subsample,
                                  const float subsample_ratio,
                                  const int snap_num);

long long io_count_black_holes_to_write(const struct space* s,
                                        const int subsample,
                                        const float subsample_ratio,
                                        const int snap_num);

long long io_count_neutrinos_to_write(const struct space* s,
                                      const int subsample,
                                      const float subsample_ratio,
                                      const int snap_num);

void io_collect_parts_to_write(const struct part* restrict parts,
                               const struct xpart* restrict xparts,
                               struct part* restrict parts_written,
                               struct xpart* restrict xparts_written,
                               const int subsample, const float subsample_ratio,
                               const int snap_num, const size_t Nparts,
                               const size_t Nparts_written);
void io_collect_sinks_to_write(const struct sink* restrict sinks,
                               struct sink* restrict sinks_written,
                               const int subsample, const float subsample_ratio,
                               const int snap_num, const size_t Nsinks,
                               const size_t Nsinks_written);
void io_collect_sparts_to_write(const struct spart* restrict sparts,
                                struct spart* restrict sparts_written,
                                const int subsample,
                                const float subsample_ratio, const int snap_num,
                                const size_t Nsparts,
                                const size_t Nsparts_written);
void io_collect_bparts_to_write(const struct bpart* restrict bparts,
                                struct bpart* restrict bparts_written,
                                const int subsample,
                                const float subsample_ratio, const int snap_num,
                                const size_t Nbparts,
                                const size_t Nbparts_written);
void io_collect_gparts_to_write(const struct gpart* restrict gparts,
                                const struct velociraptor_gpart_data* vr_data,
                                struct gpart* restrict gparts_written,
                                struct velociraptor_gpart_data* vr_data_written,
                                const int subsample,
                                const float subsample_ratio, const int snap_num,
                                const size_t Ngparts,
                                const size_t Ngparts_written, int with_stf);
void io_collect_gparts_background_to_write(
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* vr_data,
    struct gpart* restrict gparts_written,
    struct velociraptor_gpart_data* vr_data_written, const int subsample,
    const float subsample_ratio, const int snap_num, const size_t Ngparts,
    const size_t Ngparts_written, int with_stf);
void io_collect_gparts_neutrino_to_write(
    const struct gpart* restrict gparts,
    const struct velociraptor_gpart_data* vr_data,
    struct gpart* restrict gparts_written,
    struct velociraptor_gpart_data* vr_data_written, const int subsample,
    const float subsample_ratio, const int snap_num, const size_t Ngparts,
    const size_t Ngparts_written, int with_stf);

void io_prepare_dm_gparts(struct threadpool* tp, struct gpart* const gparts,
                          size_t Ndm);
void io_prepare_dm_background_gparts(struct threadpool* tp,
                                     struct gpart* const gparts, size_t Ndm);
size_t io_count_dm_background_gparts(const struct gpart* const gparts,
                                     size_t Ndm);
void io_prepare_dm_neutrino_gparts(struct threadpool* tp,
                                   struct gpart* const gparts, size_t Ndm);
size_t io_count_dm_neutrino_gparts(const struct gpart* const gparts,
                                   size_t Ndm);
void io_duplicate_hydro_gparts(struct threadpool* tp, struct part* const parts,
                               struct gpart* const gparts, size_t Ngas,
                               size_t Ndm);
void io_duplicate_stars_gparts(struct threadpool* tp,
                               struct spart* const sparts,
                               struct gpart* const gparts, size_t Nstars,
                               size_t Ndm);
void io_duplicate_sinks_gparts(struct threadpool* tp, struct sink* const sinks,
                               struct gpart* const gparts, size_t Nsinks,
                               size_t Ndm);
void io_duplicate_black_holes_gparts(struct threadpool* tp,
                                     struct bpart* const bparts,
                                     struct gpart* const gparts, size_t Nstars,
                                     size_t Ndm);

void io_prepare_output_fields(struct output_options* output_options,
                              const int with_cosmology, const int with_fof,
                              const int with_stf, int verbose);

void io_write_output_field_parameter(const char* filename, int with_cosmology,
                                     int with_fof, int with_stf);

void io_make_snapshot_subdir(const char* dirname);

void io_get_snapshot_filename(char filename[FILENAME_BUFFER_SIZE],
                              char xmf_filename[FILENAME_BUFFER_SIZE],
                              const struct output_list* output_list,
                              const int snapshots_invoke_stf,
                              const int stf_count, const int snap_count,
                              const char* default_subdir, const char* subdir,
                              const char* default_basename,
                              const char* basename);

void io_set_ids_to_one(struct gpart* gparts, const size_t Ngparts);

void io_select_hydro_fields(const struct part* const parts,
                            const struct xpart* const xparts,
                            const int with_cosmology, const int with_cooling,
                            const int with_temperature, const int with_fof,
                            const int with_stf, const int with_rt,
                            const struct engine* const e, int* const num_fields,
                            struct io_props* const list);

void io_select_dm_fields(const struct gpart* const gparts,
                         const struct velociraptor_gpart_data* gpart_group_data,
                         const int with_fof, const int with_stf,
                         const struct engine* const e, int* const num_fields,
                         struct io_props* const list);

void io_select_neutrino_fields(
    const struct gpart* const gparts,
    const struct velociraptor_gpart_data* gpart_group_data, const int with_fof,
    const int with_stf, const struct engine* const e, int* const num_fields,
    struct io_props* const list);

void io_select_sink_fields(const struct sink* const sinks,
                           const int with_cosmology, const int with_fof,
                           const int with_stf, const struct engine* const e,
                           int* const num_fields, struct io_props* const list);

void io_select_star_fields(const struct spart* const sparts,
                           const int with_cosmology, const int with_fof,
                           const int with_stf, const int with_rt,
                           const struct engine* const e, int* const num_fields,
                           struct io_props* const list);

void io_select_bh_fields(const struct bpart* const bparts,
                         const int with_cosmology, const int with_fof,
                         const int with_stf, const struct engine* const e,
                         int* const num_fields, struct io_props* const list);

#endif /* SWIFT_COMMON_IO_H */
