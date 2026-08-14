/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Stuart McAlpine (stuart.mcalpine@helsinki.fi)
 *               2024 Will J. Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_ZOOM_H
#define SWIFT_ZOOM_H

/* Config parameters. */
#include <config.h>

/* HDF5 header. */
#include <hdf5.h>

/* Local includes */
#include "part_type.h"
#include "timeline.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

/* METIS/ParMETIS headers only used when MPI is also available. */
#ifdef HAVE_PARMETIS
#include <parmetis.h>
#endif
#ifdef HAVE_METIS
#include <metis.h>
#endif

/* Avoid cyclic inclusions. */
struct swift_params;
struct space;
struct cell;
struct engine;
struct io_props;
struct partition;

/* Define a constant for the background task depth. */
#define zoom_bkg_subdepth_diff_grav_default 4
extern int zoom_bkg_subdepth_diff_grav;

/* Zoom region and cell grid initialisation functions. */
void zoom_props_init(struct swift_params *params, struct space *s,
                     const int verbose);
void zoom_region_init(struct space *s, const int regridding, const int verbose);
void zoom_get_region_dim_and_shift(struct space *s, const int verbose);
void zoom_apply_zoom_shift_to_particles(struct space *s, const int verbose);
void zoom_report_cell_properties(const struct space *s);
void zoom_dump_geometry(const struct engine *e);

/* Construct top level cells with a zoom region. */
void zoom_construct_tl_cells(struct space *s, const integertime_t ti_current,
                             int verbose);

/* Linking zoom cells to void leaves. */
void zoom_link_void_leaves(struct space *s, struct cell *c);

/* Space regridding functions. */
int zoom_need_regrid(const struct space *s, const int new_cdim[3]);
void zoom_prepare_cells(struct space *s, const int zoom_cdim[3], int verbose);
void zoom_allocate_cells(struct space *s);

/* Void cell tree construction function. */
void zoom_void_split_recursive(struct space *s, struct cell *c,
                               const short int tpid);
void zoom_void_space_split(struct space *s, int verbose);

/* Task creation functions. */
void zoom_engine_make_self_gravity_tasks(struct space *s, struct engine *e);

/* Void cell super mapper. */
void zoom_cell_set_void_super_mapper(void *map_data, int num_elements,
                                     void *extra_data);

/* Void cell gravity task creation. */
void zoom_engine_make_hierarchical_void_tasks(struct engine *e);

/* Update the void cell gravity timesteps. */
void zoom_void_timestep_collect(struct engine *e);

/* Zoom proxy creation functions. */
void zoom_engine_makeproxies(struct engine *e);

/* Zoom partitioning functions. */
void partition_zoom_grid(struct partition *initial_partition, int nr_nodes,
                         struct space *s);
void partition_zoom_vector(int nr_nodes, struct space *s);
void zoom_partition_voids(struct space *s, int nodeID);

#if defined(WITH_MPI) && (defined(HAVE_METIS) || defined(HAVE_PARMETIS))
/* Zoom metis-specific partitioning functions. */
int zoom_partition_count_vertex_edges(struct space *s, int periodic,
                                      int *cell_edge_offsets);
void zoom_partition_sizes_to_edges(struct space *s, double *counts,
                                   double *edges, const int *cell_edge_offsets);
void zoom_partition_graph_init(struct space *s, int periodic, idx_t *adjncy,
                               int *nadjcny, idx_t *xadj, int *nxadj,
                               const int *cell_edge_offsets);
#endif

/* Zoom specific IO. */
void zoom_write_metadata(hid_t root_grp, hid_t head_grp, const struct space *s);
void zoom_write_particle_counts(
    hid_t head_grp, const long long this_file[swift_type_count],
    const long long in_cells_this_file[swift_type_count],
    const long long total[swift_type_count],
    const long long in_cells_total[swift_type_count],
    const int num_fields[swift_type_count]);
void zoom_io_count_particles_in_cells(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    long long local_in_cells[swift_type_count]);
#ifdef WITH_MPI
void zoom_io_prepare_particle_layout(
    const struct engine *e, const int subsample[swift_type_count],
    const float subsample_fraction[swift_type_count],
    const long long local[swift_type_count],
    const long long total[swift_type_count],
    const long long offset[swift_type_count], MPI_Comm comm,
    long long local_in_cells[swift_type_count],
    long long total_in_cells[swift_type_count],
    long long offset_in_cells[swift_type_count],
    long long offset_outside_cells[swift_type_count]);
#endif
void zoom_io_advance_particle_pointers(struct io_props *props, size_t offset);
void zoom_io_map_virtual_particle_regions(
    hid_t h_prop, hid_t h_space, hid_t h_source_space, const char *file_name,
    const char *source_dataset_name, int dimension, hsize_t particle_count,
    hsize_t count_in_cells, hsize_t start_in_cells[2],
    hsize_t start_outside_cells[2]);
void zoom_io_write_serial_particle_regions(
    hid_t h_data, hid_t h_memspace, hid_t h_filespace, hid_t h_type,
    const void *buffer, int rank, int dimension, size_t count,
    size_t count_in_cells, long long offset_in_cells,
    long long offset_outside_cells, const char *field_name);
void zoom_unshift_pos(const struct space *s, double pos[3]);

#endif /* SWIFT_ZOOM_H */
