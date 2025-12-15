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
#include "timeline.h"

/* Avoid cyclic inclusions. */
struct swift_params;
struct space;
struct cell;
struct engine;
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

/* Zoom specific IO. */
void zoom_write_metadata(hid_t root_grp, hid_t head_grp, const struct space *s);
void zoom_unshift_pos(const struct space *s, double pos[3]);

#endif /* SWIFT_ZOOM_H */
