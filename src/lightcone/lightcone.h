/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
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

#ifndef SWIFT_LIGHTCONE_H
#define SWIFT_LIGHTCONE_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "lightcone/lightcone_map_types.h"
#include "lightcone/lightcone_particle_io.h"
#include "lightcone/lightcone_replications.h"
#include "lightcone/lightcone_shell.h"
#include "lightcone/pixel_index.h"
#include "lightcone/projected_kernel.h"
#include "parser.h"
#include "part_type.h"
#include "particle_buffer.h"
#include "threadpool.h"
#include "timeline.h"
#include "units.h"

/* Avoid cyclic inclusions */
struct cosmology;
struct engine;
struct space;

/**
 * @brief Lightcone data
 */
struct lightcone_props {

  /*! Index of this lightcone */
  int index;

  /*! Whether to write extra log messages */
  int verbose;

  /*! Which particle types we're doing */
  int use_type[swift_type_count];

  /*! Minimum redshift for particle output for each type */
  double z_min_for_type[swift_type_count];

  /*! Maximum redshift for particle output for each type */
  double z_max_for_type[swift_type_count];

  /*! Minimum a to search for lightcone crossing for each type  */
  double a_min_search_for_type[swift_type_count];

  /*! Maximum a to search for lightcone crossing for each type  */
  double a_max_search_for_type[swift_type_count];

  /*! Whether we need to do lightcone crossing checks for each type at this step
   */
  int check_type_for_crossing[swift_type_count];

  /*! Enable selective output of high redshift gas */
  int gas_filtering_enabled;

  /*! Will output all gas below this redshift */
  double min_z_for_gas_filtering;

  /*! Will output all gas after this scale factor */
  double max_a_for_gas_filtering;

  /*! At z>min_z_for_gas_filtering require gas T>min_temp_for_high_z_gas */
  double min_temp_for_filtered_gas;

  /*! At z>min_z_for_gas_filtering require gas
   * nh>min_nh_for_filtered_gas*(1+z)^4 */
  double min_nh_for_filtered_gas;

  /*! Exclude recently heated gas from xray and sz maps */
  double xray_maps_recent_AGN_injection_exclusion_time;

  /*! Don't exclude gas with temperature less than this factor times AGN_delta_T
   */
  double xray_maps_recent_AGN_min_temp_factor;

  /*! Don't exclude gas with temperature more than this factor times AGN_delta_T
   */
  double xray_maps_recent_AGN_max_temp_factor;

  /*! Output base name */
  char basename[PARSER_MAX_LINE_SIZE];

  /*! Output directory */
  char subdir[PARSER_MAX_LINE_SIZE];

  /*! Position of the observer in the simulation box */
  double observer_position[3];

  /*! Range in distance squared in which we output particles of each type */
  double r2_min_for_type[swift_type_count], r2_max_for_type[swift_type_count];

  /*! Range in expansion factor covered by particle outputs and healpix maps */
  double a_min, a_max;

  /*! Corresponding range in distance squared for a_max and a_min */
  double r2_min, r2_max;

  /*! Size of chunks in particle buffer */
  int buffer_chunk_size;

  /*! Size of chunks in HDF5 output files */
  int hdf5_chunk_size;

  /* Maximum amount of data (in megabytes) to send from any one rank when
   * updating healpix maps */
  double max_map_update_send_size_mb;

  /*! Whether to apply lossy compression */
  int particles_lossy_compression;

  /*! Lossless compression level for particles (0 to disable) */
  int particles_gzip_level;

  /*! Lossless compression level for healpix maps (0 to disable) */
  int maps_gzip_level;

  /*! Simulation box size (volume must be a cube) */
  double boxsize;

  /*! Top level cell width */
  double cell_width;

  /*! Whether list of replications exists */
  int have_replication_list;

  /*! List of periodic replications to check on this timestep */
  struct replication_list replication_list;

  /*! Number of particles written to the current file by this MPI rank */
  long long num_particles_written_to_file[swift_type_count];

  /*! Number of particles of each type which have been output on this rank */
  long long num_particles_written_this_rank[swift_type_count];

  /*! Index of the current output file for this MPI rank */
  int current_file;

  /*! Range of times used to generate the replication list */
  integertime_t ti_old, ti_current;

  /*! Expansion factors corresponding to z_min, z_max */
  double a_at_z_min, a_at_z_max;

  /*! Buffers to store particles on the lightcone */
  struct particle_buffer buffer[swift_type_count];

  /*! Will write particles to disk if buffer exceeds this size */
  int max_particles_buffered;

  /*! Whether we should make a new file on the next flush */
  int start_new_file;

  /*! Whether we have started a particle file and not finalized it yet */
  int file_needs_finalizing;

  /*! Number of pending map updates to trigger communication */
  int max_updates_buffered;

  /*! Whether to write distributed maps in MPI mode */
  int distributed_maps;

  /*! Name of the file with radii of spherical shells */
  char radius_file[PARSER_MAX_LINE_SIZE];

  /*! Healpix nside parameter */
  int nside;

  /*! Healpix pixel area */
  double pixel_area_steradians;

  /*! Number of shells */
  int nr_shells;

  /*! Array of lightcone shells */
  struct lightcone_shell *shell;

  /*! Number of healpix maps we're making for each shell */
  int nr_maps;

  /*! Types of healpix map we're making for each shell */
  struct lightcone_map_type *map_type;

  /*! Range of shells that might be updated this step */
  int shell_nr_min, shell_nr_max;

  /*! Information about each particle type contributing to the maps */
  struct lightcone_particle_type part_type[swift_type_count];

  /*! Output fields */
  struct lightcone_io_field_list particle_fields[swift_type_count];

  /*! Tabulation of projected SPH smoothing kernel */
  struct projected_kernel_table kernel_table;
};

void lightcone_init(struct lightcone_props *props, const int index,
                    const struct space *s, const struct cosmology *cosmo,
                    struct swift_params *params,
                    const struct unit_system *internal_units,
                    const struct phys_const *physical_constants,
                    const int verbose);

void lightcone_clean(struct lightcone_props *props);

void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream);

void lightcone_struct_restore(struct lightcone_props *props, FILE *stream);

void lightcone_prepare_for_step(struct lightcone_props *props,
                                const struct cosmology *cosmo,
                                const integertime_t ti_earliest_undrifted,
                                const integertime_t ti_current);

void lightcone_buffer_particle(struct lightcone_props *props,
                               const struct engine *e, const struct gpart *gp,
                               const double a_cross, const double x_cross[3]);

void lightcone_flush_particle_buffers(struct lightcone_props *props, double a,
                                      const struct unit_system *internal_units,
                                      const struct unit_system *snapshot_units,
                                      int flush_all, int end_file);

void lightcone_buffer_map_update(struct lightcone_props *props,
                                 const struct engine *e, const struct gpart *gp,
                                 const double a_cross, const double x_cross[3]);

void lightcone_flush_map_updates(struct lightcone_props *props,
                                 struct threadpool *tp);

void lightcone_dump_completed_shells(struct lightcone_props *props,
                                     struct threadpool *tp,
                                     const struct cosmology *c,
                                     const struct unit_system *internal_units,
                                     const struct unit_system *snapshot_units,
                                     const int dump_all, const int need_flush);

int lightcone_trigger_map_update(struct lightcone_props *props);

void lightcone_memory_use(struct lightcone_props *props,
                          size_t *particle_buffer_bytes,
                          size_t *map_buffer_bytes, size_t *pixel_data_bytes);

void lightcone_write_index(struct lightcone_props *props,
                           const struct unit_system *internal_units,
                           const struct unit_system *snapshot_units);

void lightcone_map_set_baseline(const struct cosmology *c,
                                struct lightcone_props *props,
                                struct lightcone_map *map);

#endif /* SWIFT_LIGHTCONE_H */
