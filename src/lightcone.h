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

#define LIGHTCONE_MAX_HEALPIX_MAPS 10
#define LIGHTCONE_MAX_SHELLS       200

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "lightcone_map.h"
#include "lightcone_replications.h"
#include "parser.h"
#include "part_type.h"
#include "particle_buffer.h"
#include "timeline.h"

/* Avoid cyclic inclusions */
struct cosmology;
struct engine;
struct space;

/**
 * @brief Lightcone data
 */
struct lightcone_props {

  /*! Whether we're doing lightcone outputs */
  int enabled;

  /*! Which particle types we're doing */
  int use_type[swift_type_count];

  /*! Output base name */
  char basename[PARSER_MAX_LINE_SIZE];

  /*! Position of the observer in the simulation box */
  double observer_position[3];

  /*! Redshift range in which we will output particles */
  double z_min_for_particles, z_max_for_particles;

  /*! Range in distance squared in which we output particles */
  double r2_min_for_particles, r2_max_for_particles;

  /*! Range in expansion factor covered by particle outputs and healpix maps */
  double a_min, a_max;
  
  /*! Corresponding range in distance squared for a_max and a_min */
  double r2_min, r2_max;
  
  /*! Size of chunks in particle buffer */
  int buffer_chunk_size;

  /*! Size of chunks in HDF5 output files */
  int hdf5_chunk_size;

  /*! Simulation box size (volume must be a cube) */
  double boxsize;

  /*! Top level cell width */
  double cell_width;

  /*! Whether list of replications exists */
  int have_replication_list;

  /*! List of periodic replications to check on this timestep */
  struct replication_list replication_list;

  /*! Total number of particles written to the lightcone by this MPI rank */
  long long tot_num_particles_written[swift_type_count];

  /*! Number of particles written to the current file by this MPI rank */
  long long num_particles_written_to_file[swift_type_count];

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

  /*! Name of the file with radii of spherical shells */
  char radius_file[PARSER_MAX_LINE_SIZE];

  /*! Healpix nside parameter */
  int nside;

  /*! Number of healpix maps we're making from this lightcone */
  int nr_maps;

  /*! Names of the healpix maps to make for this lightcone */
  char map_names[LIGHTCONE_MAX_HEALPIX_MAPS][PARSER_MAX_LINE_SIZE];

  /*! Number of shells */
  int nr_shells;

  /*! Inner radii of shells for this lightcone */
  double shell_rmin[LIGHTCONE_MAX_SHELLS];

  /*! Outer radii of shells for this lightcone */
  double shell_rmax[LIGHTCONE_MAX_SHELLS];

  /*! Minimum expansion factor for each shell */
  double shell_amin[LIGHTCONE_MAX_SHELLS];

  /*! Maximum expansion factor for each shell */
  double shell_amax[LIGHTCONE_MAX_SHELLS];

  /*! Array of pointers to lightcone maps */
  struct lightcone_map *map[LIGHTCONE_MAX_HEALPIX_MAPS][LIGHTCONE_MAX_SHELLS];

};


void lightcone_init(struct lightcone_props *props,
                    const struct space *s,
                    const struct cosmology *cosmo,
                    struct swift_params *params);

void lightcone_clean(struct lightcone_props *props);

void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream);

void lightcone_struct_restore(struct lightcone_props *props, FILE *stream);

void lightcone_init_replication_list(struct lightcone_props *props,
                                     const struct cosmology *cosmo,
                                     const integertime_t ti_old,
                                     const integertime_t ti_current,
                                     const double dt_max);

void lightcone_flush_buffers(struct lightcone_props *props,
                               int flush_all, int end_file);

#endif /* SWIFT_LIGHTCONE_H */
