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
#include "parser.h"
#include "periodic_replications.h"
#include "timeline.h"
#include "part_type.h"
#include "particle_buffer.h"


/* Avoid cyclic inclusions */
struct cosmology;
struct engine;
struct gpart;
struct space;

/**
 * @brief Lightcone data
 */
struct lightcone_props {

  /*! Whether we're doing lightcone outputs */
  int enabled;

  /*! Output base name */
  char basename[PARSER_MAX_LINE_SIZE];

  /*! Position of the observer in the simulation box */
  double observer_position[3];

  /*! Redshift range the lightcone covers */
  double z_min, z_max;

  /*! Whether we're doing a pencil beam */
  int pencil_beam;

  /*! Vector along the pencil beam */
  double view_vector[3];

  /*! Radius of the pencil beam in radians */
  double view_radius;

  /*! Simulation box size (volume must be a cube) */
  double boxsize;

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

};


void lightcone_init(struct lightcone_props *props,
                    const struct space *s,
                    const struct cosmology *cosmo,
                    struct swift_params *params,
                    const int restart);


void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream);

void lightcone_struct_restore(struct lightcone_props *props, FILE *stream);

void lightcone_init_replication_list(struct lightcone_props *props,
                                     const struct cosmology *cosmo,
                                     const integertime_t ti_old,
                                     const integertime_t ti_current,
                                     const double dt_max);

void lightcone_check_particle_crosses(struct lightcone_props *props,
                                      const struct cosmology *c, const struct gpart *gp,
                                      const double *x, const float *v_full,
                                      const double dt_drift, const integertime_t ti_old,
                                      const integertime_t ti_current,
                                      void *extra1, void *extra2);

void lightcone_flush_buffers(struct lightcone_props *props,
                               int flush_all, int end_file);

#endif /* SWIFT_LIGHTCONE_H */
