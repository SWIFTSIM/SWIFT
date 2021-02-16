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

  /*! Position of the observer in the simulation box */
  double observer_position[3];

  /*! Redshift range the lightcone covers */
  double z_min, z_max;

  /*! List of periodic replications to check on this timestep */
  struct replication_list replication_list;

};

void lightcone_init(struct lightcone_props *props, struct swift_params *params);

void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream);

void lightcone_struct_restore(struct lightcone_props *props, FILE *stream);

void lightcone_init_replication_list(struct lightcone_props *props,
                                     struct cosmology *cosmo,
                                     struct space *s);

void lightcone_check_gpart_crosses(const struct engine *e, struct gpart *gp,
                                   const double dt_drift, const integertime_t);

#endif /* SWIFT_LIGHTCONE_H */
