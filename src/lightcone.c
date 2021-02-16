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


/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>

/* This object's header. */
#include "lightcone.h"

/* Local headers */
#include "parser.h"
#include "restart.h"
#include "error.h"
#include "periodic_replications.h"
#include "cosmology.h"


void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  restart_write_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");
}


void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");
}


/**
 * @brief Initialise the properties of the lightcone code.
 *
 * @param props the #lightcone_props structure to fill.
 * @param params the parameter file parser.
 */
void lightcone_init(struct lightcone_props *props, struct swift_params *params) {
  
  /* Whether we generate lightcone output */
  props->enabled = 1;

  /* Redshift range for the lightcone */
  props->z_min = parser_get_param_double(params, "Lightcone:z_min");
  props->z_max = parser_get_param_double(params, "Lightcone:z_max");

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, "Lightcone:observer_position", 3,
                                props->observer_position);

}


void lightcone_init_replication_list(struct lightcone_props *props,
                                     struct cosmology *cosmo,
                                     struct space *s) {
  /* 
     For now, we'll check all periodic replications between z_min and z_max.
     
     TODO: on each timestep, get limits on times particles can be drifted between
     and regenerate the list of replications.
  */

  /* Get the size of the simulation box */
  double boxsize = s->dim[0];
  if(s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, 1.0/(1.0+props->z_max));
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, 1.0/(1.0+props->z_min));;
  double lightcone_boundary = 0.5*boxsize;
  if(lightcone_rmin > lightcone_rmax)
    error("Lightcone has rmin > rmax - check z_min and z_max parameters?");

  /* Determine periodic copies we need to search */
  replication_list_init(&props->replication_list, boxsize,
                        props->observer_position,
                        lightcone_rmin, lightcone_rmax,
                        lightcone_boundary);
}
