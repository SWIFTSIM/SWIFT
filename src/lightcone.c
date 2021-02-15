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


void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  if(props->replication_list)
    error("Replication list should not be allocated when dumping lightcone_props\n");

  restart_write_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");
}

void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");
  props->replication_list = NULL;
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
  props->z_min = parser_get_param_double(params, "z_min");
  props->z_max = parser_get_param_double(params, "z_max");

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, "observer_position", 3, props->observer_position);

  /* Initially have no list of periodic replications to check */
  props->replication_list = NULL;

}

