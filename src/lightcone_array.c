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
#include <string.h>

/* This object's header. */
#include "lightcone_array.h"

/* Local headers */
#include "common_io.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "lightcone.h"
#include "lightcone_particle_io.h"
#include "lightcone_replications.h"
#include "parser.h"
#include "particle_buffer.h"
#include "periodic.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"


/**
 * @brief Initialise the properties of the lightcone code.
 *
 */
void lightcone_array_init(struct lightcone_array_props *props,
                          const struct space *s,
                          const struct cosmology *cosmo,
                          struct swift_params *params,
                          const int verbose) {
  
  /* Determine number of lightcones */
  props->nr_lightcones = 0;
  for(int lightcone_nr=0; lightcone_nr<=MAX_LIGHTCONES; lightcone_nr+=1) {
    char name[PARSER_MAX_LINE_SIZE];
    snprintf(name, PARSER_MAX_LINE_SIZE, "Lightcone%d:enabled", props->nr_lightcones);
    if(parser_get_opt_param_int(params, name, 0)) {      
      props->nr_lightcones += 1;
    }
  }
  message("found %d lightcones to generate", props->nr_lightcones);
  
  /* Allocate array of lightcones */
  props->lightcone = malloc(sizeof(struct lightcone_props)*props->nr_lightcones);
  if(!props->lightcone)error("Failed to allocate lightcone array");

  /* Initialise lightcones */
  props->nr_lightcones = 0;
  for(int lightcone_nr=0; lightcone_nr<=MAX_LIGHTCONES; lightcone_nr+=1) {
    char name[PARSER_MAX_LINE_SIZE];
    snprintf(name, PARSER_MAX_LINE_SIZE, "Lightcone%d:enabled", props->nr_lightcones);
    if(parser_get_opt_param_int(params, name, 0)) {      
      snprintf(name, PARSER_MAX_LINE_SIZE, "Lightcone%d", props->nr_lightcones);
      lightcone_init(props->lightcone+props->lightcone_nr,
                     name, s, cosmo, params, verbose);
      props->nr_lightcones += 1;
    }
  }
}

void lightcone_array_clean(struct lightcone_array_props *props) {
  
  for(int i=0; i<props->nr_lightcones; i+=1)
    lightcone_clean(props->lightcone+i);
  free(props->lightcone);

}

void lightcone_array_struct_dump(const struct lightcone_array_props *props, FILE *stream) {
  
  struct lightcone_array_props tmp = *props;
  tmp.props = NULL;
  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_array_props), 1, stream,
                       "lightcone_array_props", "lightcone_array_props");
  
  for(int i=0; i<props->nr_lightcones; i+=1)
    lightcone_struct_dump(props->lightcone+i, stream);
}


void lightcone_array_struct_restore(struct lightcone_array_props *props, FILE *stream) {

  restart_read_blocks((void *)props, sizeof(struct lightcone_array_props), 1, stream,
                      NULL, "lightcone_array_props");

  props->lightcone = malloc(sizeof(struct lightcone_props)*props->nr_lightcones);
  if(!props->lightcone)error("Failed to allocate lightcone array");
  
  for(int i=0; i<props->nr_lightcones; i+=1)
    lightcone_struct_restore(props->lightcone+i, stream);

}


void lightcone_array_prepare_for_step(struct lightcone_array_props *props,
                                      const struct cosmology *cosmo,
                                      const integertime_t ti_old,
                                      const integertime_t ti_current,
                                      const double dt_max) {

  for(int i=0; i<props->nr_lightcones; i+=1)
    lightcone_prepare_for_step(props->lightcone+i, cosmo, ti_old, ti_current, dt_max);
}


int lightcone_array_trigger_map_update(struct lightcone_array_props *props) {
  
  for(int i=0; i<props->nr_lightcones; i+=1) {
    if(lightcone_trigger_map_update(props->lightcone+i))return 1;
  }
  return 0;
}
