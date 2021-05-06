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
  
  if(engine_rank==0)
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
      lightcone_init(props->lightcone+lightcone_nr,
                     name, lightcone_nr, s, cosmo, params, verbose);
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
  tmp.lightcone = NULL;
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

/**
 * @brief Flush buffers for all lightcones in the array
 *
 * Buffers are flushed if they get large or a flush is forced
 * by setting one of the input flags.
 *
 * props the #lightcone_array_props struct
 * flush_map_updates force full update of the healpix maps
 * flush_particles force output of all buffered particles
 * end_file start a new file next time particles are written out
 * dump_all_shells immediately output all remaining healpix maps
 *
 */
void lightcone_array_flush(struct lightcone_array_props *props,
                           const struct cosmology *cosmo,
                           const struct unit_system *internal_units,
                           const struct unit_system *snapshot_units,
                           int flush_map_updates, int flush_particles,
                           int end_file, int dump_all_shells) {

  /* Loop over lightcones */
  const int nr_lightcones = props->nr_lightcones;
  for(int lightcone_nr=0; lightcone_nr<nr_lightcones; lightcone_nr+=1) {

    /* Get a pointer to this lightcone */
    struct lightcone_props *lc_props = props->lightcone+lightcone_nr;

    /* Apply lightcone map updates if requested */
    if(flush_map_updates) lightcone_flush_map_updates(lc_props);

    /* Flush particle buffers if they're large or flag is set */
    lightcone_flush_particle_buffers(lc_props,
                                     internal_units, snapshot_units,
                                     flush_particles, end_file);
    
    /* Write out any completed healpix maps */
    lightcone_dump_completed_shells(lc_props, cosmo, internal_units,
                                    snapshot_units, dump_all_shells,
                                    /*need_flush=*/!flush_map_updates);
  }

}


/**
 * @brief Make a refined replication list for each lightcone
 *
 * Returns an array of struct #replication_list. Must be freed
 * with lightcone_array_free_replications().
 *
 * props the #lightcone_array_props struct
 * cell the #cell for which we're making replication lists
 *
 */
struct replication_list *lightcone_array_refine_replications(struct lightcone_array_props *props,
                                                             const struct cell *cell) {

  /* Get number of lightcones */
  const int nr_lightcones = props->nr_lightcones;

  /* Allocate a replication list for each lightcone */
  struct replication_list *lists = malloc(sizeof(struct replication_list)*nr_lightcones);

  /* Loop over lightcones */
  for(int lightcone_nr=0; lightcone_nr<nr_lightcones; lightcone_nr+=1) {

    /* Make refined replication list for this lightcone */
    struct lightcone_props *lightcone = props->lightcone+lightcone_nr;
    replication_list_subset_for_cell(&lightcone->replication_list,
                                     cell, lightcone->observer_position,
                                     lists+lightcone_nr);
  }
  
  return lists;
}


/**
 * @brief Free lists returned by lightcone_array_refine_replications
 *
 * props the #lightcone_array_props struct
 * lists the array of struct #replication_list to free
 *
 */
void lightcone_array_free_replications(struct lightcone_array_props *props,
                                       struct replication_list *lists) {
  
  /* Get number of lightcones */
  const int nr_lightcones = props->nr_lightcones;

  /* Loop over lightcones and clean replication lists */
  for(int lightcone_nr=0; lightcone_nr<nr_lightcones; lightcone_nr+=1) {
    replication_list_clean(lists+lightcone_nr);
  }

  /* Free replication list array */
  free(lists);
}
