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

#ifndef SWIFT_LIGHTCONE_ARRAY_H
#define SWIFT_LIGHTCONE_ARRAY_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "lightcone.h"
#include "lightcone_replications.h"
#include "parser.h"
#include "part_type.h"
#include "particle_buffer.h"
#include "timeline.h"

/* Avoid cyclic inclusions */
struct cosmology;
struct engine;
struct space;

#define MAX_LIGHTCONES 8

/**
 * @brief Lightcone data for multiple lightcones
 */
struct lightcone_array_props {

  /*! Number of lightcones */
  int nr_lightcones;
  
  /*! Lightcone properties */
  struct lightcone_props *lightcone;

};


void lightcone_array_init(struct lightcone_array_props *props,
                          const struct space *s,
                          const struct cosmology *cosmo,
                          struct swift_params *params,
                          const int verbose);

void lightcone_array_clean(struct lightcone_array_props *props);

void lightcone_array_struct_dump(const struct lightcone_array_props *props, FILE *stream);

void lightcone_array_struct_restore(struct lightcone_array_props *props, FILE *stream);

void lightcone_array_prepare_for_step(struct lightcone_array_props *props,
                                      const struct cosmology *cosmo,
                                      const integertime_t ti_old,
                                      const integertime_t ti_current,
                                      const double dt_max);

int lightcone_array_trigger_map_update(struct lightcone_array_props *props);

void lightcone_array_flush(struct lightcone_array_props *props,
                           const struct cosmology *cosmo,
                           const struct unit_system *internal_units,
                           const struct unit_system *snapshot_units,
                           int flush_map_updates, int flush_particles,
                           int end_file, int dump_all_shells);

struct replication_list *lightcone_array_refine_replications(struct lightcone_array_props *props,
                                                             const struct cell *cell);

void lightcone_array_free_replications(struct lightcone_array_props *props,
                                       struct replication_list *lists);

#endif /* SWIFT_LIGHTCONE_ARRAY_H */
