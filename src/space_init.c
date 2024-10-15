/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include <config.h>

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "adaptive_softening.h"
#include "black_holes.h"
#include "chemistry.h"
#include "engine.h"
#include "gravity.h"
#include "mhd.h"
#include "rt.h"
#include "sink.h"
#include "star_formation.h"
#include "stars.h"
#include "threadpool.h"
#include "tracers.h"

void space_init_parts_mapper(void *restrict map_data, int count,
                             void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;
  const struct hydro_space *restrict hs = &e->s->hs;
  const int with_cosmology = (e->policy & engine_policy_cosmology);

  size_t ind = parts - e->s->parts;
  struct xpart *restrict xparts = e->s->xparts + ind;

  for (int k = 0; k < count; k++) {
    hydro_init_part(&parts[k], hs);
    adaptive_softening_init_part(&parts[k]);
    mhd_init_part(&parts[k]);
    black_holes_init_potential(&parts[k].black_holes_data);
    chemistry_init_part(&parts[k], e->chemistry);
    rt_init_part(&parts[k]);
    rt_reset_part(&parts[k], e->cosmology);
    star_formation_init_part(&parts[k], e->star_formation);
    tracers_after_init(&parts[k], &xparts[k], e->internal_units,
                       e->physical_constants, with_cosmology, e->cosmology,
                       e->hydro_properties, e->cooling_func, e->time);
    sink_init_part(&parts[k], e->sink_properties);
  }
}

/**
 * @brief Calls the #part initialisation function on all particles in the space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_parts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_init_parts_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   s->e);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_gparts_mapper(void *restrict map_data, int count,
                              void *restrict extra_data) {

  struct gpart *gparts = (struct gpart *)map_data;
  for (int k = 0; k < count; k++) gravity_init_gpart(&gparts[k]);
}

/**
 * @brief Calls the #gpart initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_gparts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_gparts > 0)
    threadpool_map(&s->e->threadpool, space_init_gparts_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_sparts_mapper(void *restrict map_data, int scount,
                              void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  for (int k = 0; k < scount; k++) {
    stars_init_spart(&sparts[k]);
    rt_init_spart(&sparts[k]);
    rt_reset_spart(&sparts[k]);
  }
}

/**
 * @brief Calls the #spart initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_sparts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, space_init_sparts_mapper, s->sparts,
                   s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_bparts_mapper(void *restrict map_data, int bcount,
                              void *restrict extra_data) {

  struct bpart *restrict bparts = (struct bpart *)map_data;
  for (int k = 0; k < bcount; k++) black_holes_init_bpart(&bparts[k]);
}

/**
 * @brief Calls the #bpart initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_bparts(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_bparts > 0)
    threadpool_map(&s->e->threadpool, space_init_bparts_mapper, s->bparts,
                   s->nr_bparts, sizeof(struct bpart),
                   threadpool_auto_chunk_size, /*extra_data=*/NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_sinks_mapper(void *restrict map_data, int sink_count,
                             void *restrict extra_data) {

  struct sink *restrict sinks = (struct sink *)map_data;
  for (int k = 0; k < sink_count; k++) sink_init_sink(&sinks[k]);
}

/**
 * @brief Calls the #sink initialisation function on all particles in the
 * space.
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void space_init_sinks(struct space *s, int verbose) {

  const ticks tic = getticks();

  if (s->nr_sinks > 0)
    threadpool_map(&s->e->threadpool, space_init_sinks_mapper, s->sinks,
                   s->nr_sinks, sizeof(struct sink), threadpool_auto_chunk_size,
                   /*extra_data=*/NULL);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
