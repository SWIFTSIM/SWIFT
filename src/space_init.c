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
#include "part_init.h"

void space_init_parts_mapper(void *restrict map_data, int count,
                             void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;
  size_t ind = parts - e->s->parts;
  struct xpart *restrict xparts = e->s->xparts + ind;

  for (int k = 0; k < count; k++) {
    part_init(&parts[k], &xparts[k], e);
    rt_reset_part(&parts[k], e->cosmology);
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
  const struct engine *restrict e = (struct engine *)extra_data;

  for (int k = 0; k < count; k++) {
    gpart_init(&gparts[k], e);
  }
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
                   threadpool_auto_chunk_size, s->e);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_sparts_mapper(void *restrict map_data, int scount,
                              void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  for (int k = 0; k < scount; k++) {
    spart_init(&sparts[k], e);
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
                   threadpool_auto_chunk_size, s->e);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_bparts_mapper(void *restrict map_data, int bcount,
                              void *restrict extra_data) {

  struct bpart *restrict bparts = (struct bpart *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  for (int k = 0; k < bcount; k++) {
    bpart_init(&bparts[k], e);
  }
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
                   threadpool_auto_chunk_size, s->e);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_init_sinks_mapper(void *restrict map_data, int sink_count,
                             void *restrict extra_data) {

  struct sink *restrict sinks = (struct sink *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  for (int k = 0; k < sink_count; k++) {
    sink_init(&sinks[k], e);
  }
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
                   s->e);
  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
