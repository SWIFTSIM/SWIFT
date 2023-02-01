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
#include "cell.h"
#include "engine.h"

/* Some standard headers. */
#include <string.h>

/**
 * @brief Allocate memory for the extra particles used for on-the-fly creation.
 *
 * This rarely actually allocates memory. Most of the time, we convert
 * pre-allocated memory inot extra particles.
 *
 * This function also sets the extra particles' location to their top-level
 * cells. They can then be sorted into their correct memory position later on.
 *
 * @param s The current #space.
 * @param verbose Are we talkative?
 */
void space_allocate_extras(struct space *s, int verbose) {

  const int local_nodeID = s->e->nodeID;

  /* Anything to do here? (Abort if we don't want extras)*/
  if (space_extra_parts == 0 && space_extra_gparts == 0 &&
      space_extra_sparts == 0 && space_extra_bparts == 0 &&
      space_extra_sinks == 0)
    return;

  /* The top-level cells */
  const struct cell *cells = s->cells_top;
  const double half_cell_width[3] = {0.5 * cells[0].width[0],
                                     0.5 * cells[0].width[1],
                                     0.5 * cells[0].width[2]};

  /* The current number of particles (including spare ones) */
  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;
  size_t nr_bparts = s->nr_bparts;
  size_t nr_sinks = s->nr_sinks;

  /* The current number of actual particles */
  size_t nr_actual_parts = nr_parts - s->nr_extra_parts;
  size_t nr_actual_gparts = nr_gparts - s->nr_extra_gparts;
  size_t nr_actual_sparts = nr_sparts - s->nr_extra_sparts;
  size_t nr_actual_bparts = nr_bparts - s->nr_extra_bparts;
  size_t nr_actual_sinks = nr_sinks - s->nr_extra_sinks;

  /* The number of particles we allocated memory for (MPI overhead) */
  size_t size_parts = s->size_parts;
  size_t size_gparts = s->size_gparts;
  size_t size_sparts = s->size_sparts;
  size_t size_bparts = s->size_bparts;
  size_t size_sinks = s->size_sinks;

  int *local_cells = (int *)malloc(sizeof(int) * s->nr_cells);
  if (local_cells == NULL)
    error("Failed to allocate list of local top-level cells");

  /* List the local cells */
  size_t nr_local_cells = 0;
  for (int i = 0; i < s->nr_cells; ++i) {
    if (s->cells_top[i].nodeID == local_nodeID) {
      local_cells[nr_local_cells] = i;
      ++nr_local_cells;
    }
  }

  /* Number of extra particles we want for each type */
  const size_t expected_num_extra_parts = nr_local_cells * space_extra_parts;
  const size_t expected_num_extra_gparts = nr_local_cells * space_extra_gparts;
  const size_t expected_num_extra_sparts = nr_local_cells * space_extra_sparts;
  const size_t expected_num_extra_bparts = nr_local_cells * space_extra_bparts;
  const size_t expected_num_extra_sinks = nr_local_cells * space_extra_sinks;

  if (verbose) {
    message("Currently have %zd/%zd/%zd/%zd/%zd real particles.",
            nr_actual_parts, nr_actual_gparts, nr_actual_sinks,
            nr_actual_sparts, nr_actual_bparts);
    message("Currently have %zd/%zd/%zd/%zd/%zd spaces for extra particles.",
            s->nr_extra_parts, s->nr_extra_gparts, s->nr_extra_sinks,
            s->nr_extra_sparts, s->nr_extra_bparts);
    message(
        "Requesting space for future %zd/%zd/%zd/%zd/%zd "
        "part/gpart/sinks/sparts/bparts.",
        expected_num_extra_parts, expected_num_extra_gparts,
        expected_num_extra_sinks, expected_num_extra_sparts,
        expected_num_extra_bparts);
  }

  if (expected_num_extra_parts < s->nr_extra_parts)
    error(
        "Reduction in top-level cells number not handled. Please set a lower "
        "h_max or reduce the number of top level cells.");
  if (expected_num_extra_gparts < s->nr_extra_gparts)
    error(
        "Reduction in top-level cells number not handled. Please set a lower "
        "h_max or reduce the number of top level cells.");
  if (expected_num_extra_sparts < s->nr_extra_sparts)
    error(
        "Reduction in top-level cells number not handled. Please set a lower "
        "h_max or reduce the number of top level cells.");
  if (expected_num_extra_bparts < s->nr_extra_bparts)
    error(
        "Reduction in top-level cells number not handled. Please set a lower "
        "h_max or reduce the number of top level cells.");
  if (expected_num_extra_sinks < s->nr_extra_sinks)
    error(
        "Reduction in top-level cells number not handled. Please set a lower "
        "h_max or reduce the number of top level cells.");

  /* Do we have enough space for the extra gparts (i.e. we haven't used up any)
   * ? */
  if (nr_actual_gparts + expected_num_extra_gparts > nr_gparts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_gparts + expected_num_extra_gparts > size_gparts) {

      size_gparts = (nr_actual_gparts + expected_num_extra_gparts) *
                    engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating gparts array from %zd to %zd", s->size_gparts,
                size_gparts);

      /* Create more space for parts */
      struct gpart *gparts_new = NULL;
      if (swift_memalign("gparts", (void **)&gparts_new, gpart_align,
                         sizeof(struct gpart) * size_gparts) != 0)
        error("Failed to allocate new gpart data");
      memcpy(gparts_new, s->gparts, sizeof(struct gpart) * s->size_gparts);
      swift_free("gparts", s->gparts);
      s->gparts = gparts_new;

      /* Update the counter */
      s->size_gparts = size_gparts;

      /* We now need to correct all the pointers of the other particle arrays */
      part_relink_all_parts_to_gparts(gparts_new, s->nr_gparts, s->parts,
                                      s->sinks, s->sparts, s->bparts,
                                      &s->e->threadpool);
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_gparts; i < nr_actual_gparts + expected_num_extra_gparts;
         ++i) {
      bzero(&s->gparts[i], sizeof(struct gpart));
      s->gparts[i].time_bin = time_bin_not_created;
      s->gparts[i].type = swift_type_dark_matter;
      s->gparts[i].id_or_neg_offset = -1;
    }

    /* Put the spare particles in their correct cell */
    size_t local_cell_id = 0;
    int current_cell = local_cells[local_cell_id];
    int count_in_cell = 0;
#ifdef SWIFT_DEBUG_CHECKS
    size_t count_extra_gparts = 0;
#endif
    for (size_t i = 0; i < nr_actual_gparts + expected_num_extra_gparts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->gparts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->gparts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->gparts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->gparts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
#ifdef SWIFT_DEBUG_CHECKS
        count_extra_gparts++;
#endif
      }

      /* Once we have reached the number of extra gpart per cell, we move to the
       * next */
      if (count_in_cell == space_extra_gparts) {
        ++local_cell_id;

        if (local_cell_id == nr_local_cells) break;

        current_cell = local_cells[local_cell_id];
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_gparts != expected_num_extra_gparts)
      error("Constructed the wrong number of extra gparts (%zd vs. %zd)",
            count_extra_gparts, expected_num_extra_gparts);
#endif

    /* Update the counters */
    s->nr_gparts = nr_actual_gparts + expected_num_extra_gparts;
    s->nr_extra_gparts = expected_num_extra_gparts;
  }

  /* Do we have enough space for the extra parts (i.e. we haven't used up any) ?
   */
  if (nr_actual_parts + expected_num_extra_parts > nr_parts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_parts + expected_num_extra_parts > size_parts) {

      size_parts = (nr_actual_parts + expected_num_extra_parts) *
                   engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating parts array from %zd to %zd", s->size_parts,
                size_parts);

      /* Create more space for parts */
      struct part *parts_new = NULL;
      if (swift_memalign("parts", (void **)&parts_new, part_align,
                         sizeof(struct part) * size_parts) != 0)
        error("Failed to allocate new part data");
      memcpy(parts_new, s->parts, sizeof(struct part) * s->size_parts);
      swift_free("parts", s->parts);
      s->parts = parts_new;

      /* Same for xparts */
      struct xpart *xparts_new = NULL;
      if (swift_memalign("xparts", (void **)&xparts_new, xpart_align,
                         sizeof(struct xpart) * size_parts) != 0)
        error("Failed to allocate new xpart data");
      memcpy(xparts_new, s->xparts, sizeof(struct xpart) * s->size_parts);
      swift_free("xparts", s->xparts);
      s->xparts = xparts_new;

      /* Update the counter */
      s->size_parts = size_parts;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_parts; i < nr_actual_parts + expected_num_extra_parts;
         ++i) {
      bzero(&s->parts[i], sizeof(struct part));
      bzero(&s->xparts[i], sizeof(struct xpart));
      s->parts[i].time_bin = time_bin_not_created;
      s->parts[i].id = -42;
    }

    /* Put the spare particles in their correct cell */
    size_t local_cell_id = 0;
    int current_cell = local_cells[local_cell_id];
    int count_in_cell = 0;
#ifdef SWIFT_DEBUG_CHECKS
    size_t count_extra_parts = 0;
#endif
    for (size_t i = 0; i < nr_actual_parts + expected_num_extra_parts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->parts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->parts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->parts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->parts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
#ifdef SWIFT_DEBUG_CHECKS
        count_extra_parts++;
#endif
      }

      /* Once we have reached the number of extra part per cell, we move to the
       * next */
      if (count_in_cell == space_extra_parts) {
        ++local_cell_id;

        if (local_cell_id == nr_local_cells) break;

        current_cell = local_cells[local_cell_id];
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_parts != expected_num_extra_parts)
      error("Constructed the wrong number of extra parts (%zd vs. %zd)",
            count_extra_parts, expected_num_extra_parts);
#endif

    /* Update the counters */
    s->nr_parts = nr_actual_parts + expected_num_extra_parts;
    s->nr_extra_parts = expected_num_extra_parts;
  }

  /* Do we have enough space for the extra sinks (i.e. we haven't used up any)
   * ? */
  if (nr_actual_sinks + expected_num_extra_sinks > nr_sinks) {
    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_sinks + expected_num_extra_sinks > size_sinks) {

      size_sinks = (nr_actual_sinks + expected_num_extra_sinks) *
                   engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating sinks array from %zd to %zd", s->size_sinks,
                size_sinks);

      /* Create more space for parts */
      struct sink *sinks_new = NULL;
      if (swift_memalign("sinks", (void **)&sinks_new, sink_align,
                         sizeof(struct sink) * size_sinks) != 0)
        error("Failed to allocate new sink data");
      memcpy(sinks_new, s->sinks, sizeof(struct sink) * s->size_sinks);
      swift_free("sinks", s->sinks);
      s->sinks = sinks_new;

      /* Update the counter */
      s->size_sinks = size_sinks;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_sinks; i < nr_actual_sinks + expected_num_extra_sinks;
         ++i) {
      bzero(&s->sinks[i], sizeof(struct sink));
      s->sinks[i].time_bin = time_bin_not_created;
      s->sinks[i].id = -42;
    }

    /* Put the spare particles in their correct cell */
    size_t local_cell_id = 0;
    int current_cell = local_cells[local_cell_id];
    int count_in_cell = 0;
#ifdef SWIFT_DEBUG_CHECKS
    size_t count_extra_sinks = 0;
#endif
    for (size_t i = 0; i < nr_actual_sinks + expected_num_extra_sinks; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->sinks[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->sinks[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->sinks[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->sinks[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
#ifdef SWIFT_DEBUG_CHECKS
        count_extra_sinks++;
#endif
      }

      /* Once we have reached the number of extra sink per cell, we move to the
       * next */
      if (count_in_cell == space_extra_sinks) {
        ++local_cell_id;

        if (local_cell_id == nr_local_cells) break;

        current_cell = local_cells[local_cell_id];
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_sinks != expected_num_extra_sinks)
      error("Constructed the wrong number of extra sinks (%zd vs. %zd)",
            count_extra_sinks, expected_num_extra_sinks);
#endif

    /* Update the counters */
    s->nr_sinks = nr_actual_sinks + expected_num_extra_sinks;
    s->nr_extra_sinks = expected_num_extra_sinks;
  }

  /* Do we have enough space for the extra sparts (i.e. we haven't used up any)
   * ? */
  if (nr_actual_sparts + expected_num_extra_sparts > nr_sparts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_sparts + expected_num_extra_sparts > size_sparts) {

      size_sparts = (nr_actual_sparts + expected_num_extra_sparts) *
                    engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating sparts array from %zd to %zd", s->size_sparts,
                size_sparts);

      /* Create more space for parts */
      struct spart *sparts_new = NULL;
      if (swift_memalign("sparts", (void **)&sparts_new, spart_align,
                         sizeof(struct spart) * size_sparts) != 0)
        error("Failed to allocate new spart data");
      memcpy(sparts_new, s->sparts, sizeof(struct spart) * s->size_sparts);
      swift_free("sparts", s->sparts);
      s->sparts = sparts_new;

      /* Update the counter */
      s->size_sparts = size_sparts;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_sparts; i < nr_actual_sparts + expected_num_extra_sparts;
         ++i) {
      bzero(&s->sparts[i], sizeof(struct spart));
      s->sparts[i].time_bin = time_bin_not_created;
      s->sparts[i].id = -42;
    }

    /* Put the spare particles in their correct cell */
    size_t local_cell_id = 0;
    int current_cell = local_cells[local_cell_id];
    int count_in_cell = 0;
#ifdef SWIFT_DEBUG_CHECKS
    size_t count_extra_sparts = 0;
#endif
    for (size_t i = 0; i < nr_actual_sparts + expected_num_extra_sparts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->sparts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->sparts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->sparts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->sparts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
#ifdef SWIFT_DEBUG_CHECKS
        count_extra_sparts++;
#endif
      }

      /* Once we have reached the number of extra spart per cell, we move to the
       * next */
      if (count_in_cell == space_extra_sparts) {
        ++local_cell_id;

        if (local_cell_id == nr_local_cells) break;

        current_cell = local_cells[local_cell_id];
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_sparts != expected_num_extra_sparts)
      error("Constructed the wrong number of extra sparts (%zd vs. %zd)",
            count_extra_sparts, expected_num_extra_sparts);
#endif

    /* Update the counters */
    s->nr_sparts = nr_actual_sparts + expected_num_extra_sparts;
    s->nr_extra_sparts = expected_num_extra_sparts;
  }

  /* Do we have enough space for the extra bparts (i.e. we haven't used up any)
   * ? */
  if (nr_actual_bparts + expected_num_extra_bparts > nr_bparts) {

    /* Ok... need to put some more in the game */

    /* Do we need to reallocate? */
    if (nr_actual_bparts + expected_num_extra_bparts > size_bparts) {

      size_bparts = (nr_actual_bparts + expected_num_extra_bparts) *
                    engine_redistribute_alloc_margin;

      if (verbose)
        message("Re-allocating bparts array from %zd to %zd", s->size_bparts,
                size_bparts);

      /* Create more space for parts */
      struct bpart *bparts_new = NULL;
      if (swift_memalign("bparts", (void **)&bparts_new, bpart_align,
                         sizeof(struct bpart) * size_bparts) != 0)
        error("Failed to allocate new bpart data");
      memcpy(bparts_new, s->bparts, sizeof(struct bpart) * s->size_bparts);
      swift_free("bparts", s->bparts);
      s->bparts = bparts_new;

      /* Update the counter */
      s->size_bparts = size_bparts;
    }

    /* Turn some of the allocated spares into particles we can use */
    for (size_t i = nr_bparts; i < nr_actual_bparts + expected_num_extra_bparts;
         ++i) {
      bzero(&s->bparts[i], sizeof(struct bpart));
      s->bparts[i].time_bin = time_bin_not_created;
      s->bparts[i].id = -42;
    }

    /* Put the spare particles in their correct cell */
    size_t local_cell_id = 0;
    int current_cell = local_cells[local_cell_id];
    int count_in_cell = 0;
#ifdef SWIFT_DEBUG_CHECKS
    size_t count_extra_bparts = 0;
#endif
    for (size_t i = 0; i < nr_actual_bparts + expected_num_extra_bparts; ++i) {

#ifdef SWIFT_DEBUG_CHECKS
      if (current_cell == s->nr_cells)
        error("Cell counter beyond the maximal nr. cells.");
#endif

      if (s->bparts[i].time_bin == time_bin_not_created) {

        /* We want the extra particles to be at the centre of their cell */
        s->bparts[i].x[0] = cells[current_cell].loc[0] + half_cell_width[0];
        s->bparts[i].x[1] = cells[current_cell].loc[1] + half_cell_width[1];
        s->bparts[i].x[2] = cells[current_cell].loc[2] + half_cell_width[2];
        ++count_in_cell;
#ifdef SWIFT_DEBUG_CHECKS
        count_extra_bparts++;
#endif
      }

      /* Once we have reached the number of extra bpart per cell, we move to the
       * next */
      if (count_in_cell == space_extra_bparts) {
        ++local_cell_id;

        if (local_cell_id == nr_local_cells) break;

        current_cell = local_cells[local_cell_id];
        count_in_cell = 0;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (count_extra_bparts != expected_num_extra_bparts)
      error("Constructed the wrong number of extra bparts (%zd vs. %zd)",
            count_extra_bparts, expected_num_extra_bparts);
#endif

    /* Update the counters */
    s->nr_bparts = nr_actual_bparts + expected_num_extra_bparts;
    s->nr_extra_bparts = expected_num_extra_bparts;
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  if ((nr_gparts > 0 && nr_parts > 0) || (nr_gparts > 0 && nr_sparts > 0) ||
      (nr_gparts > 0 && nr_bparts > 0) || (nr_gparts > 0 && nr_sinks > 0))
    part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                      nr_parts, nr_gparts, nr_sinks, nr_sparts, nr_bparts,
                      verbose);
#endif

  /* Free the list of local cells */
  free(local_cells);
}
