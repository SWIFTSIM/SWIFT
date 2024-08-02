/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will J. Roper (w.roper@sussex.ac.uk)
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

/* Config */
#include <config.h>

/* Local includes */
#include "cell.h"
#include "engine.h"
#include "space.h"
#include "threadpool.h"
#include "zoom.h"

/**
 * @brief Recursively split a cell.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param ti_current The current time step.
 * @param tpid The thread id.
 */
void zoom_void_split_recursive(struct space *s, struct cell *c,
                               const integertime_t ti_current,
                               const short int tpid) {

  const int depth = c->depth;
  int maxdepth = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;

  /* Set the top level void cell tpid. Doing it here ensures top level void
   * cells have the same tpid as their progeny. */
  if (depth == 0) c->tpid = tpid;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we haven't found a void cell with particles. */
  if (c->subtype == cell_subtype_void && c->grav.count > 0)
    error(
        "Trying to split a Void with particles! "
        "(c->type=%s, c->subtype=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);

  /* Ensure we haven't got an empty cell. */
  if (c->subtype == cell_subtype_empty)
    error(
        "Trying to split an empty cell which "
        "shouldn't have particles! (c->type=%s, c->type=%s)",
        cellID_names[c->type], subcellID_names[c->subtype]);
#endif

  /* If the depth is too large, we have a problem and should stop. */
  if (depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting the void cells, aborting. "
        "This is most likely due to having the zoom region too deep within the "
        "background cells. Try using buffer cells or increasing "
        "region_buffer_cell_ratio (depth=%d)",
        space_cell_maxdepth, depth);
  }

  /* Construct the progeny ready to populate with particles and multipoles (if
   * doing gravity). */
  space_construct_progeny(s, c, tpid);

  for (int k = 0; k < 8; k++) {

    /* Get the progenitor */
    struct cell *cp = c->progeny[k];

#ifdef SWIFT_DEBUG_CHECKS
    if (cp->subtype == cell_subtype_void &&
        cp->width[0] <= s->zoom_props->width[0])
      error(
          "We have a zoom cell labelled as a void cell! We have gone too "
          "deep in the zoom cell tree, this could be because background "
          "cells are comprable in size to the zoom cells."
          " (cp->type=%s, cp->subtype=%s, cp->width[0]=%f, cp->depth=%d,"
          " s->zoom_props->width[0]=%f, zoom_props->zoom_cell_depth=%d)",
          cellID_names[cp->type], subcellID_names[cp->subtype], cp->width[0],
          cp->depth, s->zoom_props->width[0], s->zoom_props->zoom_cell_depth);
#endif

    /* If the next level progeny is at the zoom level then we need to
     * link the zoom cells in as the progeny of the void sub-cell. */
    if (cp->depth == s->zoom_props->zoom_cell_depth - 1) {

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that the widths are right. */
      if (fabs((cp->width[0] / 2) - s->zoom_props->width[0]) >
          (0.001 * s->zoom_props->width[0]))
        error(
            "The width of the zoom cell is not half the width of the void "
            "cell were about to link to! (cp->width[0]=%f, "
            "s->zoom_props->width[0]=%f)",
            cp->width[0] / 2, s->zoom_props->width[0]);
#endif

      zoom_link_void_leaves(s, cp, ti_current);

    } else {

      /* Recurse */
      zoom_void_split_recursive(s, cp, ti_current, tpid);

      /* Increase the depth */
      maxdepth = max(maxdepth, cp->maxdepth);

      /* Update the gravity time step properties. */
      ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
      ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);
    }
  }

  /* Update the properties of the void cell. */
  c->maxdepth = maxdepth;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;

  /* Deal with the multipole */
  if (s->with_self_gravity) {
    space_populate_multipole(c);
  }
}

/**
 * @brief Split particles between cells of the void cell hierarchy.
 *
 * This has to be done after proxy exchange to ensure we have all the
 * foreign multipoles, so is separated from all other splitting done in the
 * usual space_split.
 *
 *
 * @param s The #space.
 * @param verbose Are we talkative?
 */
void zoom_void_space_split(struct space *s, int verbose) {

#ifdef SWIFT_DEBUG_CHECKS
  /* We should never get here when not running with a zoom region. */
  if (!s->with_zoom_region) {
    error("zoom_void_space_split called when running without a zoom region.");
  }
#endif

  const ticks tic = getticks();

  /* Unpack some useful information. */
  struct cell *cells_top = s->cells_top;
  int *void_cells_top = s->zoom_props->void_cells_top;
  int nr_void_cells = s->zoom_props->nr_void_cells;
  const integertime_t ti_current = s->e->ti_current;

  /* Create the void cell trees and populate their multipoles. This is only
   * a handful of cells so no threadpool. */

  /* Loop over the void cells */
  for (int ind = 0; ind < nr_void_cells; ind++) {
    struct cell *c = &cells_top[void_cells_top[ind]];
    zoom_void_split_recursive(s, c, ti_current, /*tpid*/ 0);
  }

  if (verbose)
    message("Void cell tree and multipole construction took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure all cells are linked into the tree. */
  int notlinked = 0;
  int nr_gparts_in_zoom = 0;
  for (int k = 0; k < s->zoom_props->nr_zoom_cells; k++) {
    nr_gparts_in_zoom += s->multipoles_top[k].m_pole.num_gpart;
    if (cells_top[k].void_parent == NULL) notlinked++;
  }
  if (notlinked > 0)
    error("%d zoom cells are not linked into a void cell tree!", notlinked);

  /* Check all void cells have void children. */

  /* Compare the number of particles in the void multipole and zoom cells. */
  int nr_gparts_in_void = 0;
  for (int i = 0; i < nr_void_cells; i++)
    nr_gparts_in_void +=
        s->multipoles_top[s->zoom_props->void_cells_top[i]].m_pole.num_gpart;

  if (nr_gparts_in_void != nr_gparts_in_zoom)
    error(
        "Number of gparts is in consistent between zoom cells and "
        "void multipole (nr_gparts_in_void=%d, nr_gparts_in_zoom=%d)",
        nr_gparts_in_void, nr_gparts_in_zoom);

#endif
}
