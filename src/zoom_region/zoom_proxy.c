/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Will Roper (w.roper@sussex.ac.uk)
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

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "proxy.h"
#include "zoom.h"

/**
 * @brief Handle proxy creation for zoom cells inside void cells.
 *
 * This helper function processes the case where one or both cells are void,
 * and we need to create proxies for the actual zoom cells nested within them.
 *
 * @param e The #engine.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param proxy_type The type of proxy needed.
 * @param nodeID The current node ID.
 */
#ifdef WITH_MPI
static void zoom_get_void_cell_proxies(struct engine *e, struct cell *ci,
                                       struct cell *cj, const int proxy_type,
                                       const int nodeID) {
  const struct space *s = e->s;
  const int ci_is_void = (ci->subtype == cell_subtype_void);
  const int cj_is_void = (cj->subtype == cell_subtype_void);

  /* Both cells are void: need to find zoom cells in both */
  if (ci_is_void && cj_is_void) {
    for (int zid = 0; zid < s->zoom_props->nr_zoom_cells; zid++) {
      struct cell *zi = &s->cells_top[zid];
      if (!cell_is_inside(zi, ci)) continue;

      for (int zjd = 0; zjd < s->zoom_props->nr_zoom_cells; zjd++) {
        struct cell *zj = &s->cells_top[zjd];
        if (!cell_is_inside(zj, cj)) continue;

        /* Early abort (both same node or both foreign) */
        if ((zi->nodeID == nodeID && zj->nodeID == nodeID) ||
            (zi->nodeID != nodeID && zj->nodeID != nodeID))
          continue;

        engine_add_proxy(e, zi, zj, proxy_type);
      }
    }
  }
  /* Only ci is void: find zoom cells in ci to pair with cj */
  else if (ci_is_void) {
    for (int zid = 0; zid < s->zoom_props->nr_zoom_cells; zid++) {
      struct cell *zi = &s->cells_top[zid];
      if (!cell_is_inside(zi, ci)) continue;

      /* Only care about foreign zoom cells */
      if (zi->nodeID == nodeID) continue;

      engine_add_proxy(e, zi, cj, proxy_type);
    }
  }
  /* Only cj is void: find zoom cells in cj to pair with ci */
  else {
    for (int zjd = 0; zjd < s->zoom_props->nr_zoom_cells; zjd++) {
      struct cell *zj = &s->cells_top[zjd];
      if (!cell_is_inside(zj, cj)) continue;

      /* Only care about foreign zoom cells */
      if (zj->nodeID == nodeID) continue;

      engine_add_proxy(e, ci, zj, proxy_type);
    }
  }
}
#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies.
 *
 * This is the zoom specific form of engine_makeproxies. Unlike the uniform
 * box proxies we don't ahead of time know exactly which cells will end up
 * having direct interactions. Instead we just need to make sure all that
 * could interact do have a proxy, so we use a brute force loop over all
 * possible pairs of cells. This will account for any possible inter-cell grid
 * interactions.
 *
 * @param e The #engine.
 */
void zoom_engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const struct space *s = e->s;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Loop over all cell pairs once */
  for (int cid = s->zoom_props->bkg_cell_offset; cid < s->nr_cells; cid++) {

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Loop over the prospective neighbours */
    for (int cjd = cid + 1; cjd < s->nr_cells; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

      /* Check if we have at least one void cell */
      const int ci_is_void = (ci->subtype == cell_subtype_void);
      const int cj_is_void = (cj->subtype == cell_subtype_void);
      const int has_void = ci_is_void || cj_is_void;

      /* For non-void pairs, apply standard node checks early */
      if (!has_void) {
        /* Early abort (both same node) */
        if (ci->nodeID == nodeID && cj->nodeID == nodeID) continue;

        /* Early abort (both foreign node) */
        if (ci->nodeID != nodeID && cj->nodeID != nodeID) continue;
      }

      /* Calculate the maximum distance based on the diagonal distance of the
       * pair */
      const double ir_diag2 = ci->width[0] * ci->width[0] +
                              ci->width[1] * ci->width[1] +
                              ci->width[2] * ci->width[2];
      const double ir_diag = 0.5 * sqrt(ir_diag2);
      const double jr_diag2 = cj->width[0] * cj->width[0] +
                              cj->width[1] * cj->width[1] +
                              cj->width[2] * cj->width[2];
      const double jr_diag = 0.5 * sqrt(jr_diag2);

      /* Calculate the maximum distance between the cells */
      const double r_max = ir_diag + jr_diag;

      /* Get the proxy type */
      int proxy_type = engine_get_proxy_type(e, ci, cj, r_max);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Handle void cells: delegate to zoom cells inside */
      if (has_void) {
        zoom_get_void_cell_proxies(e, ci, cj, proxy_type, nodeID);
      }
      /* Both cells are non-void: add proxy directly */
      else {
        engine_add_proxy(e, ci, cj, proxy_type);
      }
    }
  }

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
