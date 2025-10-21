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

  /* Loop over the cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip void cells (these will never need a proxy). */
    if (ci->subtype == cell_subtype_void) continue;

    /* Loop over the prospective neighbours. */
    for (int cjd = cid + 1; cjd < s->nr_cells; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

      /* Skip void cells (these will never need a proxy). */
      if (cj->subtype == cell_subtype_void) continue;

      /* Early abort (both same node) -> Nigel is happy */
      if (ci->nodeID == nodeID && cj->nodeID == nodeID) continue;

      /* Early abort (both foreign node) -> Nigel is angry */
      if (ci->nodeID != nodeID && cj->nodeID != nodeID) continue;

      /* We might need a proxy, one cell is foreign (Like Nigel and his wife).*/

      /* Calculate the maximum distance based on the diagonal distance of the
       * pair. */
      const double ir_diag2 = ci->width[0] * ci->width[0] +
                              ci->width[1] * ci->width[1] +
                              ci->width[2] * ci->width[2];
      const double ir_diag = 0.5 * sqrt(ir_diag2);
      const double jr_diag2 = cj->width[0] * cj->width[0] +
                              cj->width[1] * cj->width[1] +
                              cj->width[2] * cj->width[2];
      const double jr_diag = 0.5 * sqrt(jr_diag2);

      /* Calculate the maximum distance between the cells. */
      const double r_max = 2 * max(ir_diag, jr_diag);

      /* Get the proxy type. We only need to do the direct check if both
       * cells are the same type. Note, the cdim is only used if
       * icdim == jcdim and we're doing a direct check. */
      int proxy_type = engine_get_proxy_type(e, ci, cj, r_max);

      if (cid == 128 && cj->type == cell_type_zoom)
        message(
            "Proxy between cell %d (node %d) and cell %d (node %d): "
            "type=%d",
            cid, ci->nodeID, cjd, cj->nodeID, proxy_type);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Ok, we need to add a proxy. */
      engine_add_proxy(e, ci, cj, proxy_type);
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
