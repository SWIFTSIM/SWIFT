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
#include "engine.h"
#include "proxy.h"
#include "zoom.h"

/**
 * @brief Create and fill the proxies.
 *
 * This is the zoom specific form of engine_makeproxies. Unlike the uniform
 * box proxies we don't ahead of time know exactly which cells will end up
 * having direct interactions. Instead we just need to make sure all that
 * could interact do have a proxy so we use a brute force loop over all
 * possible pairs of cells.
 *
 * @param e The #engine.
 */
void zoom_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = 2 * r_diag;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Loop over the cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Get the right cdim. */
    int icdim[3];
    if (cells[cid].type == cell_type_zoom) {
      icdim[0] = s->zoom_props->cdim[0];
      icdim[1] = s->zoom_props->cdim[1];
      icdim[2] = s->zoom_props->cdim[2];
    } else if (cells[cid].type == cell_type_buffer) {
      icdim[0] = s->zoom_props->buffer_cdim[0];
      icdim[1] = s->zoom_props->buffer_cdim[1];
      icdim[2] = s->zoom_props->buffer_cdim[2];
    } else {
      icdim[0] = s->cdim[0];
      icdim[1] = s->cdim[1];
      icdim[2] = s->cdim[2];
    }

    /* Void cells are never proxies. */
    if (cells[cid].type == cell_type_void) continue;

    /* Once we have finished all zoom cells (start of the array) we will no
     * longer need to consider hydro. */
    if (cells[cid].type != cell_type_zoom && with_hydro) with_hydro = 0;

    /* Loop over the prospective neighbours. */
    for (int cjd = cid + 1; cjd < s->nr_cells; cjd++) {

      /* Void cells are never proxies. */
      if (cells[cjd].type == cell_type_void) continue;

      /* Early abort (both same node) -> Nigel is happy */
      if (cells[cid].nodeID == nodeID && cells[cjd].nodeID == nodeID) continue;

      /* Early abort (both foreign node) -> Nigel is angry */
      if (cells[cid].nodeID != nodeID && cells[cjd].nodeID != nodeID) continue;

      /* We need a proxy, one cell is foreign (Like Nigel and his wife).*/

      /* Get the right cdim. */
      int jcdim[3];
      if (cells[cjd].type == cell_type_zoom) {
        jcdim[0] = s->zoom_props->cdim[0];
        jcdim[1] = s->zoom_props->cdim[1];
        jcdim[2] = s->zoom_props->cdim[2];
      } else if (cells[cjd].type == cell_type_buffer) {
        jcdim[0] = s->zoom_props->buffer_cdim[0];
        jcdim[1] = s->zoom_props->buffer_cdim[1];
        jcdim[2] = s->zoom_props->buffer_cdim[2];
      } else {
        jcdim[0] = s->cdim[0];
        jcdim[1] = s->cdim[1];
        jcdim[2] = s->cdim[2];
      }

      /* Integer indices of the cells in the top-level grid */
      const int i = cid / (icdim[1] * icdim[2]);
      const int j = (cid / icdim[2]) % icdim[1];
      const int k = cid % icdim[2];
      const int ii = cjd / (jcdim[1] * jcdim[2]);
      const int jj = (cjd / jcdim[2]) % jcdim[1];
      const int kk = cjd % jcdim[2];

      /* Get the proxy type. We only need to do the direct check if both
       * cells are the same type. Note, the cdim is only used if
       * icdim == jcdim and we're doing a direct check. */
      int proxy_type = engine_get_proxy_type(
          cells, i, j, k, ii, jj, kk, icdim, with_hydro, with_gravity, cid, cjd,
          dim, periodic, r_max, max_mesh_dist2, theta_crit,
          /*do_direct_check*/ cells[cid].type == cells[cjd].type);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Ok, we need to add a proxy. */
      engine_add_proxy(e, cells, proxies, cid, cjd, proxy_type, nodeID);
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
