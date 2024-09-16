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

#ifdef WITH_MPI
/** @brief Recursive function to make proxies.
 *
 * If a void cell is passed then it will be split until we reach the zoom
 * region, otherwise we get all the information we need to test if we need
 * a proxy and then make it if we need to.
 *
 *
 * @param e The #engine.
 * @param ci The first cell.
 * @param cj The second cell.
 * @param max_mesh_dist2 The maximum mesh distance squared.
 * @param theta_crit The critical opening angle.
 * @param with_hydro Whether hydro is enabled.
 * @param with_gravity Whether gravity is enabled.
 */
static void zoom_make_proxies_recursive(struct engine *e, struct cell *ci,
                                        struct cell *cj, double max_mesh_dist2,
                                        double theta_crit, int with_hydro,
                                        int with_gravity) {

  /* Unpack the space. */
  const struct space *s = e->s;

  /* Get the nodeID for this rank. */
  const int nodeID = e->nodeID;

  /* Get the proxies. */
  struct proxy *proxies = e->proxies;

  /* If ci or cj are void cells we need to recurse. */
  if (ci->subtype == cell_subtype_void) {
    for (int i = 0; i < s->zoom_props->nr_zoom_cells; i++) {
      zoom_make_proxies_recursive(e, &s->zoom_props->zoom_cells_top[i], cj,
                                  max_mesh_dist2, theta_crit, with_hydro,
                                  with_gravity);
    }
  } else if (cj->subtype == cell_subtype_void) {
    for (int j = 0; j < s->zoom_props->nr_zoom_cells; j++) {
      zoom_make_proxies_recursive(e, ci, &s->zoom_props->zoom_cells_top[j],
                                  max_mesh_dist2, theta_crit, with_hydro,
                                  with_gravity);
    }
  } else {

    /* Early abort (both same node) -> Nigel is happy */
    if (ci->nodeID == nodeID && cj->nodeID == nodeID) return;

    /* Early abort (both foreign node) -> Nigel is angry */
    if (ci->nodeID != nodeID && cj->nodeID != nodeID) return;

    /* We might need a proxy, one cell is foreign (Like Nigel and his wife).*/

    /* Get the right cdims and bounds for this pair. */
    int icdim[3], jcdim[3];
    double ilower_bounds[3];
    double jlower_bounds[3];
    int ioffset, joffset;
    if (ci->type == cell_type_zoom) {
      icdim[0] = s->zoom_props->cdim[0];
      icdim[1] = s->zoom_props->cdim[1];
      icdim[2] = s->zoom_props->cdim[2];
      ilower_bounds[0] = s->zoom_props->region_lower_bounds[0];
      ilower_bounds[1] = s->zoom_props->region_lower_bounds[1];
      ilower_bounds[2] = s->zoom_props->region_lower_bounds[2];
      ioffset = 0;
    } else if (ci->type == cell_type_buffer) {
      icdim[0] = s->zoom_props->buffer_cdim[0];
      icdim[1] = s->zoom_props->buffer_cdim[1];
      icdim[2] = s->zoom_props->buffer_cdim[2];
      ilower_bounds[0] = s->zoom_props->buffer_lower_bounds[0];
      ilower_bounds[1] = s->zoom_props->buffer_lower_bounds[1];
      ilower_bounds[2] = s->zoom_props->buffer_lower_bounds[2];
      ioffset = s->zoom_props->buffer_cell_offset;
    } else {
      icdim[0] = s->cdim[0];
      icdim[1] = s->cdim[1];
      icdim[2] = s->cdim[2];
      ilower_bounds[0] = 0.0;
      ilower_bounds[1] = 0.0;
      ilower_bounds[2] = 0.0;
      ioffset = s->zoom_props->bkg_cell_offset;
    }
    if (cj->type == cell_type_zoom) {
      jcdim[0] = s->zoom_props->cdim[0];
      jcdim[1] = s->zoom_props->cdim[1];
      jcdim[2] = s->zoom_props->cdim[2];
      jlower_bounds[0] = s->zoom_props->region_lower_bounds[0];
      jlower_bounds[1] = s->zoom_props->region_lower_bounds[1];
      jlower_bounds[2] = s->zoom_props->region_lower_bounds[2];
      joffset = 0;
    } else if (cj->type == cell_type_buffer) {
      jcdim[0] = s->zoom_props->buffer_cdim[0];
      jcdim[1] = s->zoom_props->buffer_cdim[1];
      jcdim[2] = s->zoom_props->buffer_cdim[2];
      jlower_bounds[0] = s->zoom_props->buffer_lower_bounds[0];
      jlower_bounds[1] = s->zoom_props->buffer_lower_bounds[1];
      jlower_bounds[2] = s->zoom_props->buffer_lower_bounds[2];
      joffset = s->zoom_props->buffer_cell_offset;
    } else {
      jcdim[0] = s->cdim[0];
      jcdim[1] = s->cdim[1];
      jcdim[2] = s->cdim[2];
      jlower_bounds[0] = 0.0;
      jlower_bounds[1] = 0.0;
      jlower_bounds[2] = 0.0;
      joffset = s->zoom_props->bkg_cell_offset;
    }

    /* Get the i, j, k coordinates for each cell. */
    const int i = (ci->loc[0] - ilower_bounds[0]) / ci->width[0];
    const int j = (ci->loc[1] - ilower_bounds[1]) / ci->width[1];
    const int k = (ci->loc[2] - ilower_bounds[2]) / ci->width[2];
    const int ii = (cj->loc[0] - jlower_bounds[0]) / cj->width[0];
    const int jj = (cj->loc[1] - jlower_bounds[1]) / cj->width[1];
    const int kk = (cj->loc[2] - jlower_bounds[2]) / cj->width[2];

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
    const double r_max = ir_diag + jr_diag;

    /* Get the indices for each cell. */
    const int cid = cell_getid_offset(icdim, ioffset, i, j, k);
    const int cjd = cell_getid_offset(jcdim, joffset, ii, jj, kk);

    /* Get the proxy type. We only need to do the direct check if both
     * cells are the same type. Note, the cdim is only used if
     * icdim == jcdim and we're doing a direct check. */
    int proxy_type = engine_get_proxy_type(
        s->cells_top, i, j, k, ii, jj, kk, icdim,
        (ci->type == cell_type_zoom && cj->type == cell_type_zoom) ? with_hydro
                                                                   : 0,
        with_gravity, cid, cjd, s->dim,
        (ci->type != cell_type_bkg && cj->type != cell_type_bkg) ? 0
                                                                 : s->periodic,
        r_max, max_mesh_dist2, theta_crit,
        /*do_direct_check*/ ci->type == cj->type);

    /* Abort if not in range at all */
    if (proxy_type == proxy_cell_type_none) return;

    /* Ok, we need to add a proxy. */
    engine_add_proxy(e, s->cells_top, proxies, cid, cjd, proxy_type, e->nodeID);
  }
}
#endif

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
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;

  /* Get some info about the physics */
  int with_hydro = (e->policy & engine_policy_hydro);
  int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Loop over the cells */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Loop over the prospective neighbours. */
    for (int cjd = cid + 1; cjd < s->nr_cells; cjd++) {

      /* Test these cells. */
      zoom_make_proxies_recursive(e, &cells[cid], &cells[cjd], max_mesh_dist2,
                                  theta_crit, with_hydro, with_gravity);
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
