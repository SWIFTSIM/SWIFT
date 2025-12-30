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
 * @brief Structure to hold a cell pair and the proxy type.
 */
struct cell_type_pair {
  struct cell *ci, *cj;
  int type;
};

#ifdef WITH_MPI
/**
 * @brief Get proxies for void cell pairs by checking nested zoom cells.
 *
 * This function handles cases where one or both cells are void (background
 * cells), requiring us to check the nested zoom cells within them.
 *
 * @param e The #engine.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param proxy_type The type of proxy needed.
 * @param nodeID The current node ID.
 */
static void zoom_get_void_cell_proxies(struct engine *e, struct cell *ci,
                                       struct cell *cj, const int proxy_type,
                                       const int nodeID) {

  const struct space *s = e->s;

  /* Which are void? */
  const int ci_is_void = (ci->subtype == cell_subtype_void);
  const int cj_is_void = (cj->subtype == cell_subtype_void);

  /* Both are void - need to check zoom cells in both */
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
  /* Only ci is void - find zoom cells in ci to pair with cj */
  else if (ci_is_void) {
    for (int zid = 0; zid < s->zoom_props->nr_zoom_cells; zid++) {
      struct cell *zi = &s->cells_top[zid];
      if (!cell_is_inside(zi, ci)) continue;

      /* Early abort (both same node or both foreign) */
      if ((zi->nodeID == nodeID && cj->nodeID == nodeID) ||
          (zi->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      engine_add_proxy(e, zi, cj, proxy_type);
    }
  }
  /* Only cj is void - find zoom cells in cj to pair with ci */
  else if (cj_is_void) {
    for (int zjd = 0; zjd < s->zoom_props->nr_zoom_cells; zjd++) {
      struct cell *zj = &s->cells_top[zjd];
      if (!cell_is_inside(zj, cj)) continue;

      /* Early abort (both same node or both foreign) */
      if ((ci->nodeID == nodeID && zj->nodeID == nodeID) ||
          (ci->nodeID != nodeID && zj->nodeID != nodeID))
        continue;

      engine_add_proxy(e, ci, zj, proxy_type);
    }
  }
}
#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies.
 *
 * This is the zoom specific form of engine_makeproxies. We use structured grid
 * walking over the background cells, with special handling for void cells that
 * checks nested zoom cells.
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
  struct cell *cells = s->cells_top;

  const int bkg_cdim[3] = {s->zoom_props->bkg_cdim[0],
                           s->zoom_props->bkg_cdim[1],
                           s->zoom_props->bkg_cdim[2]};
  const int periodic = s->periodic;

  /* Get some info about the physics */
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;

  /* Prepare the proxies and the proxy index */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Compute maximum distance for background cell proxy interactions */
  const double bkg_cell_width[3] = {
      cells[s->zoom_props->bkg_cell_offset].width[0],
      cells[s->zoom_props->bkg_cell_offset].width[1],
      cells[s->zoom_props->bkg_cell_offset].width[2]};

  const double bkg_r_diag2 = bkg_cell_width[0] * bkg_cell_width[0] +
                             bkg_cell_width[1] * bkg_cell_width[1] +
                             bkg_cell_width[2] * bkg_cell_width[2];
  const double bkg_r_diag = 0.5 * sqrt(bkg_r_diag2);
  const double bkg_r_max = 2 * bkg_r_diag;

  /* Compute how many background cells away we need to walk */
  int bkg_delta_cells = 1; /* hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * bkg_r_max / theta_crit;
    bkg_delta_cells =
        (int)(distance / cells[s->zoom_props->bkg_cell_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int bkg_delta_m = bkg_delta_cells;
  int bkg_delta_p = bkg_delta_cells;

  /* Special case where every background cell is in range of every other one */
  if (bkg_delta_cells >= bkg_cdim[0] / 2) {
    bkg_delta_m = bkg_cdim[0] / 2;
    bkg_delta_p = bkg_cdim[0] / 2;
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d background cells away (delta_m=%d "
        "delta_p=%d)",
        bkg_delta_cells, bkg_delta_m, bkg_delta_p);

  /* Loop over each background cell using structured grid walk */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the cell ID */
        const int cid = cell_getid_offset(
            bkg_cdim, s->zoom_props->bkg_cell_offset, i, j, k);

        struct cell *ci = &cells[cid];

        /* Loop over all its neighbours in range */
        for (int ii = -bkg_delta_m; ii <= bkg_delta_p; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= bkg_cdim[0])) continue;
          iii = (iii + bkg_cdim[0]) % bkg_cdim[0];

          for (int jj = -bkg_delta_m; jj <= bkg_delta_p; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= bkg_cdim[1])) continue;
            jjj = (jjj + bkg_cdim[1]) % bkg_cdim[1];

            for (int kk = -bkg_delta_m; kk <= bkg_delta_p; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= bkg_cdim[2])) continue;
              kkk = (kkk + bkg_cdim[2]) % bkg_cdim[2];

              /* Get the cell ID */
              const int cjd = cell_getid_offset(
                  bkg_cdim, s->zoom_props->bkg_cell_offset, iii, jjj, kkk);

              /* Early abort */
              if (cid >= cjd) continue;

              struct cell *cj = &cells[cjd];

              /* Check if we have void cells */
              const int has_void = (ci->subtype == cell_subtype_void ||
                                    cj->subtype == cell_subtype_void);

              /* For non-void pairs, apply standard node checks early */
              if (!has_void) {
                /* Early abort (both same node or both foreign) */
                if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                    (ci->nodeID != nodeID && cj->nodeID != nodeID))
                  continue;
              }

              /* Calculate maximum distance for this pair */
              const double ir_diag2 = ci->width[0] * ci->width[0] +
                                      ci->width[1] * ci->width[1] +
                                      ci->width[2] * ci->width[2];
              const double ir_diag = 0.5 * sqrt(ir_diag2);
              const double jr_diag2 = cj->width[0] * cj->width[0] +
                                      cj->width[1] * cj->width[1] +
                                      cj->width[2] * cj->width[2];
              const double jr_diag = 0.5 * sqrt(jr_diag2);
              const double pair_r_max = ir_diag + jr_diag;

              /* Get the proxy type */
              int proxy_type = engine_get_proxy_type(e, ci, cj, pair_r_max);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* If any void cells are involved, use special handling */
              if (has_void) {
                zoom_get_void_cell_proxies(e, ci, cj, proxy_type, nodeID);
              } else {
                /* Regular cells - add proxy directly */
                engine_add_proxy(e, ci, cj, proxy_type);
              }
            }
          }
        }
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Validate the proxies */

  /* Loop over the proxies and add the send tasks, which also generates the
   * cell tags for super-cells. */
  int max_num_send_cells = 0;
  for (int pid = 0; pid < e->nr_proxies; pid++)
    max_num_send_cells += e->proxies[pid].nr_cells_out;
  struct cell_type_pair *send_cell_type_pairs = NULL;
  if ((send_cell_type_pairs = (struct cell_type_pair *)malloc(
           sizeof(struct cell_type_pair) * max_num_send_cells)) == NULL)
    error("Failed to allocate temporary cell pointer list.");
  int num_send_cells = 0;

  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    struct proxy *p = &e->proxies[pid];

    for (int k = 0; k < p->nr_cells_out; k++) {
      send_cell_type_pairs[num_send_cells].ci = p->cells_out[k];
      send_cell_type_pairs[num_send_cells].cj = p->cells_in[k];
      send_cell_type_pairs[num_send_cells++].type = p->cells_out_type[k];
    }
  }

  /* Ensure we have valid cells */
  for (int n = 0; n < num_send_cells; n++) {
    struct cell *ci = send_cell_type_pairs[n].ci;
    struct cell *cj = send_cell_type_pairs[n].cj;
    if (ci == NULL) {
      error(
          "ci is NULL. cj=%p (%s/%s, cj->depth=%d, cj->count=%d, "
          "cj->nodeID=%d)",
          (void *)cj, cellID_names[cj->type], subcellID_names[cj->subtype],
          cj->depth, cj->grav.count, cj->nodeID);
    }
    if (cj == NULL) {
      error(
          "cj is NULL. ci=%p (%s/%s, ci->depth=%d, ci->count=%d, "
          "ci->nodeID=%d)",
          (void *)ci, cellID_names[ci->type], subcellID_names[ci->subtype],
          ci->depth, ci->grav.count, ci->nodeID);
    }
  }
  free(send_cell_type_pairs);
#endif

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
