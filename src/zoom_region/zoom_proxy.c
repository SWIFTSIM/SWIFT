/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will Roper (w.roper@sussex.ac.uk)
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

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "proxy.h"
#include "zoom.h"

#ifdef WITH_MPI

/**
 * @brief Apply gravity proxy criteria to a pair of cells.
 *
 * Note that cells must have equal widths, i.e. recursion is done in step.
 *
 * We only consider gravity in this function because bkg <-> bkg, and zoom <->
 * bkg pairs can only ever interact through gravity and there is therefore no
 * point in checking them for hydro pairs. Hydro considerations are done by in
 * zoom_get_gravity_proxy_type.
 *
 * @param e The #engine.
 * @param proxy_ci First cell used for geometric proxy criteria.
 * @param proxy_cj Second cell used for geometric proxy criteria.
 *
 * @return #proxy_cell_type_gravity if a proxy is required, otherwise zero.
 */
static int zoom_get_gravity_proxy_type(const struct engine *e,
                                       const struct cell *proxy_ci,
                                       const struct cell *proxy_cj) {

  const struct space *s = e->s;
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist2 = e->mesh->r_cut_max * e->mesh->r_cut_max;

  /* Nothing to do if gravity is not enabled. */
  if (!with_gravity) return proxy_cell_type_none;

  /* Derive the geometric criterion for a proxy. */
  const double r_diag2 = proxy_ci->width[0] * proxy_ci->width[0] +
                         proxy_ci->width[1] * proxy_ci->width[1] +
                         proxy_ci->width[2] * proxy_ci->width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);
  const double r_max = 2.0 * r_diag;

  /* Minimal distance between any two points in the cells. */
  const double min_dist_CoM2 =
      cell_min_dist2(proxy_ci, proxy_cj, s->periodic, s->dim);

  /* No gravity proxy is required beyond the mesh cutoff. */
  if (s->periodic && min_dist_CoM2 >= max_mesh_dist2) {
    return proxy_cell_type_none;
  }

  /* Check whether a gravity proxy is required. */
  if (!(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2)) {
    return proxy_cell_type_gravity;
  }

  /* We don't need a proxy */
  return proxy_cell_type_none;
}

/**
 * @brief Are cells direct zoom neighbours?
 *
 * @param e The #engine.
 * @param ci First zoom cell.
 * @param cj Second zoom cell.
 *
 * @return 1 if cells are adjacent and both zoom cells.
 */
static int zoom_is_direct_zoom_pair(const struct engine *e,
                                    const struct cell *ci,
                                    const struct cell *cj) {

  /* If these aren't both zoom cells its an immediate return. */
  if (ci->type != cell_type_zoom || cj->type != cell_type_zoom) {
    return 0
  }

  /* Get the zoom properties. */
  const struct zoom_region_properties *zp = e->s->zoom_props;

  /* Get the ijk integer cell positions */
  const int i =
      (int)((ci->loc[0] - zp->region_lower_bounds[0]) * zp->iwidth[0]);
  const int j =
      (int)((ci->loc[1] - zp->region_lower_bounds[1]) * zp->iwidth[1]);
  const int k =
      (int)((ci->loc[2] - zp->region_lower_bounds[2]) * zp->iwidth[2]);
  const int ii =
      (int)((cj->loc[0] - zp->region_lower_bounds[0]) * zp->iwidth[0]);
  const int jj =
      (int)((cj->loc[1] - zp->region_lower_bounds[1]) * zp->iwidth[1]);
  const int kk =
      (int)((cj->loc[2] - zp->region_lower_bounds[2]) * zp->iwidth[2]);

  /* Are the cells adjacent? */
  if (abs(ai - bi) <= 1 && abs(aj - bj) <= 1 && abs(ak - bk) <= 1) {
    return 1;
  } else {
    return 0;
  }
}

/**
 * @brief Apply proxy criteria to a pair of cells.
 *
 * This function checks first for a gravity proxy and then checks if cells are
 * adjacent and are thus hydro proxies (if running with hydro)
 *
 * @param e The #engine.
 * @param ci First zoom cell.
 * @param cj Second zoom cell.
 *
 * @return The proxy type.
 */
static int zoom_get_zoom_proxy_type(const struct engine *e,
                                    const struct cell *ci,
                                    const struct cell *cj) {

  /* Direct zoom neighbours? */
  int is_direct_neighbour = zoom_is_direct_zoom_pair(e, ci, cj);

  /* Do we have a gravity interaction between these pairs? */
  int proxy_type = zoom_get_gravity_proxy_type(e, ci, cj);

  /* Do we also have a hydro interaction? i.e. are these cells adjacent? */
  if ((e->policy & engine_policy_hydro) && is_direct_neighbour) {
    proxy_type |= (int)proxy_cell_type_hydro;
  }

  return proxy_type;
}

/**
 * @brief Create a proxy for a pair of cells if its needed.
 *
 */
static int zoom_make_proxy_for_pair(struct engine *e, const struct cell *ci,
                                    const struct cell *cj) {

  /* What node are we on? */
  const int nodeID = e->nodeID;

  /* Skip entirely local and entirely foreign pairs. */
  if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
      (ci->nodeID != nodeID && cj->nodeID != nodeID)) {
    return;
  }

  /* Get the proxy type for this pair. */
  const int proxy_type =
      engine_get_proxy_type(e, ci, i, j, k, cj, iii, jjj, kkk, bkg_r_max);

  /* No proxy needed? */
  if (proxy_type == proxy_cell_type_none) return;

  /* Ok, we need one. Add the proxy */
  engine_add_proxy(e, ci, cj, proxy_type);
}

/**
 * @brief Create proxies based on geometric pair recursion.
 *
 * This function will recurse through the pair of cells until we git the zoom
 * top level depth. Once at the zoom top level depth we will check what proxies
 * we need.
 *
 * Note that this could lead to us creating proxies for cells that may never
 * exist but we don't have enough information yet.
 *
 * @param e The #engine.
 * @param s The #space.
 * @param loc_i Lower corner of first region.
 * @param loc_j Lower corner of second region.
 * @param width Width of both regions.
 * @param depth Current recursion depth below the void cells.
 */
static void zoom_make_proxies_pair_recursive(
    struct engine *e, const struct space *s, const double loc_i[3],
    const double loc_j[3], const double width[3], const int depth) {

  const int zoom_depth = e->s->zoom_props->zoom_cell_depth;

  /* The local nodeID */
  const int nodeID = e->nodeID;

  /* At the zoom level we make the proxy (if necessary) */
  if (depth == zoom_depth) {
    /* Get the cells */
    struct cell *ci = &s->cells_top[cell_getid_from_pos(
        s, loc_i[0] + 0.5 * width[0], loc_i[1] + 0.5 * width[1],
        loc_i[2] + 0.5 * width[2];)];
    struct cell *cj = &s->cells_top[cell_getid_from_pos(
        s, loc_j[0] + 0.5 * width[0], loc_j[1] + 0.5 * width[1],
        loc_j[2] + 0.5 * width[2];)];

    /* Make the proxies (if we need them) */
    zoom_make_proxy_for_pair(e, ci, cj);
  }

  /* Recurse (geometrically) to the next level. */
  const double sub_width[3] = {0.5 * width[0], 0.5 * width[1], 0.5 * width[2]};
  for (int i = 0; i < 8; i++) {
    double sub_loc_i[3] = {loc_i[0] + ((i & 4) ? sub_width[0] : 0.0),
                           loc_i[1] + ((i & 2) ? sub_width[1] : 0.0),
                           loc_i[2] + ((i & 1) ? sub_width[2] : 0.0)};
    for (int j = 0; j < 8; j++) {
      double sub_loc_j[3] = {loc_j[0] + ((j & 4) ? sub_width[0] : 0.0),
                             loc_j[1] + ((j & 2) ? sub_width[1] : 0.0),
                             loc_j[2] + ((j & 1) ? sub_width[2] : 0.0)};
      zoom_make_proxies_pair_recursive(e, s, sub_loc_i, sub_loc_j, sub_width,
                                       depth + 1);
    }
  }
}

/**
 * @brief Recursively create zoom cells proxies nested in a single void cell.
 *
 * Recurse through the void hierarchy until we hit the zoom top level. Once at
 * the top level we will check what proxies we need.
 *
 * @param e The #engine.
 * @param s The #space.
 * @param loc Lower-left corner of the current geometric region.
 * @param width Width of the current geometric region.
 * @param depth Recursion depth below the top-level void.
 */
static void zoom_make_proxies_self_recursive(struct engine *e,
                                             const struct space *s,
                                             const double loc[3],
                                             const double width[3],
                                             const int depth) {

  const int zoom_depth = e->s->zoom_props->zoom_cell_depth;

  /* Once at the zoom depth there is no more proxies to make below it. */
  if (depth == zoom_depth) return;

  /* Create the geometry for all the progeny. */
  double sub_loc[8][3];
  double sub_width[3] = {0.5 * width[0], 0.5 * width[1], 0.5 * width[2]};
  for (int k = 0; k < 8; k++) {
    sub_loc[k][0] = loc[0] + ((k & 4) ? sub_width[0] : 0.0);
    sub_loc[k][1] = loc[1] + ((k & 2) ? sub_width[1] : 0.0);
    sub_loc[k][2] = loc[2] + ((k & 1) ? sub_width[2] : 0.0);
  }

  /* Keep recursing to handle zoom top level pairs inside the cell. */
  for (int k = 0; k < 8; k++) {
    zoom_make_void_self_proxies(e, s, sub_loc[k], sub_width, depth + 1);
  }

  /* Recurse into every pair of cells. */
  for (int a = 0; a < 8; a++) {
    for (int b = a + 1; b < 8; b++) {
      zoom_make_proxies_pair_recursive(e, s, sub_loc[a], sub_loc[b], sub_width,
                                       depth + 1);
    }
  }
}

#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies for a zoom simulation.
 *
 * Background <-> Background pairs behave the same as a normal uniform
 * simulation. Zoom <-> Zoom pairs are reached via recursion in void cells and
 * similarly any Background <-> Zoom pair that will exist is handled through
 * recursion in the void cell to find foreign zoom cells in the neighbouring
 * background cells.
 *
 * Note that this is more complex because of need to match the recursion
 * inherent in task creation. Since we don't have the full cell tree we have
 * to "recurse" based on geometry down to the zoom top level cells. We never
 * need to go lower than zoom top level depth since this is where the
 * partition is done.
 *
 * @param e The #engine.
 */
void zoom_engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  const ticks tic = getticks();

  /* Unpack useful information */
  const struct space *s = e->s;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;
  const int bkg_cdim[3] = {s->zoom_props->bkg_cdim[0],
                           s->zoom_props->bkg_cdim[1],
                           s->zoom_props->bkg_cdim[2]};
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
  const int periodic = s->periodic;
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_distance = e->mesh->r_cut_max;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Background cell scale. */
  const double bkg_cell_width[3] = {cells[bkg_offset].width[0],
                                    cells[bkg_offset].width[1],
                                    cells[bkg_offset].width[2]};
  const double bkg_r_diag2 = bkg_cell_width[0] * bkg_cell_width[0] +
                             bkg_cell_width[1] * bkg_cell_width[1] +
                             bkg_cell_width[2] * bkg_cell_width[2];
  const double bkg_r_diag = 0.5 * sqrt(bkg_r_diag2);
  const double bkg_r_max = 2 * bkg_r_diag;

  /* Background stencil derived from the gravity opening angle because we
   * don't have access to the MAC quantites yet. */
  int bkg_delta_cells = 1;
  if (with_gravity) {
    const double distance = 2. * bkg_r_max / theta_crit;

    /* If the mesh distance is smaller then use that instead. */
    if (periodic && distance < max_distance) {
      distance = max_distance;
    }

    /* Convert to a number of cells */
    bkg_delta_cells = (int)(distance / cells[bkg_offset].dmin) + 1;
  }

  /* Convert to upper and lower bounds */
  int bkg_delta_m = bkg_delta_cells;
  int bkg_delta_p = bkg_delta_cells;
  if (bkg_delta_cells >= bkg_cdim[0] / 2) {
    bkg_delta_m = bkg_cdim[0] / 2;
    bkg_delta_p = bkg_cdim[0] / 2;
  }

  if (e->verbose) {
    message(
        "Looking for proxies up to %d background cells away "
        "(delta_m=%d delta_p=%d)",
        bkg_delta_cells, bkg_delta_m, bkg_delta_p);
  }

  /* Loop over the background top-cell grid. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        /* Get the cell */
        const int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);
        struct cell *ci = &cells[cid];

        /* Is ci a void cell? */
        const int ci_is_void = (ci->subtype == cell_subtype_void);

        /* If the cell is a void cell, recurse inside it to handle zoom pairs
         * within it */
        if (ci_is_void) {
          zoom_make_void_self_proxies(e, s, ci->loc, ci->width, /*depth=*/0);
        }

        /* Pair walk over the stencil. */
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

              /* Convert ijk indices to neighbouring cell index. */
              const int cjd =
                  cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);

              /* Each unordered pair handled once. */
              if (cid >= cjd) continue;

              /* Get the cell */
              struct cell *cj = &cells[cjd];

              /* Is cj a void cell? */
              const int cj_is_void = (cj->subtype == cell_subtype_void);

              /* Do we have an normal pair? */
              if (!ci_is_void && !cj_is_void) {

                /* Bkg <-> Bkg pair: handle at the top level just like normal */
                zoom_make_proxy_for_pair(e, ci, cj);

              } else {
                /* Otherwise, recurse until we hit the zoom top level (either in
                 * both for a void <-> void pair or on the void side in a void
                 * <-> bkg pair */
                zoom_make_proxies_pair_recursive(e, s, cj->loc, cj->width,
                                                 /*depth=*/0, cj, ci);
              }
            }
          }
        }
      }
    }
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
