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
 * @brief Apply proxy criteria to a pair of equal-width cells.
 *
 * @param e The #engine.
 * @param ci First cell whose proxy is being constructed.
 * @param cj Second cell whose proxy is being constructed.
 * @param proxy_ci First cell used for geometric proxy criteria.
 * @param proxy_cj Second cell used for geometric proxy criteria.
 * @param is_direct_neighbour Whether the cells are direct zoom neighbours.
 * @return Bit mask describing required proxy data.
 */
static int zoom_get_proxy_type(const struct engine *e, const struct cell *ci,
                               const struct cell *cj,
                               const struct cell *proxy_ci,
                               const struct cell *proxy_cj,
                               const int is_direct_neighbour) {

  const struct space *s = e->s;
  const int with_hydro_policy = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist2 = e->mesh->r_cut_max * e->mesh->r_cut_max;
  int proxy_type = 0;

  /* Hydro is only possible between direct zoom neighbours. */
  const int both_zoom =
      (ci->type == cell_type_zoom && cj->type == cell_type_zoom);
  if (with_hydro_policy && both_zoom && is_direct_neighbour) {
    proxy_type |= (int)proxy_cell_type_hydro;
  }

  if (!with_gravity) return proxy_type;

  if (is_direct_neighbour) {
    proxy_type |= (int)proxy_cell_type_gravity;
    return proxy_type;
  }

  const double r_diag2 = proxy_ci->width[0] * proxy_ci->width[0] +
                         proxy_ci->width[1] * proxy_ci->width[1] +
                         proxy_ci->width[2] * proxy_ci->width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);
  const double r_max = 2.0 * r_diag;
  const double min_dist_CoM2 =
      cell_min_dist2(proxy_ci, proxy_cj, s->periodic, s->dim);

  if (s->periodic) {
    if ((min_dist_CoM2 < max_mesh_dist2) &&
        !(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2))
      proxy_type |= (int)proxy_cell_type_gravity;
  } else {
    if (!(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2))
      proxy_type |= (int)proxy_cell_type_gravity;
  }
  return proxy_type;
}

/**
 * @brief Get zoom cell corresponding to a geometric region centre.
 *
 * @param s The #space.
 * @param loc Lower corner of the region.
 * @param width Width of the region.
 * @return The zoom top-level cell in the region.
 */
static struct cell *zoom_cell_from_region(const struct space *s,
                                          const double loc[3],
                                          const double width[3]) {
  const double x = loc[0] + 0.5 * width[0];
  const double y = loc[1] + 0.5 * width[1];
  const double z = loc[2] + 0.5 * width[2];
  return &s->cells_top[cell_getid_from_pos(s, x, y, z)];
}

/**
 * @brief Recursively create proxies for zoom-zoom pairs inside two void cells.
 *
 * Both regions are subdivided in lockstep until zoom top-level cells are
 * reached. No particle or cell-tree information is available here, so every
 * geometric pair is considered.
 *
 * @param e The #engine.
 * @param s The #space.
 * @param loc_i Lower corner of first region.
 * @param loc_j Lower corner of second region.
 * @param width Width of both regions.
 * @param depth Current recursion depth below the void cells.
 * @param zoom_depth Depth at which zoom top-level cells are reached.
 */
static void zoom_make_zoom_pair_proxies_recursive(
    struct engine *e, const struct space *s, const double loc_i[3],
    const double loc_j[3], const double width[3], const int depth,
    const int zoom_depth) {

  const int nodeID = e->nodeID;

  if (depth == zoom_depth) {
    struct cell *ci = zoom_cell_from_region(s, loc_i, width);
    struct cell *cj = zoom_cell_from_region(s, loc_j, width);

#ifdef SWIFT_DEBUG_CHECKS
    if (ci->subtype == cell_subtype_void || cj->subtype == cell_subtype_void)
      error("Zoom recursion reached a void cell (ci=%s, cj=%s)",
            subcellID_names[ci->subtype], subcellID_names[cj->subtype]);
#endif

    if (ci == cj) return;
    if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
        (ci->nodeID != nodeID && cj->nodeID != nodeID))
      return;

    int is_direct_neighbour = 0;
    const struct zoom_region_properties *zp = s->zoom_props;
    const int ai =
        (int)((ci->loc[0] - zp->region_lower_bounds[0]) * zp->iwidth[0]);
    const int aj =
        (int)((ci->loc[1] - zp->region_lower_bounds[1]) * zp->iwidth[1]);
    const int ak =
        (int)((ci->loc[2] - zp->region_lower_bounds[2]) * zp->iwidth[2]);
    const int bi =
        (int)((cj->loc[0] - zp->region_lower_bounds[0]) * zp->iwidth[0]);
    const int bj =
        (int)((cj->loc[1] - zp->region_lower_bounds[1]) * zp->iwidth[1]);
    const int bk =
        (int)((cj->loc[2] - zp->region_lower_bounds[2]) * zp->iwidth[2]);
    if (abs(ai - bi) <= 1 && abs(aj - bj) <= 1 && abs(ak - bk) <= 1)
      is_direct_neighbour = 1;

    const int proxy_type =
        zoom_get_proxy_type(e, ci, cj, ci, cj, is_direct_neighbour);
    if (proxy_type == proxy_cell_type_none) return;

    engine_add_proxy(e, ci, cj, proxy_type);
    return;
  }

  const double sub_width[3] = {0.5 * width[0], 0.5 * width[1],
                               0.5 * width[2]};
  for (int i = 0; i < 8; i++) {
    double sub_loc_i[3] = {loc_i[0] + ((i & 4) ? sub_width[0] : 0.0),
                           loc_i[1] + ((i & 2) ? sub_width[1] : 0.0),
                           loc_i[2] + ((i & 1) ? sub_width[2] : 0.0)};
    for (int j = 0; j < 8; j++) {
      double sub_loc_j[3] = {loc_j[0] + ((j & 4) ? sub_width[0] : 0.0),
                             loc_j[1] + ((j & 2) ? sub_width[1] : 0.0),
                             loc_j[2] + ((j & 1) ? sub_width[2] : 0.0)};
      zoom_make_zoom_pair_proxies_recursive(e, s, sub_loc_i, sub_loc_j,
                                            sub_width, depth + 1, zoom_depth);
    }
  }
}

/**
 * @brief Recursively create proxies for zoom cells inside one void cell.
 *
 * Recurse through the 8 octants of a top-level void cell, generating pair
 * recursions between every pair of octants (and self recursions on each
 * octant). At the leaves this produces zoom<->zoom proxy decisions between
 * zoom top cells contained in the same top-level void.
 *
 * @param e         The #engine.
 * @param s         The #space.
 * @param loc       Lower-left corner of the current geometric region.
 * @param width     Width of the current geometric region.
 * @param depth     Recursion depth below the top-level void.
 * @param zoom_depth Depth at which the leaves are zoom top cells.
 */
static void zoom_make_void_self_proxies(struct engine *e, const struct space *s,
                                        const double loc[3],
                                        const double width[3], const int depth,
                                        const int zoom_depth) {

  /* At zoom depth the "self" of a single zoom cell carries no proxy work. */
  if (depth == zoom_depth) return;

  /* Build the 8 octants. */
  double sub_loc[8][3];
  double sub_width[3] = {0.5 * width[0], 0.5 * width[1], 0.5 * width[2]};
  for (int k = 0; k < 8; k++) {
    sub_loc[k][0] = loc[0] + ((k & 4) ? sub_width[0] : 0.0);
    sub_loc[k][1] = loc[1] + ((k & 2) ? sub_width[1] : 0.0);
    sub_loc[k][2] = loc[2] + ((k & 1) ? sub_width[2] : 0.0);
  }

  /* Self recursion on each octant. */
  for (int k = 0; k < 8; k++) {
    zoom_make_void_self_proxies(e, s, sub_loc[k], sub_width, depth + 1,
                                zoom_depth);
  }

  /* Pair recursion on every pair of octants. */
  for (int a = 0; a < 8; a++) {
    for (int b = a + 1; b < 8; b++) {
      zoom_make_zoom_pair_proxies_recursive(
          e, s, sub_loc[a], sub_loc[b], sub_width, depth + 1, zoom_depth);
    }
  }
}

/**
 * @brief Recursively create proxies between zoom cells and one background cell.
 *
 * @param e The #engine.
 * @param s The #space.
 * @param loc Lower corner of the enclosing void region.
 * @param width Width of the enclosing void region.
 * @param depth Current recursion depth below the void cell.
 * @param void_cell The top-level void cell used for geometric criteria.
 * @param bkg_cell The top-level background cell receiving the proxy.
 * @param zoom_depth Depth at which the zoom top-level cells are reached.
 */
static void zoom_make_zoom_bkg_proxies_recursive(
    struct engine *e, const struct space *s, const double loc[3],
    const double width[3], const int depth, const struct cell *void_cell,
    struct cell *bkg_cell, const int zoom_depth) {

  if (depth == zoom_depth) {
    struct cell *zoom_cell = zoom_cell_from_region(s, loc, width);
    if ((zoom_cell->nodeID == e->nodeID && bkg_cell->nodeID == e->nodeID) ||
        (zoom_cell->nodeID != e->nodeID && bkg_cell->nodeID != e->nodeID))
      return;

    const int proxy_type = zoom_get_proxy_type(
        e, zoom_cell, bkg_cell, void_cell, bkg_cell, 0);
    if (proxy_type != proxy_cell_type_none)
      engine_add_proxy(e, zoom_cell, bkg_cell, proxy_type);
    return;
  }

  const double sub_width[3] = {0.5 * width[0], 0.5 * width[1],
                               0.5 * width[2]};
  for (int i = 0; i < 8; i++) {
    const double sub_loc[3] = {
        loc[0] + ((i & 4) ? sub_width[0] : 0.0),
        loc[1] + ((i & 2) ? sub_width[1] : 0.0),
        loc[2] + ((i & 1) ? sub_width[2] : 0.0)};
    zoom_make_zoom_bkg_proxies_recursive(e, s, sub_loc, sub_width, depth + 1,
                                         void_cell, bkg_cell, zoom_depth);
  }
}

#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies (zoom variant).
 *
 * Top-level loop iterates background cells with the gravity-opening-angle
 * stencil and handles background-background pairs directly. Void cells are
 * expanded geometrically only as far as zoom top-level cells to construct
 * zoom-zoom and zoom-background proxies.
 *
 * @param e The #engine.
 */
void zoom_engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  const ticks tic = getticks();

  const struct space *s = e->s;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;
  const int bkg_cdim[3] = {s->zoom_props->bkg_cdim[0],
                           s->zoom_props->bkg_cdim[1],
                           s->zoom_props->bkg_cdim[2]};
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
  const int periodic = s->periodic;
  const int zoom_depth = s->zoom_props->zoom_cell_depth;

  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;

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

  /* Background stencil derived from the gravity opening angle. */
  int bkg_delta_cells = 1;
  if (with_gravity) {
    const double distance = 2. * bkg_r_max / theta_crit;
    bkg_delta_cells = (int)(distance / cells[bkg_offset].dmin) + 1;
  }
  int bkg_delta_m = bkg_delta_cells;
  int bkg_delta_p = bkg_delta_cells;
  if (bkg_delta_cells >= bkg_cdim[0] / 2) {
    bkg_delta_m = bkg_cdim[0] / 2;
    bkg_delta_p = bkg_cdim[0] / 2;
  }

  if (e->verbose)
    message(
        "Looking for proxies up to %d background cells away "
        "(delta_m=%d delta_p=%d)",
        bkg_delta_cells, bkg_delta_m, bkg_delta_p);

  /* Loop over the background top-cell grid. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        const int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);
        struct cell *ci = &cells[cid];

        /* Construct zoom-zoom proxies within this void cell. */
        if (ci->subtype == cell_subtype_void) {
          zoom_make_void_self_proxies(e, s, ci->loc, ci->width, /*depth=*/0,
                                      zoom_depth);
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

              const int cjd =
                  cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);

              /* Each unordered pair handled once. */
              if (cid >= cjd) continue;

              struct cell *cj = &cells[cjd];

              const int ci_void = (ci->subtype == cell_subtype_void);
              const int cj_void = (cj->subtype == cell_subtype_void);

              if (!ci_void && !cj_void) {

                /* bkg <-> bkg pair: handle at the top level. */

                /* Skip entirely local and entirely foreign pairs. */
                if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                    (ci->nodeID != nodeID && cj->nodeID != nodeID))
                  continue;

                const int proxy_type = engine_get_proxy_type(
                    e, ci, i, j, k, cj, iii, jjj, kkk, bkg_r_max);
                if (proxy_type == proxy_cell_type_none) continue;
                engine_add_proxy(e, ci, cj, proxy_type);

              } else if (ci_void && cj_void) {

                /* Construct zoom-zoom proxies between two void cells. */
                zoom_make_zoom_pair_proxies_recursive(
                    e, s, ci->loc, cj->loc, ci->width, /*depth=*/0,
                    zoom_depth);

              } else if (ci_void) {

                /* Construct zoom-background proxies from the first void. */
                zoom_make_zoom_bkg_proxies_recursive(
                    e, s, ci->loc, ci->width, /*depth=*/0, ci, cj,
                    zoom_depth);

              } else {

                /* Construct zoom-background proxies from the second void. */
                zoom_make_zoom_bkg_proxies_recursive(
                    e, s, cj->loc, cj->width, /*depth=*/0, cj, ci,
                    zoom_depth);
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
