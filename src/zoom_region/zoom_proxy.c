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
 * @brief The square of the minimal distance between two equal sized regions.
 *
 * This is cell_min_dist2 but using raw geometry rather than to a pair of
 * cells. This is because we do not have the cell trees to recurse through here.
 *
 * Both regions must have the same width, exactly as cell_min_dist2 demands of
 * its cells.
 *
 * @param loc_i Lower corner of the first region.
 * @param loc_j Lower corner of the second region.
 * @param width Width of both regions.
 * @param periodic Are we using periodic BCs?
 * @param dim The dimensions of the simulation volume.
 *
 * @return The square of the minimal distance between the two regions.
 */
static double zoom_region_min_dist2(const double loc_i[3],
                                    const double loc_j[3],
                                    const double width[3], const int periodic,
                                    const double dim[3]) {

  /* Prepare the distance. */
  double dist2 = 0.0;

  /* Loop over the three dimensions. */
  for (int k = 0; k < 3; k++) {

    const double i_min = loc_i[k];
    const double i_max = loc_i[k] + width[k];
    const double j_min = loc_j[k];
    const double j_max = loc_j[k] + width[k];

    /* Get the minimal distance in this dimension, taking periodicity into
     * account (the minimum distance between any two corners of the two
     * "cells"). */
    double dx;
    if (periodic) {
      dx = min4(fabs(nearest(i_min - j_min, dim[k])),
                fabs(nearest(i_min - j_max, dim[k])),
                fabs(nearest(i_max - j_min, dim[k])),
                fabs(nearest(i_max - j_max, dim[k])));
    } else {
      dx = min4(fabs(i_min - j_min), fabs(i_min - j_max), fabs(i_max - j_min),
                fabs(i_max - j_max));
    }

    dist2 += dx * dx;
  }

  return dist2;
}

/**
 * @brief Apply the proxy criteria to a pair of regions.
 *
 * Note that regions must have equal widths, i.e. recursion is done in step.
 *
 * Only zoom <-> zoom pairs can ever be hydro pairs.
 *
 * @param e The #engine.
 * @param loc_i Lower corner of the first region.
 * @param loc_j Lower corner of the second region.
 * @param width Width of both regions.
 * @param type_i Type of the cell the first region resolves to.
 * @param type_j Type of the cell the second region resolves to.
 *
 * @return The proxy type these regions require.
 */
static int zoom_get_proxy_type(const struct engine *e, const double loc_i[3],
                               const double loc_j[3], const double width[3],
                               const int type_i, const int type_j) {

  const struct space *s = e->s;
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const int with_hydro = (e->policy & engine_policy_hydro);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist2 = e->mesh->r_cut_max * e->mesh->r_cut_max;

  /* Derive the geometric criterion for a proxy. */
  const double r_diag2 =
      width[0] * width[0] + width[1] * width[1] + width[2] * width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);
  const double r_max = 2.0 * r_diag;

  /* Minimal distance between any two points in the regions. */
  const double min_dist_CoM2 =
      zoom_region_min_dist2(loc_i, loc_j, width, s->periodic, s->dim);

  int proxy_type = proxy_cell_type_none;

  /* Do we need a gravity proxy? Nothing is needed beyond the mesh cutoff,
   * where the truncated forces are zero. */
  if (with_gravity && !(s->periodic && min_dist_CoM2 >= max_mesh_dist2) &&
      !(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2)) {
    proxy_type |= (int)proxy_cell_type_gravity;
  }

  /* Do we need a hydro proxy? Only a pair of adjacent zoom cells is ever a
   * hydro pair. Regions on the same grid either touch or lie at least a full
   * width apart, so half a width is a safe discriminator against the rounding
   * in the corners we accumulate as we recurse. */
  const double tol = 0.5 * min3(width[0], width[1], width[2]);
  if (with_hydro && type_i == cell_type_zoom && type_j == cell_type_zoom &&
      min_dist_CoM2 < tol * tol) {
    proxy_type |= (int)proxy_cell_type_hydro;
  }

  return proxy_type;
}

/**
 * @brief Create a proxy for a pair of cells if its needed.
 *
 * @param e The #engine.
 * @param ci First cell.
 * @param cj Second cell.
 * @param proxy_type The verdict for the regions these cells were reached
 * through.
 */
static void zoom_make_proxy_for_pair(struct engine *e, struct cell *ci,
                                     struct cell *cj, const int proxy_type) {

  /* No proxy needed? */
  if (proxy_type == proxy_cell_type_none) return;

  /* What node are we on? */
  const int nodeID = e->nodeID;

  /* Skip entirely local and entirely foreign pairs. */
  if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
      (ci->nodeID != nodeID && cj->nodeID != nodeID)) {
    return;
  }

  /* Ok, we need one. Add the proxy */
  engine_add_proxy(e, ci, cj, proxy_type);
}

/**
 * @brief Create proxies based on geometric pair recursion.
 *
 * This function will recurse through the pair of cells until we reach the zoom
 * top level depth. Once at the zoom top level depth we will check what proxies
 * we need.
 *
 * Note that this could lead to us creating proxies for cells that may never
 * exist but we don't have enough information yet to avoid this.
 *
 * We don't have to recurse the whole way down. If a pair is found not not need
 * a proxy at a given level then all levels below it will also not need a proxy.
 *
 * @param e The #engine.
 * @param s The #space.
 * @param loc_i Lower corner of first region.
 * @param loc_j Lower corner of second region.
 * @param width Width of both regions.
 * @param type_i Type of the cells the first region resolves to at the zoom
 * level, i.e. #cell_type_zoom if we are descending inside a void cell.
 * @param type_j As type_i, for the second region.
 * @param depth Current recursion depth below the void cells.
 * @param proxy_loc_i Initial geometry for side i, used for mixed pairs.
 * @param proxy_loc_j Initial geometry for side j, used for mixed pairs.
 * @param proxy_width Initial width, used for mixed pairs.
 */
static void zoom_make_proxies_pair_recursive(
    struct engine *e, const struct space *s, const double loc_i[3],
    const double loc_j[3], const double width[3], const int type_i,
    const int type_j, const int depth, const double proxy_loc_i[3],
    const double proxy_loc_j[3], const double proxy_width[3]) {

  const int zoom_depth = e->s->zoom_props->zoom_cell_depth;

  /* Out of range? Then so is everything below us. */
  const int mixed_pair =
      (type_i == cell_type_zoom) != (type_j == cell_type_zoom);
  const int proxy_type = zoom_get_proxy_type(
      e, mixed_pair ? proxy_loc_i : loc_i, mixed_pair ? proxy_loc_j : loc_j,
      mixed_pair ? proxy_width : width, type_i, type_j);
  if (proxy_type == proxy_cell_type_none) return;

  /* At the zoom level we make the proxy (if necessary) */
  if (depth == zoom_depth) {
    /* Get the cells */
    double xi = loc_i[0] + 0.5 * width[0];
    double yi = loc_i[1] + 0.5 * width[1];
    double zi = loc_i[2] + 0.5 * width[2];
    struct cell *ci = &s->cells_top[cell_getid_from_pos(s, xi, yi, zi)];
    double xj = loc_j[0] + 0.5 * width[0];
    double yj = loc_j[1] + 0.5 * width[1];
    double zj = loc_j[2] + 0.5 * width[2];
    struct cell *cj = &s->cells_top[cell_getid_from_pos(s, xj, yj, zj)];

    /* Make the proxies (if we need them), passing the proxy type we found
     * above to avoid the duplicate proxy type calculation. */
    zoom_make_proxy_for_pair(e, ci, cj, proxy_type);
    return;
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
                                       type_i, type_j, depth + 1, proxy_loc_i,
                                       proxy_loc_j, proxy_width);
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
    zoom_make_proxies_self_recursive(e, s, sub_loc[k], sub_width, depth + 1);
  }

  /* Recurse into every pair of cells. */
  for (int a = 0; a < 8; a++) {
    for (int b = a + 1; b < 8; b++) {
      zoom_make_proxies_pair_recursive(e, s, sub_loc[a], sub_loc[b], sub_width,
                                       cell_type_zoom, cell_type_zoom,
                                       depth + 1, sub_loc[a], sub_loc[b],
                                       sub_width);
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
    double distance = 2. * bkg_r_max / theta_crit;

    /* If the mesh distance is smaller then use that instead. */
    if (periodic && max_distance < distance) {
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
          zoom_make_proxies_self_recursive(e, s, ci->loc, ci->width,
                                           /*depth=*/0);
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
                zoom_make_proxy_for_pair(
                    e, ci, cj,
                    zoom_get_proxy_type(e, ci->loc, cj->loc, ci->width,
                                        ci->type, cj->type));

              } else {
                /* Otherwise, recurse until we hit the zoom top level (either in
                 * both for a void <-> void pair or on the void side in a void
                 * <-> bkg pair */
                /* A void side resolves to zoom cells at the zoom level, a
                 * background side stays the cell it already is. */
                zoom_make_proxies_pair_recursive(
                    e, s, ci->loc, cj->loc, ci->width,
                    ci_is_void ? cell_type_zoom : ci->type,
                    cj_is_void ? cell_type_zoom : cj->type, /*depth=*/0,
                    ci->loc, cj->loc, ci->width);
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
