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
 * @brief Axis-aligned box-to-box minimum squared distance with periodicity.
 *
 * Both boxes are described by their lower corner and width. Periodic wrapping
 * is applied per axis when @c periodic is non-zero. This tolerates unequal
 * widths between the two boxes (unlike @c cell_min_dist2 under debug checks).
 */
static double zoom_proxy_box_min_dist2(const double loc_i[3],
                                       const double width_i[3],
                                       const double loc_j[3],
                                       const double width_j[3],
                                       const int periodic,
                                       const double dim[3]) {
  double r2 = 0.0;
  for (int d = 0; d < 3; d++) {
    const double ai_lo = loc_i[d];
    const double ai_hi = loc_i[d] + width_i[d];
    const double aj_lo = loc_j[d];
    const double aj_hi = loc_j[d] + width_j[d];
    /* Signed gap on this axis (positive = no overlap). */
    double gap;
    if (ai_hi <= aj_lo)
      gap = aj_lo - ai_hi;
    else if (aj_hi <= ai_lo)
      gap = ai_lo - aj_hi;
    else
      gap = 0.0;
    if (periodic) {
      /* Try the wrapped case: the smallest distance between any two points
       * accounting for box wrapping. */
      const double L = dim[d];
      const double wrapped_gap = L - (width_i[d] + width_j[d]) -
                                 ((ai_hi <= aj_lo) ? (aj_lo - ai_hi)
                                  : (aj_hi <= ai_lo) ? (ai_lo - aj_hi)
                                                     : 0.0);
      if (wrapped_gap > 0.0 && wrapped_gap < gap) gap = wrapped_gap;
    }
    r2 += gap * gap;
  }
  return r2;
}

/**
 * @brief Apply the gravity proxy decision at a recursion leaf with cells of
 *        possibly unequal widths.
 *
 * Mirrors @c engine_get_proxy_type's gravity branch (theta_crit MAC test and
 * periodic mesh-cutoff gate) but uses our hand-rolled box-to-box min-distance
 * so unequal widths are handled cleanly (the equal-width assertion in
 * @c cell_min_dist2 would otherwise trip in debug builds).
 *
 * Hydro: enabled only when both leaves are zoom cells AND they are direct
 * neighbours in the zoom grid (matches splittask hydro-pair eligibility and
 * @c engine_get_proxy_type's zoom-only-hydro rule).
 */
static int zoom_proxy_get_leaf_proxy_type(const struct engine *e,
                                          const struct cell *ci,
                                          const struct cell *cj,
                                          const int is_direct_neighbour) {
  const struct space *s = e->s;
  const int with_hydro_policy = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist2 = e->mesh->r_cut_max * e->mesh->r_cut_max;
  int proxy_type = 0;

  /* Hydro: only zoom<->zoom direct neighbours are eligible. */
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

  /* r_max from the actual leaf cells (max diagonal). */
  const double diag2_i = ci->width[0] * ci->width[0] +
                         ci->width[1] * ci->width[1] +
                         ci->width[2] * ci->width[2];
  const double diag2_j = cj->width[0] * cj->width[0] +
                         cj->width[1] * cj->width[1] +
                         cj->width[2] * cj->width[2];
  const double r_diag = 0.5 * sqrt(diag2_i > diag2_j ? diag2_i : diag2_j);
  const double r_max = 2.0 * r_diag;

  const double min_dist_CoM2 =
      zoom_proxy_box_min_dist2(ci->loc, ci->width, cj->loc, cj->width,
                               s->periodic, s->dim);

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
 * @brief Geometric pair recursion mirror of zoom splittask void pair.
 *
 * One side of the pair (the "void side(s)") is described purely geometrically
 * by (loc, width); we DO NOT walk real cell tree pointers because at proxy
 * creation time we have not yet exchanged neighbour cell trees. The other
 * side may be either a real top-level non-void cell (passed via @c top_other,
 * not subdivided), or another void side being recursed in lock-step.
 *
 * No proxy logic is performed in this function except at the leaf level,
 * where one or both geometric sides have reached the zoom depth. At the
 * leaf:
 *   - If both sides are void-recursing, they have both reached zoom depth
 *     and we resolve real zoom top cells from positions.
 *   - If only one side is void-recursing (the other is a fixed bkg top
 *     cell), the recursing side reaches zoom depth and we resolve its real
 *     zoom top cell.
 *
 * @param e          The #engine.
 * @param s          The #space.
 * @param loc_i      Lower-left corner of side i (geometric region).
 * @param width_i    Width of side i.
 * @param depth_i    Recursion depth on side i below the top-level void
 *                   (0 == top-level void). Only meaningful when side i is
 *                   void-recursing.
 * @param recurse_i  1 if side i is descending through a void; 0 if side i
 *                   is a fixed real top-level cell (held in @c top_i).
 * @param top_i      Top-level cell pointer for side i. Always set; for a
 *                   void-recursing side this is the top-level void cell at
 *                   the start of the recursion (used only for foreign-pair
 *                   short-circuiting at the top level).
 * @param loc_j ... top_j  Same as above for side j.
 * @param zoom_depth Depth at which a void's progeny are zoom top cells.
 */
static void zoom_proxy_void_pair_recursive(
    struct engine *e, const struct space *s, const double loc_i[3],
    const double width_i[3], const int depth_i, const int recurse_i,
    struct cell *top_i, const double loc_j[3], const double width_j[3],
    const int depth_j, const int recurse_j, struct cell *top_j,
    const int zoom_depth) {

  const int nodeID = e->nodeID;

  /* Have we reached the leaf level on every recursing side? A non-recursing
   * side is always at its leaf (it is the top-level real cell). */
  const int leaf_i = !recurse_i || depth_i == zoom_depth;
  const int leaf_j = !recurse_j || depth_j == zoom_depth;

  if (leaf_i && leaf_j) {

    /* Resolve the real top cell on each recursing side from its geometric
     * centre. The non-recursing side already has its real top cell. */
    struct cell *ci = top_i;
    struct cell *cj = top_j;

    if (recurse_i) {
      const double cx = loc_i[0] + 0.5 * width_i[0];
      const double cy = loc_i[1] + 0.5 * width_i[1];
      const double cz = loc_i[2] + 0.5 * width_i[2];
      ci = &s->cells_top[cell_getid_from_pos(s, cx, cy, cz)];
    }
    if (recurse_j) {
      const double cx = loc_j[0] + 0.5 * width_j[0];
      const double cy = loc_j[1] + 0.5 * width_j[1];
      const double cz = loc_j[2] + 0.5 * width_j[2];
      cj = &s->cells_top[cell_getid_from_pos(s, cx, cy, cz)];
    }

#ifdef SWIFT_DEBUG_CHECKS
    if (ci->subtype == cell_subtype_void || cj->subtype == cell_subtype_void)
      error("Recursion bottomed out on a void cell (ci=%s, cj=%s)",
            subcellID_names[ci->subtype], subcellID_names[cj->subtype]);
#endif

    /* Skip self pairs in the symmetric case (e.g. a void self-pair recursion
     * that happens to land on the same zoom cell from both sides). */
    if (ci == cj) return;

    /* Skip entirely local and entirely foreign pairs. */
    if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
        (ci->nodeID != nodeID && cj->nodeID != nodeID))
      return;

    /* Decide and add the proxy. r_max and the min-dist are computed from the
     * actual leaf cells inside the helper, so unequal widths (zoom<->bkg)
     * are handled correctly. is_direct_neighbour is determined here from
     * the actual leaf cells: zoom<->zoom adjacency in the zoom grid (for
     * hydro eligibility), otherwise 0 (the gravity MAC test handles it). */
    int is_direct_neighbour = 0;
    if (ci->type == cell_type_zoom && cj->type == cell_type_zoom) {
      const struct zoom_region_properties *zp = s->zoom_props;
      const int ai = (int)((ci->loc[0] - zp->region_lower_bounds[0]) *
                           zp->iwidth[0]);
      const int aj = (int)((ci->loc[1] - zp->region_lower_bounds[1]) *
                           zp->iwidth[1]);
      const int ak = (int)((ci->loc[2] - zp->region_lower_bounds[2]) *
                           zp->iwidth[2]);
      const int bi = (int)((cj->loc[0] - zp->region_lower_bounds[0]) *
                           zp->iwidth[0]);
      const int bj = (int)((cj->loc[1] - zp->region_lower_bounds[1]) *
                           zp->iwidth[1]);
      const int bk = (int)((cj->loc[2] - zp->region_lower_bounds[2]) *
                           zp->iwidth[2]);
      if (abs(ai - bi) <= 1 && abs(aj - bj) <= 1 && abs(ak - bk) <= 1)
        is_direct_neighbour = 1;
    }
    const int proxy_type =
        zoom_proxy_get_leaf_proxy_type(e, ci, cj, is_direct_neighbour);
    if (proxy_type == proxy_cell_type_none) return;

    engine_add_proxy(e, ci, cj, proxy_type);
    return;
  }

  /* Not at the leaf yet on at least one side: subdivide the recursing
   * side(s) and recurse. */

  /* Iterate progeny on side i (or treat side i as a single "child" if it is
   * not recursing). */
  const int ni = recurse_i ? 8 : 1;
  const int nj = recurse_j ? 8 : 1;

  for (int ki = 0; ki < ni; ki++) {

    double sub_loc_i[3];
    double sub_width_i[3];
    int sub_depth_i;
    if (recurse_i && !leaf_i) {
      sub_width_i[0] = 0.5 * width_i[0];
      sub_width_i[1] = 0.5 * width_i[1];
      sub_width_i[2] = 0.5 * width_i[2];
      sub_loc_i[0] = loc_i[0] + ((ki & 4) ? sub_width_i[0] : 0.0);
      sub_loc_i[1] = loc_i[1] + ((ki & 2) ? sub_width_i[1] : 0.0);
      sub_loc_i[2] = loc_i[2] + ((ki & 1) ? sub_width_i[2] : 0.0);
      sub_depth_i = depth_i + 1;
    } else {
      sub_loc_i[0] = loc_i[0];
      sub_loc_i[1] = loc_i[1];
      sub_loc_i[2] = loc_i[2];
      sub_width_i[0] = width_i[0];
      sub_width_i[1] = width_i[1];
      sub_width_i[2] = width_i[2];
      sub_depth_i = depth_i;
    }

    for (int kj = 0; kj < nj; kj++) {

      double sub_loc_j[3];
      double sub_width_j[3];
      int sub_depth_j;
      if (recurse_j && !leaf_j) {
        sub_width_j[0] = 0.5 * width_j[0];
        sub_width_j[1] = 0.5 * width_j[1];
        sub_width_j[2] = 0.5 * width_j[2];
        sub_loc_j[0] = loc_j[0] + ((kj & 4) ? sub_width_j[0] : 0.0);
        sub_loc_j[1] = loc_j[1] + ((kj & 2) ? sub_width_j[1] : 0.0);
        sub_loc_j[2] = loc_j[2] + ((kj & 1) ? sub_width_j[2] : 0.0);
        sub_depth_j = depth_j + 1;
      } else {
        sub_loc_j[0] = loc_j[0];
        sub_loc_j[1] = loc_j[1];
        sub_loc_j[2] = loc_j[2];
        sub_width_j[0] = width_j[0];
        sub_width_j[1] = width_j[1];
        sub_width_j[2] = width_j[2];
        sub_depth_j = depth_j;
      }

      zoom_proxy_void_pair_recursive(e, s, sub_loc_i, sub_width_i, sub_depth_i,
                                     recurse_i, top_i, sub_loc_j, sub_width_j,
                                     sub_depth_j, recurse_j, top_j, zoom_depth);
    }
  }
}

/**
 * @brief Geometric self recursion mirror of zoom splittask void self.
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
 * @param top       The top-level void cell (carried for context).
 * @param zoom_depth Depth at which the leaves are zoom top cells.
 */
static void zoom_proxy_void_self_recursive(struct engine *e,
                                           const struct space *s,
                                           const double loc[3],
                                           const double width[3],
                                           const int depth, struct cell *top,
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
    zoom_proxy_void_self_recursive(e, s, sub_loc[k], sub_width, depth + 1, top,
                                   zoom_depth);
  }

  /* Pair recursion on every pair of octants. */
  for (int a = 0; a < 8; a++) {
    for (int b = a + 1; b < 8; b++) {
      zoom_proxy_void_pair_recursive(
          e, s, sub_loc[a], sub_width, depth + 1, /*recurse_i=*/1, top,
          sub_loc[b], sub_width, depth + 1, /*recurse_j=*/1, top, zoom_depth);
    }
  }
}

#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies (zoom variant).
 *
 * Top-level loop iterates background top cells with the gravity-opening-angle
 * stencil and:
 *   - bkg<->bkg pairs: standard top-level proxy decision (with index-stride
 *     is_direct_neighbour test and bkg-scale r_max).
 *   - any pair involving a void cell (one or both sides): hand off to the
 *     geometric pair recursion (zoom_proxy_void_pair_recursive) which mirrors
 *     the splittask void pair recursion. Proxy decisions for these are made
 *     only at the leaves (zoom top cells).
 *   - top-level void self pairs: hand off to the geometric self recursion
 *     (zoom_proxy_void_self_recursive) which generates zoom<->zoom proxies
 *     within the same top-level void.
 *
 * @param e The #engine.
 */
void zoom_engine_makeproxies(struct engine *e) {

#ifdef WITH_MPI
  const ticks tic = getticks();

  const struct space *s = e->s;
  const int nodeID = e->nodeID;
  struct cell *cells = s->cells_top;
  struct cell *zoom_cells = s->zoom_props->zoom_cells_top;
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
    message("Looking for proxies up to %d background cells away "
            "(delta_m=%d delta_p=%d)",
            bkg_delta_cells, bkg_delta_m, bkg_delta_p);

  /* Loop over the background top-cell grid. */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        const int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);
        struct cell *ci = &cells[cid];

        /* Top-level void self: recurse to find zoom<->zoom pairs inside it. */
        if (ci->subtype == cell_subtype_void) {
          zoom_proxy_void_self_recursive(e, s, ci->loc, ci->width, /*depth=*/0,
                                         ci, zoom_depth);
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

                const int is_direct_neighbour =
                    ((abs(i - iii) <= 1 || abs(i - iii - bkg_cdim[0]) <= 1 ||
                      abs(i - iii + bkg_cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - bkg_cdim[1]) <= 1 ||
                      abs(j - jjj + bkg_cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - bkg_cdim[2]) <= 1 ||
                      abs(k - kkk + bkg_cdim[2]) <= 1));
                const int proxy_type = engine_get_proxy_type(
                    e, ci, cj, is_direct_neighbour, bkg_r_max);
                if (proxy_type == proxy_cell_type_none) continue;
                engine_add_proxy(e, ci, cj, proxy_type);

              } else {

                /* Anything involving a void cell: recurse geometrically. The
                 * recursion makes the proxy decision only at the leaves. */
                zoom_proxy_void_pair_recursive(
                    e, s, ci->loc, ci->width, /*depth_i=*/0,
                    /*recurse_i=*/ci_void, ci, cj->loc, cj->width,
                    /*depth_j=*/0, /*recurse_j=*/cj_void, cj, zoom_depth);
              }
            }
          }
        }
      }
    }
  }

  /* Silence unused warning if we never look up zoom_cells directly. */
  (void)zoom_cells;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
