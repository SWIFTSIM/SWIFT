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
#include "error.h"
#include "proxy.h"
#include "zoom.h"

#ifdef WITH_MPI

/**
 * @brief Compute proxy type from a pair of background cells.
 *
 * Proxy type evaluation uses geometry helpers that assume same-size cells, so
 * we always evaluate this on background cells and then propagate the result to
 * zoom/background proxy pairs.
 */
static INLINE int zoom_get_proxy_type_from_bkg_pair(struct engine *e,
                                                    const struct cell *ci,
                                                    const struct cell *cj) {

  const double ir_diag2 = ci->width[0] * ci->width[0] +
                          ci->width[1] * ci->width[1] +
                          ci->width[2] * ci->width[2];
  const double ir_diag = 0.5 * sqrt(ir_diag2);
  const double jr_diag2 = cj->width[0] * cj->width[0] +
                          cj->width[1] * cj->width[1] +
                          cj->width[2] * cj->width[2];
  const double jr_diag = 0.5 * sqrt(jr_diag2);
  const double pair_r_max = ir_diag + jr_diag;

  return engine_get_proxy_type(e, ci, cj, pair_r_max);
}

/**
 * @brief Compute proxy type from a pair of zoom cells.
 */
static INLINE int zoom_get_proxy_type_from_zoom_pair(struct engine *e,
                                                     const struct cell *ci,
                                                     const struct cell *cj) {

  const double ir_diag2 = ci->width[0] * ci->width[0] +
                          ci->width[1] * ci->width[1] +
                          ci->width[2] * ci->width[2];
  const double ir_diag = 0.5 * sqrt(ir_diag2);
  const double jr_diag2 = cj->width[0] * cj->width[0] +
                          cj->width[1] * cj->width[1] +
                          cj->width[2] * cj->width[2];
  const double jr_diag = 0.5 * sqrt(jr_diag2);
  const double pair_r_max = ir_diag + jr_diag;

  return engine_get_proxy_type(e, ci, cj, pair_r_max);
}
#endif /* WITH_MPI */

/**
 * @brief Create and fill the proxies.
 *
 * This is the zoom specific form of engine_makeproxies.
 *
 * The loops follow the same structure as zoom_partition_graph_init:
 * - Walk over zoom cells, building zoom-zoom and zoom-background proxies.
 * - Then walk over background-background neighbours.
 *
 * The only difference with the partition graph build is that we walk over a
 * wider stencil (delta_m/delta_p) set by the gravity opening angle.
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
  struct cell *zoom_cells = s->zoom_props->zoom_cells_top;

  const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                            s->zoom_props->cdim[2]};
  const int bkg_cdim[3] = {s->zoom_props->bkg_cdim[0],
                           s->zoom_props->bkg_cdim[1],
                           s->zoom_props->bkg_cdim[2]};
  const int bkg_offset = s->zoom_props->bkg_cell_offset;
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

  /* Compute maximum distance for zoom cell proxy interactions */
  const double zoom_cell_width[3] = {
      zoom_cells[0].width[0], zoom_cells[0].width[1], zoom_cells[0].width[2]};

  const double zoom_r_diag2 = zoom_cell_width[0] * zoom_cell_width[0] +
                              zoom_cell_width[1] * zoom_cell_width[1] +
                              zoom_cell_width[2] * zoom_cell_width[2];
  const double zoom_r_diag = 0.5 * sqrt(zoom_r_diag2);
  const double zoom_r_max = 2 * zoom_r_diag;

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

  /* Compute how many zoom/background cells away we need to walk */
  int zoom_delta_cells = 1; /* hydro case */
  int bkg_delta_cells = 1;  /* hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double zoom_distance = 2. * zoom_r_max / theta_crit;
    const double distance = 2. * bkg_r_max / theta_crit;
    zoom_delta_cells = (int)(zoom_distance / zoom_cells[0].dmin) + 1;
    bkg_delta_cells =
        (int)(distance / cells[s->zoom_props->bkg_cell_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int zoom_delta_m = zoom_delta_cells;
  int zoom_delta_p = zoom_delta_cells;
  int bkg_delta_m = bkg_delta_cells;
  int bkg_delta_p = bkg_delta_cells;

  if (zoom_delta_cells > zoom_cdim[0]) {
    zoom_delta_m = zoom_cdim[0];
    zoom_delta_p = zoom_cdim[0];
  }

  /* Special case where every background cell is in range of every other one */
  if (bkg_delta_cells >= bkg_cdim[0] / 2) {
    bkg_delta_m = bkg_cdim[0] / 2;
    bkg_delta_p = bkg_cdim[0] / 2;
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d zoom cells away (delta_m=%d delta_p=%d)"
        " and %d background cells away (delta_m=%d delta_p=%d)",
        zoom_delta_cells, zoom_delta_m, zoom_delta_p, bkg_delta_cells,
        bkg_delta_m, bkg_delta_p);

  /*  Zoom -> zoom and zoom -> background neighbours */
  for (int i = 0; i < zoom_cdim[0]; i++) {
    for (int j = 0; j < zoom_cdim[1]; j++) {
      for (int k = 0; k < zoom_cdim[2]; k++) {

        const int zid = cell_getid(zoom_cdim, i, j, k);
        struct cell *zi = &zoom_cells[zid];

        /* Zoom -> zoom neighbours. Zoom cells are never periodic. */
        for (int ii = -zoom_delta_m; ii <= zoom_delta_p; ii++) {
          const int iii = i + ii;
          if (iii < 0 || iii >= zoom_cdim[0]) continue;

          for (int jj = -zoom_delta_m; jj <= zoom_delta_p; jj++) {
            const int jjj = j + jj;
            if (jjj < 0 || jjj >= zoom_cdim[1]) continue;

            for (int kk = -zoom_delta_m; kk <= zoom_delta_p; kk++) {
              const int kkk = k + kk;
              if (kkk < 0 || kkk >= zoom_cdim[2]) continue;

              const int zjd = cell_getid(zoom_cdim, iii, jjj, kkk);
              if (zid >= zjd) continue;

              struct cell *zj = &zoom_cells[zjd];

              if ((zi->nodeID == nodeID && zj->nodeID == nodeID) ||
                  (zi->nodeID != nodeID && zj->nodeID != nodeID))
                continue;

              const int proxy_type =
                  zoom_get_proxy_type_from_zoom_pair(e, zi, zj);
              if (proxy_type == proxy_cell_type_none) continue;

              engine_add_proxy(e, zi, zj, proxy_type);
            }
          }
        }

        /* Find this zoom cell's enclosing background (void) cell. */
        const int bkg_i = (zi->loc[0] + 0.5 * zi->width[0]) * s->iwidth[0];
        const int bkg_j = (zi->loc[1] + 0.5 * zi->width[1]) * s->iwidth[1];
        const int bkg_k = (zi->loc[2] + 0.5 * zi->width[2]) * s->iwidth[2];

        const int ci_void_id =
            cell_getid_offset(bkg_cdim, bkg_offset, bkg_i, bkg_j, bkg_k);
        struct cell *ci_void = &cells[ci_void_id];

        for (int ii = -bkg_delta_m; ii <= bkg_delta_p; ii++) {
          int iii = bkg_i + ii;
          if (!periodic && (iii < 0 || iii >= bkg_cdim[0])) continue;
          iii = (iii + bkg_cdim[0]) % bkg_cdim[0];

          for (int jj = -bkg_delta_m; jj <= bkg_delta_p; jj++) {
            int jjj = bkg_j + jj;
            if (!periodic && (jjj < 0 || jjj >= bkg_cdim[1])) continue;
            jjj = (jjj + bkg_cdim[1]) % bkg_cdim[1];

            for (int kk = -bkg_delta_m; kk <= bkg_delta_p; kk++) {
              int kkk = bkg_k + kk;
              if (!periodic && (kkk < 0 || kkk >= bkg_cdim[2])) continue;
              kkk = (kkk + bkg_cdim[2]) % bkg_cdim[2];

              const int cjd =
                  cell_getid_offset(bkg_cdim, bkg_offset, iii, jjj, kkk);
              struct cell *cj = &cells[cjd];

              const int proxy_type =
                  zoom_get_proxy_type_from_bkg_pair(e, ci_void, cj);
              if (proxy_type == proxy_cell_type_none) continue;

              if (cj->subtype == cell_subtype_void) continue;

              if ((zi->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (zi->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              engine_add_proxy(e, zi, cj, proxy_type);
            }
          }
        }
      }
    }
  }

  /* Background -> background neighbours */
  for (int i = 0; i < bkg_cdim[0]; i++) {
    for (int j = 0; j < bkg_cdim[1]; j++) {
      for (int k = 0; k < bkg_cdim[2]; k++) {

        const int cid = cell_getid_offset(bkg_cdim, bkg_offset, i, j, k);
        struct cell *ci = &cells[cid];

        if (ci->subtype == cell_subtype_void) continue;

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
              if (cid >= cjd) continue;

              struct cell *cj = &cells[cjd];

              if (cj->subtype == cell_subtype_void) continue;

              if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              const int proxy_type =
                  zoom_get_proxy_type_from_bkg_pair(e, ci, cj);
              if (proxy_type == proxy_cell_type_none) continue;

              engine_add_proxy(e, ci, cj, proxy_type);
            }
          }
        }
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
