/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#include "../config.h"
#include "cell.h"
#include "engine.h"
#include "gravity_properties.h"
#include "partition.h"
#include "proxy.h"
#include "space.h"
#include "zoom_region.h"

#include <float.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

#ifdef WITH_MPI
/**
 * @brief What type of proxy do we need between these cells?
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param is_adjacent Do there cells touch?
 * @param r_max The minimum separation of ci and cj.
 */
int find_proxy_type(struct cell *ci, struct cell *cj, struct engine *e,
                    int is_adjacent, double r_max, const double *dim,
                    const int periodic) {

  /* Gravity information */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  int proxy_type = 0;

  /* In the hydro case, only care about direct neighbours of the same
   * type. */
  if (with_hydro && ci->type == zoom && ci->type == cj->type) {

    /* Check for direct neighbours without periodic BC */
    if (is_adjacent) proxy_type |= (int)proxy_cell_type_hydro;
  }

  /* In the gravity case, check distances using the MAC. */
  if (with_gravity) {

    /* We don't have multipoles yet (or their CoMs) so we will
       have to cook up something based on cell locations only. We
       hence need a lower limit on the distance that the CoMs in
       those cells could have and an upper limit on the distance
       of the furthest particle in the multipole from its CoM.
       We then can decide whether we are too close for an M2L
       interaction and hence require a proxy as this pair of cells
       cannot rely on just an M2L calculation. */

    /* Minimal distance between any two points in the cells. */
    const double min_dist_CoM2 = cell_min_dist2(ci, cj, periodic, dim);

    /* Are we beyond the distance where the truncated forces are 0
     * but not too far such that M2L can be used? */
    if (periodic) {

      if ((min_dist_CoM2 < max_mesh_dist2) &&
          !(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2))
        proxy_type |= (int)proxy_cell_type_gravity;

    } else {

      if (!(4. * r_max * r_max < theta_crit * theta_crit * min_dist_CoM2)) {
        proxy_type |= (int)proxy_cell_type_gravity;
      }
    }
  }

  return proxy_type;
}

/**
 * @brief What type of proxy do we need between these cells?
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param e The #engine.
 * @param proxies The proxies themselves.
 * @param nodeID What rank is this?
 * @param proxy_type What sort of proxy is this?
 */
void add_proxy(struct cell *ci, struct cell *cj, struct engine *e,
               struct proxy *proxies, const int nodeID, const int proxy_type) {

  /* Add to proxies? */
  if (ci->nodeID == nodeID && cj->nodeID != nodeID) {

    /* Do we already have a relationship with this node? */
    int proxy_id = e->proxy_ind[cj->nodeID];
    if (proxy_id < 0) {
      if (e->nr_proxies == engine_maxproxies)
        error("Maximum number of proxies exceeded.");

      /* Ok, start a new proxy for this pair of nodes */
      proxy_init(&proxies[e->nr_proxies], nodeID, cj->nodeID);

      /* Store the information */
      e->proxy_ind[cj->nodeID] = e->nr_proxies;
      proxy_id = e->nr_proxies;
      e->nr_proxies += 1;

      /* Check the maximal proxy limit */
      if ((size_t)proxy_id > 8 * sizeof(long long))
        error(
            "Created more than %zd proxies. cell.mpi.sendto will "
            "overflow.",
            8 * sizeof(long long));
    }

    /* Add the cell to the proxy */
    proxy_addcell_in(&proxies[proxy_id], cj, proxy_type);
    proxy_addcell_out(&proxies[proxy_id], ci, proxy_type);

    /* Store info about where to send the cell */
    ci->mpi.sendto |= (1ULL << proxy_id);
  }

  /* Same for the symmetric case? */
  if (cj->nodeID == nodeID && ci->nodeID != nodeID) {

    /* Do we already have a relationship with this node? */
    int proxy_id = e->proxy_ind[ci->nodeID];
    if (proxy_id < 0) {
      if (e->nr_proxies == engine_maxproxies)
        error("Maximum number of proxies exceeded.");

      /* Ok, start a new proxy for this pair of nodes */
      proxy_init(&proxies[e->nr_proxies], e->nodeID, ci->nodeID);

      /* Store the information */
      e->proxy_ind[ci->nodeID] = e->nr_proxies;
      proxy_id = e->nr_proxies;
      e->nr_proxies += 1;

      /* Check the maximal proxy limit */
      if ((size_t)proxy_id > 8 * sizeof(long long))
        error(
            "Created more than %zd proxies. cell.mpi.sendto will "
            "overflow.",
            8 * sizeof(long long));
    }

    /* Add the cell to the proxy */
    proxy_addcell_in(&proxies[proxy_id], ci, proxy_type);
    proxy_addcell_out(&proxies[proxy_id], cj, proxy_type);

    /* Store info about where to send the cell */
    cj->mpi.sendto |= (1ULL << proxy_id);
  }
}

/**
 * @brief Get the zoom cell proxies for leaves of this void cell.
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param e The #engine.
 * @param proxies The proxies themselves.
 * @param nodeID What rank is this?
 * @param proxy_type What sort of proxy is this?
 */
void get_void_proxy(struct cell *c, struct cell *void_c, struct space *s,
                    struct proxy *proxies, const int nodeID, double void_rmax,
                    double zoom_rmax) {

  /* How many zoom cells? */
  int nr_zoom_cells = s->zoom_props->nr_zoom_cells;

  /* Loop over zoom cells. */
  for (int zoom_cjd = 0; zoom_cjd < nr_zoom_cells; zoom_cjd++) {

    /* Get the cell. */
    struct cell *zoom_cj = &s->cells_top[zoom_cjd];

    /* Avoid completely local and foreign pairs */
    if ((c->nodeID == nodeID && zoom_cj->nodeID == nodeID) ||
        (c->nodeID != nodeID && zoom_cj->nodeID != nodeID))
      continue;

    /* What type of proxy do we need?
     * (proxy_cell_type_none if no proxy needed). */
    int proxy_type =
        find_proxy_type(zoom_cj, c, s->e, 0 /*is_adjacent*/,
                        void_rmax + zoom_rmax, s->dim, s->periodic);

    /* Abort if not in range at all */
    if (proxy_type == proxy_cell_type_none) return;

    /* Make the proxies. */
    add_proxy(zoom_cj, c, s->e, proxies, nodeID, proxy_type);
  }
}
#endif

/**
 * @brief Create and fill the proxies including the zoom region.
 *
 * This replaces the function in engine_proxy when running with a zoom region.
 *
 * @param e The #engine.
 */
void engine_makeproxies_with_zoom_region(struct engine *e) {
#ifdef WITH_MPI

  /* Let's time this */
  const ticks tic = getticks();

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL)
    if ((e->proxy_ind = (int *)malloc(sizeof(int) * e->nr_nodes)) == NULL)
      error("Failed to allocate proxy index.");
  for (int k = 0; k < e->nr_nodes; k++) e->proxy_ind[k] = -1;
  e->nr_proxies = 0;

  /* Useful local information */
  struct space *s = e->s;
  const int nodeID = e->nodeID;

  /* Handle on the cells and proxies */
  struct cell *cells = s->cells_top;
  struct proxy *proxies = e->proxies;

  /* Gravity information */
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;

  /* Some info about the domain */
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;

  /* Set up some width and distance variables. */
  double r_diag2, r_diag, r_max_zoom, r_max_buff, r_max_bkg;

  /* Calculate r_max for each level. */

  /* Distance between centre of the cell and corners */
  r_diag2 = ((cells[0].width[0] * cells[0].width[0]) +
             (cells[0].width[1] * cells[0].width[1]) +
             (cells[0].width[2] * cells[0].width[2]));
  r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  r_max_zoom = r_diag;

  /* Distance between centre of the cell and corners */
  r_diag2 = ((cells[bkg_offset].width[0] * cells[bkg_offset].width[0]) +
             (cells[bkg_offset].width[1] * cells[bkg_offset].width[1]) +
             (cells[bkg_offset].width[2] * cells[bkg_offset].width[2]));
  r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  r_max_bkg = r_diag;

  /* Do we have buffer cells? */
  if (s->zoom_props->with_buffer_cells) {

    /* Distance between centre of the cell and corners */
    r_diag2 = ((cells[buff_offset].width[0] * cells[buff_offset].width[0]) +
               (cells[buff_offset].width[1] * cells[buff_offset].width[1]) +
               (cells[buff_offset].width[2] * cells[buff_offset].width[2]));
    r_diag = 0.5 * sqrt(r_diag2);

    /* Maximal distance from shifted CoM to any corner */
    r_max_buff = r_diag;
  } else {
    r_max_buff = 0;
  }

  /* Define an array of r_maxs which can be indexed by cell type. */
  const double r_maxs[3] = {r_max_bkg, r_max_zoom, r_max_buff};

  /* Define an array of cdims which can be indexed by cell type. */
  const int cdims[3][3] = {s->cdim, s->zoom_props->cdim,
                           s->zoom_props->buffer_cdim};
  /* Loop over all cells.
   * Note, this is the dumb and wasteful way to do this but simpler for
   * debugging purposes! */
  for (int cid = 0; cid < s->nr_cells; cid++) {

    /* Get the cell. */
    struct cell *ci = &cells[cid];

    /* Skip void or empty cells. */
    if (ci->subtype == void_cell || ci->subtype == empty) continue;

    /* Get the cdim for this cell type. */
    const int cdim[3] = cdims[ci->type];

    /* Integer indices of the cell in the top-level grid */
    const int i = cid / (cdim[1] * cdim[2]);
    const int j = (cid / cdim[2]) % cdim[1];
    const int k = cid % cdim[2];

    /* Loop over all other cells we haven't checked against this one yet. */
    for (int cjd = cid + 1; cjd < s->nr_cells; cjd++) {

      /* Get the other cell. */
      struct cell *cj = &cells[cjd];

      /* Skip void cells and empty cells. */
      if (cj->subtype == void_cell || cj->subtype == empty) continue;

      /* Avoid completely local and foreign pairs */
      if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
          (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      /* Get the cdim for this cell type. */
      const int cdjm[3] = cdims[cj->type];

      /* Integer indices of the cell in the top-level grid */
      const int ii = cjd / (cdjm[1] * cdjm[2]);
      const int jj = (cjd / cdjm[2]) % cdjm[1];
      const int kk = cjd % cdjm[2];

      /* Get the r_max for this combination of cells */
      const double r_max = r_maxs[ci->type] + r_maxs[cj > type];

      /* Are these cells adjacent? */
      int is_adjacent = 0;
      if (ci->type == cj->type) {
        is_adjacent = (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
      }

      /* What type of proxy do we need?
       * (proxy_cell_type_none if no proxy needed). */
      int proxy_type =
          find_proxy_type(ci, cj, e, is_adjacent, r_max, dim, periodic);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) continue;

      /* Make the proxies. */
      add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
    }
  }

  if (e->verbose) {
    message("Number of proxies: %d", e->nr_proxies);
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
