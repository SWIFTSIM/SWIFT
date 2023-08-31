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
 * @param i integer coordinate of ci.
 * @param j integer coordinate of ci.
 * @param k integer coordinate of ci.
 * @param ii integer coordinate of cj.
 * @param jj integer coordinate of cj.
 * @param kk integer coordinate of cj.
 * @param r_max The minimum separation of ci and cj.
 */
int find_proxy_type(struct cell *ci, struct cell *cj, struct engine *e,
                    int i, int j, int k, int ii, int jj, int kk,
                    double r_max, const double *dim, const int periodic) {

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
    if (abs(i - ii) <= 1 && abs(j - jj) <= 1 && abs(k - kk) <= 1)
      proxy_type |= (int)proxy_cell_type_hydro;
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
      proxy_init(&proxies[e->nr_proxies], nodeID,
                 cj->nodeID);
          
      /* Store the information */
      e->proxy_ind[cj->nodeID] = e->nr_proxies;
      proxy_id = e->nr_proxies;
      e->nr_proxies += 1;
          
      /* Check the maximal proxy limit */
      if ((size_t)proxy_id > 8 * sizeof(long long))
        error("Created more than %zd proxies. cell.mpi.sendto will "
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
      proxy_init(&proxies[e->nr_proxies], e->nodeID,
                 ci->nodeID);
      
      /* Store the information */
      e->proxy_ind[ci->nodeID] = e->nr_proxies;
      proxy_id = e->nr_proxies;
      e->nr_proxies += 1;
      
      /* Check the maximal proxy limit */
      if ((size_t)proxy_id > 8 * sizeof(long long))
        error("Created more than %zd proxies. cell.mpi.sendto will "
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
 * @brief What type of proxy do we need between these cells?
 *
 * @param ci First #cell.
 * @param cj Neighbour #cell.
 * @param e The #engine.
 * @param proxies The proxies themselves.
 * @param nodeID What rank is this?
 * @param proxy_type What sort of proxy is this?
 */
void get_void_proxy(struct cell *ci, struct cell *cj, struct engine *e,
                    struct proxy *proxies, const int nodeID,
                    double rmax_i, double rmax_j) {

  /* Avoid completely local and foreign pairs */
  if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
      (ci->nodeID != nodeID && (cj->nodeID != nodeID && cj->nodeID >= 0)))
    return;

  /* If all leaves of a void cell are on the same node make a proxy for the
   * whole void cell. */
  if (cj->nodeID >= 0 && cj->nodeID < 0) {

      /* What type of proxy do we need?
       * (proxy_cell_type_none if no proxy needed). */
      int proxy_type =
        find_proxy_type(ci, cj, e, 0, 0, 0, 10, 10, 10,
                        2 * rmax_i, e->s->dim, e->s->periodic);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) return;

      /* Make the proxies. */
      add_proxy(ci, cj, e, proxies, nodeID, proxy_type);

  } else {

    /* Ok, we have to check the zoom cells inside this void cell but don't have
     * the full tree yet. */

    /* How many zoom cells? */
    int nr_zoom_cells = e->s->zoom_props->nr_zoom_cells;

    /* Loop over zoom cells. */
    for (int zoom_cjd = 0; zoom_cjd < nr_zoom_cells; zoom_cjd++) {

      /* Get the cell. */
      struct cell *zoom_cj = &e->s->cells_top[zoom_cjd];

      /* Avoid completely local and foreign pairs */
      if ((ci->nodeID == nodeID && zoom_cj->nodeID == nodeID) ||
          (ci->nodeID != nodeID && zoom_cj->nodeID != nodeID))
        continue;

      /* Is this zoom cell inside cj? */
      double zoom_loc[3] = {
      zoom_cj->loc[0] + (zoom_cj->width[0] / 2),
      zoom_cj->loc[1] + (zoom_cj->width[1] / 2),
      zoom_cj->loc[2] + (zoom_cj->width[2] / 2)
    };
      if (!((zoom_loc[0] >= cj->loc[0]
             && zoom_loc[0] < (cj->loc[0] + cj->width[0])) &&
            (zoom_loc[1] >= cj->loc[1]
             && zoom_loc[1] < (cj->loc[1] + cj->width[1])) &&
            (zoom_loc[2] >= cj->loc[2]
             && zoom_loc[2] < (cj->loc[2] + cj->width[2]))))
        continue;

      /* What type of proxy do we need?
       * (proxy_cell_type_none if no proxy needed). */
      int proxy_type =
        find_proxy_type(ci, zoom_cj, e, 0, 0, 0, 10, 10, 10,
                        rmax_i + rmax_j, e->s->dim, e->s->periodic);

      /* Abort if not in range at all */
      if (proxy_type == proxy_cell_type_none) return;

      /* Make the proxies. */
      add_proxy(ci, zoom_cj, e, proxies, nodeID, proxy_type);

    }
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

  /* Set up cell offsets. */
  const int bkg_offset = s->zoom_props->bkg_cell_offset ;
  const int buff_offset = s->zoom_props->buffer_cell_offset;

  /* Set up some width and distance variables. */
  double r_diag2, r_diag, r_max_zoom, r_max_buff, r_max_bkg;

  /* Declare looping variables. */
  int delta_cells, delta_m, delta_p, *cdim;

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


  /* ======================= Star with buffer cells ======================= */

  if (s->zoom_props->with_buffer_cells) {

    /* Get cdim. */
    cdim = s->zoom_props->buffer_cdim;

    /* Compute how many cells away we need to walk */
    delta_cells = 1; /*hydro case */

    /* Gravity needs to take the opening angle into account */
    if (with_gravity) {
      const double distance = 2. * r_max_buff / theta_crit;
      delta_cells = (int)(distance / cells[buff_offset].dmin) + 1;
    }

    /* Turn this into upper and lower bounds for loops */
    delta_m = delta_cells;
    delta_p = delta_cells;

    /* Special case where every cell is in range of every other one */
    if (delta_cells > cdim[0]) {
      delta_m = cdim[0];
      delta_p = cdim[0];
    }

    /* Let's be verbose about this choice */
    if (e->verbose)
      message(
          "Looking for proxies up to %d top-level buffer cells away "
          "(delta_m=%d delta_p=%d)", delta_cells, delta_m, delta_p);

    /* Loop over each cell in the space. */
    for (int i = 0; i < cdim[0]; i++) {
      for (int j = 0; j < cdim[1]; j++) {
        for (int k = 0; k < cdim[2]; k++) {

          /* Get the cell ID. */
          const int cid = cell_getid(cdim, i, j, k) + buff_offset;

          /* Get the cell. */
          struct cell *ci = &cells[cid];

          /* Skip void cells. */
          if (ci->subtype == void_cell) continue;

          /* Loop over all its neighbours in range. */
          for (int ii = -delta_m; ii <= delta_p; ii++) {
            int iii = i + ii;
            if (iii < 0 || iii >= cdim[0]) continue;
            for (int jj = -delta_m; jj <= delta_p; jj++) {
              int jjj = j + jj;
              if (jjj < 0 || jjj >= cdim[1]) continue;
              for (int kk = -delta_m; kk <= delta_p; kk++) {
                int kkk = k + kk;
                if (kkk < 0 || kkk >= cdim[2]) continue;

                /* Get the cell ID. */
                const int cjd = cell_getid(cdim, iii, jjj, kkk) + buff_offset;

                /* Get the cell. */
                struct cell *cj = &cells[cjd];

                /* Handle void cells. */
                if (ci->subtype == void_cell || cj->subtype == void_cell) {
                  
                  get_void_proxy(ci, cj, e, proxies, nodeID,
                                 r_max_buff, r_max_zoom);
                  continue;
                }

                /* Early abort  */
                if (cid >= cjd) continue;

                /* Avoid completely local and foreign pairs */
                if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                    (ci->nodeID != nodeID && cj->nodeID != nodeID))
                  continue;

                /* What type of proxy do we need?
                 * (proxy_cell_type_none if no proxy needed). */
                int proxy_type  = find_proxy_type(ci, cj, e, i, j, k, iii, jjj,
                                                  kkk, 2 * r_max_buff, dim,
                                                  periodic);

                /* Abort if not in range at all */
                if (proxy_type == proxy_cell_type_none) continue;

                /* Make the proxies. */
                add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
              }
            }
          } 
        }
      }
    }
  }

  /* ======================= Now background cells ======================= */

  /* Get cdim. */
  cdim = s->cdim;

  /* Compute how many cells away we need to walk */
  delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max_bkg / theta_crit;
    delta_cells = (int)(distance / cells[bkg_offset].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  delta_m = delta_cells;
  delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
    }
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d top-level background cells away "
        "(delta_m=%d delta_p=%d)", delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k) + bkg_offset;

        /* Get the cell. */
        struct cell *ci = &cells[cid];

        /* Skip void and empty cells. */
        if (ci->subtype == void_cell || ci->subtype == empty) continue;
        
        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (!periodic && (iii < 0 || iii >= cdim[0])) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (!periodic && (jjj < 0 || jjj >= cdim[1])) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (!periodic && (kkk < 0 || kkk >= cdim[2])) continue;
              kkk = (kkk + cdim[2]) % cdim[2];

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk) + bkg_offset;
              
              /* Get the cell. */
              struct cell *cj = &cells[cjd];

              /* Handle void cells. */
              if (cj->subtype == void_cell) {

                get_void_proxy(ci, cj, e, proxies, nodeID,
                               r_max_bkg, r_max_zoom);
                continue;
              }

              /* Handle empty cells. */
              else if (cj->subtype == empty) {

                /* Loop over zoom cells. */
                for (int buff_cjd = buff_offset; buff_cjd < s->nr_cells;
                     buff_cjd++) {

                  /* Get the cell. */
                  struct cell *buff_cj = &cells[buff_cjd];

                  /* Skip void cells. */
                  if (buff_cj->subtype == void_cell)
                    continue;

                  /* Avoid completely local and foreign pairs */
                  if ((ci->nodeID == nodeID && buff_cj->nodeID == nodeID) ||
                      (ci->nodeID != nodeID && buff_cj->nodeID != nodeID))
                    continue;

                  /* What type of proxy do we need?
                   * (proxy_cell_type_none if no proxy needed). */
                  int proxy_type  = find_proxy_type(ci, buff_cj, e, i, j, k,
                                                    iii, jjj, kkk,
                                                    r_max_buff + r_max_bkg,
                                                    dim, periodic);

                  /* Abort if not in range at all */
                  if (proxy_type == proxy_cell_type_none) continue;

                  /* Make the proxies. */
                  add_proxy(ci, buff_cj, e, proxies, nodeID, proxy_type);
                }

                continue;
              }

              /* Early abort  */
              if (cid >= cjd) continue;

              /* Avoid completely local and foreign pairs */
              if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              /* What type of proxy do we need?
               * (proxy_cell_type_none if no proxy needed). */
              int proxy_type  = find_proxy_type(ci, cj, e, i, j, k,
                                                iii, jjj, kkk,
                                                2 * r_max_bkg, dim, periodic);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;
              
              /* Make the proxies. */
              add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
            }
          }
        }
            }
    }
  }

  /* ======================= Finally, zoom cells ======================= */

  /* Get cdim. */
  cdim = s->zoom_props->cdim;

  /* Compute how many cells away we need to walk */
  delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max_zoom / theta_crit;
    delta_cells = (int)(distance / cells[0].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  delta_m = delta_cells;
  delta_p = delta_cells;

  /* Special case where every cell is in range of every other one */
  if (delta_cells > cdim[0]) {
    delta_m = cdim[0];
    delta_p = cdim[0];
  }

  /* Let's be verbose about this choice */
  if (e->verbose)
    message(
        "Looking for proxies up to %d top-level zoom cells away (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Get the cell. */
        struct cell *ci = &cells[cid];

        /* Loop over all its neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;

              /* Get the cell ID. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Early abort  */
              if (cid >= cjd) continue;

              /* Get the cell. */
              struct cell *cj = &cells[cjd];

              /* Avoid completely local and foreign pairs */
              if ((ci->nodeID == nodeID && cj->nodeID == nodeID) ||
                  (ci->nodeID != nodeID && cj->nodeID != nodeID))
                continue;

              /* What type of proxy do we need?
               * (proxy_cell_type_none if no proxy needed). */
              int proxy_type  = find_proxy_type(ci, cj, e, i, j, k,
                                                iii, jjj, kkk,
                                                2 * r_max_zoom, dim, periodic);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Make the proxies. */
              add_proxy(ci, cj, e, proxies, nodeID, proxy_type);
            }
          }
        }
      }
    }
  }

  /* Be clear about the time */
  if (e->verbose) {
    message("Number of proxies: %d", e->nr_proxies);
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

