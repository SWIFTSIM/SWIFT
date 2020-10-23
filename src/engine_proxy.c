/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "proxy.h"

/**
 * @brief Create and fill the proxies.
 *
 * @param e The #engine.
 */
void engine_makeproxies(struct engine *e) {

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
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = s->periodic;
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double theta_crit_inv = 1. / e->gravity_properties->theta_crit;
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

  /* Compute how many cells away we need to walk */
  int delta_cells = 1; /*hydro case */

  /* Gravity needs to take the opening angle into account */
  if (with_gravity) {
    const double distance = 2. * r_max * theta_crit_inv;
    delta_cells = (int)(distance / cells[0].dmin) + 1;
  }

  /* Turn this into upper and lower bounds for loops */
  int delta_m = delta_cells;
  int delta_p = delta_cells;

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
        "Looking for proxies up to %d top-level cells away (delta_m=%d "
        "delta_p=%d)",
        delta_cells, delta_m, delta_p);

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Loop over all its neighbours neighbours in range. */
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
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Early abort  */
              if (cid >= cjd) continue;

              /* Early abort (both same node) */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID == nodeID)
                continue;

              /* Early abort (both foreign node) */
              if (cells[cid].nodeID != nodeID && cells[cjd].nodeID != nodeID)
                continue;

              int proxy_type = 0;

              /* In the hydro case, only care about direct neighbours */
              if (with_hydro) {

                // MATTHIEU: to do: Write a better expression for the
                // non-periodic case.

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1)))
                  proxy_type |= (int)proxy_cell_type_hydro;
              }

              /* In the gravity case, check distances using the MAC. */
              if (with_gravity) {

                /* First just add the direct neighbours. Then look for
                   some further out if the opening angle demands it */

                /* This is super-ugly but checks for direct neighbours */
                /* with periodic BC */
                if (((abs(i - iii) <= 1 || abs(i - iii - cdim[0]) <= 1 ||
                      abs(i - iii + cdim[0]) <= 1) &&
                     (abs(j - jjj) <= 1 || abs(j - jjj - cdim[1]) <= 1 ||
                      abs(j - jjj + cdim[1]) <= 1) &&
                     (abs(k - kkk) <= 1 || abs(k - kkk - cdim[2]) <= 1 ||
                      abs(k - kkk + cdim[2]) <= 1))) {

                  proxy_type |= (int)proxy_cell_type_gravity;
                } else {

                  /* We don't have multipoles yet (or their CoMs) so we will
                     have to cook up something based on cell locations only. We
                     hence need a lower limit on the distance that the CoMs in
                     those cells could have and an upper limit on the distance
                     of the furthest particle in the multipole from its CoM.
                     We then can decide whether we are too close for an M2L
                     interaction and hence require a proxy as this pair of cells
                     cannot rely on just an M2L calculation. */

                  /* Minimal distance between any two points in the cells */
                  const double min_dist_CoM2 = cell_min_dist2_same_size(
                      &cells[cid], &cells[cjd], periodic, dim);

                  /* Are we beyond the distance where the truncated forces are 0
                   * but not too far such that M2L can be used? */
                  if (periodic) {

                    if ((min_dist_CoM2 < max_mesh_dist2) &&
                        !(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2))
                      proxy_type |= (int)proxy_cell_type_gravity;

                  } else {

                    if (!(4. * r_max * r_max <
                          theta_crit * theta_crit * min_dist_CoM2)) {
                      proxy_type |= (int)proxy_cell_type_gravity;
                    }
                  }
                }
              }

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Add to proxies? */
              if (cells[cid].nodeID == nodeID && cells[cjd].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cells[cjd].nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cjd].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cjd].nodeID] = e->nr_proxies;
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
                proxy_addcell_in(&proxies[proxy_id], &cells[cjd], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cid], proxy_type);

                /* Store info about where to send the cell */
                cells[cid].mpi.sendto |= (1ULL << proxy_id);
              }

              /* Same for the symmetric case? */
              if (cells[cjd].nodeID == nodeID && cells[cid].nodeID != nodeID) {

                /* Do we already have a relationship with this node? */
                int proxy_id = e->proxy_ind[cells[cid].nodeID];
                if (proxy_id < 0) {
                  if (e->nr_proxies == engine_maxproxies)
                    error("Maximum number of proxies exceeded.");

                  /* Ok, start a new proxy for this pair of nodes */
                  proxy_init(&proxies[e->nr_proxies], e->nodeID,
                             cells[cid].nodeID);

                  /* Store the information */
                  e->proxy_ind[cells[cid].nodeID] = e->nr_proxies;
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
                proxy_addcell_in(&proxies[proxy_id], &cells[cid], proxy_type);
                proxy_addcell_out(&proxies[proxy_id], &cells[cjd], proxy_type);

                /* Store info about where to send the cell */
                cells[cjd].mpi.sendto |= (1ULL << proxy_id);
              }
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
