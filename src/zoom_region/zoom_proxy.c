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
#include "engine.h"
#include "proxy.h"
#include "zoom.h"

/**
 * @brief Create and fill the zoom proxies.
 *
 * Handles:
 *  - zoom -> zoom.
 *  - zoom -> buffer (if using buffer cells).
 *  - zoom -> bkg (if not using buffer cells).
 *
 * @param e The #engine.
 */
void zoom_make_zoom_proxies(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->zoom_props->zoom_cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const int cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                       s->zoom_props->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = 0; /* Zoom cells are never periodic. */
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
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
    const double distance = 2. * r_max / theta_crit;
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

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid_offset(cdim, /*offset*/ 0, i, j, k);

        /* Loop over all its neighbours neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;
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

              /* Get the proxy type. */
              int proxy_type = engine_get_proxy_type(
                  cells, i, j, k, iii, jjj, kkk, cdim, with_hydro, with_gravity,
                  cid, cjd, dim, periodic, r_max, max_mesh_dist2, theta_crit,
                  /*do_direct_check*/ 1);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Ok, we need to add a proxy. */
              engine_add_proxy(e, cells, proxies, cid, cjd, proxy_type, nodeID);
            }
          }
        }

        /* We also need a proxy for all neighbouring cells. */
        for (int cjd = 0; cjd < s->zoom_props->nr_neighbour_cells; cjd++) {

          /* Get the cell. */
          struct cell *cj =
              &s->zoom_props
                   ->cells_top[s->zoom_props->neighbour_cells_top[cjd]];

          /* Early abort (both same node) */
          if (cells[cid].nodeID == nodeID && cj->nodeID == nodeID) continue;

          /* Early abort (both foreign node) */
          if (cells[cid].nodeID != nodeID && cj->nodeID != nodeID) continue;

          /* Ok, we need to add a proxy. */
          engine_add_proxy(e, cells, proxies, cid, cjd, proxy_cell_type_gravity,
                           nodeID);
        }
      }
    }
  }

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#endif
}

/**
 * @brief Create and fill the background proxies.
 *
 * Handles:
 *  - bkg -> bkg.
 *  - bkg -> buffer (if using buffer cells).
 *
 * @param e The #engine.
 */
void zoom_make_bkg_proxies(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->zoom_props->bkg_cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = 0; /* Zoom cells are never periodic. */
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  double r_diag2 = cell_width[0] * cell_width[0] +
                   cell_width[1] * cell_width[1] +
                   cell_width[2] * cell_width[2];
  double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  double r_max = 2 * r_diag;

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
    const double distance = 2. * r_max / theta_crit;
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

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid =
            cell_getid_offset(cdim, s->zoom_props->bkg_cell_offset, i, j, k);

        /* Loop over all its neighbours neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;
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

              /* Get the proxy type. */
              int proxy_type = engine_get_proxy_type(
                  cells, i, j, k, iii, jjj, kkk, cdim, with_hydro, with_gravity,
                  cid, cjd, dim, periodic, r_max, max_mesh_dist2, theta_crit,
                  /*do_direct_check*/ 1);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Ok, we need to add a proxy. */
              engine_add_proxy(e, cells, proxies, cid, cjd, proxy_type, nodeID);
            }
          }
        }

        /* We also need a proxy for all buffer cells (if they exist). To do this
         * we need the combined r_max. */
        r_diag2 = s->zoom_props->buffer_cells_top[0].width[0] *
                      s->zoom_props->buffer_cells_top[0].width[0] +
                  s->zoom_props->buffer_cells_top[1].width[1] *
                      s->zoom_props->buffer_cells_top[1].width[1] +
                  s->zoom_props->buffer_cells_top[2].width[2] *
                      s->zoom_props->buffer_cells_top[2].width[2];
        r_diag = 0.5 * sqrt(r_diag2);
        double buff_bkg_r_max = (r_max / 2) + r_diag;

        for (int cjd = 0; cjd < s->zoom_props->nr_buffer_cells; cjd++) {

          /* Get the cell. */
          struct cell *cj = &s->zoom_props->buffer_cells_top[cjd];

          /* Early abort (both same node) */
          if (cells[cid].nodeID == nodeID && cj->nodeID == nodeID) continue;

          /* Early abort (both foreign node) */
          if (cells[cid].nodeID != nodeID && cj->nodeID != nodeID) continue;

          /* Get the proxy type. */
          int proxy_type = engine_get_proxy_type(
              cells, i, j, k, iii, jjj, kkk, cdim, with_hydro, with_gravity,
              cid, cjd, dim, periodic, r_max, max_mesh_dist2, theta_crit,
              /*do_direct_check*/ 0);

          /* Abort if not in range at all */
          if (proxy_type == proxy_cell_type_none) continue;

          /* Ok, we need to add a proxy. */
          engine_add_proxy(e, cells, proxies, cid, cjd, proxy_type, nodeID);
        }
      }
    }
  }

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#endif
}

/**
 * @brief Create and fill the background proxies.
 *
 * Handles:
 *  - buffer - > buffer.
 *
 * All others are handled by the other functions already.
 *
 * @param e The #engine.
 */
void zoom_make_buffer_proxies(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  /* Useful local information */
  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  /* Handle on the cells and proxies */
  struct cell *cells = s->zoom_props->bkg_cells_top;
  struct proxy *proxies = e->proxies;

  /* Some info about the domain */
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int periodic = 0; /* Zoom cells are never periodic. */
  const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
                                cells[0].width[2]};

  /* Get some info about the physics */
  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_gravity = (e->policy & engine_policy_self_gravity);
  const double theta_crit = e->gravity_properties->theta_crit;
  const double max_mesh_dist = e->mesh->r_cut_max;
  const double max_mesh_dist2 = max_mesh_dist * max_mesh_dist;

  /* Distance between centre of the cell and corners */
  double r_diag2 = cell_width[0] * cell_width[0] +
                   cell_width[1] * cell_width[1] +
                   cell_width[2] * cell_width[2];
  double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from shifted CoM to any corner */
  double r_max = 2 * r_diag;

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
    const double distance = 2. * r_max / theta_crit;
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

  /* Loop over each cell in the space. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid =
            cell_getid_offset(cdim, s->zoom_props->bkg_cell_offset, i, j, k);

        /* Loop over all its neighbours neighbours in range. */
        for (int ii = -delta_m; ii <= delta_p; ii++) {
          int iii = i + ii;
          if (iii < 0 || iii >= cdim[0]) continue;
          iii = (iii + cdim[0]) % cdim[0];
          for (int jj = -delta_m; jj <= delta_p; jj++) {
            int jjj = j + jj;
            if (jjj < 0 || jjj >= cdim[1]) continue;
            jjj = (jjj + cdim[1]) % cdim[1];
            for (int kk = -delta_m; kk <= delta_p; kk++) {
              int kkk = k + kk;
              if (kkk < 0 || kkk >= cdim[2]) continue;
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

              /* Get the proxy type. */
              int proxy_type = engine_get_proxy_type(
                  cells, i, j, k, iii, jjj, kkk, cdim, with_hydro, with_gravity,
                  cid, cjd, dim, periodic, r_max, max_mesh_dist2, theta_crit,
                  /*do_direct_check*/ 1);

              /* Abort if not in range at all */
              if (proxy_type == proxy_cell_type_none) continue;

              /* Ok, we need to add a proxy. */
              engine_add_proxy(e, cells, proxies, cid, cjd, proxy_type, nodeID);
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
#endif
}

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
  /* Do the zoom cells. */
  zoom_make_zoom_proxies(e);

  /* Do the background cells. */
  zoom_make_bkg_proxies(e);

  /* Do the buffer cells (if we have them). */
  if (e->s->zoom_props->with_buffer_cells) zoom_make_buffer_proxies(e);
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
