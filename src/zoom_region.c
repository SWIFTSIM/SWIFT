/* Config parameters. */
#include "../config.h"

#include <float.h>

#include "engine.h"
#include "space.h"
#include "cell.h"
#include "proxy.h"

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

void zoom_region_init(struct swift_params *params, struct space *s) {

  /* Initialise the zoom region. */
  s->with_zoom_region =
      parser_get_opt_param_int(params, "ZoomRegion:enable", 0);

  if (s->with_zoom_region) {
    s->zoom_props = (struct zoom_region_properties *)malloc(
        sizeof(struct zoom_region_properties));
    if (s->zoom_props == NULL)
      error("Error allocating memory for the zoom parameters.");

    s->zoom_props->mask_parttype =
        parser_get_param_int(params, "ZoomRegion:mask_parttype");
    s->zoom_props->max_size =
        parser_get_opt_param_double(params, "ZoomRegion:max_size", FLT_MAX);
    s->zoom_props->boost_factor =
        parser_get_opt_param_double(params, "ZoomRegion:boost_factor", 1.0);
  }
}

void construct_zoom_region(struct space *s, int verbose) {

  double new_zoom_boundary[6] = {1e20, -1e20, 1e20, -1e20, 1e20, -1e20};
  const size_t nr_gparts = s->nr_gparts;
  double mtot = 0.0;
  double com[3] = {0.0, 0.0, 0.0};

  for (size_t k = 0; k < nr_gparts; k++) {
    if (s->gparts[k].type != s->zoom_props->mask_parttype) continue;

    if (s->gparts[k].x[0] < new_zoom_boundary[0])
      new_zoom_boundary[0] = s->gparts[k].x[0];
    if (s->gparts[k].x[0] > new_zoom_boundary[1])
      new_zoom_boundary[1] = s->gparts[k].x[0];
    if (s->gparts[k].x[1] < new_zoom_boundary[2])
      new_zoom_boundary[2] = s->gparts[k].x[1];
    if (s->gparts[k].x[1] > new_zoom_boundary[3])
      new_zoom_boundary[3] = s->gparts[k].x[1];
    if (s->gparts[k].x[2] < new_zoom_boundary[4])
      new_zoom_boundary[4] = s->gparts[k].x[2];
    if (s->gparts[k].x[2] > new_zoom_boundary[5])
      new_zoom_boundary[5] = s->gparts[k].x[2];

    mtot += s->gparts[k].mass;
    com[0] += s->gparts[k].x[0] * s->gparts[k].mass;
    com[1] += s->gparts[k].x[1] * s->gparts[k].mass;
    com[2] += s->gparts[k].x[2] * s->gparts[k].mass;
  }
#ifdef WITH_MPI
  /* Make sure all nodes agree on the zoom region. */
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[0], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[1], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[2], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[3], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);   
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[4], 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &new_zoom_boundary[5], 1, MPI_DOUBLE, MPI_MAX,
                MPI_COMM_WORLD);

  /* CoM. */
  MPI_Allreduce(MPI_IN_PLACE, com, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  const double imass = 1.0 / mtot;
  com[0] *= imass;
  com[1] *= imass;
  com[2] *= imass;

  for (int k = 0; k < 3; k++) s->zoom_props->loc[k] = com[k];

  s->zoom_props->dim[0] = (new_zoom_boundary[1] - new_zoom_boundary[0]) *
                          s->zoom_props->boost_factor;
  s->zoom_props->dim[1] = (new_zoom_boundary[3] - new_zoom_boundary[2]) *
                          s->zoom_props->boost_factor;
  s->zoom_props->dim[2] = (new_zoom_boundary[5] - new_zoom_boundary[4]) *
                          s->zoom_props->boost_factor;

  for (int k = 0; k < 3; k++)
    if (s->zoom_props->dim[k] > s->zoom_props->max_size)
      s->zoom_props->dim[k] = s->zoom_props->max_size;

  if (verbose)
    message("zoom loc: [%f %f %f] zoom dim: [%f %f %f]", s->zoom_props->loc[0],
            s->zoom_props->loc[1], s->zoom_props->loc[2],
            s->zoom_props->dim[0], s->zoom_props->dim[1],
            s->zoom_props->dim[2]);
}

void check_zoom_region(const struct space *s, const int verbose) {
  if (verbose) {
    message("top level zoom cell dimensions [ %i %i %i ].", s->zoom_props->cdim[0],
			s->zoom_props->cdim[1], s->zoom_props->cdim[2]);
    message(
        "zoom region centered on [ %.3f %.3f %.3f ] with dimensions [ %.3f "
        "%.3f %.3f ]",
        s->zoom_props->loc[0], s->zoom_props->loc[1], s->zoom_props->loc[2],
        s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  }
}

int cell_getid_zoom(const int cdim[3], const double x, const double y, const double z,
	const struct zoom_region_properties *zoom_props,
	const int i, const int j, const int k) {

  /* Properties of the zoom region. */
  const int zoom_cell_offset = zoom_props->tl_cell_offset;
  const double zoom_region_bounds[6] = {
      zoom_props->region_bounds[0], zoom_props->region_bounds[1],
      zoom_props->region_bounds[2], zoom_props->region_bounds[3],
      zoom_props->region_bounds[4], zoom_props->region_bounds[5]};
  const double ih_x_zoom = zoom_props->iwidth[0];
  const double ih_y_zoom = zoom_props->iwidth[1];
  const double ih_z_zoom = zoom_props->iwidth[2];
  int cell_id;

  /* Are the passed coordinates within the zoom region? */
  if (x >= zoom_region_bounds[0] && x <= zoom_region_bounds[1] &&
      y >= zoom_region_bounds[2] && y <= zoom_region_bounds[3] &&
      z >= zoom_region_bounds[4] && z <= zoom_region_bounds[5]) {

    const int zoom_index =
        cell_getid(cdim, (x - zoom_region_bounds[0]) * ih_x_zoom,
                   (y - zoom_region_bounds[2]) * ih_y_zoom,
                   (z - zoom_region_bounds[4]) * ih_z_zoom);
    cell_id = zoom_cell_offset + zoom_index;
#ifdef SWIFT_DEBUG_CHECKS
  if (zoom_index < 0 || zoom_index >= cdim[0]*cdim[1]*cdim[2])
    error("zoom_index out of range %i (%f %f %f)", cell_id, x, y, z);
#endif
  } else {
    cell_id = cell_getid(cdim, i, j, k);
#ifdef SWIFT_DEBUG_CHECKS
  if (cell_id < 0 || cell_id >= cdim[0]*cdim[1]*cdim[2])
    error("cell_id out of range %i (i:%i j:%i k:%i) (x:%f y:%f z:%f)",
            cell_id, i, j, k, x, y, z);
#endif
  }

  return cell_id;
}

void makeproxies_between_top_levels(struct engine *e) {

#ifdef WITH_MPI
  /* Let's time this */
  const ticks tic = getticks();

  const int nodeID = e->nodeID;
  const struct space *s = e->s;

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const int periodic = s->periodic; 
  struct cell *cells = s->cells_top;
  const int nr_cells = cdim[0] * cdim[1] * cdim[2];
  const int nr_zoom_cells = s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2];

  const int delta_m = 5;
  const int delta_p = 5;

  int neighbour_count = 0;
  int *neighbour_cells = NULL;
  if ((neighbour_cells = (int *)malloc(sizeof(int) * nr_cells)) == NULL)
            error("Failed to neighbour cells array.");
  bzero(neighbour_cells, sizeof(int) * nr_cells);

  /* Loop over each cell in the space to find the neighbouring top level cells
   * surrounding the zoom region. */
  for (int i = 0; i < cdim[0]; i++) {
    for (int j = 0; j < cdim[1]; j++) {
      for (int k = 0; k < cdim[2]; k++) {

        /* Get the cell ID. */
        const int cid = cell_getid(cdim, i, j, k);

        /* Only interested in cells hosting zoom top level cells. */
        if (cells[cid].tl_cell_type != void_tl_cell) continue;

        /* Loop over all its direct neighbours. */
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

              /* Get the cell ID of the neighbour. */
              const int cjd = cell_getid(cdim, iii, jjj, kkk);

              /* Only interested in neighbouring cells that don't host zoom top level cells. */
              //if (cells[cjd].tl_cell_type == void_tl_cell) continue;

              /* Record that we've found a neighbour. */
              neighbour_cells[cjd] = 1;
              neighbour_count++;
            }
          }
        }
      }
    }
  }

  message("%i total neighbours found around the zoom region.", neighbour_count);
  struct proxy *proxies = e->proxies;
//
//  /* Some info about the domain */
//  //const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
//  //const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
//  const int periodic = s->periodic;
//  //const double cell_width[3] = {cells[0].width[0], cells[0].width[1],
//  //                              cells[0].width[2]};
//
  /* Get some info about the physics */
//  //const int with_hydro = (e->policy & engine_policy_hydro);
  //const int with_gravity = (e->policy & engine_policy_self_gravity);
  int proxy_type = 0;
  proxy_type |= (int)proxy_cell_type_gravity;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL) error("Why hasn't proxy_ind been allocated yet?");

  /* Loop over each zoom top level cell in the space. */
  for (int i = 0; i < nr_zoom_cells; i++) {

    /* Get the zoom top level cell. */
    struct cell *ci = &cells[i+s->zoom_props->tl_cell_offset];

    /* Loop over all non top level zoom cells. */
    for (int j = 0; j < nr_cells; j++) {

      if (neighbour_cells[j] == 0) continue;

      /* Get the cell. */
      struct cell *cj = &cells[j];

      /* Early abort (both same node) */
      if (ci->nodeID == nodeID && cj->nodeID == nodeID)
        continue;

      /* Early abort (both foreign node) */
      if (ci->nodeID != nodeID && cj->nodeID != nodeID)
        continue;

      /* Add to proxies? */
      if (ci->nodeID == nodeID && cj->nodeID != nodeID) {
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[cj->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");

          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     cj->nodeID);

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
          proxy_init(&proxies[e->nr_proxies], e->nodeID,
                     ci->nodeID);

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
  }

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

