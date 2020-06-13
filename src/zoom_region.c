/* Config parameters. */
#include "../config.h"

#include <float.h>

#include "cell.h"
#include "engine.h"
#include "proxy.h"
#include "space.h"

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
            s->zoom_props->loc[1], s->zoom_props->loc[2], s->zoom_props->dim[0],
            s->zoom_props->dim[1], s->zoom_props->dim[2]);
}

void check_zoom_region(const struct space *s, const int verbose) {
  if (verbose) {
    message("top level zoom cell dimensions [ %i %i %i ].",
            s->zoom_props->cdim[0], s->zoom_props->cdim[1],
            s->zoom_props->cdim[2]);
    message(
        "zoom region centered on [ %.3f %.3f %.3f ] with dimensions [ %.3f "
        "%.3f %.3f ]",
        s->zoom_props->loc[0], s->zoom_props->loc[1], s->zoom_props->loc[2],
        s->zoom_props->dim[0], s->zoom_props->dim[1], s->zoom_props->dim[2]);
  }
}

int cell_getid_zoom(const int cdim[3], const double x, const double y,
                    const double z,
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
  if (x > zoom_region_bounds[0] && x < zoom_region_bounds[1] &&
      y > zoom_region_bounds[2] && y < zoom_region_bounds[3] &&
      z > zoom_region_bounds[4] && z < zoom_region_bounds[5]) {

    const int zoom_index =
        cell_getid(cdim, (x - zoom_region_bounds[0]) * ih_x_zoom,
                   (y - zoom_region_bounds[2]) * ih_y_zoom,
                   (z - zoom_region_bounds[4]) * ih_z_zoom);
    cell_id = zoom_cell_offset + zoom_index;
#ifdef SWIFT_DEBUG_CHECKS
    if (zoom_index < 0 || zoom_index >= cdim[0] * cdim[1] * cdim[2])
      error("zoom_index out of range %i (%f %f %f)", cell_id, x, y, z);
#endif
  } else {
    cell_id = cell_getid(cdim, i, j, k);
#ifdef SWIFT_DEBUG_CHECKS
    if (cell_id < 0 || cell_id >= cdim[0] * cdim[1] * cdim[2])
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
  // const int periodic = s->periodic;
  struct cell *cells = s->cells_top;
  const int nr_cells = cdim[0] * cdim[1] * cdim[2];
  const int nr_zoom_cells =
      s->zoom_props->cdim[0] * s->zoom_props->cdim[1] * s->zoom_props->cdim[2];

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
  // const int with_gravity = (e->policy & engine_policy_self_gravity);
  int proxy_type = (int)proxy_cell_type_gravity;

  /* Prepare the proxies and the proxy index. */
  if (e->proxy_ind == NULL) error("Why hasn't proxy_ind been allocated yet?");

  /* Loop over each zoom top level cell in the space. */
  for (int i = 0; i < nr_zoom_cells; i++) {

    /* Get the zoom top level cell. */
    struct cell *ci = &cells[i + s->zoom_props->tl_cell_offset];

    /* Loop over all non top level zoom cells. */
    for (int j = 0; j < nr_cells; j++) {

      /* Get the cell. */
      struct cell *cj = &cells[j];

      if (cj->tl_cell_type != tl_cell_neighbour) continue;

      /* Early abort (both same node) */
      if (ci->nodeID == nodeID && cj->nodeID == nodeID) continue;

      /* Early abort (both foreign node) */
      if (ci->nodeID != nodeID && cj->nodeID != nodeID) continue;

      /* Add to proxies? */
      if (ci->nodeID == nodeID && cj->nodeID != nodeID) {
        /* Do we already have a relationship with this node? */
        int proxy_id = e->proxy_ind[cj->nodeID];
        if (proxy_id < 0) {
          if (e->nr_proxies == engine_maxproxies)
            error("Maximum number of proxies exceeded.");

          /* Ok, start a new proxy for this pair of nodes */
          proxy_init(&proxies[e->nr_proxies], e->nodeID, cj->nodeID);

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
  }

  /* Be clear about the time */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

double cell_min_dist2_diff_size(const struct cell *restrict ci,
                                const struct cell *restrict cj,
                                const int periodic, const double dim[3]) {

  const double cix = ci->loc[0] + ci->width[0] / 2.;
  const double ciy = ci->loc[1] + ci->width[1] / 2.;
  const double ciz = ci->loc[2] + ci->width[2] / 2.;
  const double ci_diag2 = ci->width[0] / 2. * ci->width[0] / 2. +
                          ci->width[1] / 2. * ci->width[1] / 2. +
                          ci->width[2] / 2. * ci->width[2] / 2.;

  const double cjx = cj->loc[0] + cj->width[0] / 2.;
  const double cjy = cj->loc[1] + cj->width[1] / 2.;
  const double cjz = cj->loc[2] + cj->width[2] / 2.;
  const double cj_diag2 = cj->width[0] / 2. * cj->width[0] / 2. +
                          cj->width[1] / 2. * cj->width[1] / 2. +
                          cj->width[2] / 2. * cj->width[2] / 2.;

  if (periodic) {

    const double dx = nearest(cix - cjx, dim[0]);
    const double dy = nearest(ciy - cjy, dim[1]);
    const double dz = nearest(ciz - cjz, dim[2]);

    const double dr2 = (dx * dx + dy * dy + dz * dz) - ci_diag2 - cj_diag2;
    return max(0.0, dr2);

  } else {
    const double dx = cix - cjx;
    const double dy = ciy - cjy;
    const double dz = ciz - cjz;

    const double dr2 = (dx * dx + dy * dy + dz * dz) - ci_diag2 - cj_diag2;
    return max(0.0, dr2);
  }
}

double cell_min_dist2(const struct cell *restrict ci,
                      const struct cell *restrict cj, const int periodic,
                      const double dim[3]) {

  double dist2;

  if (ci->tl_cell_type <= 1 && cj->tl_cell_type <= 1) {
    dist2 = cell_min_dist2_same_size(ci, cj, periodic, dim);
  } else if (ci->tl_cell_type == zoom_tl_cell &&
             cj->tl_cell_type == zoom_tl_cell) {
    dist2 = cell_min_dist2_same_size(ci, cj, 0, dim);
  } else {
    dist2 = cell_min_dist2_diff_size(ci, cj, periodic, dim);
  }

  return dist2;
}

void engine_make_self_gravity_tasks_mapper_between_toplevels(void *map_data,
                                                             int num_elements,
                                                             void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  struct cell *cells = s->cells_top;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;
  const int nr_cells = s->nr_cells;

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    if (cid < s->zoom_props->tl_cell_offset) continue;

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    for (int cjd = 0; cjd < nr_cells; cjd++) {

      /* Get the cell */
      struct cell *cj = &cells[cjd];

      if (cj->tl_cell_type != tl_cell_neighbour) continue;

      /* Avoid duplicates, empty cells and completely foreign pairs */
      if (cj->grav.count == 0 || (ci->nodeID != nodeID && cj->nodeID != nodeID))
        continue;

      /* Recover the multipole information */
      const struct gravity_tensors *multi_i = ci->grav.multipole;
      const struct gravity_tensors *multi_j = cj->grav.multipole;

      if (multi_i == NULL && ci->nodeID != nodeID)
        error("Multipole of ci was not exchanged properly via the proxies");
      if (multi_j == NULL && cj->nodeID != nodeID)
        error("Multipole of cj was not exchanged properly via the proxies");

      /* Minimal distance between any pair of particles */
      const double min_radius2 =
          cell_min_dist2_diff_size(ci, cj, periodic, dim);

      /* Are we beyond the distance where the truncated forces are 0 ?*/
      if (periodic && min_radius2 > max_distance2) continue;

      /* Are the cells too close for a MM interaction ? */
      if (!cell_can_use_pair_mm_rebuild(ci, cj, e, s, s->periodic)) {

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI
        /* Let's cross-check that we had a proxy for that cell */
        if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[cj->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", cj->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == cj) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cjd);
        } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

          /* Find the proxy for this node */
          const int proxy_id = e->proxy_ind[ci->nodeID];
          if (proxy_id < 0)
            error("No proxy exists for that foreign node %d!", ci->nodeID);

          const struct proxy *p = &e->proxies[proxy_id];

          /* Check whether the cell exists in the proxy */
          int n = 0;
          for (; n < p->nr_cells_in; n++)
            if (p->cells_in[n] == ci) {
              break;
            }
          if (n == p->nr_cells_in)
            error(
                "Cell %d not found in the proxy but trying to construct "
                "grav task!",
                cid);
        }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */

        // message("%i adding task", nodeID);
        /* Ok, we need to add a direct pair calculation */
        scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0, ci,
                          cj);
      }
    }
  }
}

void engine_make_self_gravity_tasks_mapper_zoom(void *map_data,
                                                int num_elements,
                                                void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  struct scheduler *sched = &e->sched;
  const int nodeID = e->nodeID;
  const int periodic = 0;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  struct cell *cells = s->cells_top;
  const double theta_crit_inv = e->gravity_properties->theta_crit_inv;
  // const double theta_crit = e->gravity_properties->theta_crit;
  const double max_distance = e->mesh->r_cut_max;
  const double max_distance2 = max_distance * max_distance;
  const int offset = s->zoom_props->tl_cell_offset;
  const double cell_width[3] = {cells[offset].width[0], cells[offset].width[1],
                                cells[offset].width[2]};
  const double dmin = cells[offset].dmin;

  /* Distance between centre of the cell and corners */
  const double r_diag2 = cell_width[0] * cell_width[0] +
                         cell_width[1] * cell_width[1] +
                         cell_width[2] * cell_width[2];
  const double r_diag = 0.5 * sqrt(r_diag2);

  /* Maximal distance from a shifted CoM to centre of cell */
  const double delta_CoM = engine_max_proxy_centre_frac * r_diag;

  /* Maximal distance from shifted CoM to any corner */
  const double r_max = r_diag + 2. * delta_CoM;

  /* Compute how many cells away we need to walk */
  const double distance = 2. * r_max * theta_crit_inv;
  const int delta = (int)(distance / dmin) + 1;
  int delta_m = delta;
  int delta_p = delta;

  /* Special case where every cell is in range of every other one */
  if (delta >= cdim[0] / 2) {
    if (cdim[0] % 2 == 0) {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2 - 1;
    } else {
      delta_m = cdim[0] / 2;
      delta_p = cdim[0] / 2;
    }
  }

  /* Loop through the elements, which are just byte offsets from NULL. */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Get the cell index. */
    const int cid = (size_t)(map_data) + ind;

    if (cid < offset) continue;

    const int new_cid = cid - offset;

    /* Integer indices of the cell in the top-level grid */
    const int i = new_cid / (cdim[1] * cdim[2]);
    const int j = (new_cid / cdim[2]) % cdim[1];
    const int k = new_cid % cdim[2];

    /* Get the cell */
    struct cell *ci = &cells[cid];

    /* Skip cells without gravity particles */
    if (ci->grav.count == 0) continue;

    /* If the cell is local build a self-interaction */
    if (ci->nodeID == nodeID) {
      scheduler_addtask(sched, task_type_self, task_subtype_grav, 0, 0, ci,
                        NULL);
    }

    /* Loop over every other cell within (Manhattan) range delta */
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

          /* Get the cell */
          const int cjd = cell_getid(cdim, iii, jjj, kkk);
          const int new_cjd = cjd + offset;

          struct cell *cj = &cells[new_cjd];

          /* Avoid duplicates, empty cells and completely foreign pairs */
          if (new_cid >= cjd || cj->grav.count == 0 ||
              (ci->nodeID != nodeID && cj->nodeID != nodeID))
            continue;

          /* Recover the multipole information */
          const struct gravity_tensors *multi_i = ci->grav.multipole;
          const struct gravity_tensors *multi_j = cj->grav.multipole;

          if (multi_i == NULL && ci->nodeID != nodeID)
            error("Multipole of ci was not exchanged properly via the proxies");
          if (multi_j == NULL && cj->nodeID != nodeID)
            error("Multipole of cj was not exchanged properly via the proxies");

          /* Minimal distance between any pair of particles */
          const double min_radius2 =
              cell_min_dist2_same_size(ci, cj, periodic, dim);

          /* Are we beyond the distance where the truncated forces are 0 ?*/
          if (min_radius2 > max_distance2) continue;

          /* Are the cells too close for a MM interaction ? */
          if (!cell_can_use_pair_mm_rebuild(ci, cj, e, s, periodic)) {
#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI

            /* Let's cross-check that we had a proxy for that cell */
            if (ci->nodeID == nodeID && cj->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[cj->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", cj->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == cj) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cjd);
            } else if (cj->nodeID == nodeID && ci->nodeID != engine_rank) {

              /* Find the proxy for this node */
              const int proxy_id = e->proxy_ind[ci->nodeID];
              if (proxy_id < 0)
                error("No proxy exists for that foreign node %d!", ci->nodeID);

              const struct proxy *p = &e->proxies[proxy_id];

              /* Check whether the cell exists in the proxy */
              int n = 0;
              for (; n < p->nr_cells_in; n++)
                if (p->cells_in[n] == ci) {
                  break;
                }
              if (n == p->nr_cells_in)
                error(
                    "Cell %d not found in the proxy but trying to construct "
                    "grav task!",
                    cid);
            }
#endif /* WITH_MPI */
#endif /* SWIFT_DEBUG_CHECKS */

            /* Ok, we need to add a direct pair calculation */
            scheduler_addtask(sched, task_type_pair, task_subtype_grav, 0, 0,
                              ci, cj);
          }
        }
      }
    }
  }
}

void find_neighbouring_cells(struct space *s, const int verbose) {

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const int periodic = s->periodic;
  struct cell *cells = s->cells_top;

  const int delta_m = 1;
  const int delta_p = 1;

  int neighbour_count = 0;

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

              if (cells[cjd].tl_cell_type == tl_cell) {

                /* Record that we've found a neighbour. */
                cells[cjd].tl_cell_type = tl_cell_neighbour;
                neighbour_count++;
              }
            }
          }
        }
      }
    }
  }

  if (verbose) message("%i neighbour cells found", neighbour_count);
}
