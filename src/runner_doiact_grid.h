//
// Created by yuyttenh on 23/03/22.
//

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_H

/* Local headers. */
#ifdef MOVING_MESH
#include "active.h"
#include "cell.h"
#ifdef SHADOWSWIFT_BVH
#include "shadowswift/bvh.h"
#endif
#include "shadowswift/delaunay.h"
#include "shadowswift/shadowswift.h"
#include "timers.h"

__attribute__((always_inline)) INLINE static void runner_build_grid(
    struct runner *r, struct cell *c, int timer) {
  const struct engine *restrict e = r->e;

  TIMER_TIC;

  /* Before doing anyting, reset the grid specific timers for this cell */
  for (int i = 0; i < grid_timers_count; i++) {
    c->grid.extra_info.timers[i] = 0;
  }
  ticks grid_tic;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
#ifdef SHADOWSWIFT_FIX_PARTICLES
  if (c->grid.voronoi != NULL) return;
#endif

  /* Is the cell active and local? */
  if (!cell_is_active_hydro(c, e) || c->nodeID != e->nodeID)
    error("Running construction task for inactive cell!");

  /* Are we on the construction level? */
  if (c->grid.construction_level != c)
    error("Trying to build grid, but not on construction level!");

  /* Check that cells are drifted to the current time. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  /* Now add ghost particles (i.e. particles from neighbouring cells and/or
   * inactive particles) */

  /* Init the list of unconverged particles (Initially just the active
   * particles) that have to be updated and their search radii. */
  int *pid_unconverged = NULL;
  if ((pid_unconverged =
           (int *)malloc(c->hydro.count * sizeof(*pid_unconverged))) == NULL)
    error("Can't allocate memory for pid_unconverged.");
  double *search_radii = NULL;
  if ((search_radii =
           (double *)malloc(c->hydro.count * sizeof(*search_radii))) == NULL)
    error("Can't allocate memory for search radii.");
  /* Init the list of ghost candidate particles (initially just the inactive
   * particles) */
  int *pid_ghost_candidate = NULL;
  if ((pid_ghost_candidate = (int *)malloc(
           c->hydro.count * sizeof(*pid_ghost_candidate))) == NULL)
    error("Can't allocate memory for pid_ghost_candidate.");

  struct part *parts = c->hydro.parts;
  int count_unconverged = 0, count_ghost = 0;
  float r_max = 0.f;
  for (int i = 0; i < c->hydro.count; i++) {
#ifdef SHADOWSWIFT_FIX_PARTICLES
    /* Just add all the particles at every rebuild (since we will only be
     * building the tesselation at those rebuilds...) */
    pid_unconverged[count_unconverged] = i;
    search_radii[count_unconverged] = 0.;
    r_max = fmaxf(r_max, parts[i].geometry.search_radius);
    ++count_unconverged;
    // Nothing left to do
    continue;
#endif
    if (part_is_active(&parts[i], e)) {
      pid_unconverged[count_unconverged] = i;
      search_radii[count_unconverged] = 0.;
      r_max = fmaxf(r_max, parts[i].geometry.search_radius);
      ++count_unconverged;
    } else {
      pid_ghost_candidate[count_ghost] = i;
      count_ghost++;
      /* Invalidate delaunay vertex of inactive particles */
      parts[i].geometry.delaunay_vertex = -1;
    }
  }

  /* Now that we now the number of active particles, malloc the delaunay */
  struct delaunay *d = delaunay_malloc(c->loc, c->width, count_unconverged);

#ifdef SHADOWSWIFT_BVH
  /* Malloc the bvh */
  grid_tic = getticks();
  struct flat_bvh *bvh = flat_bvh_malloc(count_unconverged);
  /* Build bvh of unconverged particles */
  flat_bvh_populate(bvh, parts, pid_unconverged, count_unconverged);
  c->grid.extra_info.timers[bvh_construction] += getticks() - grid_tic;
#else
  struct flat_bvh *bvh = NULL;
#endif

  /* First add all the active particles to the delaunay tesselation */
  cell_add_local_parts_grid(d, c, parts, bvh, pid_unconverged,
                            count_unconverged);

  /* Now add ghost particles (i.e. particles from neighbouring cells and/or
   * inactive particles) until all active particles have converged */
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  float r_max_unconverged = r_max, r_max_active = 0.f;
  int redo, num_reruns;
  for (num_reruns = 0; count_unconverged > 0 && num_reruns < max_smoothing_iter;
       num_reruns++) {

    /* Add ghost particles from this cell */
    cell_add_ghost_parts_grid_self(d, c, e, parts, bvh, pid_ghost_candidate,
                                   count_ghost, pid_unconverged,
                                   count_unconverged);

    /* Add ghost particles from neighbouring cells */
    for (struct link *l = c->grid.sync_in; l != NULL; l = l->next) {
      if (l->t->type == task_type_self) continue;
      struct cell *c_in = l->t->cj;

      /* Add ghost particles from cj */
      cell_add_ghost_parts_grid_pair(d, c, c_in, e, parts, bvh, pid_unconverged,
                                     r_max_unconverged, count_unconverged);
    }

    if (!e->s->periodic) {
      /* Add boundary particles */
      cell_add_boundary_parts_grid(d, c, parts, pid_unconverged,
                                   count_unconverged);
    }

    /* Check if particles have converged */
    delaunay_get_search_radii(d, parts, pid_unconverged, count_unconverged,
                              search_radii);
    redo = 0;
    for (int i = 0; i < count_unconverged; i++) {
      /* Get a direct pointer on the part. */
      struct part *p = &parts[pid_unconverged[i]];

      float r_new = (float)search_radii[i];
      if (r_new >= p->geometry.search_radius) {
        /* Un-converged particle */
        p->geometry.search_radius =
            fminf(1.01f * r_new, 1.2f * p->geometry.search_radius);
        r_max_unconverged = fmaxf(r_max_unconverged, p->geometry.search_radius);
        pid_unconverged[redo] = pid_unconverged[i];
        redo += 1;
      } else {
        /* Particle has converged. Add a small buffer zone to compensate for
         * particle movement in the next iteration. */
        p->geometry.search_radius = 1.05f * r_new;
      }
      r_max_active = fmaxf(r_max_active, p->geometry.search_radius);
    }
    count_unconverged = redo;
    if (r_max_active > c->dmin) {
      error("Particle search radii grew larger than cell dimensions!");
    }

#ifdef SHADOWSWIFT_BVH
    /* rebuild bvh of unconverged particles */
    grid_tic = getticks();
    flat_bvh_populate(bvh, parts, pid_unconverged, count_unconverged);
    c->grid.extra_info.timers[bvh_rebuild] += getticks() - grid_tic;
#endif
  }

  if (count_unconverged) {
    warning(
        "Search radius failed to converge for the following gas "
        "particles:");
    for (int i = 0; i < count_unconverged; i++) {
      struct part *p = &parts[pid_unconverged[i]];
      warning("ID: %lld, search radius: %g", p->id, p->geometry.search_radius);
    }

    error("Search radius failed to converge on %i particles.",
          count_unconverged);
  }

  /* Finally build the voronoi grid */
  int n_cells = d->vertex_end - d->vertex_start;
  grid_tic = getticks();
  if (c->grid.voronoi == NULL) {
    c->grid.voronoi = voronoi_malloc(n_cells, c->width[0]);
  } else {
    voronoi_reset(c->grid.voronoi, n_cells, c->width[0]);
  }
  voronoi_build(c->grid.voronoi, d, parts);
  c->grid.extra_info.timers[voronoi_construction] += getticks() - grid_tic;

  /* Be clean */
  delaunay_destroy(d);
#ifdef SHADOWSWIFT_BVH
  flat_bvh_destroy(bvh);
#endif
  free(pid_unconverged);
  free(pid_ghost_candidate);
  free(search_radii);

  if (timer) TIMER_TOC(timer_do_grid_construction);
}

#else
/* No grid construction needed */

__attribute__((always_inline)) INLINE static void runner_build_grid(
    struct runner *r, struct cell *c, int timer) {}

#endif

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_H
