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

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Is the cell active and local? */
  if (!cell_is_active_hydro(c, e) || c->nodeID != e->nodeID)
    error("Running construction task for inactive cell!");

  /* Are we on the construction level? */
  if (c->grid.construction_level != c)
    error("Trying to build grid, but not on construction level!");

  /* Check that cells are drifted to the current time. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  struct part *parts = c->hydro.parts;
  struct delaunay *d = delaunay_malloc(c->loc, c->width, c->hydro.count);

  /* Now add ghost particles (i.e. particles from neighbouring cells and/or
   * inactive particles) */

  /* Init the list of active particles that have to be updated and their
   * search radii. */
  int *pid = NULL;
  if ((pid = (int *)malloc(c->hydro.count * sizeof(*pid))) == NULL)
    error("Can't allocate memory for pid.");
  double *search_radii = NULL;
  if ((search_radii =
           (double *)malloc(c->hydro.count * sizeof(*search_radii))) == NULL)
    error("Can't allocate memory for search radii.");

  int count = 0;
  float h_max = 0.f;
  for (int i = 0; i < c->hydro.count; i++)
    if (part_is_active(&parts[i], e)) {
      pid[count] = i;
      search_radii[count] = 0.;
      h_max = fmaxf(h_max, parts[i].h);
      ++count;
    }

  /* First add all the active particles to the delaunay tesselation */
  cell_add_local_parts_grid(d, c, parts, pid, count);

  /* Now add ghost particles (i.e. particles from neighbouring cells and/or
   * inactive particles) until all active particles have converged */
#ifdef SHADOWSWIFT_BVH
  struct flat_bvh *bvh = flat_bvh_malloc(count);
#else
  struct flat_bvh *bvh = NULL;
#endif
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  float h_max_unconverged = h_max, h_max_active = 0.f;
  int redo, num_reruns;
  for (num_reruns = 0; count > 0 && num_reruns < max_smoothing_iter;
       num_reruns++) {

#ifdef SHADOWSWIFT_BVH
    /* Build bvh of unconverged particles */
    flat_bvh_populate(bvh, parts, pid, count);
#endif

    /* Add ghost particles from this cell */
    cell_add_ghost_parts_grid_self(d, c, e, parts, bvh, pid, count);

    /* Add ghost particles from neighbouring cells */
    for (struct link *l = c->grid.sync_in; l != NULL; l = l->next) {
      if (l->t->type == task_type_self) continue;
      struct cell *c_in = l->t->cj;

      /* Add ghost particles from cj */
      cell_add_ghost_parts_grid_pair(d, c, c_in, e, parts, bvh, pid,
                                     h_max_unconverged, count);
    }

    if (!e->s->periodic) {
      /* Add boundary particles */
      cell_add_boundary_parts_grid(d, c, parts, pid, count);
    }

    /* Check if particles have converged */
    delaunay_get_search_radii(d, parts, pid, count, search_radii);
    redo = 0;
    h_max_unconverged = 0.f;
    for (int i = 0; i < count; i++) {
      /* Get a direct pointer on the part. */
      struct part *p = &parts[pid[i]];

      float r_new = (float)search_radii[i];
      if (r_new >= p->h) {
        /* Un-converged particle */
        p->h = fminf(1.01f * r_new, 1.2f * p->h);
        h_max_unconverged = fmaxf(h_max_unconverged, p->h);
        pid[redo] = pid[i];
        redo += 1;
      } else {
        /* Particle has converged. Add a small buffer zone to compensate for
         * particle movement in the next iteration. */
        p->h = 1.05f * r_new;
      }
      h_max_active = fmaxf(h_max_active, p->h);
    }
    count = redo;
    if (h_max_active > c->dmin) {
      error("Particle search radii grew larger than cell dimensions!");
    }
  }

  if (count) {
    warning(
        "Search radius failed to converge for the following gas "
        "particles:");
    for (int i = 0; i < count; i++) {
      struct part *p = &parts[pid[i]];
      warning("ID: %lld, search radius: %g", p->id, p->h);
    }

    error("Search radius failed to converge on %i particles.", count);
  }

  /* Finally build the voronoi grid */
  if (c->grid.voronoi == NULL) {
    c->grid.voronoi = voronoi_malloc(c->hydro.count, c->width[0]);
  } else {
    voronoi_reset(c->grid.voronoi, c->hydro.count, c->width[0]);
  }
  voronoi_build(c->grid.voronoi, d, parts);

  /* Be clean */
  delaunay_destroy(d);
#ifdef SHADOWSWIFT_BVH
  flat_bvh_destroy(bvh);
#endif
  free(pid);
  free(search_radii);

  if (timer) TIMER_TOC(timer_do_grid_construction);
}

#else
/* No grid construction needed */

__attribute__((always_inline)) INLINE static void runner_build_grid(
    struct runner *r, struct cell *c, int timer) {}

#endif

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_H
