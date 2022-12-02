//
// Created by yuyttenh on 30/11/22.
//

#include "active.h"
#include "cell.h"
#include "engine.h"
#include "shadowswift/delaunay.h"
#include "shadowswift/voronoi.h"
#include "space_getsid.h"

/**
 * @brief Add ghost particles from #cell cj to a #delaunay tesselation.
 *
 * The BVH is assumed to be constructed from particles from cell ci.
 *
 * @param c The #cell containing the particles for which we are building the
 * local deltess
 * @param e The #engine
 * @param deltess The local #delaunay we are constructing
 * @param pid The local #parts we are constructing the local grid for
 * @param search_radii The current search radii of the #parts
 * @param count The number of parts in the #pid and #search_radii arrays.
 */
__attribute__((always_inline)) INLINE static void
cell_grid_add_boundary_particles(const struct cell *c, const struct engine *e,
                                 struct delaunay *deltess, const int *pid,
                                 const double *search_radii, int count) {
  /* Anything to do here? */
  if (c->grid.sid_is_inside_face_mask == 0b111111111111111111111111111ul)
    return;

  /* Loop over unconverged parts*/
  for (int i = 0; i < count; i++) {
    /* Retrieve particle */
    int p_idx = pid[i];
    struct part *p = &c->hydro.parts[p_idx];
    double radius = search_radii[i];

    /* Do we need to add a mirror of this particle as a boundary particle? */
    for (int sid = 0; sid < 27; sid++) {
      /* Inside face? */
      if (c->grid.sid_is_inside_face_mask & 1 << sid) continue;

      /* Calculate reflected coordinates of particle */
      double reflected_x[3];
      cell_reflect_coordinates(c, p->x, sid, reflected_x);

      /* Compute the pairwise distance. */
      double dx[3] = {p->x[0] - reflected_x[0], p->x[1] - reflected_x[1],
                      p->x[2] - reflected_x[2]};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < radius * radius) {
        delaunay_add_new_vertex(deltess, reflected_x[0], reflected_x[1],
                                reflected_x[2], sid, p_idx, -1, 1);
      }
    }
  }
}

/**
 * @brief Add ghost particles from #cell cj to a #delaunay tesselation.
 *
 * The BVH is assumed to be constructed from particles from cell ci.
 *
 * @param ci The #cell containing the particles for which we are building the
 * local deltess
 * @param cj The #cell containing the potential ghost particles (may be NULL!)
 * @param e The #engine
 * @param deltess The #delaunay struct to update
 * @param bvh The #BVH used for neighbour finding
 */
__attribute__((always_inline)) INLINE static void cell_grid_add_ghost_particles(
    struct cell *ci, struct cell *cj, const struct engine *e,
    struct delaunay *deltess, const struct BVH *bvh) {

  if (cj == NULL) {
    /* add local ghost particles */

    for (int i = 0; i < ci->hydro.count; i++) {
      /* Retrieve particle */
      const struct part *pi = &ci->hydro.parts[i];
      const double pix = pi->x[0];
      const double piy = pi->x[1];
      const double piz = pi->x[2];

      int j = bvh_hit(bvh, ci->hydro.parts, pix, piy, piz);
      if (j >= 0) {
        delaunay_add_new_vertex(deltess, pix, piy, piz, 13, i, -1, 0);
      }
    }
  } else {
    /* Ghost particles from other cell, i.e. we can use sorts */

    /* some usefull variables */
    int count_j = cj->hydro.count;
    const struct part *parts_j = cj->hydro.parts;

    /* Get sid and shift*/
    double shift[3] = {0.0, 0.0, 0.0};
    struct cell *ci_temp = ci;
    struct cell *cj_temp = cj;
    int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
    int flipped = ci != ci_temp;
#ifdef SWIFT_DEBUG_CHECKS
    /* Has cj been sorted? */
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
      error("Interacting unsorted cells.");
#endif
    struct sort_entry *sort_j = cell_get_hydro_sorts(cj, sid);

    /* ci on the right? */
    if (flipped) {
      /* Loop over the neighbouring particles parts_j until they are definitely
       * too far to be a candidate ghost particle. */
      double d_min = bvh_get_d(bvh, sid);
      for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d > d_min; pjd--) {
        /* Get a pointer to the jth particle. */
        int pj_idx = sort_j[pjd].i;
        const struct part *pj = &parts_j[pj_idx];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Shift pj so that it is in the frame of ci (with cj on the left) */
        const double pjx = pj->x[0] - shift[0];
        const double pjy = pj->x[1] - shift[1];
        const double pjz = pj->x[2] - shift[2];

        /* Check for a hit using the bvh */
        int j = bvh_hit(bvh, ci->hydro.parts, pjx, pjy, pjz);
        if (j >= 0) {
          delaunay_add_new_vertex(deltess, pjx, pjy, pjz, sid, pj_idx, -1, 0);
        }
      }
    }
    /* ci on the left! */
    else {
      /* Loop over the neighbouring particles parts_j until they are definitely
       * too far to be a candidate ghost particle. */
      sid = 26 - sid;
      double d_max = bvh_get_d(bvh, sid);
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < d_max; pjd++) {
        /* Get a pointer to the jth particle. */
        int pj_idx = sort_j[pjd].i;
        const struct part *pj = &parts_j[pj_idx];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Shift pj so that it is in the frame of ci (with cj on the left) */
        const double pjx = pj->x[0] + shift[0];
        const double pjy = pj->x[1] + shift[1];
        const double pjz = pj->x[2] + shift[2];

        /* Check for a hit using the bvh */
        int j = bvh_hit(bvh, ci->hydro.parts, pjx, pjy, pjz);
        if (j >= 0) {
          delaunay_add_new_vertex(deltess, pjx, pjy, pjz, sid, pj_idx, -1, 0);
        }
      }
    } /* cj on the right? */
  }   /* cj == NULL? */
}

/**
 * @brief Construct the local voronoi grid around a given set of particles.
 *
 * This method constructs the voronoi cells of the direct neighbours of the
 * #part pi, which may or may not be part of #cell c.
 *
 * @param c The #cell for which we want to construct the local #voronoi
 * @param e The #engine of the current node.
 * @param pid The indices of the #parts to construct the local #voronoi around
 * @param count The number of #parts to construct the local #voronoi around
 * @param vortess (return) Pointer to the empty voronoi struct we want to
 * build
 * @param parts_out (return) Empty buffer of #parts of size #count to store
 * new geometries in.
 */
__attribute__((always_inline)) INLINE void cell_grid_construct_local_voronoi(
    struct cell *c, const struct engine *e, const int *pid, const int count,
    struct voronoi *vortess, struct part *parts_out) {

  // Get some useful constants
  int periodic = e->s->periodic;

  // 1. add the local parts to the deltess
  struct delaunay *deltess = delaunay_malloc(c->loc, c->width, count);

  for (int i = 0; i < count; i++) {
    int pidx = pid[i];
    struct part *p = &c->hydro.parts[pidx];
    delaunay_add_local_vertex(deltess, i, p->x[0], p->x[1], p->x[2], -1);
  }

  // 2. add particles from c and it's neighbours to local deltess until the
  // search radius of local parts has converged

  // Get initial search radii
  double *search_radii = malloc(count * sizeof(*search_radii));
  int *idx_unconverged = malloc(count * sizeof(*idx_unconverged));
  int *pid_unconverged = malloc(count * sizeof(*pid_unconverged));
  int count_unconverged = count;
  for (int i = 0; i < count; i++) {
    search_radii[i] = c->hydro.parts[pid[i]].h;
    idx_unconverged[i] = i;
    pid_unconverged[i] = pid[i];
  }

  // add ghost particles until convergence
  struct BVH bvh;
  for (int num_iter = 0; count_unconverged > 0 && num_iter < 1000; num_iter++) {

    /* Build bvh from unconverged particles (initially all particles) */
    bvh_populate(&bvh, c->hydro.parts, search_radii, pid_unconverged,
                 count_unconverged);

    if (!periodic) {
      /* add boundary particles */
      cell_grid_add_boundary_particles(c, e, deltess, pid_unconverged,
                                       search_radii, count_unconverged);
    }

    /* We are already at the construction level, so we can run through this
     * cell's grid construction interactions directly. */
    for (struct link *l = c->grid.construction; l != NULL; l = l->next) {

      /* Skip cells that are not linked to c in the ci slot. */
      if (l->t->ci != c) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (l->t->ti_run < e->ti_current)
        error("Construction task should have been run.");
#endif

      /* Self-interaction? */
      if (l->t->type == task_type_self) {
        cell_grid_add_ghost_particles(c, NULL, e, deltess, &bvh);
      }
      /* Otherwise, pair interaction? */
      else if (l->t->type == task_type_pair) {
        cell_grid_add_ghost_particles(c, l->t->cj, e, deltess, &bvh);
        /* We only support pair and self construction tasks */
      } else {
        error("Unsupported interaction encountered!");
      }
    }

    // recompute search radii and check convergence
    int redo = 0;
    for (int i = 0; i < count_unconverged; i++) {
      int idx = idx_unconverged[i];
      double search_radius = delaunay_get_search_radius(deltess, idx);
      if (search_radius > search_radii[i]) {
        // unconverged
        search_radii[redo] = 1.2 * search_radius;
        idx_unconverged[redo] = idx;
        pid_unconverged[redo] = pid_unconverged[i];
        redo++;
      }
    }
    count_unconverged = redo;

    /* This bvh is no longer valid and will be rebuilt in the next iteration */
    bvh_clear(&bvh);
  }
  if (count_unconverged) {
    error("Failed to converge in local delaunay construction for %i particles!",
          count_unconverged);
  }

  // 3. Construct voronoi cells for the local parts.
  // TODO get rid of part_is_active
  int *part_is_active = malloc(count * sizeof(*part_is_active));
  for (int i = 0; i < count; i++) part_is_active[i] = 1;
  voronoi_build(vortess, deltess, parts_out, part_is_active, count);

  // be clean
  delaunay_destroy(deltess);
  free(search_radii);
  free(idx_unconverged);
  free(pid_unconverged);
  free(part_is_active);
}

void cell_grid_doself_apoptosis(struct cell *c, struct engine *e, int p_idx) {
  error("UNIMPLEMENTED");
}

void cell_grid_dopair_apoptosis(struct cell *restrict ci,
                                struct cell *restrict cj, struct engine *e,
                                int p_idx) {
  error("UNIMPLEMENTED");
}