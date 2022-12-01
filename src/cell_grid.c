//
// Created by yuyttenh on 30/11/22.
//

#include "cell.h"
#include "engine.h"
#include "shadowswift/delaunay.h"
#include "shadowswift/voronoi.h"
#include "space_getsid.h"

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
 * @param vortess (return) Pointer to the voronoi struct we want to build
 * @param parts_out (return) Empty buffer of #parts of size #count to store new
 * geometries in.
 */
__attribute__((always_inline)) INLINE void cell_construct_local_voronoi(
    struct cell *c, const struct engine *e, const int *pid,
    const int count, struct voronoi *vortess, struct part *parts_out) {

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
//      error("unimplemented");
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
        /* add local ghost particles */

        for (int i = 0; i < c->hydro.count; i++) {
          /* Retrieve particle */
          const struct part *pi = &c->hydro.parts[i];
          const double pix = pi->x[0];
          const double piy = pi->x[1];
          const double piz = pi->x[2];

          int j = bvh_hit(&bvh, c->hydro.parts, pix, piy, piz);
          if (j >= 0) {
            delaunay_add_new_vertex(deltess, pix, piy, piz, 13, i, -1, 0);
          }
        }
      }
      /* Otherwise, pair interaction? */
      else if (l->t->type == task_type_pair) {
        /* add remote ghost particles */
        struct cell *cj = l->t->cj;

        /* Get the sort ID. */
        double shift[3] = {0.0, 0.0, 0.0};
        struct cell *ci_temp = c;
        struct cell *cj_temp = cj;
        int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
        int flipped = c != ci_temp;
        if (!flipped) sid = 26 - sid;

        for (int i = 0; i < cj->hydro.count; i++) {
          /* Retrieve particle */
          const struct part *pi = &cj->hydro.parts[i];
          double pix, piy, piz;
          // TODO: clean this up + use sortlists properly?
          if (flipped) {
            pix = pi->x[0] - shift[0];
            piy = pi->x[1] - shift[1];
            piz = pi->x[2] - shift[2];
          } else {
            pix = pi->x[0] + shift[0];
            piy = pi->x[1] + shift[1];
            piz = pi->x[2] + shift[2];
          }

          int j = bvh_hit(&bvh, c->hydro.parts, pix, piy, piz);
          if (j >= 0) {
            delaunay_add_new_vertex(deltess, pix, piy, piz, sid, i, -1, 0);
          }
        }
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