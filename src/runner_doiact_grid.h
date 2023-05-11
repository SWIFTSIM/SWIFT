//
// Created by yuyttenh on 23/03/22.
//

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_H

/* Local headers. */
#include "part.h"
#ifdef MOVING_MESH
#include "active.h"
#include "cell.h"
#include "shadowswift/bvh.h"
#include "shadowswift/delaunay.h"
#include "shadowswift/shadowswift.h"
#include "space_getsid.h"
#include "timers.h"

#ifdef SHADOWSWIFT_HILBERT_ORDERING
#include "shadowswift/hilbert.h"

/*! @brief Calculate the hilbert keys of the vertices
 *
 * @param c Cell containing the vertices
 */
__attribute__((always_inline)) INLINE static void get_hilbert_keys(
    const struct cell *restrict c, unsigned long *restrict keys) {
  /* TODO only calculate (and later sort) the keys for the active particles! */
  for (int i = 0; i < c->hydro.count; i++) {
#if defined(HYDRO_DIMENSION_1D)
    /* In 1D, we just use the mantissa of the rescaled x coordinate to sort the
     * particles. */
    /* Get rescaled x coordinate within [1, 2) */
    double x = c->hydro.parts[i].x[0];
    double x_scaled = (x - c->loc[0] + c->width[0]) / (3. * c->width[0]) + 1;
    keys[i] = delaunay_double_to_int(x_scaled);
#elif defined(HYDRO_DIMENSION_2D)
    float dx_max = c->hydro.dx_max_part;
    unsigned long bits[2];
    int nbits = 32;
    double max_width = max(c->width[0], c->width[1]) + 2 * dx_max;
    bits[0] = (c->hydro.parts[i].x[0] + dx_max - c->loc[0]) / max_width *
              (1ul << (nbits - 1));
    bits[1] = (c->hydro.parts[i].x[1] + dx_max - c->loc[1]) / max_width *
              (1ul << (nbits - 1));
    keys[i] = hilbert_get_key_2d(bits, nbits);
#elif defined(HYDRO_DIMENSION_3D)
    float dx_max = c->hydro.dx_max_part;
    unsigned long bits[3];
    int nbits = 21;
    double max_width = max3(c->width[0], c->width[1], c->width[2]) + 2 * dx_max;
    bits[0] = (c->hydro.parts[i].x[0] + dx_max - c->loc[0]) / max_width *
              (1ul << (nbits - 1));
    bits[1] = (c->hydro.parts[i].x[1] + dx_max - c->loc[1]) / max_width *
              (1ul << (nbits - 1));
    bits[2] = (c->hydro.parts[i].x[2] + dx_max - c->loc[2]) / max_width *
              (1ul << (nbits - 1));
    keys[i] = hilbert_get_key_3d(bits, nbits);
#endif
  }
}
#endif

__attribute__((always_inline)) INLINE static void runner_build_bvh(
    struct runner *r, struct cell *c, int timer) {
#ifdef SHADOWSWIFT_BVH
  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grid.construction_level != on_construction_level)
    error("Trying to build bvh, but not on construction level!");
#endif

  const struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  int count = c->hydro.count;

  /* Allocate array with pids of active particles */
  int *pid_active = malloc(count * sizeof(int));
  if (pid_active == NULL)
    error("Failed to allocate memory for active particle pids!");
  int count_active = 0;
  for (int i = 0; i < count; i++) {
    if (part_is_active(&parts[i], e)) {
      pid_active[count_active] = i;
      count_active++;
    }
  }

  /* Construct the bvh */
  struct BVH *bvh = malloc(sizeof(*bvh));
  bvh_populate(bvh, parts, pid_active, count_active, count);
  c->grid.bvh = bvh;

  /* Be clean */
  free(pid_active);

  if (timer) TIMER_TOC(timer_bvh);
#endif
}

__attribute((always_inline)) INLINE static void
runner_dopair_grid_construction_naive(
    struct cell *restrict ci, const struct cell *restrict cj,
    const struct engine *e, const struct sort_entry *restrict sort_i,
    const struct sort_entry *restrict sort_j, int flipped, const double *shift,
    const double rshift, const double hi_max, const double dx_max, int sid) {

  /* Get some useful values. */
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  if (flipped) { /* cj on the left */

    const double di_min = sort_i[0].d - hi_max - dx_max + rshift;

    /* Loop over parts in cj */
    for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d > di_min; pjd--) {

      /* Recover pj */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the left) */
      const double pjx = pj->x[0] - shift[0];
      const double pjy = pj->x[1] - shift[1];
      const double pjz = pj->x[2] - shift[2];

      delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx, -1,
                              0);
      /* Update delaunay flags to signal that the particle was added for
       * this sid */
      pj->geometry.delaunay_flags |= 1 << sid;
    }
  } else {
    const double di_max = sort_i[count_i - 1].d + hi_max + dx_max - rshift;

    /* Loop over parts in cj */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di_max; pjd++) {

      /* Recover pj */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the right) */
      const double pjx = pj->x[0] + shift[0];
      const double pjy = pj->x[1] + shift[1];
      const double pjz = pj->x[2] + shift[2];

      delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx, -1,
                              0);
      /* Update delaunay flags to signal that the particle was added for
       * this sid */
      pj->geometry.delaunay_flags |= 1 << sid;
    }
  }
}

__attribute((always_inline)) INLINE static void runner_dopair_grid_construction(
    struct cell *restrict ci, const struct cell *restrict cj,
    const struct engine *e, const struct sort_entry *restrict sort_i,
    const struct sort_entry *restrict sort_j,
    const struct sort_entry *restrict sort_active_i, int count_active,
    int flipped, const double *shift, const double rshift, const double hi_max,
    const double dx_max, int sid) {

#ifdef SHADOWSWIFT_BVH
  /* Get some useful values. */
  struct part *restrict parts_i = ci->hydro.parts;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  if (flipped) {
    /* Loop over the parts in cj (on the left) */
    for (int pjd = count_j - 1; pjd >= 0; pjd--) {
      /* Recover pj */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip particles that are already added */
      if (pj->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the left) */
      const double pjx = pj->x[0] - shift[0];
      const double pjy = pj->x[1] - shift[1];
      const double pjz = pj->x[2] - shift[2];

      /* If pj is no longer contained in the bbox of all the active particles,
       * we are done here. */
      if (!bbox_contains(&ci->grid.bvh->bbox, pjx, pjy, pjz)) break;

      /* Check if there is an active particle of ci that contains pj in its
       * search radius */
      int pi_idx = bvh_hit(ci->grid.bvh, parts_i, pjx, pjy, pjz);

      /* Hit or miss? */
      if (pi_idx >= 0) {
        delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                pi_idx, 0);
        /* Update delaunay flags to signal that the particle was added for
         * this sid */
        pj->geometry.delaunay_flags |= 1 << sid;
      }
    }
  } else {
    /* Loop over the parts in cj (on the right) */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Recover pj */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[sort_j[pjd].i];

      /* Skip particles that are already added */
      if (pj->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the right) */
      const double pjx = pj->x[0] + shift[0];
      const double pjy = pj->x[1] + shift[1];
      const double pjz = pj->x[2] + shift[2];

      /* If pj is no longer contained in the bbox of all the active particles,
       * we are done here. */
      if (!bbox_contains(&ci->grid.bvh->bbox, pjx, pjy, pjz)) break;

      /* Check if there is an active particle of ci that contains pj in its
       * search radius */
      int pi_idx = bvh_hit(ci->grid.bvh, parts_i, pjx, pjy, pjz);

      /* Hit or miss? */
      if (pi_idx >= 0) {
        delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                pi_idx, 0);
        /* Update delaunay flags to signal that the particle was added for
         * this sid */
        pj->geometry.delaunay_flags |= 1 << sid;
      }
    }
  }

#else

  /* Get some useful values. */
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  if (flipped) {
    /* ci on the right */

    /* Loop over the sorted active parts in ci (on the right) */
    for (int i = 0; i < count_active; i++) {

      /* Get a hold of pi. */
      struct sort_entry sort_pi = sort_active_i[i];
      struct part *restrict pi = &parts_i[sort_pi.i];
      const float ri = pi->h;

      const double dj_min = sort_pi.d - ri - dx_max + rshift;

      /* Loop over the parts in cj (on the left) */
      for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d > dj_min; pjd--) {

        /* Recover pj */
        int pj_idx = sort_j[pjd].i;
        struct part *restrict pj = &parts_j[pj_idx];

        /* Skip particles that are already added */
        if (pj->geometry.delaunay_flags & 1 << sid) continue;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Shift pj so that it is in the frame of ci (with cj on the left) */
        const double pjx = pj->x[0] - shift[0];
        const double pjy = pj->x[1] - shift[1];
        const double pjz = pj->x[2] - shift[2];

        /* Compute the pairwise distance. */
        double dx[3] = {pi->x[0] - pjx, pi->x[1] - pjy, pi->x[2] - pjz};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < ri * ri) {
          delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                  sort_pi.i, 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          pj->geometry.delaunay_flags |= 1 << sid;
        }
      }
    }
  } else {
    /* ci on the left */

    /* Loop over the sorted active parts in ci (on the left) */
    for (int i = 0; i < count_active; i++) {

      /* Get a hold of pi. */
      struct sort_entry sort_pi = sort_active_i[i];
      struct part *restrict pi = &parts_i[sort_pi.i];
      const float ri = pi->h;

      const double dj_max = sort_pi.d + ri + dx_max - rshift;

      /* Loop over the parts in cj (on the right) */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < dj_max; pjd++) {

        /* Recover pj */
        int pj_idx = sort_j[pjd].i;
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip particles that are already added */
        if (pj->geometry.delaunay_flags & 1 << sid) continue;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Shift pj so that it is in the frame of ci (with cj on the right) */
        const double pjx = pj->x[0] + shift[0];
        const double pjy = pj->x[1] + shift[1];
        const double pjz = pj->x[2] + shift[2];

        /* Compute the pairwise distance. */
        double dx[3] = {pi->x[0] - pjx, pi->x[1] - pjy, pi->x[2] - pjz};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < ri * ri) {
          delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                  sort_pi.i, 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          pj->geometry.delaunay_flags |= 1 << sid;
        }
      }
    }
  } /* Flipped? */
#endif
}

__attribute__((always_inline)) INLINE static void
runner_dopair_branch_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct cell *restrict cj) {

  TIMER_TIC;

  const struct engine *restrict e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  assert(ci->hydro.count != 0 && cj->hydro.count != 0);
#endif

  /* Is the cell active and local? */
  assert((cell_is_active_hydro(ci, e) && ci->nodeID == e->nodeID));

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  if (ci == cj) error("Interacting cell with itself!");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  struct cell *ci_temp = ci;
  struct cell *cj_temp = cj;
  int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
  const int flipped = ci != ci_temp;

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

  /* Make sure there is no voronoi allocated. */
  if (ci->grid.voronoi != NULL) {
    voronoi_destroy(ci->grid.voronoi);
    ci->grid.voronoi = NULL;
  }

  /* Delaunay already allocated? */
  if (ci->grid.delaunay == NULL) {
    ci->grid.delaunay = delaunay_malloc(ci->loc, ci->width, ci->hydro.count);
    ci->grid.ti_old = e->ti_current;
  } else {
    /* Check if reset is needed */
    if (ci->grid.ti_old < e->ti_current) {
      delaunay_reset(ci->grid.delaunay, ci->loc, ci->width, ci->hydro.count);
      ci->grid.ti_old = e->ti_current;
    }
  }

  /* We are good to go!*/

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Correct sid if needed */
  if (!flipped) sid = 26 - sid;

  /* Mark cell face as inside of simulation volume */
  ci->grid.delaunay->sid_is_inside_face_mask |= 1 << sid;

  /* Get some other useful values. */
  const float hi_max = ci->hydro.h_max_active;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_GRID
  /* Do a naive pair interaction */
  runner_dopair_grid_construction_naive(ci, cj, e, sort_i, sort_j, flipped,
                                        shift, rshift, hi_max, dx_max, sid);
#else
  /* Do a smart pair interaction. (we only construct voronoi cells of active
   * particles). */

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;

  int count_active = 0;
  struct sort_entry *restrict sort_active_i =
      (struct sort_entry *)malloc(count_i * sizeof(struct sort_entry));
  if (flipped) { /* ci on the right */

    const double di_max = sort_j[count_j - 1].d + hi_max + dx_max - rshift;

    /* Loop over parts in ci that could interact with parts in cj */
    for (int pid = 0; pid < count_i && sort_i[pid].d < di_max; pid++) {

      /* Found active particle? */
      if (part_is_active(&parts_i[sort_i[pid].i], e)) {
        sort_active_i[count_active] = sort_i[pid];
        count_active++;
      }
    }
  } else {
    const double di_min = sort_j[0].d - hi_max - dx_max + rshift;

    /* Loop over parts in ci that could interact with parts in cj */
    for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > di_min; pid--) {

      /* Found active particle? */
      if (part_is_active(&parts_i[sort_i[pid].i], e)) {
        sort_active_i[count_active] = sort_i[pid];
        count_active++;
      }
    }
  }

  runner_dopair_grid_construction(ci, cj, e, sort_i, sort_j, sort_active_i,
                                  count_active, flipped, shift, rshift, hi_max,
                                  dx_max, sid);

  /* Be clean */
  free(sort_active_i);
#endif

  TIMER_TOC(timer_dopair_grid_construction);
}

__attribute__((always_inline)) INLINE static void
runner_doself_grid_construction_naive(struct cell *restrict c) {

  /* Get useful variables */
  int count = c->hydro.count;
  const struct part *restrict parts = c->hydro.parts;

  /* Loop over the parts in c. */
  for (int i = 0; i < count; i++) {
#ifdef SHADOWSWIFT_HILBERT_ORDERING
    int pid = c->grid.hilbert_r_sort[i];
#else
    int pid = i;
#endif
    /* Get a pointer to the idx-th particle. */
    const struct part *restrict pi = &parts[pid];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Add all particles to the delaunay tesselation */
    delaunay_add_local_vertex(c->grid.delaunay, pid, pix, piy, piz, -1);
  }
}

__attribute__((always_inline)) INLINE static void
runner_doself_grid_construction(const struct engine *restrict e,
                                struct cell *restrict c) {

  /* Get useful variables */
  int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

#ifndef SHADOWSWIFT_BVH
  int *pid_active = malloc(c->hydro.count * sizeof(int));
  if (pid_active == NULL) error("Allocation of pid_active array failed.");
  int count_active = 0;
#endif
  int *pid_inactive = malloc(c->hydro.count * sizeof(int));
  if (pid_inactive == NULL) error("Allocation of pid_inactive array failed.");
  int count_inactive = 0;

  /* Loop over the parts in c. */
  for (int i = 0; i < count; i++) {
#ifdef SHADOWSWIFT_HILBERT_ORDERING
    int pid = c->grid.hilbert_r_sort[i];
#else
    int pid = i;
#endif
    /* Get a pointer to the idx-th particle. */
    struct part *restrict pi = &parts[pid];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* Add inactive particles only if they fall in the search radius of an
     * active particle */
    if (part_is_active(pi, e)) {
      delaunay_add_local_vertex(c->grid.delaunay, pid, pix, piy, piz, -1);
      /* Update delaunay flags to signal that the particle was added for
       * the self interaction */
      pi->geometry.delaunay_flags |= 1 << 13;

#ifndef SHADOWSWIFT_BVH
      /* Add it to the array of active particles. */
      pid_active[count_active] = pid;
      count_active++;
#endif
    } else {
      /* Add it to the array of inactive particles. */
      pid_inactive[count_inactive] = pid;
      count_inactive++;
    }
  }

  if (count_inactive == 0) {
    /* Be clean */
#ifndef SHADOWSWIFT_BVH
    free(pid_active);
#endif
    free(pid_inactive);
    /* We are done here */
    return;
  }

#ifdef SHADOWSWIFT_BVH
  /* Loop through the inactive particles and check if they are contained in
   * the search radius of an active particle using the bvh */
  for (int i = 0; i < count_inactive; i++) {
    /* Get a pointer to the i-th inactive particle. */
    const int pid = pid_inactive[i];
    struct part *restrict pi = &parts[pid];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    int pjd = bvh_hit(c->grid.bvh, parts, pix, piy, piz);
    /* Hit or miss? */
    if (pjd >= 0) {
      delaunay_add_local_vertex(c->grid.delaunay, pid, pix, piy, piz, pjd);
      /* Update delaunay flags to signal that the particle was added for
       * the self interaction */
      pi->geometry.delaunay_flags |= 1 << 13;
    }
  }
#else
  for (int i = 0; i < count_inactive; i++) {
    /* Get a pointer to the i-th inactive particle. */
    const int pid = pid_inactive[i];
    struct part *restrict pi = &parts[pid];
    const double pix = pi->x[0];
    const double piy = pi->x[1];
    const double piz = pi->x[2];

    /* pi is inactive, check if there is an active pj containing pi in its
     * search radius. */
    for (int j = 0; j < count_active; j++) {
      /* Get a pointer to the j-th active particle. */
      const int pjd = pid_active[j];
      const struct part *restrict pj = &parts[pjd];
      const double rj = pj->h;

      /* Compute pairwise distance */
      const double dx[3] = {pj->x[0] - pix, pj->x[1] - piy, pj->x[2] - piz};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < rj * rj) {
        delaunay_add_local_vertex(c->grid.delaunay, pid, pix, piy, piz, pjd);
        /* Update delaunay flags to signal that the particle was added for
         * the self interaction */
        pi->geometry.delaunay_flags |= 1 << 13;
        break;
      }
    }
  }
#endif

  /* Be clean */
#ifndef SHADOWSWIFT_BVH
  free(pid_active);
#endif
  free(pid_inactive);
}

__attribute__((always_inline)) INLINE static void
runner_doself_branch_grid_construction(struct runner *restrict r,
                                       struct cell *restrict c) {

  TIMER_TIC;

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Is the cell active and local? */
  assert((cell_is_active_hydro(c, e) && c->nodeID == e->nodeID));

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  /* Make sure there is no voronoi allocated. */
  if (c->grid.voronoi != NULL) {
    voronoi_destroy(c->grid.voronoi);
    c->grid.voronoi = NULL;
  }

  /* Delaunay already allocated? */
  if (c->grid.delaunay == NULL) {
    c->grid.delaunay = delaunay_malloc(c->loc, c->width, c->hydro.count);
    c->grid.ti_old = e->ti_current;
  } else {
    /* Check if rebuild is needed */
    if (c->grid.ti_old < e->ti_current) {
      delaunay_reset(c->grid.delaunay, c->loc, c->width, c->hydro.count);
      c->grid.ti_old = e->ti_current;
    }
  }

  /* We are good to go!*/

#ifdef SHADOWSWIFT_HILBERT_ORDERING
  const int count = c->hydro.count;

  /* Calculate hilbert keys + sort */
  unsigned long *hilbert_keys =
      (unsigned long *)malloc(count * sizeof(unsigned long));
  get_hilbert_keys(c, hilbert_keys);

  c->grid.hilbert_r_sort = (int *)malloc(count * sizeof(int));
  for (int i = 0; i < count; i++) {
    c->grid.hilbert_r_sort[i] = i;
  }
  qsort_r(c->grid.hilbert_r_sort, count, sizeof(int), sort_h_comp,
          hilbert_keys);
#endif

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_GRID
  /* Do a naive self interaction */
  runner_doself_grid_construction_naive(c);
#else
  /* Do a smart self interaction. (we only construct voronoi cells of active
   * particles). */
  runner_doself_grid_construction(e, c);
#endif

#ifdef SHADOWSWIFT_HILBERT_ORDERING
  /* Be clean */
  free(hilbert_keys);
  free(c->grid.hilbert_r_sort);
  c->grid.hilbert_r_sort = NULL;
#endif

  TIMER_TOC(timer_doself_grid_construction);
}

__attribute__((always_inline)) INLINE static void
runner_dopair_subset_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *restrict ind, float r_max,
                                       int count, struct cell *restrict cj) {
  const struct engine *restrict e = r->e;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  if (!cell_is_active_hydro(ci, e) && cj->nodeID == e->nodeID)
    error("Running construction task for inactive cell!");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  struct cell *ci_temp = ci;
  struct cell *cj_temp = cj;
  int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

  /* Get the cutoff shift. */
  double rshift = 0.0f;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Useful variables*/
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  const int flipped = ci != ci_temp;
  if (flipped) {
    /* ci on the right */

    /* Get the minimal position of any particle of ci along the sorting axis. */
    const double di_min = sort_i[0].d - dx_max - r_max + rshift;

    /* Loop over the neighbouring particles parts_j until they are definitely
     * too far to be a candidate ghost particle. */
    for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d > di_min; pjd--) {

      /* Get a pointer to the jth particle. */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip particles that were already added */
      if (pj->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the left) */
      const double pjx = pj->x[0] - shift[0];
      const double pjy = pj->x[1] - shift[1];
      const double pjz = pj->x[2] - shift[2];

      /* Loop over all the unconverged particles in parts_i and check if pj
       * falls within the new search radius, but outside the old search radius
       * of the unconverged particle pi. */
      for (int pid = 0; pid < count; pid++) {

        /* Get a hold of the ith part in ci. */
        struct part *restrict pi = &parts_i[ind[pid]];
        const double ri = pi->h;

#ifdef SWIFT_DEBUG_CHECKS
        if (!part_is_active(pi, e)) {
          error(
              "Encountered inactive unconverged particle in ghost "
              "construction "
              "task!");
        }
#endif

        /* Compute the pairwise distance. */
        const double dx[3] = {pi->x[0] - pjx, pi->x[1] - pjy, pi->x[2] - pjz};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < ri * ri) {
          delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                  ind[pid], 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          atomic_or(&pj->geometry.delaunay_flags, 1 << sid);
          break;
        }
      } /* Loop over unconverged particles in ci */
    }   /* Loop over particles in cj */
  } else {
    /* ci on the left */

    /* Correct sid if the cells have not been flipped */
    sid = 26 - sid;

    /* Get the maximal position of any particle of ci along the sorting axis. */
    const double di_max =
        sort_i[ci->hydro.count - 1].d + dx_max + r_max - rshift;

    /* Loop over the neighbouring particles parts_j until they are definitely
     * too far to be a candidate ghost particle. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di_max; pjd++) {

      /* Get a pointer to the jth particle. */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip particles that were already added */
      if (pj->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the right) */
      const double pjx = pj->x[0] + shift[0];
      const double pjy = pj->x[1] + shift[1];
      const double pjz = pj->x[2] + shift[2];

      /* Loop over all the unconverged particles in parts_i and check if pj
       * falls within the new search radius, but outside the old search radius
       * of the unconverged particle pi. */
      for (int pid = 0; pid < count; pid++) {

        /* Get a hold of the ith part in ci. */
        struct part *restrict pi = &parts_i[ind[pid]];
        const double ri = pi->h;

#ifdef SWIFT_DEBUG_CHECKS
        if (!part_is_active(pi, e)) {
          error(
              "Encountered inactive unconverged particle in ghost "
              "construction "
              "task!");
        }
#endif

        /* Compute the pairwise distance. */
        const double dx[3] = {pi->x[0] - pjx, pi->x[1] - pjy, pi->x[2] - pjz};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < ri * ri) {
          delaunay_add_new_vertex(ci->grid.delaunay, pjx, pjy, pjz, sid, pj_idx,
                                  ind[pid], 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          atomic_or(&pj->geometry.delaunay_flags, 1 << sid);
          break;
        }
      } /* Loop over unconverged particles in ci */
    }   /* Loop over particles in cj */
  }     /* Flipped? */
}

__attribute__((always_inline)) INLINE static void
runner_doself_subset_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *restrict ind, int count) {

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_GRID
  /* Nothing to do here, all particles already added */
  return;
#else
  struct engine *restrict e = r->e;

  /* Loop over all inactive particles in ci */
  for (int pjd = 0; pjd < ci->hydro.count; pjd++) {

    /* Retrieve particle */
    struct part *restrict pj = &parts_i[pjd];
    const double pjx = pj->x[0];
    const double pjy = pj->x[1];
    const double pjz = pj->x[2];

    /* Skip already added particles */
    if (pj->geometry.delaunay_flags & 1 << 13) continue;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pj, e)) continue;

    /* Loop over all unconverged particles */
    for (int i = 0; i < count; i++) {

      /* Retrieve particle */
      const int pid = ind[i];
      struct part *restrict pi = &parts_i[pid];
      const double ri = pi->h;

      /* Compute pairwise distance */
      const double dx[3] = {pi->x[0] - pjx, pi->x[1] - pjy, pi->x[2] - pjz};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < ri * ri) {
        delaunay_add_local_vertex(ci->grid.delaunay, pjd, pjx, pjy, pjz, pid);
        /* Update delaunay flags to signal that the particle was added for
         * the self interaction */
        atomic_or(&pj->geometry.delaunay_flags, 1 << 13);
        break;
      }
    }
  }
#endif
}

__attribute__((always_inline)) INLINE static void
grid_construction_get_cell_corner_to_compare(const struct cell *c, int sid,
                                             double *x_out) {
  const double cell_loc[3] = {c->loc[0], c->loc[1], c->loc[2]};
  const double cell_width[3] = {c->width[0], c->width[1], c->width[2]};

  for (int i = 0; i < 3; i++) {
    if (sortlist_shift_vector[sid][i] > 0) {
      x_out[i] = cell_loc[i] + cell_width[i];
    } else {
      x_out[i] = cell_loc[i];
    }
  }
}

__attribute__((always_inline)) INLINE static void
runner_add_boundary_particles_grid_construction(struct runner *restrict r,
                                                struct cell *restrict c,
                                                double r_max) {

  struct engine *e = r->e;

  /* Anything to do here? */
  if (e->s->periodic) return;
  if (!cell_is_active_hydro(c, e)) return;

  int count = c->hydro.count;
  struct part *parts = c->hydro.parts;
  for (int sid = 0; sid < 27; sid++) {
    /* Do we need to add periodic boundary particles for this cell? */
    if (c->grid.delaunay->sid_is_inside_face_mask & 1 << sid) continue;

    /* Pick the correct corner of the cell to compare the sorted positions
     * along the axis corresponding to the sid with. */
    double cell_corner[3];
    grid_construction_get_cell_corner_to_compare(c, sid, cell_corner);
    /* Calculate the position of the cell_corner along the sid axis */
    int sortlist_id = sortlistID[sid];
    double cell_corner_d = runner_shift[sortlist_id][0] * cell_corner[0] +
                           runner_shift[sortlist_id][1] * cell_corner[1] +
                           runner_shift[sortlist_id][2] * cell_corner[2];

    /* Get the sort entries of the particles for the sid axis */
    struct sort_entry *restrict sort = NULL;
    if (c->hydro.sort_allocated & (1 << sortlist_id))
      sort = cell_get_hydro_sorts(c, sortlist_id);

    /* Declare variables */
    struct part *p;
    double reflected_x[3];

    /* Do we need to flip the sorting direction? */
    if (runner_flip[sid]) {
      /* c on the right, reflection on the left */

      /* Loop over the sorted parts in c (on the right) */
      for (int i = 0; i < count; i++) {

        int p_idx;
        /* Can we use the sorts?
         * TODO: better to make sure that the sorts are activated...
         * Should be always the case for active cells unless the entire
         * simulation is one big cell... */
        if (sort != NULL) {
          if (sort[i].d > cell_corner_d + r_max) break;
          p_idx = sort[i].i;
        } else {
          p_idx = i;
        }

        /* Get a hold of the particle. */
        p = &parts[p_idx];

        /* Skip inactive particles */
        if (!part_is_active(p, e)) continue;

        /* Skip inhibited particles. */
        if (part_is_inhibited(p, e)) continue;

        /* Calculate reflected coordinates of particle */
        cell_reflect_coordinates(c, p->x, sid, reflected_x);

        /* Compute the pairwise distance. */
        double dx[3] = {p->x[0] - reflected_x[0], p->x[1] - reflected_x[1],
                        p->x[2] - reflected_x[2]};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        double radius = p->h;
        if (r2 < radius * radius) {
          delaunay_add_new_vertex(c->grid.delaunay, reflected_x[0],
                                  reflected_x[1], reflected_x[2], sid, p_idx,
                                  p_idx, 1);
        }
      }
    } else {
      /* c on the left, reflection on the right */
      for (int i = count - 1; i >= 0; i--) {

        int p_idx;
        /* Can we use the sorts? */
        if (sort != NULL) {
          if (sort[i].d < cell_corner_d - r_max) break;
          p_idx = sort[i].i;
        } else {
          p_idx = i;
        }

        /* Get a hold of the particle. */
        p = &parts[p_idx];

        /* Skip inactive particles */
        if (!part_is_active(p, e)) continue;

        /* Skip inhibited particles. */
        if (part_is_inhibited(p, e)) continue;

        /* Calculate reflected coordinates of particle */
        cell_reflect_coordinates(c, p->x, sid, reflected_x);

        /* Compute the pairwise distance. */
        double dx[3] = {p->x[0] - reflected_x[0], p->x[1] - reflected_x[1],
                        p->x[2] - reflected_x[2]};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        double radius = p->h;
        if (r2 < radius * radius) {
          /* We add boundary particles under sid 27 (for now) */
          delaunay_add_new_vertex(c->grid.delaunay, reflected_x[0],
                                  reflected_x[1], reflected_x[2], sid, p_idx,
                                  p_idx, 1);
        }
      }
    }
  }
}

__attribute__((always_inline)) INLINE static void
runner_add_boundary_particles_subset_grid_construction(
    struct runner *restrict r, struct cell *restrict c,
    struct part *restrict parts, const int *restrict ind, int count) {

  struct engine *e = r->e;

  /* Anything to do here? */
  if (e->s->periodic) return;
  if (c->grid.delaunay->sid_is_inside_face_mask ==
      0b111111111111111111111111111ul)
    return;

  /* Loop over unconverged parts*/
  for (int i = 0; i < count; i++) {
    /* Retrieve particle */
    int p_idx = ind[i];
    struct part *p = &parts[p_idx];

    /* Do we need to add a mirror of this particle as a boundary particle? */
    for (int sid = 0; sid < 27; sid++) {
      /* Inside face? */
      if (c->grid.delaunay->sid_is_inside_face_mask & 1 << sid) continue;

      /* Calculate reflected coordinates of particle */
      double reflected_x[3];
      cell_reflect_coordinates(c, p->x, sid, reflected_x);

      /* Compute the pairwise distance. */
      double dx[3] = {p->x[0] - reflected_x[0], p->x[1] - reflected_x[1],
                      p->x[2] - reflected_x[2]};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      double radius = p->h;
      if (r2 < radius * radius) {
        delaunay_add_new_vertex(c->grid.delaunay, reflected_x[0],
                                reflected_x[1], reflected_x[2], sid, p_idx,
                                p_idx, 1);
      }
    }
  }
}

__attribute__((always_inline)) INLINE static void runner_build_grid(
    struct runner *r, struct cell *c, int timer) {
  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Is the cell active and local? */
  if (!cell_is_active_hydro(c, e) || c->nodeID != e->nodeID)
    error("Running construction task for inactive cell!");

  /* Are we on the construction level? */
  if (c->grid.construction_level != on_construction_level)
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
  int *mask_active;
  if ((mask_active = (int *)malloc(c->hydro.count * sizeof(*mask_active))) ==
      NULL) {
    error("Can't allocate memory for mask_active");
  }
  int count = 0;
  float h_max = 0.f;
  for (int i = 0; i < c->hydro.count; i++)
    if (part_is_active(&parts[i], e)) {
      pid[count] = i;
      search_radii[count] = 0.;
      mask_active[i] = 1;
      h_max = fmaxf(h_max, parts[i].h);
      ++count;
    } else {
      mask_active[i] = 0;
    }

  /* First add all the active particles to the delaunay tesselation */
  cell_add_local_parts_grid(d, c, parts, pid, count);

  /* Now add ghost particles (i.e. particles from neighbouring cells and/or
   * inactive particles) until all active particles have converged */
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  float h_max_unconverged = h_max, h_max_active = 0.f;
  int redo, num_reruns;
  for (num_reruns = 0; count > 0 && num_reruns < max_smoothing_iter;
       num_reruns++) {

    /* Build bvh of unconverged particles */
    struct BVH bvh;
#ifdef SHADOWSWIFT_BVH
    bvh_populate(&bvh, parts, pid, count, count);
#endif

    /* Add ghost particles from this cell */
    cell_add_ghost_parts_grid_self(d, c, e, parts, &bvh, pid, count);

    /* Add ghost particles from neighbouring cells */
    for (struct link *l = c->grid.pair_sync_in; l != NULL; l = l->next) {
      struct cell *c_in = l->t->cj;

      /* Add ghost particles from cj */
      cell_add_ghost_parts_grid_pair(d, c, c_in, e, parts, &bvh, pid,
                                     h_max_unconverged, count);
    }

    if (!e->s->periodic) {
      /* Add boundary particles */
      cell_add_boundary_parts_grid(d, c, parts, pid, count);
    }

    /* Check if particles have converged */
    delaunay_get_search_radii(d, pid, count, search_radii);
    redo = 0;
    h_max_unconverged = 0.f;
    for (int i = 0; i < count; i++) {
      /* Get a direct pointer on the part. */
      struct part *p = &parts[pid[i]];

      float r_new = (float)search_radii[i];
      if (r_new >= p->h) {
        /* Un-converged particle */
        p->h *= 1.2f;
        h_max_unconverged = fmaxf(h_max_unconverged, p->h);
        pid[redo] = pid[i];
        redo += 1;
      } else {
        /* Particle has converged. Add a small buffer zone to compensate for
         * particle movement in the next iteration. */
        p->h = 1.1f * r_new;
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
  voronoi_build(c->grid.voronoi, d, parts, mask_active, c->hydro.count);

  /* Be clean */
  delaunay_destroy(d);
  free(pid);
  free(search_radii);
  free(mask_active);
}

#else
/* No grid construction needed */

__attribute__((always_inline)) INLINE static void
runner_dopair_branch_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct cell *restrict cj) {}

__attribute__((always_inline)) INLINE static void
runner_doself_branch_grid_construction(struct runner *restrict r,
                                       struct cell *restrict c) {}

__attribute__((always_inline)) INLINE static void
runner_dopair_subset_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *restrict ind,
                                       const double *restrict r_prev,
                                       double r_max, int count,
                                       struct cell *restrict cj) {}

__attribute__((always_inline)) INLINE static void
runner_doself_subset_grid_construction(struct runner *restrict r,
                                       struct cell *restrict ci,
                                       struct part *restrict parts_i,
                                       const int *restrict ind,
                                       const double *restrict r_prev,
                                       int count) {}
__attribute__((always_inline)) INLINE static void runner_build_bvh(
    struct runner *r, struct cell *c, int timer) {}

__attribute__((always_inline)) INLINE static void runner_build_grid(
    struct runner *r, struct cell *c, int timer) {}

#endif

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_H
