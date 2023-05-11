//
// Created by yuyttenh on 11/05/23.
//

#ifndef SWIFTSIM_SHADOWSWIFT_H
#define SWIFTSIM_SHADOWSWIFT_H

#ifdef SHADOWSWIFT_HILBERT_ORDERING
#include "hilbert.h"

/*! @brief Calculate the hilbert keys of the vertices
 *
 * @param c Cell containing the vertices
 */
__attribute__((always_inline)) INLINE static void get_hilbert_keys_active(
    const struct cell *c, const struct part *restrict parts, const int *pid,
    int count, unsigned long *restrict keys) {
  for (int i = 0; i < count; i++) {
    float dx_max = c->hydro.dx_max_part;
#if defined(HYDRO_DIMENSION_1D)
    /* In 1D, we just use the mantissa of the rescaled x coordinate to sort the
     * particles. */
    /* Get rescaled x coordinate within [1, 2) */
    double x = parts[pid[i]].x[0];
    double x_scaled =
        (x + dx_max - c->loc[0]) / (c->width[0] + 2. * dx_max) + 1.;
    keys[i] = delaunay_double_to_int(x_scaled);
#elif defined(HYDRO_DIMENSION_2D)
    unsigned long bits[2];
    int nbits = 32;
    double max_width = max(c->width[0], c->width[1]) + 2 * dx_max;
    bits[0] = (parts[pid[i]].x[0] + dx_max - c->loc[0]) / max_width *
              (1ul << (nbits - 1));
    bits[1] = (parts[pid[i]].x[1] + dx_max - c->loc[1]) / max_width *
              (1ul << (nbits - 1));
    keys[i] = hilbert_get_key_2d(bits, nbits);
#elif defined(HYDRO_DIMENSION_3D)
    unsigned long bits[3];
    int nbits = 21;
    double max_width = max3(c->width[0], c->width[1], c->width[2]) + 2 * dx_max;
    bits[0] = (parts[pid[i]].x[0] + dx_max - c->loc[0]) / max_width *
              (1ul << (nbits - 1));
    bits[1] = (parts[pid[i]].x[1] + dx_max - c->loc[1]) / max_width *
              (1ul << (nbits - 1));
    bits[2] = (parts[pid[i]].x[2] + dx_max - c->loc[2]) / max_width *
              (1ul << (nbits - 1));
    keys[i] = hilbert_get_key_3d(bits, nbits);
#endif
  }
}
#endif

__attribute__((always_inline)) INLINE static void cell_add_local_parts_grid(
    struct delaunay *d, struct cell *c, struct part *parts, const int *pid,
    int count) {
#ifdef SHADOWSWIFT_HILBERT_ORDERING
  /* Calculate hilbert keys of active particles + sort */
  unsigned long *hilbert_keys =
      (unsigned long *)malloc(count * sizeof(*hilbert_keys));
  get_hilbert_keys_active(c, parts, pid, count, hilbert_keys);

  int *hilbert_r_sort = (int *)malloc(count * sizeof(*hilbert_r_sort));
  for (int i = 0; i < count; i++) {
    hilbert_r_sort[i] = i;
  }
  qsort_r(hilbert_r_sort, count, sizeof(int), sort_h_comp, hilbert_keys);
#endif

  /* Add the active particles to the delaunay tesselation */
  /* Loop over the parts in c. */
  for (int i = 0; i < count; i++) {
#ifdef SHADOWSWIFT_HILBERT_ORDERING
    int p_idx = pid[hilbert_r_sort[i]];
#else
    int p_idx = pid[i];
#endif
    /* Get a pointer to the idx-th particle. */
    struct part *restrict p = &parts[p_idx];
    const double p_x = p->x[0];
    const double p_y = p->x[1];
    const double p_z = p->x[2];

    /* Add all particles to the delaunay tesselation */
    delaunay_add_local_vertex(d, p_idx, p_x, p_y, p_z, -1);
    /* Update delaunay flags to signal that the particle was added for
     * the self interaction */
    atomic_or(&p->geometry.delaunay_flags, 1 << 13);
  }

#ifdef SHADOWSWIFT_HILBERT_ORDERING
  /* Be clean */
  free(hilbert_keys);
  free(hilbert_r_sort);
#endif
}

__attribute__((always_inline)) INLINE static void
cell_add_ghost_parts_grid_self(struct delaunay *d, struct cell *restrict c,
                               const struct engine *e,
                               struct part *restrict parts,
                               const struct BVH *bvh, const int *restrict pid,
                               int count) {

  /* Loop over all inactive particles in ci */
  for (int i = 0; i < c->hydro.count; i++) {

    /* Retrieve particle */
    struct part *restrict p = &parts[i];

    /* Skip already added particles (this also skips active particles) */
    if (p->geometry.delaunay_flags & 1 << 13) continue;
    /* Skip inhibited particles. */
    if (part_is_inhibited(p, e)) continue;

    const double p_x = p->x[0];
    const double p_y = p->x[1];
    const double p_z = p->x[2];

#ifdef SHADOWSWIFT_BVH
    /* Find a bvh hit (if any) for this part */
    int ngb_id = bvh_hit(bvh, parts, p_x, p_y, p_z);
    if (ngb_id >= 0) {
      delaunay_add_local_vertex(d, i, p_x, p_y, p_z, ngb_id);
      /* Update delaunay flags to signal that the particle was added for
       * the self interaction */
      atomic_or(&p->geometry.delaunay_flags, 1 << 13);
    }
#else
    /* Loop over all unconverged particles to find a neighbour (if any) */
    for (int j = 0; j < count; j++) {

      /* Retrieve particle */
      const int ngb_id = pid[j];
      struct part *restrict ngb = &parts[ngb_id];
      const double r = ngb->h;

      /* Compute pairwise distance */
      const double dx[3] = {ngb->x[0] - p_x, ngb->x[1] - p_y, ngb->x[2] - p_z};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < r * r) {
        delaunay_add_local_vertex(d, i, p_x, p_y, p_z, ngb_id);
        /* Update delaunay flags to signal that the particle was added for
         * the self interaction */
        atomic_or(&ngb->geometry.delaunay_flags, 1 << 13);
        break;
      }
    }
#endif
  }
}

__attribute__((always_inline)) INLINE static void
cell_add_ghost_parts_grid_pair(struct delaunay *d, struct cell *c,
                               struct cell *c_in, const struct engine *e,
                               struct part *parts, const struct BVH *bvh,
                               const int *pid, float r_max, int count) {

  const int count_in = c_in->hydro.count;
  struct part *restrict parts_in = c_in->hydro.parts;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  struct cell *ci_temp = c;
  struct cell *cj_temp = c_in;
  int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
  const int flipped = c != ci_temp;

#ifdef SHADOWSWIFT_BVH
  /* Pick out the sorted lists of the incoming particles */
  const struct sort_entry *restrict sort_in = cell_get_hydro_sorts(c_in, sid);

  /* Get the cut-off distance from the bvh */
  double bvh_width[3];
  bvh_get_width(bvh, bvh_width);
  double cut_off = sort_get_cell_min_dist(sid, bvh->bbox.anchor, bvh_width);
  if (!flipped) {
    for (int k = 0; k < 3; k++)
      cut_off -=
          runner_shift[sid][k] * sortlist_shift_vector[sid][k] * bvh_width[k];
  }

  /* Get the cutoff shift. */
  double rshift = 0.0f;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];
  double dx_max = c_in->hydro.dx_max_sort;

  if (flipped) {
    /* c_in on the left */
    double d_min = cut_off + rshift;
    for (int i = count_in - 1; i >= 0 && sort_in[i].d + dx_max > d_min; i--) {
      /* Retrieve the i-th particle */
      const int p_idx = sort_in[i].i;
      struct part *p = &parts_in[p_idx];

      /* Skip particles that were already added */
      if (p->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(p, e)) continue;

      /* Shift p so that it is in the frame of c (with c_in on the left) */
      const double p_x = p->x[0] - shift[0];
      const double p_y = p->x[1] - shift[1];
      const double p_z = p->x[2] - shift[2];

      /* Find a bvh_hit if any. */
      int ngb_id = bvh_hit(bvh, parts, p_x, p_y, p_z);
      if (ngb_id >= 0) {
        delaunay_add_new_vertex(d, p_x, p_y, p_z, sid, p_idx, ngb_id, 0);
        /* Update delaunay flags to signal that the particle was added for
         * the self interaction */
        atomic_or(&p->geometry.delaunay_flags, 1 << sid);
      }
    }
  } else {
    /* c_in on the right */
    sid = 26 - sid;
    double d_max = cut_off - rshift;
    for (int i = 0; i < count_in && sort_in[i].d - dx_max < d_max; i++) {
      /* Retrieve the i-th particle */
      const int p_idx = sort_in[i].i;
      struct part *p = &parts_in[p_idx];

      /* Skip particles that were already added */
      if (p->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(p, e)) continue;

      /* Shift p so that it is in the frame of c (with c_in on the left) */
      const double p_x = p->x[0] + shift[0];
      const double p_y = p->x[1] + shift[1];
      const double p_z = p->x[2] + shift[2];

      /* Find a bvh_hit if any. */
      int ngb_id = bvh_hit(bvh, parts, p_x, p_y, p_z);
      if (ngb_id >= 0) {
        delaunay_add_new_vertex(d, p_x, p_y, p_z, sid, p_idx, ngb_id, 0);
        /* Update delaunay flags to signal that the particle was added for
         * the self interaction */
        atomic_or(&p->geometry.delaunay_flags, 1 << sid);
      }
    }
  }
#else
  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort = cell_get_hydro_sorts(c, sid);
  const struct sort_entry *restrict sort_in = cell_get_hydro_sorts(c_in, sid);

  /* Get the cutoff shift. */
  double rshift = 0.0f;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Useful variables*/
  const float dx_max = (c->hydro.dx_max_sort + c_in->hydro.dx_max_sort);

  if (flipped) {
    /* c on the right */

    /* Get the minimal position of any particle of ci along the sorting axis. */
    const double d_min = sort[0].d - dx_max - r_max + rshift;

    /* Loop over the neighbouring particles parts_j until they are definitely
     * too far to be a candidate ghost particle. */
    for (int i = count_in - 1; i >= 0 && sort_in[i].d > d_min; i--) {

      /* Get a pointer to the jth particle. */
      int p_in_idx = sort_in[i].i;
      struct part *restrict p_in = &parts_in[p_in_idx];

      /* Skip particles that were already added */
      if (p_in->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(p_in, e)) continue;

      /* Shift p_in so that it is in the frame of ci (with cj on the left) */
      const double p_in_x = p_in->x[0] - shift[0];
      const double p_in_y = p_in->x[1] - shift[1];
      const double p_in_z = p_in->x[2] - shift[2];

      /* Loop over all the unconverged particles in parts and check if p_in
       * falls within the new search radius, but outside the old search radius
       * of the unconverged particle p. */
      for (int j = 0; j < count; j++) {

        /* Get a hold of the ith part in c. */
        struct part *restrict p = &parts[pid[j]];
        const double r = p->h;

#ifdef SWIFT_DEBUG_CHECKS
        if (!part_is_active(p, e)) {
          error(
              "Encountered inactive unconverged particle in ghost "
              "construction "
              "task!");
        }
#endif

        /* Compute the pairwise distance. */
        const double dx[3] = {p->x[0] - p_in_x, p->x[1] - p_in_y,
                              p->x[2] - p_in_z};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < r * r) {
          delaunay_add_new_vertex(d, p_in_x, p_in_y, p_in_z, sid, p_in_idx,
                                  pid[j], 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          atomic_or(&p_in->geometry.delaunay_flags, 1 << sid);
          break;
        }
      } /* Loop over unconverged particles in ci */
    }   /* Loop over particles in cj */
  } else {
    /* ci on the left */

    /* Correct sid if the cells have not been flipped */
    sid = 26 - sid;

    /* Get the maximal position of any particle of c along the sorting axis. */
    const double d_max = sort[c->hydro.count - 1].d + dx_max + r_max - rshift;

    /* Loop over the neighbouring particles parts_in until they are definitely
     * too far to be a candidate ghost particle. */
    for (int i = 0; i < count_in && sort_in[i].d < d_max; i++) {

      /* Get a pointer to the jth particle. */
      int p_in_id = sort_in[i].i;
      struct part *restrict p_in = &parts_in[p_in_id];

      /* Skip particles that were already added */
      if (p_in->geometry.delaunay_flags & 1 << sid) continue;

      /* Skip inhibited particles. */
      if (part_is_inhibited(p_in, e)) continue;

      /* Shift pj so that it is in the frame of ci (with cj on the right) */
      const double p_in_x = p_in->x[0] + shift[0];
      const double p_in_y = p_in->x[1] + shift[1];
      const double p_in_z = p_in->x[2] + shift[2];

      /* Loop over all the unconverged particles in parts_i and check if pj
       * falls within the new search radius, but outside the old search radius
       * of the unconverged particle pi. */
      for (int j = 0; j < count; j++) {

        /* Get a hold of the ith part in ci. */
        struct part *restrict p = &parts[pid[j]];
        const double r = p->h;

#ifdef SWIFT_DEBUG_CHECKS
        if (!part_is_active(p, e)) {
          error(
              "Encountered inactive unconverged particle in ghost "
              "construction "
              "task!");
        }
#endif

        /* Compute the pairwise distance. */
        const double dx[3] = {p->x[0] - p_in_x, p->x[1] - p_in_y,
                              p->x[2] - p_in_z};
        const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < r * r) {
          delaunay_add_new_vertex(d, p_in_x, p_in_y, p_in_z, sid, p_in_id,
                                  pid[j], 0);
          /* Update delaunay flags to signal that the particle was added for
           * this sid */
          atomic_or(&p_in->geometry.delaunay_flags, 1 << sid);
          break;
        }
      } /* Loop over unconverged particles in ci */
    }   /* Loop over particles in cj */
  }     /* Flipped? */
#endif
}

__attribute__((always_inline)) INLINE static void cell_add_boundary_parts_grid(
    struct delaunay *d, struct cell *restrict c, struct part *restrict parts,
    const int *restrict pid, int count) {

  /* Anything to do here? */
  if (d->sid_is_inside_face_mask == 0b111111111111111111111111111ul) return;

  /* Loop over unconverged parts*/
  for (int i = 0; i < count; i++) {
    /* Retrieve particle */
    int p_idx = pid[i];
    struct part *p = &parts[p_idx];

    /* Do we need to add a mirror of this particle as a boundary particle? */
    for (int sid = 0; sid < 27; sid++) {
      /* Inside face? */
      if (d->sid_is_inside_face_mask & 1 << sid) continue;

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
        delaunay_add_new_vertex(d, reflected_x[0], reflected_x[1],
                                reflected_x[2], sid, p_idx, p_idx, 1);
      }
    }
  }
}

#endif  // SWIFTSIM_SHADOWSWIFT_H
