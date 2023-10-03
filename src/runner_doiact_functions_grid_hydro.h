//
// Created by yuyttenh on 11/04/22.
//

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_grid_hydro.h"
#include "swift.h"

#ifdef MOVING_MESH

/** @brief Do the interaction between an interior and exterior boundary particle
 * if necessary.
 *
 * This function applies to boundary particles inside the simulation (obstacles)
 * and not the global boundary conditions of the simulation. Faces between an
 * interior (only connected to other boundary particles) and exterior boundary
 * particle are treated as reflective. This function assumes that boundary
 * particles are either stationary or move together as a rigid body.
 *
 * @param part_left The left particle to consider
 * @param part_right The right particle to consider
 * @param left_is_active Whether to consider the left particle (we only treat
 * local active boundary particles).
 * @param right_is_active Whether to consider the right particle.
 * @param centroid The centroid of the face between both particles
 * @param surface_area The surface area of the face between both particles.
 * @param shift The shift to apply to the right particle to bring it into the
 * frame of the left
 * @return 1 When we found a pair of boundary particles, 0 otherwise. */
__attribute__((always_inline)) INLINE static int IACT_BOUNDARY_PARTICLES(
    struct part *part_left, struct part *part_right, int left_is_active,
    int right_is_active, const double *centroid, double surface_area,
    const double *shift) {
#ifdef SWIFT_BOUNDARY_PARTICLES
  /* Interaction between an interior and exterior boundary particle? */
  if (part_left->id < SWIFT_BOUNDARY_PARTICLES &&
      part_left->id >= space_boundary_parts_interior &&
      part_right->id < space_boundary_parts_interior) {
    /* Anything to do here? */
    if (!left_is_active) return 1;

    /* Create placeholder particle to apply boundary conditions on */
    struct part p_boundary = *part_right;
#if FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE
    /* Use an exact riemann solver to solve for reflective boundaries */
    runner_iact_boundary_reflective_flux_exchange(part_left, &p_boundary,
                                                  surface_area, centroid);
#else
    /* Normal functions for other interactions */
    runner_reflect_primitives(&p_boundary, part_left, centroid);
    IACT(part_left, &p_boundary, centroid, surface_area, shift, 0);
#endif
    /* Nothing left to do for this pair*/
    return 1;
  } else if (part_right->id < SWIFT_BOUNDARY_PARTICLES &&
             part_right->id >= space_boundary_parts_interior &&
             part_left->id < space_boundary_parts_interior) {
    /* Anything to do here? */
    if (!right_is_active) return 1;

    /* Create placeholder particle to apply boundary conditions on */
    struct part p_boundary = *part_left;
#if FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE
    runner_iact_boundary_reflective_flux_exchange(part_right, &p_boundary,
                                                  surface_area, centroid);
#else
    /* The pair needs to be flipped around */
    /* The midpoint from the reference frame of the right particle */
    double midpoint[3] = {centroid[0] - shift[0], centroid[1] - shift[1],
                          centroid[2] - shift[2]};
    /* The reversed shift */
    double r_shift[3] = {-shift[0], -shift[1], -shift[2]};
    runner_reflect_primitives(&p_boundary, part_right, centroid);
    IACT(part_right, &p_boundary, midpoint, surface_area, r_shift, 0);
#endif
    /* Nothing left to do for this pair*/
    return 1;
  } else if (part_left->id < space_boundary_parts_interior &&
             part_right->id < space_boundary_parts_interior) {
    /* No flux exchange between two interior boundary particles */
    return 1;
  }
#endif
  /* No interaction between boundary particles took place */
  return 0;
}

/*! @brief Method to do the flux exchange over the faces between a pair of
 * cells.
 *
 * There are 2 modes this function can operate in. In the first mode (0), the
 * interaction (e.g. flux exchange) is carried out for all faces of active
 * voronoi cells of ci. In the second mode (1), the interaction is carried out
 * for all faces of active voronoi cells of ci iff the neighbouring particle in
 * cj is inactive. That way we are sure that any faces between particles of ci
 * and cj are treated exactly once.
 *
 * @param r The current runner
 * @param ci The first cell of the pair. We use this cells voronoi tesseltation
 * @param cj The second cell of the pair.
 * @param mode Flag indicating in which mode we operate (0 or 1).
 * */
void DOPAIR(struct runner *restrict r, struct cell *ci, struct cell *cj,
            int mode) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (mode != 0 && mode != 1) error("Unknown mode for flux exchange: %d", mode);
#endif

  struct engine *e = r->e;

  /* Recurse? If the cells have been flipped in the branch function, ci might
   * be above its construction level. */
  if (ci->grid.construction_level == NULL) {
#ifdef SWIFT_DEBUG_CHECKS
    if (cell_is_active_hydro(cj, e) && mode == 0)
      error(
          "Hydro interaction for cell above construction level, indicating the "
          "pair has been flipped, but cj is active!");
#endif
    /* Retrieve SID and shift */
    struct cell *ci_old = ci;
    double shift[3];
    int sid = space_getsid(e->s, &ci, &cj, shift);
    int flipped = ci != ci_old;

    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      struct cell *ci_sub = ci->progeny[csp->pairs[k].pid];
      struct cell *cj_sub = cj->progeny[csp->pairs[k].pjd];
      if (ci_sub == NULL || cj_sub == NULL) continue;

      if (flipped) {
        if (cell_is_active_hydro(cj_sub, e)) DOPAIR(r, cj_sub, ci_sub, mode);
      } else {
        if (cell_is_active_hydro(ci_sub, e)) DOPAIR(r, ci_sub, cj_sub, mode);
      }
    }
    return;
  }

#ifdef SWIFT_DEBUG_CHECKS
  assert(cell_is_active_hydro(ci, e));
  assert(ci->grid.voronoi != NULL);
#endif

  /* Retrieve SID and shift */
  struct cell *ci_temp = ci;
  struct cell *cj_temp = cj;
  double shift[3];
  int sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);

  /* Correct shift and sid if flipped */
  if (ci != ci_temp) {
    shift[0] = -shift[0];
    shift[1] = -shift[1];
    shift[2] = -shift[2];
  } else {
    sid = 26 - sid;
  }

  /* loop over voronoi faces between ci and cj. */
  int face_count;
  const struct voronoi_pair *faces =
      voronoi_get_sid_faces(ci->grid.voronoi, sid, &face_count);
  for (int i = 0; i < face_count; ++i) {
    const struct voronoi_pair *pair = &faces[i];

    /* Retrieve the particles */
    struct part *part_left = &ci->hydro.parts[pair->left_idx];
    struct part *part_right = &cj->hydro.parts[pair->right_idx];

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(part_left, e))
      error("Encountered face of Voronoi cell of inactive particle!");
#endif

    /* Anything to do here? If the mode is 0, we always proceed, else we only
     * treat faces between an active particle of ci and an inactive particle of
     * cj. The other faces should already have been treated in another function
     * call with ci and cj (or their parents) reversed. */
    if (mode == 0 || !part_is_active(part_right, e)) {
      int ci_local = ci->nodeID == e->nodeID;
      int cj_local = cj->nodeID == e->nodeID;
      int left_active = part_is_active(part_left, e);
      int right_active = part_is_active(part_right, e);

      if (IACT_BOUNDARY_PARTICLES(
              part_left, part_right, left_active && ci_local,
              right_active && cj_local, pair->midpoint, pair->surface_area,
              shift)) {
        /* Face between boundary particle has been treated, nothing left to do*/
        continue;
      }

#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
      /* Flux exchange always symmetric */
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift, 1);
#else
      /* Only do the gradient calculations for local active particles */
      if (ci_local && left_active && cj_local && right_active) {
        IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift,
             1);
      } else if (ci_local && left_active) {
        IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift,
             0);
      } else if (cj_local && right_active) {
        /* The pair needs to be flipped around */
        /* The midpoint from the reference frame of the right particle */
        double midpoint[3] = {pair->midpoint[0] - shift[0],
                              pair->midpoint[1] - shift[1],
                              pair->midpoint[2] - shift[2]};
        /* The reversed shift */
        double r_shift[3] = {-shift[0], -shift[1], -shift[2]};
        IACT(part_right, part_left, midpoint, pair->surface_area, r_shift, 0);
      }
#endif
    }
  } /* loop over voronoi faces between ci and cj */

  TIMER_TOC(TIMER_DOPAIR);
}

void DOPAIR_BRANCH(struct runner *restrict r, struct cell *ci,
                   struct cell *cj) {

  /* Swap ci and cj if ci is above its construction level */
  if (ci->grid.construction_level == NULL) {
    struct cell *tmp = ci;
    ci = cj;
    cj = tmp;
  }
#ifdef SWIFT_DEBUG_CHECKS
  if (ci->grid.construction_level == NULL)
    error("ci should be at the construction level!");
#endif

  const struct engine *e = r->e;

  int ci_active = cell_is_active_hydro(ci, e);
  int cj_active = cell_is_active_hydro(cj, e);

  if (!ci_active) {
    if (!cj_active)
      error("Hydro interaction activated between two inactive cells!");

    /* Exchange flux from cj to ci only (in mode 0, but does not actually
     * matter in this case). */
    DOPAIR(r, cj, ci, 0);

  } else { /* ci_active */

    /* First do flux exchange between ci an cj in mode 0. */
    DOPAIR(r, ci, cj, 0);

    /* If cj is also active, we might have missed some faces between an active
     * particle of cj and an inactive particle of ci (not present in voronoi
     * tesselation of ci). */
    if (cj_active) {
      /* Also do flux exchange between cj and ci, but this time in mode 1, to
       * treat the remaining faces between an active particle of cj and inactive
       * particle of ci (those are not present in the voronoi tesselation of
       * ci. */
      DOPAIR(r, cj, ci, 1);
    }
  }
}

void DOPAIR_BOUNDARY(struct runner *restrict r, struct cell *restrict c) {

  const struct engine *e = r->e;

  /* Recurse? */
  if (c->grid.construction_level == NULL) {
    for (int i = 0; i < 8; i++) {
      struct cell *cp = c->progeny[i];
      if (cp != NULL && cell_is_active_hydro(cp, e)) DOPAIR_BOUNDARY(r, cp);
    }
    return;
  }

  /* Loop over boundary faces to apply hydro interaction */
  int face_count;
  const struct voronoi_pair *faces =
      voronoi_get_boundary_faces(c->grid.voronoi, &face_count);
  for (int idx = 0; idx < face_count; idx++) {

    /* Extract pair */
    const struct voronoi_pair *pair = &faces[idx];
    int sid = pair->sid;
    int part_idx = pair->left_idx;

    /* Extract particle */
    struct part *p = &c->hydro.parts[part_idx];

#ifdef SWIFT_BOUNDARY_PARTICLES
    if (p->id < space_boundary_parts_interior) {
      /* No flux exchange between interior boundary particles */
      continue;
    }
#endif

    /* Reconstruct boundary particle */
    struct part p_boundary = *p;
    cell_reflect_coordinates(c, p->x, sid, p_boundary.x);

    /* Make sure boundary faces do not move perpendicular to the boundary. */
    for (int i = 0; i < 3; i++)
      if (sortlist_shift_vector[sid][i] != 0) {
        p_boundary.v_full[i] *= -1.f;
        /* Also fix centroid of reflected particle */
        p_boundary.geometry.centroid[i] *= -1.f;
      }

    /* Interact with boundary particle */
    IACT_BOUNDARY(p, &p_boundary, pair->midpoint, pair->surface_area,
                  &r->e->s->hs);
  }
}

void DOSELF(struct runner *restrict r, struct cell *restrict c) {

  TIMER_TIC;

  struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  assert(c->grid.voronoi != NULL);
  if (c->nodeID != e->nodeID)
    error("Activated self hydro task for non-local cell!");
#endif

  double shift[3] = {0., 0., 0.};

  /* Loop over local pairs */
  int face_count;
  const struct voronoi_pair *faces =
      voronoi_get_local_faces(c->grid.voronoi, &face_count);
  for (int i = 0; i < face_count; ++i) {
    const struct voronoi_pair *pair = &faces[i];

    /* Retrieve particles */
    struct part *part_left = &c->hydro.parts[pair->left_idx];
    struct part *part_right = &c->hydro.parts[pair->right_idx];

    const int left_is_active = part_is_active(part_left, e);
    const int right_is_active = part_is_active(part_right, e);

    if (IACT_BOUNDARY_PARTICLES(part_left, part_right, left_is_active,
                                        right_is_active, pair->midpoint,
                                        pair->surface_area, shift)) {
      /* Face between boundary particle has been treated, nothing left to do */
      continue;
    }

#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    /* Flux exchange always symmetric */
    IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift, 1);
#else
    if (left_is_active && right_is_active) {
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift, 1);
    } else if (left_is_active) {
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift, 0);
    } else if (right_is_active) {
      IACT(part_right, part_left, pair->midpoint, pair->surface_area, shift, 0);
    }
#endif
  }

  TIMER_TOC(TIMER_DOSELF);
}

#else

void DOPAIR_BRANCH(struct runner *restrict r, struct cell *ci,
                   struct cell *cj) {}

void DOPAIR_BOUNDARY(struct runner *restrict r, struct cell *restrict c) {}

void DOSELF(struct runner *restrict r, struct cell *restrict c) {}

#endif
