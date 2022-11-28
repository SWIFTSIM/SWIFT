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

/*! @brief Method to do the flux exchange over the faces between a pair of
 * cells.
 *
 * There are 2 modes this function can operate in. In the first mode (0), flux
 * exchange is carried out for all faces of active voronoi cells of ci. In
 * the second mode (1), flux exchange is carried out for all faces of active
 * voronoi cells of ci iff the neighbouring particle is inactive.
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
  if (ci->grid.construction_level == above_construction_level) {
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
  struct voronoi *vortess = ci->grid.voronoi;
  for (int i = 0; i < vortess->pair_index[sid]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[sid][i];

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
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
      /* Only do flux exchange (always symmetric) between normal particles */
      if (!part_do_apoptosis(part_left, e) &&
          !part_do_apoptosis(part_right, e)) {
        IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift,
             1);
      } else {
        /* One of the particles has to be killed */
        if (part_do_apoptosis(part_left, e)) {
          /* add portion of conserved quantities of part_left to part_right as
           * fluxes */
          runner_iact_apoptosis(part_left, part_right, pair->midpoint,
                                pair->surface_area);
        } else {
          /* add portion of conserved quantities of part_right to part_left as
           * fluxes */
          runner_iact_apoptosis(part_right, part_left, pair->midpoint,
                                pair->surface_area);
        }
      }
#else
      int ci_local = ci->nodeID == e->nodeID;
      int cj_local = cj->nodeID == e->nodeID;
      int left_active = part_is_active(part_left, e);
      int right_active = part_is_active(part_right, e);

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
  if (c->grid.construction_level == above_construction_level) {
    for (int i = 0; i < 8; i++) {
      struct cell *cp = c->progeny[i];
      if (cp != NULL && cell_is_active_hydro(cp, e)) DOPAIR_BOUNDARY(r, cp);
    }
    return;
  }

  /* Loop over boundary faces to apply hydro interaction */
  struct voronoi *vortess = c->grid.voronoi;
  double shift[3] = {0., 0., 0.};
  for (int idx = 0; idx < vortess->pair_index[27]; idx++) {

    /* Extract pair */
    struct voronoi_pair *pair = &vortess->pairs[27][idx];
    int sid = pair->sid;
    int part_idx = pair->left_idx;

    /* Extract particle */
    struct part *p = &c->hydro.parts[part_idx];

    /* Reconstruct boundary particle */
    struct part p_boundary = *p;
    cell_reflect_coordinates(c, p->x, sid, &p_boundary.x[0]);

    /* Make sure boundary faces do not move perpendicular to the boundary. */
    for (int i = 0; i < 3; i++)
      if (sortlist_shift_vector[sid][i] != 0) p_boundary.v_full[i] *= -1;

#if (SHADOWSWIFT_BC == VACUUM_BC)
    /* Set all primitive quantities and gradients of the vacuum particle to 0 */
    p_boundary.rho = 0.;
    p_boundary.v[0] = 0.;
    p_boundary.v[1] = 0.;
    p_boundary.v[2] = 0.;
    p_boundary.P = 0.;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    for (int i = 0; i < 3; i++) {
      p_boundary.gradients.rho[i] = 0.;
      p_boundary.gradients.v[0][i] = 0.;
      p_boundary.gradients.v[1][i] = 0.;
      p_boundary.gradients.v[2][i] = 0.;
      p_boundary.gradients.P[i] = 0.;
    }
#endif
#elif (SHADOWSWIFT_BC == OPEN_BC)
        /* We Treat inflow BC the same as open BC */
        /* Here we just flip the gradients for the flipped axis to ensure that
         * the extrapolated quantities on both sides of the face are the same
         * during the flux exchange */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    for (int i = 0; i < 3; i++) {
      if (sortlist_shift_vector[sid][i] != 0) {
        p_boundary.gradients.rho[i] *= -1;
        p_boundary.gradients.v[0][i] *= -1;
        p_boundary.gradients.v[1][i] *= -1;
        p_boundary.gradients.v[2][i] *= -1;
        p_boundary.gradients.P[i] *= -1;
      }
    }
#endif
#elif SHADOWSWIFT_BC == REFLECTIVE_BC
    /* Here we need to flip the fluid velocities and the gradients
     * perpendicular to the boundary (except the gradient of the flipped
     * velocity components). */
    for (int i = 0; i < 3; i++) {
      if (sortlist_shift_vector[sid][i] != 0) {
        /* Reflect the velocity along this axis */
        p_boundary.v[i] = -p_boundary.v[i];

#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
        p_boundary.gradients.rho[i] *= -1;
        p_boundary.gradients.v[(i + 1) % 3][i] *= -1;
        p_boundary.gradients.v[(i + 2) % 3][i] *= -1;
        p_boundary.gradients.P[i] *= -1;

        /* ... and the limiter info of the velocity */
        float temp = p_boundary.limiter.v[i][0];
        p_boundary.limiter.v[i][0] = -p_boundary.limiter.v[i][1];
        p_boundary.limiter.v[i][1] = -temp;
#endif
      }
    }
#elif (SHADOWSWIFT_BC == INFLOW_BC)
    /* Treat the boundary as open, except for the vertical boundaries, where we
     * set the correct inflow properties */
    struct hydro_space *hs = &r->e->s->hs;
    if (sid == 4 || sid == 22) {
      p_boundary.rho = hs->density;
      p_boundary.v[0] = sid == 4 ? hs->velocity : -hs->velocity;
      p_boundary.v[1] = 0.f;
      p_boundary.v[2] = 0.f;
      p_boundary.P = hs->pressure;
    }
    /* Set gradients to 0 for inflow boundary particles and flip them for the
     * rest (similar to open boundary conditions) */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    if (sid == 4 || sid == 22) {
      for (int i = 0; i < 3; i++) {
        p_boundary.gradients.rho[i] = 0.;
        p_boundary.gradients.v[0][i] = 0.;
        p_boundary.gradients.v[1][i] = 0.;
        p_boundary.gradients.v[2][i] = 0.;
        p_boundary.gradients.P[i] = 0.;
      }
    } else {
      for (int i = 0; i < 3; i++) {
        if (sortlist_shift_vector[sid][i] != 0) {
          p_boundary.gradients.rho[i] *= -1;
          p_boundary.gradients.v[0][i] *= -1;
          p_boundary.gradients.v[1][i] *= -1;
          p_boundary.gradients.v[2][i] *= -1;
          p_boundary.gradients.P[i] *= -1;
        }
      }
    }
#endif
#elif (SHADOWSWIFT_BC == RADIAL_INFLOW_BC)
    /* Set a radial inflow velocity everywhere and gradients to 0 */
    struct space *s = r->e->s;
    struct hydro_space *hs = &s->hs;
    double dx[3] = {
        0.5 * s->dim[0] - p_boundary.x[0],
        0.5 * s->dim[1] - p_boundary.x[1],
        0.5 * s->dim[2] - p_boundary.x[2],
    };
#ifdef HYDRO_DIMENSION_2D
    dx[2] = 0.;
#endif
    double norm = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    p_boundary.rho = hs->density;
    p_boundary.v[0] = hs->velocity * dx[0] / norm;
    p_boundary.v[1] = hs->velocity * dx[1] / norm;
    p_boundary.v[2] = hs->velocity * dx[2] / norm;
    p_boundary.P = hs->pressure;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    for (int i = 0; i < 3; i++) {
      p_boundary.gradients.rho[i] = 0.;
      p_boundary.gradients.v[0][i] = 0.;
      p_boundary.gradients.v[1][i] = 0.;
      p_boundary.gradients.v[2][i] = 0.;
      p_boundary.gradients.P[i] = 0.;
    }
#endif
#else
    error("Unknown boundary condition for non periodic run!");
#endif
    /* Now interact the real particle with the reconstructed boundary particle.
     */
    IACT(p, &p_boundary, pair->midpoint, pair->surface_area, shift, 0);
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

  /* Retrieve voronoi grid */
  struct voronoi *vortess = c->grid.voronoi;

  /* Loop over local pairs (sid 13) */
  for (int i = 0; i < vortess->pair_index[13]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[13][i];

    /* Retrieve particles */
    struct part *part_left = &c->hydro.parts[pair->left_idx];
    struct part *part_right = &c->hydro.parts[pair->right_idx];

#if (FUNCTION_TASK_LOOP == TASK_LOOP_FLUX_EXCHANGE)
    /* Only do flux exchange (always symmetric) between normal particles */
    if (!part_do_apoptosis(part_left, e) && !part_do_apoptosis(part_right, e)) {
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift, 1);
    } else {
      /* One of the particles has to be killed */
      if (part_do_apoptosis(part_left, e)) {
        /* add portion of conserved quantities of part_left to part_right as
         * fluxes */
        runner_iact_apoptosis(part_left, part_right, pair->midpoint,
                              pair->surface_area);
      } else {
        /* add portion of conserved quantities of part_right to part_left as
         * fluxes */
        runner_iact_apoptosis(part_right, part_left, pair->midpoint,
                              pair->surface_area);
      }
    }
#else
    const int left_is_active = part_is_active(part_left, e);
    const int right_is_active = part_is_active(part_right, e);

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
