//
// Created by yuyttenh on 11/04/22.
//

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_grid_hydro.h"

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
          "Flux exchange for cell above construction level, indicating the "
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
    sid = 26 - sid;
    shift[0] = -shift[0];
    shift[1] = -shift[1];
    shift[2] = -shift[2];
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
     * call. */
    if (mode == 0 || !part_is_active(part_right, e)) {
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift);
    }
  } /* loop over voronoi faces between ci and cj */

  TIMER_TOC(TIMER_DOPAIR);
}

void DOPAIR_BRANCH(
    struct runner *restrict r, struct cell *ci, struct cell *cj) {
  struct engine *e = r->e;

  int ci_active = cell_is_active_hydro(ci, e);
  int cj_active = cell_is_active_hydro(cj, e);

  if (!ci_active) {
    if (!cj_active)
      error("Flux exchange activated between two inactive cells!");

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

void DOSELF(
    struct runner *restrict r, struct cell *restrict c) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  assert(c->grid.voronoi != NULL);
#endif

  struct engine *e = r->e;

  double shift[3] = {0., 0., 0.};

  /* Retrieve voronoi grid */
  struct voronoi *vortess = c->grid.voronoi;

  /* Loop over local pairs (sid 13) */
  for (int i = 0; i < vortess->pair_index[13]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[13][i];

    /* Retrieve particles */
    struct part *part_left = &c->hydro.parts[pair->left_idx];
    struct part *part_right = &c->hydro.parts[pair->right_idx];

    if (part_is_active(part_left, e)) {
      IACT(part_left, part_right, pair->midpoint, pair->surface_area, shift);
    } else if (part_is_active(part_right, e)) {
      IACT(part_right, part_left, pair->midpoint, pair->surface_area, shift);
    }
  }

  TIMER_TOC(TIMER_DOSELF);
}
