//
// Created by yuyttenh on 29/03/22.
//

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "space_getsid.h"
#include "timers.h"

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H

__attribute__((always_inline)) INLINE static void
runner_dopair_grid_flux_exchange(struct runner *restrict r,
                                 struct cell *restrict ci,
                                 struct cell *restrict cj) {
  struct engine *e = r->e;

  /* anything to do here? */
  int ci_active = cell_is_active_hydro(ci, e);
  int cj_active = cell_is_active_hydro(cj, e);

  assert(ci_active);
#ifdef SWIFT_DEBUG_CHECKS
  assert(ci->grid.voronoi != NULL);
  assert(cj->grid.voronoi != NULL);
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

  struct voronoi *vortess = ci->grid.voronoi;

  double inverse_shift[3] = {-shift[0], -shift[1], -shift[2]};

  /* loop over voronoi faces between ci and cj */
  for (int i = 0; i < vortess->pair_index[sid]; ++i) {
    struct voronoi_pair *pair = &vortess->pairs[sid][i];

    /* Retrieve the particles */
    struct part *part_left = &ci->hydro.parts[pair->left_idx];
    struct part *part_right = &cj->hydro.parts[pair->right_idx];

    /* Anything to do here? */
    int left_active = part_is_active(part_left, e);
    int right_active = cj_active && part_is_active(part_right, e);
    if (left_active) {
      runner_iact_flux(part_left, part_right, pair->midpoint,
                       pair->surface_area, shift);
    } else if (right_active) {
      runner_iact_flux(part_left, part_right, pair->midpoint,
                       pair->surface_area, inverse_shift);
    }
  } /* loop over voronoi faces between ci and cj */
}

__attribute__((always_inline)) INLINE static void
runner_doself_grid_flux_exchange(struct runner *restrict r,
                                 struct cell *restrict c) {
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
      runner_iact_flux(part_left, part_right, pair->midpoint,
                       pair->surface_area, shift);
    } else if (part_is_active(part_right, e)) {
      runner_iact_flux(part_left, part_right, pair->midpoint,
                       pair->surface_area, shift);
    }
  }
}

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_HYDRO_H
