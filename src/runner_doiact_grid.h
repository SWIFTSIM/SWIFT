//
// Created by yuyttenh on 23/03/22.
//

#ifndef SWIFTSIM_RUNNER_DOIACT_GRID_H
#define SWIFTSIM_RUNNER_DOIACT_GRID_H

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "shadowswift/delaunay.h"
#include "space_getsid.h"
#include "timers.h"

__attribute__((always_inline)) INLINE static void runner_dopair_grid_construction(struct runner *restrict r,
                                     struct cell *restrict ci,
                                     struct cell *restrict cj) {
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

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

  /* Delaunay already allocated? */
  if (ci->grid.delaunay == NULL) {
    ci->grid.delaunay = delaunay_malloc(ci->loc, ci->width, ci->hydro.count);
  }

  /* We are good to go!*/

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

  /* Correct sid if the cells have been flipped */
  if (ci != ci_temp) sid = 26 - sid;

  /* Get some other useful values. */
  const double hi_max = ci->hydro.h_max - rshift;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Mark cell face as inside of simulation volume */
  ci->grid.delaunay->sid_is_inside_face[sid] |= 1;

  /* Loop over the parts in ci. */
  for (int pid = count_i - 1;
       pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[sort_i[pid].i];
    const float hi = pi->h;

    /* Skip inactive particles */
    if (!part_is_active(pi, e)) {
      /* TODO what should we do here?
       * For the moment, also build grid for inactive particles... */
      continue;
    }

    /* Is there anything we need to interact with ? */
    //      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max -
    //      rshift;
    const double di = sort_i[pid].d + hi + dx_max - rshift;
    if (di < dj_min) continue;

    /* Get some additional information about pi */
    //      const float hig2 = hi * hi * kernel_gamma2;
    const float hig2 = hi * hi;
    const double pix = pi->x[0] - (cj->loc[0] + shift[0]);
    const double piy = pi->x[1] - (cj->loc[1] + shift[1]);
    const double piz = pi->x[2] - (cj->loc[2] + shift[2]);

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

      /* Recover pj */
      int pj_idx = sort_j[pjd].i;
      struct part *pj = &parts_j[pj_idx];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0] - cj->loc[0];
      const double pjy = pj->x[1] - cj->loc[1];
      const double pjz = pj->x[2] - cj->loc[2];

      /* Compute the pairwise distance. */
      double dx[3] = {pix - pjx, piy - pjy, piz - pjz};
      const double r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hig2) {
        delaunay_add_new_vertex(ci->grid.delaunay, pj->x[0] + shift[0],
                                pj->x[1] + shift[1], pj->x[2] + shift[2], sid,
                                pj_idx);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

__attribute__((always_inline)) INLINE static void runner_doself_grid_construction(struct runner *restrict r,
                                     struct cell *restrict c) {
  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Is the cell active and local? */
  assert((cell_is_active_hydro(c, e) && c->nodeID == e->nodeID));

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  /* Delaunay already allocated? */
  if (c->grid.delaunay == NULL) {
    c->grid.delaunay = delaunay_malloc(c->loc, c->width, c->hydro.count);
  }

  /* We are good to go!*/

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

#ifdef SHADOWFAX_HILBERT_ORDERING
  /* Update hilbert keys + sort */
  cell_update_hilbert_keys(c);
  for (int i = 0; i < count; i++) {
    c->hydro.hilbert_r_sort[i] = i;
  }
  qsort_r(c->hydro.hilbert_r_sort, count, sizeof(int), sort_h_comp,
          c->hydro.hilbert_keys);
#endif

  /* Loop over the parts in c. */
  for (int i = 0; i < count; i++) {
#ifdef SHADOWFAX_HILBERT_ORDERING
    int idx = c->hydro.hilbert_r_sort[i];
#else
    int idx = i;
#endif
    /* Get a pointer to the idx-th particle. */
    struct part *restrict p = &parts[idx];
    /* TODO skip inactive particles */
    delaunay_add_local_vertex(c->grid.delaunay, idx, p->x[0], p->x[1], p->x[2]);
  }
}

__attribute__((always_inline)) INLINE static void runner_dopair_subset_grid_construction(struct runner *restrict r,
                                            struct cell *restrict ci,
                                            struct part *restrict parts_i,
                                            int *restrict ind,
                                            double *restrict h_prev, int count,
                                            struct cell *restrict cj) {
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
  const struct sort_entry *sort_j = cell_get_hydro_sorts(cj, sid);
  const float dxj = cj->hydro.dx_max_sort;


  /* Loop over the parts_i. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[ind[pid]];
    const double pix = pi->x[0] - (shift[0]);
    const double piy = pi->x[1] - (shift[1]);
    const double piz = pi->x[2] - (shift[2]);
    const float hi = pi->h;
    const float hig2 = hi * hi;
    const double hi_prev = h_prev[ind[pid]];
    const double hi_prev2 = hi_prev * hi_prev;
    const double di = hi + dxj + pix * runner_shift[sid][0] +
                      piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

      /* Get a pointer to the jth particle. */
      int pj_idx = sort_j[pjd].i;
      struct part *restrict pj = &parts_j[pj_idx];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                     (float)(piz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Hit or miss? */
      if (r2 < hig2 && r2 >= hi_prev2) {
        delaunay_add_new_vertex(ci->grid.delaunay, pj->x[0] + shift[0],
                                pj->x[1] + shift[1], pj->x[2] + shift[2],
                                26 - sid, pj_idx);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

__attribute__((always_inline)) INLINE static void runner_doself_subset_grid_construction(struct runner *restrict r,
                                            struct cell *restrict ci,
                                            struct part *restrict parts_i,
                                            int *restrict ind,
                                            double *restrict h_prev,
                                            int count) {
  /* TODO */
}

#endif  // SWIFTSIM_RUNNER_DOIACT_GRID_H
