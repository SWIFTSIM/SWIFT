/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_limiter.h"

/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj, const int limit_min_h,
                   const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(ci, e) && !cell_is_starting_hydro(cj, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->dmin != cj->dmin) error("Cells of different size!");
#endif

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_starting(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      const char depth_j = pj->depth_h;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = part_is_starting(pj, e);

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      const int doi = pi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = pj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

      /* Hit or miss? */
      if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
      }
      if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
      }

    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1_NAIVE(struct runner *r, const struct cell *c,
                   const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count = c->hydro.count;
  struct part *parts = c->hydro.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_starting(pi, e);
    const char depth_i = pi->depth_h;
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = part_is_starting(pj, e);
      const char depth_j = pj->depth_h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = pi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = pj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        IACT(r2, dx, hi, hj, pi, pj, a, H);
      } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
      } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions between a cell pair (non-symmetric).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1(struct runner *r, const struct cell *restrict ci,
             const struct cell *restrict cj, const int limit_min_h,
             const int limit_max_h, const int sid, const double shift[3]) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;

  TIMER_TIC;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
  /* Some constants used to checks that the parts are in the right frame */
  const float shift_threshold_x =
      2. * ci->width[0] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2. * ci->width[1] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2. * ci->width[2] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

  /* Get the limits in h (if any). Note ci and cj are the same size */
#ifdef SWIFT_DEBUG_CHECKS
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
#endif
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;

  /* Get some other useful values. */
  const double hi_max =
      min(h_max, ci->hydro.h_max_active) * kernel_gamma - rshift;
  const double hj_max = min(h_max, cj->hydro.h_max_active) * kernel_gamma;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  if (cell_is_starting_hydro(ci, e)) {

    /* Loop over the *active* parts in ci that are within range (on the axis)
       of any particle in cj. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      const char depth_i = pi->depth_h;
      const float hi = pi->h;

      /* Skip inactive particles */
      if (!part_is_starting(pi, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hi > ci->hydro.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_i < min_depth) continue;
      if (depth_i > max_depth) continue;

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        }
      } /* loop over the parts in cj. */
    } /* loop over the parts in ci. */
  } /* Cell ci is active */

  if (cell_is_starting_hydro(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *pj = &parts_j[sort_j[pjd].i];
      const char depth_j = pj->depth_h;
      const float hj = pj->h;

      /* Skip inactive particles */
      if (!part_is_starting(pj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hj > cj->hydro.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_j < min_depth) continue;
      if (depth_j > max_depth) continue;

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = pj->x[0] - cj->loc[0];
      const float pjy = pj->x[1] - cj->loc[1];
      const float pjz = pj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        }
      } /* loop over the parts in ci. */
    } /* loop over the parts in cj. */
  } /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(ci, e) && !cell_is_starting_hydro(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID.
   * Note: this may swap the ci and cj pointers!! */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells (ci).");

  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells (cj).");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR1_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
#else
  DOPAIR1(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#endif
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1(struct runner *r, const struct cell *c, const int limit_min_h,
             const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  struct part *parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Set up a list of the particles for which we want to compute interactions */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++) {
    const struct part *p = &parts[k];
    const char depth = p->depth_h;
    if (part_is_starting(&parts[k], e) && (depth >= min_depth) &&
        (depth <= max_depth)) {
      indt[countdt] = k;
      countdt += 1;
    }
  }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over *all* the particles (i.e. the ones to update and not to update).
   *
   * Note the additional condition to make the loop abort if all the active
   * particles have been processed. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];
    const char depth_i = pi->depth_h;

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    /* Get the particle position and (square of) search radius. */
    const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle active and in the range of h we care about? */
    const int update_i = part_is_starting(pi, e) && (depth_i >= min_depth) &&
                         (depth_i <= max_depth);

    /* If false then it can only act as a neighbour of others */
    if (!update_i) {

      /* Loop over the particles we want to update. */
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. (by construction pi != pj) */
        struct part *restrict pj = &parts[indt[pjd]];

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        const float dx[3] = {(float)(pjx[0] - pix[0]), (float)(pjx[1] - pix[1]),
                             (float)(pjx[2] - pix[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
        }
      } /* loop over all the particles we want to update. */
    }

    /* Otherwise, interact with all candidates. */
    else {

      /* We caught a live one!
       * Move the start of the list of active ones by one slot as it will have
       * been fully processed after the following loop so no need to consider it
       * in the previous loop any more. */
      firstdt += 1;

      /* Loop over *all* the particles (i.e. the ones to update and not to
       * update) but starting from where we are in the overall list. */
      for (int pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle (by construction pi != pj). */
        struct part *restrict pj = &parts[pjd];
        const char depth_j = pj->depth_h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* This particle's (square of) search radius. */
        const float hj = pj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double pjx[3] = {pj->x[0], pj->x[1], pj->x[2]};
        float dx[3] = {(float)(pix[0] - pjx[0]), (float)(pix[1] - pjx[1]),
                       (float)(pix[2] - pjx[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Decide which of the two particles to update */

        /* We know pi is active and in the right range of h
         * -> Only check the distance to pj */
        const int doi = (r2 < hig2);

        /* We know nothing about pj
         * -> Check whether it is active
         * -> Check whether it is in the right range of h
         * -> Check the distance to pi */
        const int doj = (part_is_starting(pj, e)) && (depth_j >= min_depth) &&
                        (depth_j <= max_depth) && (r2 < hjg2);

        /* Hit or miss? */
        if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif
          /* Update both pi and pj */

          IACT(r2, dx, hi, hj, pi, pj, a, H);
        } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update only pi */

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
        } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update only pj */

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);

        } /* Hit or miss */
      } /* loop over all other particles. */
    } /* pi is active */
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Determine which version of DOSELF1 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH(struct runner *r, const struct cell *c,
                    const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    if (!limit_max_h && c->hydro.h_max_active * kernel_gamma > c->dmin)
      error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");
#endif

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_NAIVE(r, c, limit_min_h, limit_max_h);
#else
  DOSELF1(r, c, limit_min_h, limit_max_h);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int recurse_below_h_max, const int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Should we even bother? */
  /* const int do_i = cell_get_flag(ci, cell_flag_do_hydro_limiter); */
  /* const int do_j = cell_get_flag(cj, cell_flag_do_hydro_limiter); */
  /* const int do_sub_i = cell_get_flag(ci, cell_flag_do_hydro_sub_limiter); */
  /* const int do_sub_j = cell_get_flag(cj, cell_flag_do_hydro_sub_limiter); */

  /* if (!do_i && !do_j && !do_sub_i && !do_sub_j) return; */
  if (!cell_is_starting_hydro(ci, e) && !cell_is_starting_hydro(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->hydro.count < space_recurse_size_pair_hydro ||
      !cj->split || cj->hydro.count < space_recurse_size_pair_hydro) {

    /* Do any of the cells need to be sorted first?
     * Since h_max might have changed, we may not have sorted at this level */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
      /* Bert: RT probably broken here! */
      runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                           /*rt_request=*/0, /*clock=*/0);
    }

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH(r, ci, cj, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Both ci and cj are split */

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max && (!cell_can_recurse_in_subpair_hydro_task(ci) ||
                                 !cell_can_recurse_in_subpair_hydro_task(cj))) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Do any of the cells need to be sorted first?
       * Since h_max might have changed, we may not have sorted at this level */
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }
      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
        /* Bert: RT probably broken here! */
        runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*rt_request=*/0, /*clock=*/0);
      }

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR1_BRANCH(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1(r, ci->progeny[pid], cj->progeny[pjd], recurse_below_h_max,
                    /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1(struct runner *r, struct cell *c, int recurse_below_h_max,
                 const int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  /* const int do_i = cell_get_flag(c, cell_flag_do_hydro_limiter); */
  /* const int do_sub_i = cell_get_flag(c, cell_flag_do_hydro_sub_limiter); */

  /* if (!do_i && !do_sub_i) return; */
  if (!cell_is_starting_hydro(c, r->e)) return;
  if (c->hydro.count == 0) return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->hydro.count < space_recurse_size_self_hydro) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH(r, c, /*limit_h_min=*/0,
                   /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself_hydro_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1(r, c->progeny[k], recurse_below_h_max, /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1(r, c->progeny[k], c->progeny[j], recurse_below_h_max,
                        /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}
