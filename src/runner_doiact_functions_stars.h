/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

#include "runner_doiact_stars.h"

/**
 * @brief Calculate the number density of #part around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_STARS(struct runner *r, struct cell *c, const int limit_min_h,
                   const int limit_max_h) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int count = c->hydro.count;
  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  struct xpart *restrict xparts = c->hydro.xparts;
#endif

  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->dmin * 0.5 * (1. / kernel_gamma) : 0.;
  const float h_max = limit_max_h ? c->dmin * (1. / kernel_gamma) : FLT_MAX;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts[sid];
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Skip inactive particles */
    if (!feedback_is_active(si, e->time, cosmo, with_cosmology)) continue;

    /* Skip particles not in the range of h we care about */
    if (hi >= h_max) continue;
    if (hi < h_min) continue;

    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
      struct xpart *restrict xpj = &xparts[pjd];
#endif
      const float hj = pj->h;

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, si, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, si, pj, NULL, cosmo,
                                            e->feedback_props, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                          e->feedback_props, ti_current);
#endif
      }
    } /* loop over the parts in ci. */
  }   /* loop over the sparts in ci. */

  TIMER_TOC(TIMER_DOSELF_STARS);
}

/**
 * @brief Calculate the number density of cj #part around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#else
  if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->hydro.count == 0 || ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  struct xpart *restrict xparts_j = cj->hydro.xparts;
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts_i[sid];

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Skip inactive particles */
    if (!feedback_is_active(si, e->time, cosmo, with_cosmology)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
      struct xpart *restrict xpj = &xparts_j[pjd];
#endif
      const float hj = pj->h;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, si, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, si, pj, NULL, cosmo,
                                            e->feedback_props, ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, si, pj, xpj, cosmo,
                                          e->feedback_props, ti_current);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DO_SYM_PAIR1_STARS(struct runner *r, struct cell *restrict ci,
                        struct cell *restrict cj, const int limit_min_h,
                        const int limit_max_h, const int sid,
                        const double shift[3]) {

  TIMER_TIC;

  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_stars = (ci->nodeID == e->nodeID) && (ci->stars.count != 0) &&
                          (cj->hydro.count != 0) && cell_is_active_stars(ci, e);
  const int do_cj_stars = (cj->nodeID == e->nodeID) && (cj->stars.count != 0) &&
                          (ci->hydro.count != 0) && cell_is_active_stars(cj, e);
#else
  /* here we are updating the hydro -> switch ci, cj for local */
  const int do_ci_stars = (cj->nodeID == e->nodeID) && (ci->stars.count != 0) &&
                          (cj->hydro.count != 0) && cell_is_active_stars(ci, e);
  const int do_cj_stars = (ci->nodeID == e->nodeID) && (cj->stars.count != 0) &&
                          (ci->hydro.count != 0) && cell_is_active_stars(cj, e);
#endif

  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->dmin * 0.5 * (1. / kernel_gamma) : 0.;
  const float h_max = limit_max_h ? ci->dmin * (1. / kernel_gamma) : FLT_MAX;

  if (do_ci_stars) {

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
    const struct sort_entry *restrict sort_i = cell_get_stars_sorts(ci, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->stars.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hi_max = min(h_max, ci->stars.h_max) * kernel_gamma - rshift;
    const int count_i = ci->stars.count;
    const int count_j = cj->hydro.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct part *restrict parts_j = cj->hydro.parts;
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
    struct xpart *restrict xparts_j = cj->hydro.xparts;
#endif
    const double dj_min = sort_j[0].d;
    const float dx_max = (ci->stars.dx_max_sort + cj->hydro.dx_max_sort);
    const float hydro_dx_max_rshift = cj->hydro.dx_max_sort - rshift;

    /* Loop over the sparts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct spart *restrict spi = &sparts_i[sort_i[pid].i];
      const float hi = spi->h;

      /* Skip inactive particles */
      if (!spart_is_active(spi, e)) continue;

      /* Skip inactive particles */
      if (!feedback_is_active(spi, e->time, cosmo, with_cosmology)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hi > ci->stars.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (hi >= h_max) continue;
      if (hi < h_min) continue;

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].i + hi * kernel_gamma + hydro_dx_max_rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = spi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = spi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = spi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];
#if (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        struct xpart *xpj = &xparts_j[sort_j[pjd].i];
#endif

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
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
        if (spi->ti_drift != e->ti_current)
          error("Particle spi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, NULL,
                                              cosmo, e->feedback_props,
                                              ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj, cosmo,
                                            e->feedback_props, ti_current);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* do_ci_stars */

  if (do_cj_stars) {
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
    const struct sort_entry *restrict sort_j = cell_get_stars_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
    const float shift_threshold_x =
        2. * ci->width[0] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_y =
        2. * ci->width[1] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
    const float shift_threshold_z =
        2. * ci->width[2] +
        2. * max(ci->hydro.dx_max_part, cj->stars.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hj_max = min(h_max, cj->hydro.h_max) * kernel_gamma;
    const int count_i = ci->hydro.count;
    const int count_j = cj->stars.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct xpart *restrict xparts_i = ci->hydro.xparts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float dx_max = (ci->hydro.dx_max_sort + cj->stars.dx_max_sort);
    const float hydro_dx_max_rshift = ci->hydro.dx_max_sort - rshift;

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct spart *spj = &sparts_j[sort_j[pjd].i];
      const float hj = spj->h;

      /* Skip inactive particles */
      if (!spart_is_active(spj, e)) continue;

      /* Skip inactive particles */
      if (!feedback_is_active(spj, e->time, cosmo, with_cosmology)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hj > cj->stars.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (hj >= h_max) continue;
      if (hj < h_min) continue;

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - hydro_dx_max_rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = spj->x[0] - cj->loc[0];
      const float pjy = spj->x[1] - cj->loc[1];
      const float pjz = spj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];
        struct xpart *xpi = &xparts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
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
        if (spj->ti_drift != e->ti_current)
          error("Particle spj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

          IACT_STARS(r2, dx, hj, hi, spj, pi, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hj, hi, spj, pi, xpi,
                                              cosmo, e->feedback_props,
                                              ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          runner_iact_nonsym_feedback_apply(r2, dx, hj, hi, spj, pi, xpi, cosmo,
                                            e->feedback_props, ti_current);
#endif
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR_STARS);
}

void DOPAIR1_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {

  TIMER_TIC;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_stars = ci->nodeID == r->e->nodeID;
  const int do_cj_stars = cj->nodeID == r->e->nodeID;
#else
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_stars = cj->nodeID == r->e->nodeID;
  const int do_cj_stars = ci->nodeID == r->e->nodeID;
#endif
  if (do_ci_stars && ci->stars.count != 0 && cj->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS_NAIVE(r, ci, cj);
  if (do_cj_stars && cj->stars.count != 0 && ci->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS_NAIVE(r, cj, ci);

  TIMER_TOC(TIMER_DOPAIR_STARS);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param flipped Flag to check whether the cells have been flipped or not.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_STARS(struct runner *r, struct cell *restrict ci,
                          struct spart *restrict sparts_i, int *restrict ind,
                          int scount, struct cell *restrict cj, const int sid,
                          const int flipped, const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
  const float dxj = cj->hydro.dx_max_sort;

  /* Sparts are on the left? */
  if (!flipped) {

    /* Loop over the sparts_i. */
    for (int pid = 0; pid < scount; pid++) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[ind[pid]];
      const double pix = spi->x[0] - (shift[0]);
      const double piy = spi->x[1] - (shift[1]);
      const double piz = spi->x[2] - (shift[2]);
      const float hi = spi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];
        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, NULL,
                                              cosmo, e->feedback_props,
                                              e->ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          error("No subset feedback iact functions do (or should) exist!");
          /* runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj,
           * cosmo, ti_current); */
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the sparts in ci. */
  }

  /* Sparts are on the right. */
  else {

    /* Loop over the sparts_i. */
    for (int pid = 0; pid < scount; pid++) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[ind[pid]];
      const double pix = spi->x[0] - (shift[0]);
      const double piy = spi->x[1] - (shift[1]);
      const double piz = spi->x[2] - (shift[2]);
      const float hi = spi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];
        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (spi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, NULL,
                                              cosmo, e->feedback_props,
                                              e->ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
          error("No subset feedback iact functions do (or should) exist!");
          /* runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj,
           * cosmo, ti_current); */
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the sparts in ci. */
  }
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_STARS_NAIVE(struct runner *r, struct cell *restrict ci,
                                struct spart *restrict sparts_i,
                                int *restrict ind, int scount,
                                struct cell *restrict cj, const double *shift) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < scount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *restrict spi = &sparts_i[ind[pid]];

    const double pix = spi->x[0] - (shift[0]);
    const double piy = spi->x[1] - (shift[1]);
    const double piz = spi->x[2] - (shift[2]);
    const float hi = spi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                     (float)(piz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, hj, spi, pj, NULL,
                                            cosmo, e->feedback_props,
                                            e->ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        error("No subset feedback iact functions do (or should) exist! .");
        /* runner_iact_nonsym_feedback_apply(r2, dx, hi, hj, spi, pj, xpj,
         * cosmo, ti_current); */
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts The #spart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_STARS(struct runner *r, struct cell *restrict ci,
                          struct spart *restrict sparts, int *restrict ind,
                          int scount) {
#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;

  /* Early abort? */
  if (count_i == 0) return;

  /* Loop over the parts in ci. */
  for (int spid = 0; spid < scount; spid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *spi = &sparts[ind[spid]];
    const float spix[3] = {(float)(spi->x[0] - ci->loc[0]),
                           (float)(spi->x[1] - ci->loc[1]),
                           (float)(spi->x[2] - ci->loc[2])};
    const float hi = spi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                            (float)(pj->x[1] - ci->loc[1]),
                            (float)(pj->x[2] - ci->loc[2])};
      float dx[3] = {spix[0] - pjx[0], spix[1] - pjx[1], spix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, pj->h, spi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_feedback_density(r2, dx, hi, pj->h, spi, pj, NULL,
                                            cosmo, e->feedback_props,
                                            e->ti_current);
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
        error("No subset feedback iact functions do (or should) exist!");
        /* runner_iact_nonsym_feedback_apply(r2, dx, hi, pj->h, spi, pj, xpj, */
        /*                                   cosmo, e, ti_current); */
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Determine which version of DOSELF1_SUBSET_STARS needs to be called
 * depending on the optimisation level.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts The #spart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts,
                                 const int *restrict ind, const int scount) {

  //  DOSELF1_SUBSET_STARS(r, ci, sparts, ind, scount);
}

/**
 * @brief Determine which version of DOPAIR1_SUBSET_STARS needs to be called
 * depending on the orientation of the cells or whether DOPAIR1_SUBSET_STARS
 * needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #spart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR1_SUBSET_BRANCH_STARS(struct runner *r, struct cell *restrict ci,
                                 struct spart *restrict sparts_i,
                                 const int *restrict ind, const int scount,
                                 struct cell *restrict cj) {

#if 0
  
  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->hydro.count == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_STARS
  DOPAIR1_SUBSET_STARS_NAIVE(r, ci, sparts_i, ind, scount, cj, shift);
#else
  /* Get the sorting index. */
  int sid = 0;
  for (int k = 0; k < 3; k++)
    sid = 3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0)
                         ? 0
                         : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1);

  /* Switch the cells around? */
  const int flipped = runner_flip[sid];
  sid = sortlistID[sid];

  /* Has the cell cj been sorted? */
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

  DOPAIR1_SUBSET_STARS(r, ci, sparts_i, ind, scount, cj, sid, flipped, shift);
#endif

#endif
}

/**
 * @brief Determine which version of DOSELF1_STARS needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_STARS(struct runner *r, struct cell *c,
                          const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->stars.count == 0 || c->hydro.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

  /* Did we mess up the recursion? */
  if (!limit_max_h && c->stars.h_max_active * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e) || !cell_are_spart_drifted(c, e))
    error("Interacting undrifted cell.");

  DOSELF1_STARS(r, c, limit_min_h, limit_max_h);
}

#define RUNNER_CHECK_SORT(TYPE, PART, cj, ci, sid)                          \
  ({                                                                        \
    const struct sort_entry *restrict sort_j =                              \
        cell_get_##TYPE##_sorts(cj, sid);                                   \
                                                                            \
    for (int pjd = 0; pjd < cj->TYPE.count; pjd++) {                        \
      const struct PART *p = &cj->TYPE.parts[sort_j[pjd].i];                \
      if (PART##_is_inhibited(p, e)) continue;                              \
                                                                            \
      const float d = p->x[0] * runner_shift[sid][0] +                      \
                      p->x[1] * runner_shift[sid][1] +                      \
                      p->x[2] * runner_shift[sid][2];                       \
      if ((fabsf(d - sort_j[pjd].d) - cj->TYPE.dx_max_sort) >               \
              1.0e-4 * max(fabsf(d), cj->TYPE.dx_max_sort_old) &&           \
          (fabsf(d - sort_j[pjd].d) - cj->TYPE.dx_max_sort) >               \
              cj->width[0] * 1.0e-10)                                       \
        error(                                                              \
            "particle shift diff exceeds dx_max_sort in cell cj. "          \
            "cj->nodeID=%d "                                                \
            "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->" #TYPE                \
            ".dx_max_sort=%e "                                              \
            "cj->" #TYPE                                                    \
            ".dx_max_sort_old=%e, cellID=%lld super->cellID=%lld"           \
            "cj->depth=%d cj->maxdepth=%d",                                 \
            cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->TYPE.dx_max_sort, \
            cj->TYPE.dx_max_sort_old, cj->cellID, cj->hydro.super->cellID,  \
            cj->depth, cj->maxdepth);                                       \
    }                                                                       \
  })

/**
 * @brief Determine which version of DOPAIR1_STARS needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_STARS needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                          const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  const int do_ci = ci->stars.count != 0 && cj->hydro.count != 0 &&
                    cell_is_active_stars(ci, e);
  const int do_cj = cj->stars.count != 0 && ci->hydro.count != 0 &&
                    cell_is_active_stars(cj, e);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_spart_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_part_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_ci && (!(ci->stars.sorted & (1 << sid)) ||
                ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_ci && (!(cj->hydro.sorted & (1 << sid)) ||
                cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

  if (do_cj && (!(ci->hydro.sorted & (1 << sid)) ||
                ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_cj && (!(cj->stars.sorted & (1 << sid)) ||
                cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_STARS
  DOPAIR1_STARS_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
#else
  DO_SYM_PAIR1_STARS(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
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
void DOSUB_PAIR1_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                       int recurse_below_h_max, const int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  message("sid=%d", sid);

  /* What kind of pair are we doing here? */
  const int do_ci = ci->stars.count != 0 && cj->hydro.count != 0 &&
                    cell_is_active_stars(ci, e);
  const int do_cj = cj->stars.count != 0 && ci->hydro.count != 0 &&
                    cell_is_active_stars(cj, e);

  /* Should we even bother? */
  if (!do_ci && !do_cj) return;

  if (!ci->split || ci->stars.count < space_recurse_size_pair_stars ||
      !cj->split || cj->stars.count < space_recurse_size_pair_stars) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH_STARS(r, ci, cj, /*limit_h_min=*/0,
                         /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Should we change the recursion regime because we encountered a large
         particle? */
      if (!recurse_below_h_max && (!cell_can_recurse_in_pair_stars_task1(ci) ||
                                   !cell_can_recurse_in_pair_stars_task1(cj)))
        recurse_below_h_max = 1;

      /* message("Multi-level PAIR! ci->count=%d cj->count=%d", ci->hydro.count,
       */
      /* 	      cj->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR1_BRANCH_STARS(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_STARS(r, ci->progeny[pid], cj->progeny[pjd],
                          recurse_below_h_max,
                          /*gettimer=*/0);
      }
    }
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_STARS);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_STARS(struct runner *r, struct cell *c,
                       int recurse_below_h_max, const int gettimer) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (c->hydro.count == 0 || c->stars.count == 0 ||
      !cell_is_active_stars(c, r->e))
    return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->stars.count < space_recurse_size_self_stars) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH_STARS(r, c, /*limit_h_min=*/0,
                         /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_self_stars_task1(c))
      recurse_below_h_max = 1;

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH_STARS(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1_STARS(r, c->progeny[k], recurse_below_h_max,
                          /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_STARS(r, c->progeny[k], c->progeny[j],
                              recurse_below_h_max,
                              /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_STARS);
}

struct cell *FIND_SUB_STARS(const struct cell *const ci,
                            const struct spart *const sparts, const int *ind) {

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&sparts[ind[0]] >= &ci->progeny[k]->stars.parts[0] &&
            &sparts[ind[0]] <
                &ci->progeny[k]->stars.parts[ci->progeny[k]->stars.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  return sub;
}

void DOSUB_PAIR_SUBSET_STARS(struct runner *r, struct cell *ci,
                             struct spart *sparts, const int *ind,
                             const int scount, struct cell *cj,
                             const int gettimer) {
  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (ci->stars.count == 0 || cj->hydro.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse? */
  if (ci->split && cell_can_recurse_in_pair_stars_task1(ci) && cj->split &&
      cell_can_recurse_in_pair_stars_task1(cj)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_STARS(ci, sparts, ind);

    /* Get the type of pair and flip ci/cj if needed. */
    double shift[3];
    const int sid = space_getsid(s, &ci, &cj, shift);

    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
        DOSUB_PAIR_SUBSET_STARS(r, ci->progeny[pid], sparts, ind, scount,
                                cj->progeny[pjd], /*gettimer=*/0);
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
        DOSUB_PAIR_SUBSET_STARS(r, cj->progeny[pjd], sparts, ind, scount,
                                ci->progeny[pid], /*gettimer=*/0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_stars(ci, e) && cj->hydro.count > 0) {

    /* Do any of the cells need to be drifted first? */
    if (cell_is_active_stars(ci, e)) {
      if (!cell_are_spart_drifted(ci, e)) error("Cell should be drifted!");
      if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");
    }

    DOPAIR1_SUBSET_BRANCH_STARS(r, ci, sparts, ind, scount, cj);
  }
}

void DOSUB_SELF_SUBSET_STARS(struct runner *r, struct cell *ci,
                             struct spart *sparts, const int *ind,
                             const int scount, const int gettimer) {

  const struct engine *e = r->e;

  /* Should we even bother? */
  if (ci->hydro.count == 0 || ci->stars.count == 0) return;
  if (!cell_is_active_hydro(ci, e)) return;

  /* Recurse? */
  if (ci->split && cell_can_recurse_in_self_stars_task1(ci)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_STARS(ci, sparts, ind);

    /* Loop over all progeny. */
    DOSUB_SELF_SUBSET_STARS(r, sub, sparts, ind, scount, /*gettimer=*/0);
    for (int j = 0; j < 8; j++)
      if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
        DOSUB_PAIR_SUBSET_STARS(r, sub, sparts, ind, scount, ci->progeny[j],
                                /*gettimer=*/0);
  }

  /* Otherwise, compute self-interaction. */
  else
    DOSELF1_SUBSET_BRANCH_STARS(r, ci, sparts, ind, scount);
}
