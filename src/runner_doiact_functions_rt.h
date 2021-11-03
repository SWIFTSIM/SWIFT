/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#include "rt.h"
#include "rt_active.h"
#include "runner_doiact_rt.h"

/**
 * @brief Function for self-type interaction between stars and hydro particles
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_RT(struct runner *r, const struct cell *c, const int limit_min_h,
                const int limit_max_h) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->stars.count == 0) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int count = c->hydro.count;
  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the sparts in cell */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in c. */
    struct spart *si = &sparts[sid];
    const char depth_i = si->depth_h;
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!rt_is_spart_active_in_loop(si, e)) continue;

    /* Skip particles not in the range of h we care about */
    if (depth_i > max_depth) continue;
    if (depth_i < min_depth) continue;

    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cell */
    for (int pid = 0; pid < count; pid++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pid];
      const float hj = pj->h;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Skip inactive particles. */
      if (!rt_is_part_active_in_loop(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");

      if (r2 < hig2 && e->rt_props->hydro_controlled_injection) {

        const integertime_t ti_current = e->ti_current;
        const integertime_t pti_end =
            get_integer_time_end(ti_current, pj->time_bin);
        const integertime_t sti_end =
            get_integer_time_end(ti_current, si->time_bin);

        if (sti_end < ti_current)
          error(
              "s-particle in an impossible time-zone! sp->ti_end=%lld "
              "e->ti_current=%lld",
              sti_end, ti_current);
        if (pti_end < ti_current)
          error(
              "particle in an impossible time-zone! p->ti_end=%lld "
              "e->ti_current=%lld",
              pti_end, ti_current);
        if (pti_end > sti_end)
          message(
              "WARNING: Got star that whose time step ends before the "
              "interacting "
              "hydro particle's time step ends. This needs to be dealt with.");
      }
#endif

      if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_RT(r2, dx, hi, hj, si, pj, a, H);
      }
    } /* loop over the parts in ci. */
  }   /* loop over the sparts in ci. */

  TIMER_TOC(TIMER_DOSELF_RT);
}

/**
 * @brief Function for non-symmetric pair-type interaction between stars
 *        and hydro particles. Will interact star particles of cell i
 *        with hydro particles of cell j.
 *
 *
 * @param r runner task
 * @param ci the first cell, where we take star particles from
 * @param cj the second cell, where we take hydro particles from
 */
void DOPAIR1_NONSYM_RT_NAIVE(struct runner *r, const struct cell *restrict ci,
                             const struct cell *restrict cj,
                             const int limit_min_h, const int limit_max_h) {

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
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

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts_i[sid];
    const char depth_i = si->depth_h;

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!rt_is_spart_active_in_loop(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

#ifdef SWIFT_DEBUG_CHECKS
    if (hi > ci->stars.h_max_active)
      error("Particle has h larger than h_max_active");
#endif

    /* Skip particles not in the range of h we care about */
    if (depth_i > max_depth) continue;
    if (depth_i < min_depth) continue;

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      const float hj = pj->h;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Skip inactive particles. */
      if (!rt_is_part_active_in_loop(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef RT_DEBUG
      if (r2 < hig2 && e->rt_props->hydro_controlled_injection) {

        const integertime_t ti_current = e->ti_current;
        const integertime_t pti_end =
            get_integer_time_end(ti_current, pj->time_bin);
        const integertime_t sti_end =
            get_integer_time_end(ti_current, si->time_bin);

        if (sti_end < ti_current)
          error(
              "s-particle in an impossible time-zone! sp->ti_end=%lld "
              "e->ti_current=%lld",
              sti_end, ti_current);
        if (pti_end < ti_current)
          error(
              "particle in an impossible time-zone! p->ti_end=%lld "
              "e->ti_current=%lld",
              pti_end, ti_current);
        if (pti_end > sti_end)
          message(
              "WARNING: Got star that whose time step ends before the "
              "interacting "
              "hydro particle's time step ends. This needs to be dealt with.");
      }
#endif

      if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_RT(r2, dx, hi, hj, si, pj, a, H);
      }

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Function for pair-type interaction between stars
 *        and hydro particles. Will interact hydro particles of cell i
 *        with star particles of cell j.
 *
 * @param r #runner task
 * @param ci the first #cell
 * @param cj the second #cell
 * @param timer 1 if the time is to be recorded.
 */
void DO_SYM_PAIR1_RT(struct runner *r, const struct cell *restrict ci,
                     const struct cell *restrict cj, const int limit_min_h,
                     const int limit_max_h, const int sid,
                     const double shift[3]) {

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  const int do_ci_stars =
      (ci->nodeID == e->nodeID) && rt_should_iact_cell_pair(ci, cj, e);
  const int do_cj_stars =
      (cj->nodeID == e->nodeID) && rt_should_iact_cell_pair(cj, ci, e);

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->dmin != cj->dmin) error("Cells of different size!");
#endif

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
#endif
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;

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
    const double hi_max =
        min(h_max, ci->stars.h_max_active) * kernel_gamma - rshift;
    const int count_i = ci->stars.count;
    const int count_j = cj->hydro.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct part *restrict parts_j = cj->hydro.parts;
    const double dj_min = sort_j[0].d;
    const float dx_max = (ci->stars.dx_max_sort + cj->hydro.dx_max_sort);

    /* TODO: maybe change the order of the loops for better performance? */

    /* Loop over the *active* sparts in ci that are within range (on the axis)
       of any particle in cj. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[sort_i[pid].i];
      const char depth_i = spi->depth_h;
      const float hi = spi->h;

      /* Skip inhibited particles */
      if (spart_is_inhibited(spi, e)) continue;

      /* Skip inactive particles */
      if (!rt_is_spart_active_in_loop(spi, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hi > ci->stars.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_i > max_depth) continue;
      if (depth_i < min_depth) continue;

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
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

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Skip inactive particles. */
        if (!rt_is_part_active_in_loop(pj, e)) continue;

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
        if (spi->ti_drift != e->ti_current)
          error("Particle spi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");

        if (r2 < hig2 && e->rt_props->hydro_controlled_injection) {

          const integertime_t ti_current = e->ti_current;
          const integertime_t pti_end =
              get_integer_time_end(ti_current, pj->time_bin);
          const integertime_t sti_end =
              get_integer_time_end(ti_current, spi->time_bin);

          if (sti_end < ti_current)
            error(
                "s-particle in an impossible time-zone! sp->ti_end=%lld "
                "e->ti_current=%lld",
                sti_end, ti_current);
          if (pti_end < ti_current)
            error(
                "particle in an impossible time-zone! p->ti_end=%lld "
                "e->ti_current=%lld",
                pti_end, ti_current);
          if (pti_end > sti_end)
            message(
                "WARNING: Got star that whose time step ends before the "
                "interacting "
                "hydro particle's time step ends. This needs to be dealt "
                "with.");
        }

#endif

        /* Hit or miss? */
        if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_RT(r2, dx, hi, hj, spi, pj, a, H);
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
    const double hj_max = min(h_max, cj->stars.h_max_active) * kernel_gamma;
    const int count_i = ci->hydro.count;
    const int count_j = cj->stars.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float dx_max = (ci->hydro.dx_max_sort + cj->stars.dx_max_sort);

    /* TODO: maybe change the order of the loops for better performance? */
    /* Loop over the *active* sparts in cj that are within range (on the axis)
       of any particle in ci. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth spart in cj. */
      struct spart *spj = &sparts_j[sort_j[pjd].i];
      const char depth_j = spj->depth_h;
      const float hj = spj->h;

      /* Skip inhibited particles */
      if (spart_is_inhibited(spj, e)) continue;

      /* Skip inactive particles */
      if (!rt_is_spart_active_in_loop(spj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (hj > cj->stars.h_max_active)
        error("Particle has h larger than h_max_active");
#endif

      /* Skip particles not in the range of h we care about */
      if (depth_j > max_depth) continue;
      if (depth_j < min_depth) continue;

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
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

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        /* Skip inactive particles. */
        if (!rt_is_part_active_in_loop(pi, e)) continue;

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
        if (spj->ti_drift != e->ti_current)
          error("Particle spj not drifted to current time");

        if (r2 < hjg2 && e->rt_props->hydro_controlled_injection) {

          const integertime_t ti_current = e->ti_current;
          const integertime_t pti_end =
              get_integer_time_end(ti_current, pi->time_bin);
          const integertime_t sti_end =
              get_integer_time_end(ti_current, spj->time_bin);

          if (sti_end < ti_current)
            error(
                "s-particle in an impossible time-zone! sp->ti_end=%lld "
                "e->ti_current=%lld",
                sti_end, ti_current);
          if (pti_end < ti_current)
            error(
                "particle in an impossible time-zone! p->ti_end=%lld "
                "e->ti_current=%lld",
                pti_end, ti_current);
          if (pti_end > sti_end)
            message(
                "WARNING: Got star that whose time step ends before the "
                "interacting "
                "hydro particle's time step ends. This needs to be dealt "
                "with.");
        }
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_RT(r2, dx, hj, hi, spj, pi, a, H);
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR_RT);
}

/**
 * @brief Function for pair-type interaction between stars
 *        and hydro particles. Will interact hydro particles of cell i
 *        with star particles of cell j.
 *
 * @param r #runner task
 * @param ci the first #cell
 * @param cj the second #cell
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_RT_NAIVE(struct runner *r, const struct cell *restrict ci,
                      const struct cell *restrict cj, const int limit_min_h,
                      const int limit_max_h) {

  TIMER_TIC;

  const int do_stars_in_ci =
      (cj->nodeID == r->e->nodeID) && rt_should_iact_cell_pair(ci, cj, r->e);
  if (do_stars_in_ci)
    DOPAIR1_NONSYM_RT_NAIVE(r, ci, cj, limit_min_h, limit_max_h);

  const int do_stars_in_cj =
      (ci->nodeID == r->e->nodeID) && rt_should_iact_cell_pair(cj, ci, r->e);
  if (do_stars_in_cj)
    DOPAIR1_NONSYM_RT_NAIVE(r, cj, ci, limit_min_h, limit_max_h);

  TIMER_TOC(TIMER_DOPAIR_RT);
}

/**
 * @brief Determine which version of DOSELF1_RT needs to be called
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_BRANCH_RT(struct runner *r, const struct cell *c,
                       const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

#ifdef RT_DEBUG
  /* Before an early exit, loop over all parts and sparts in this cell
   * and mark that we checked these particles */

  if (c->hydro.count > 0) {
    struct part *restrict parts = c->hydro.parts;

    /* Loop over the parts in cell */
    for (int pid = 0; pid < c->hydro.count; pid++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pid];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Skip inactive particles. */
      if (!rt_is_part_active_in_loop(pj, e)) continue;

      rt_debugging_check_injection_part(pj, e->rt_props);
    }
  }

  if (c->stars.count > 0) {
    struct spart *restrict sparts = c->stars.parts;

    /* Loop over the parts in cell */
    for (int sid = 0; sid < c->stars.count; sid++) {

      /* Get a pointer to the ith spart. */
      struct spart *restrict si = &sparts[sid];

      /* Skip inhibited particles. */
      if (spart_is_inhibited(si, e)) continue;

      /* Skip inactive particles */
      if (!rt_is_spart_active_in_loop(si, e)) continue;

      rt_debugging_check_injection_spart(si, e->rt_props);
    }
  }
#endif

  /* Anything to do here? */
  if (!rt_should_iact_cell(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->stars.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  if (!limit_max_h && c->stars.h_max_active * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");
#endif

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e))
    error("Interacting undrifted cell (hydro).");
  if (!cell_are_spart_drifted(c, e))
    error("Interacting undrifted cell (stars).");

  DOSELF1_RT(r, c, limit_min_h, limit_max_h);
}

/**
 * @brief Determine which version of DOPAIR1_RT needs to be called
 *
 * @param r #runner
 * @param ci The first #cell
 * @param cj The second #cell
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_BRANCH_RT(struct runner *r, struct cell *ci, struct cell *cj,
                       const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  const int do_ci =
      (cj->nodeID == r->e->nodeID) && rt_should_iact_cell_pair(ci, cj, e);
  const int do_cj =
      (ci->nodeID == r->e->nodeID) && rt_should_iact_cell_pair(cj, ci, e);

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
    error("Interacting unsorted cells (ci stars).");

  if (do_ci && (!(cj->hydro.sorted & (1 << sid)) ||
                cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells (cj hydro).");

  /* Have the cells been sorted? */
  if (do_cj && (!(ci->hydro.sorted & (1 << sid)) ||
                ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells. (ci hydro)");

  if (do_cj && (!(cj->stars.sorted & (1 << sid)) ||
                cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells. (cj stars)");

#ifdef SWIFT_USE_NAIVE_INTERACTIONS_RT
  DOPAIR1_RT_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
#else
  DO_SYM_PAIR1_RT(r, ci, cj, limit_min_h, limit_max_h, sid, shift);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_RT(struct runner *r, struct cell *c, int recurse_below_h_max,
                    const int timer) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (!rt_should_iact_cell(c, r->e)) {

#ifdef RT_DEBUG
    /* Before an early exit, loop over all parts and sparts in this cell
     * and mark that we checked these particles */

    const struct engine *e = r->e;

    if (c->hydro.count > 0) {
      struct part *restrict parts = c->hydro.parts;
      const int count = c->hydro.count;

      /* Loop over the parts in cell */
      for (int pid = 0; pid < count; pid++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pid];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Skip inactive particles. */
        if (!rt_is_part_active_in_loop(pj, e)) continue;

        rt_debugging_check_injection_part(pj, e->rt_props);
      }
    }

    if (c->stars.count > 0) {
      struct spart *restrict sparts = c->stars.parts;
      const int scount = c->stars.count;

      /* Loop over the parts in cell */
      for (int sid = 0; sid < scount; sid++) {

        /* Get a pointer to the ith spart. */
        struct spart *restrict si = &sparts[sid];

        /* Skip inhibited particles. */
        if (spart_is_inhibited(si, e)) continue;

        /* Skip inactive particles */
        if (!rt_is_spart_active_in_loop(si, e)) continue;

        rt_debugging_check_injection_spart(si, e->rt_props);
      }
    }
#endif

    /* exit early if there is nothing to do */
    return;
  }

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->stars.count < space_recurse_size_self_stars) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH_RT(r, c, /*limit_h_min=*/0,
                      /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself_stars_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH_RT(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1_RT(r, c->progeny[k], recurse_below_h_max,
                       /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_RT(r, c->progeny[k], c->progeny[j], recurse_below_h_max,
                           /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_SELF_RT);
}

/**
 * @brief Compute grouped sub-cell interactions for pair tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj,
                    int recurse_below_h_max, const int timer) {
  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Should we even bother? */
  const int do_ci = rt_should_iact_cell_pair(ci, cj, e);
  const int do_cj = rt_should_iact_cell_pair(cj, ci, e);
  if (!do_ci && !do_cj) return;

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->stars.count < space_recurse_size_pair_stars ||
      !cj->split || cj->stars.count < space_recurse_size_pair_stars) {

    /* Do any of the cells need to be sorted first?
     * Since h_max might have changed, we may not have sorted at this level */
    if (do_ci) {
      if (!(ci->stars.sorted & (1 << sid)) ||
          ci->stars.dx_max_sort_old > ci->dmin * space_maxreldx) {
        runner_do_stars_sort(r, ci, (1 << sid), 0, 0);
      }
      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
        runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*clock=*/0);
      }
    }
    if (do_cj) {
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
        runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                             /*clock=*/0);
      }
      if (!(cj->stars.sorted & (1 << sid)) ||
          cj->stars.dx_max_sort_old > cj->dmin * space_maxreldx) {
        runner_do_stars_sort(r, cj, (1 << sid), 0, 0);
      }
    }

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH_RT(r, ci, cj, /*limit_h_min=*/0,
                      /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Both ci and cj are split */

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max && (!cell_can_recurse_in_subpair_stars_task(ci) ||
                                 !cell_can_recurse_in_subpair_stars_task(cj))) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* Do any of the cells need to be sorted first?
       * Since h_max might have changed, we may not have sorted at this level */
      if (do_ci) {
        if (!(ci->stars.sorted & (1 << sid)) ||
            ci->stars.dx_max_sort_old > ci->dmin * space_maxreldx) {
          runner_do_stars_sort(r, ci, (1 << sid), 0, 0);
        }
        if (!(cj->hydro.sorted & (1 << sid)) ||
            cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx) {
          runner_do_hydro_sort(r, cj, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                               /*clock=*/0);
        }
      }
      if (do_cj) {
        if (!(ci->hydro.sorted & (1 << sid)) ||
            ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx) {
          runner_do_hydro_sort(r, ci, (1 << sid), /*cleanup=*/0, /*lock=*/1,
                               /*clock=*/0);
        }
        if (!(cj->stars.sorted & (1 << sid)) ||
            cj->stars.dx_max_sort_old > cj->dmin * space_maxreldx) {
          runner_do_stars_sort(r, cj, (1 << sid), 0, 0);
        }
      }

      /* message("Multi-level PAIR! ci->count=%d cj->count=%d", ci->hydro.count,
       */
      /* 	      cj->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR1_BRANCH_RT(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_RT(r, ci->progeny[pid], cj->progeny[pjd],
                       recurse_below_h_max,
                       /*gettimer=*/0);
      }
    }
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_PAIR_RT);
}
