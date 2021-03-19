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
void DOSELF1_RT(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Cosmological terms */
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->stars.count == 0) return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Drift the cell to the current timestep if needed. */
  if (!cell_are_spart_drifted(c, r->e))
    error("Interacting undrifted cell (spart): %lld", c->cellID);
  if (!cell_are_part_drifted(c, r->e))
    error("Interacting undrifted cell (part): %lld", c->cellID);
#endif

  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;

  const int scount = c->stars.count;
  const int count = c->hydro.count;

  /* Loop over the sparts in cell */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in c. */
    struct spart *restrict si = &sparts[sid];

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!rt_is_spart_active_in_loop(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
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
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
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

      if (r2 < hig2) IACT_RT(r2, dx, hi, hj, si, pj, a, H);
    }
  }

  if (timer) TIMER_TOC(TIMER_DOSELF_RT);
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
void DOPAIR1_NONSYM_RT_NAIVE(struct runner *r, struct cell *ci,
                             struct cell *cj) {

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Cosmological terms */
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;

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

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!rt_is_spart_active_in_loop(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

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
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
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

      if (r2 < hig2) IACT_RT(r2, dx, hi, hj, si, pj, a, H);

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
void DO_SYM_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj,
                     const int sid, const double *shift) {

  TIMER_TIC;

  const struct engine *e = r->e;

  /* Cosmological terms */
  const struct cosmology *cosmo = e->cosmology;
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  const int do_ci_stars =
      (ci->nodeID == e->nodeID) && rt_should_do_cell_pair(ci, cj, e);
  const int do_cj_stars =
      (cj->nodeID == e->nodeID) && rt_should_do_cell_pair(cj, ci, e);

  if (do_ci_stars) {

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
    const struct sort_entry *restrict sort_i = cell_get_stars_sorts(ci, sid);

    /* Get some other useful values. */
    const double hi_max = ci->stars.h_max * kernel_gamma - rshift;
    const int count_i = ci->stars.count;
    const int count_j = cj->hydro.count;
    struct spart *restrict sparts_i = ci->stars.parts;
    struct part *restrict parts_j = cj->hydro.parts;
    const double dj_min = sort_j[0].d;
    const float dx_max = (ci->stars.dx_max_sort + cj->hydro.dx_max_sort);
    const float hydro_dx_max_rshift = cj->hydro.dx_max_sort - rshift;

    /* TODO: maybe change the order of the loops for better performance? */
    /* Loop over the sparts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict spi = &sparts_i[sort_i[pid].i];
      const float hi = spi->h;

      /* Skip inhibited particles */
      if (spart_is_inhibited(spi, e)) continue;

      /* Skip inactive particles */
      if (!rt_is_spart_active_in_loop(spi, e)) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spi->x[0], spi->x[1], spi->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double di = dist + hi * kernel_gamma + hydro_dx_max_rshift;
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
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
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
          IACT_RT(r2, dx, hi, hj, spi, pj, a, H);
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* do_ci_stars */

  if (do_cj_stars) {
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_hydro_sorts(ci, sid);
    const struct sort_entry *restrict sort_j = cell_get_stars_sorts(cj, sid);

    /* Get some other useful values. */
    const double hj_max = cj->stars.h_max * kernel_gamma;
    const int count_i = ci->hydro.count;
    const int count_j = cj->stars.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct spart *restrict sparts_j = cj->stars.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const float dx_max = (ci->hydro.dx_max_sort + cj->stars.dx_max_sort);
    const float hydro_dx_max_rshift = ci->hydro.dx_max_sort - rshift;

    /* TODO: maybe change the order of the loops for better performance? */
    /* Loop over the sparts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth spart in cj. */
      struct spart *spj = &sparts_j[sort_j[pjd].i];
      const float hj = spj->h;

      /* Skip inhibited particles */
      if (spart_is_inhibited(spj, e)) continue;

      /* Skip inactive particles */
      if (!rt_is_spart_active_in_loop(spj, e)) continue;

      /* Compute distance from the other cell. */
      const double px[3] = {spj->x[0], spj->x[1], spj->x[2]};
      float dist = px[0] * runner_shift[sid][0] + px[1] * runner_shift[sid][1] +
                   px[2] * runner_shift[sid][2];

      /* Is there anything we need to interact with ? */
      const double dj = dist - hj * kernel_gamma - hydro_dx_max_rshift;
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
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
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
void DOPAIR1_RT_NAIVE(struct runner *r, struct cell *ci, struct cell *cj,
                      int timer) {

  TIMER_TIC;
  const struct engine *restrict e = r->e;

  const int do_stars_in_ci =
      (cj->nodeID == r->e->nodeID) && rt_should_do_cell_pair(ci, cj, e);
  if (do_stars_in_ci) DOPAIR1_NONSYM_RT_NAIVE(r, ci, cj);

  const int do_stars_in_cj =
      (ci->nodeID == r->e->nodeID) && rt_should_do_cell_pair(cj, ci, e);
  if (do_stars_in_cj) DOPAIR1_NONSYM_RT_NAIVE(r, cj, ci);

  if (timer) TIMER_TOC(TIMER_DOPAIR_RT);
}

/**
 * @brief Determine which version of DOSELF1_RT needs to be called
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_BRANCH_RT(struct runner *r, struct cell *c, int timer) {

#ifdef RT_DEBUG
  /* Before an early exit, loop over all parts and sparts in this cell
   * and mark that we checked these particles */

  const struct engine *e = r->e;

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

  /* early exit? */
  if (!rt_should_do_cell(c, r->e)) return;
  DOSELF1_RT(r, c, timer);
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
                       int timer) {

  const struct engine *restrict e = r->e;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  const int do_stars_ci =
      (cj->nodeID == r->e->nodeID) && rt_should_do_cell_pair(ci, cj, e);
  const int do_stars_cj =
      (ci->nodeID == r->e->nodeID) && rt_should_do_cell_pair(cj, ci, e);

  /* Anything to do here? */
  if (!do_stars_ci && !do_stars_cj) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (do_stars_ci) {

    /* Check that cells are drifted. */
    if (!cell_are_spart_drifted(ci, e))
      error("Interacting undrifted stars in cells i=%12lld, j=%12lld",
            ci->cellID, cj->cellID);
    if (!cell_are_part_drifted(cj, e))
      error("Interacting undrifted hydro in cells i=%12lld, j=%12lld",
            ci->cellID, cj->cellID);

    /* Have the cells been sorted? */
    if (!(ci->stars.sorted & (1 << sid)) ||
        ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin)
      error("Interacting unsorted cells: %lld %d", ci->cellID, sid);
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
      error("Interacting unsorted cells: %lld %lld", ci->cellID, cj->cellID);
  }

  if (do_stars_cj) {

    /* Check that cells are drifted. */
    if (!cell_are_spart_drifted(cj, e)) {
      error("Interacting undrifted stars in cells j=%12lld, i=%12lld",
            cj->cellID, ci->cellID);
    }
    if (!cell_are_part_drifted(ci, e))
      error("Interacting undrifted hydro in cells j=%12lld, i=%12lld",
            cj->cellID, ci->cellID);

    /* Have the cells been sorted? */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
      error("Interacting unsorted cells: %lld %d", cj->cellID, sid);

    if ((!(cj->stars.sorted & (1 << sid)) ||
         cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
      error("Interacting unsorted cells: %lld %lld", ci->cellID, cj->cellID);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  if (do_stars_ci || do_stars_cj) {
#ifdef SWIFT_USE_NAIVE_INTERACTIONS_RT
    DOPAIR1_RT_NAIVE(r, ci, cj, 1);
#else
    DO_SYM_PAIR1_RT(r, ci, cj, sid, shift);
#endif
  }
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_RT(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (!rt_should_do_cell(c, r->e)) {

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

  /* Recurse? */
  if (cell_can_recurse_in_self_stars_task(c)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1_RT(r, c->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (c->progeny[j] != NULL)
            DOSUB_PAIR1_RT(r, c->progeny[k], c->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {
    DOSELF1_BRANCH_RT(r, c, 0);
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
                    int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = rt_should_do_cell_pair(ci, cj, e);
  const int should_do_cj = rt_should_do_cell_pair(cj, ci, e);
  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_stars_task(ci, cj) &&
      cell_can_recurse_in_pair_stars_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_RT(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    /* do full checks again, space_getsid() might swap ci/cj pointers */
    const int do_ci_stars =
        (cj->nodeID == e->nodeID) && rt_should_do_cell_pair(ci, cj, e);
    const int do_cj_stars =
        (ci->nodeID == e->nodeID) && rt_should_do_cell_pair(cj, ci, e);

    if (do_ci_stars || do_cj_stars) DOPAIR1_BRANCH_RT(r, ci, cj, 0);
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_PAIR_RT);
}
