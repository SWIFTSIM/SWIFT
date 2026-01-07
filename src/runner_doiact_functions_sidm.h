/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2025 Katy Proctor (katy.proctor@fysik.su.se)
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

#include "runner_doiact_sidm.h"
/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1_SIDM(struct runner *r, const struct cell *c, const int limit_min_h,
                  const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  struct sipart *siparts = c->sidm.parts;
  const int sicount = c->sidm.count;

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
  int sicountdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     sicount * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < sicount; k++) {
    const struct sipart *sip = &siparts[k];
    const char depth = sip->depth_h;
    if ((sipart_is_active(sip, e)) && (depth >= min_depth) &&
        (depth <= max_depth)) {
      indt[sicountdt] = k;
      sicountdt += 1;
    }
  }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over *all* the particles (i.e. the ones to update and not to update).
   *
   * Note the additional condition to make the loop abort if all the active
   * particles have been processed. */
  for (int siid = 0; siid < sicount && firstdt < sicountdt; siid++) {

    /* Get a pointer to the ith particle. */
    struct sipart *restrict sipi = &siparts[siid];
    const char depth_i = sipi->depth_h;

    /* Skip inhibited particles. */
    if (sipart_is_inhibited(sipi, e)) continue;

    /* Get the particle position and (square of) search radius. */
    const double sipix[3] = {sipi->x[0], sipi->x[1], sipi->x[2]};
    const float hi = sipi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle in the range of h we care about? */
    const int update_i = (sipart_is_active(sipi, e)) &&
                         (depth_i >= min_depth) && (depth_i <= max_depth);

    /* If false then it can only act as a neighbour of others */
    if (!update_i) {

      /* Loop over the particles we want to update. */
      for (int sipjd = firstdt; sipjd < sicountdt; sipjd++) {

        /* Get a pointer to the jth particle. (by construction sipi != sipj) */
        struct sipart *restrict sipj = &siparts[indt[sipjd]];

        /* Skip inhibited particles. */
        if (sipart_is_inhibited(sipj, e)) continue;

        /* This particle's (square of) search radius. */
        const float hj = sipj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (sipi->ti_drift != e->ti_current)
          error("Particle sipi not drifted to current time");
        if (sipj->ti_drift != e->ti_current)
          error("Particle sipj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double sipjx[3] = {sipj->x[0], sipj->x[1], sipj->x[2]};
        const float dx[3] = {(float)(sipjx[0] - sipix[0]),
                             (float)(sipjx[1] - sipix[1]),
                             (float)(sipjx[2] - sipix[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Hit or miss? */
        if (r2 < hjg2) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          IACT_NONSYM_SIDM(r2, dx, hj, hi, sipj, sipi, a, H);
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
      for (int sijd = siid + 1; sijd < sicount; sijd++) {

        /* Get a pointer to the jth particle (by construction pi != pj). */
        struct sipart *restrict sipj = &siparts[sijd];
        const char depth_j = sipj->depth_h;

        /* Skip inhibited particles. */
        if (sipart_is_inhibited(sipj, e)) continue;

        /* This particle's (square of) search radius. */
        const float hj = sipj->h;
        const float hjg2 = hj * hj * kernel_gamma2;

#if defined(SWIFT_DEBUG_CHECKS)
        /* Check that particles have been drifted to the current time */
        if (sipi->ti_drift != e->ti_current)
          error("Particle sipi not drifted to current time");
        if (sipj->ti_drift != e->ti_current)
          error("Particle sipj not drifted to current time");
#endif

        /* Compute the (square of) pairwise distance. */
        const double sipjx[3] = {sipj->x[0], sipj->x[1], sipj->x[2]};
        float dx[3] = {(float)(sipix[0] - sipjx[0]),
                       (float)(sipix[1] - sipjx[1]),
                       (float)(sipix[2] - sipjx[2])};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

        /* Decide which of the two particles to update */

        /* We know pi is active and in the right range of h
         * -> Only check the distance to pj */
        const int doi = (r2 < hig2);

        /* We know nothing about pj
         * -> Check whether it is in the right range of h
         * -> Check the distance to pi */
        const int doj = (sipart_is_active(sipj, e)) && (depth_j >= min_depth) &&
                        (depth_j <= max_depth) && (r2 < hjg2);

        /* Hit or miss? */
        if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif
          /* Update both sipi and sipj */

          IACT_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);
        } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hi < h_min || hi >= h_max)
            error("Inappropriate h for this level!");
#endif
          /* Update only sipi */

          IACT_NONSYM_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);
        } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
          if (hj < h_min || hj >= h_max)
            error("Inappropriate h for this level!");
#endif

          /* Update only sipj */

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];
          IACT_NONSYM_SIDM(r2, dx, hj, hi, sipj, sipi, a, H);
        } /* Hit or miss */
      } /* loop over all other particles. */
    } /* pi is active */
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF_SIDM);
}

/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR1_SIDM_NAIVE(struct runner *r, const struct cell *restrict ci,
                        const struct cell *restrict cj, const int limit_min_h,
                        const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_sidm(ci, e) && !cell_is_active_sidm(cj, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int sicount_i = ci->sidm.count;
  const int sicount_j = cj->sidm.count;
  struct sipart *restrict siparts_i = ci->sidm.parts;
  struct sipart *restrict siparts_j = cj->sidm.parts;

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

  /* Loop over the siparts in ci. */
  for (int siid = 0; siid < sicount_i; siid++) {

    /* Get a hold of the ith part in ci. */
    struct sipart *restrict sipi = &siparts_i[siid];

    /* Skip inhibited particles. */
    if (sipart_is_inhibited(sipi, e)) continue;

    const int sipi_active = sipart_is_active(sipi, e);
    const char depth_i = sipi->depth_h;
    const float hi = sipi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float sipix[3] = {(float)(sipi->x[0] - (cj->loc[0] + shift[0])),
                            (float)(sipi->x[1] - (cj->loc[1] + shift[1])),
                            (float)(sipi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int sijd = 0; sijd < sicount_j; sijd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sijd];

      /* Skip inhibited particles. */
      if (sipart_is_inhibited(sipj, e)) continue;

      const int sipj_active = sipart_is_active(sipj, e);
      const char depth_j = sipj->depth_h;
      const float hj = sipj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - cj->loc[0]),
                              (float)(sipj->x[1] - cj->loc[1]),
                              (float)(sipj->x[2] - cj->loc[2])};
      float dx[3] = {sipix[0] - sipjx[0], sipix[1] - sipjx[1],
                     sipix[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (sipi->ti_drift != e->ti_current)
        error("Particle sipi not drifted to current time");
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      const int doi = sipi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = sipj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

      /* Hit or miss? */
      if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);
      }
      if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM_SIDM(r2, dx, hj, hi, sipj, sipi, a, H);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR_SIDM);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1_SIDM_NAIVE(struct runner *r, const struct cell *c,
                        const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_sidm(c, e)) return;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int sicount = c->sidm.count;
  struct sipart *siparts = c->sidm.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the parts in ci. */
  for (int siid = 0; siid < sicount; siid++) {

    /* Get a hold of the ith part in ci. */
    struct sipart *restrict sipi = &siparts[siid];

    /* Skip inhibited particles. */
    if (sipart_is_inhibited(sipi, e)) continue;

    const int sipi_active = sipart_is_active(sipi, e);
    const char depth_i = sipi->depth_h;
    const float hi = sipi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float sipix[3] = {(float)(sipi->x[0] - c->loc[0]),
                            (float)(sipi->x[1] - c->loc[1]),
                            (float)(sipi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int sijd = siid + 1; sijd < sicount; sijd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts[sijd];

      /* Skip inhibited particles. */
      if (sipart_is_inhibited(sipj, e)) continue;

      const float hj = sipj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int sipj_active = sipart_is_active(sipj, e);
      const char depth_j = sipj->depth_h;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - c->loc[0]),
                              (float)(sipj->x[1] - c->loc[1]),
                              (float)(sipj->x[2] - c->loc[2])};
      float dx[3] = {sipix[0] - sipjx[0], sipix[1] - sipjx[1],
                     sipix[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = sipi_active && (r2 < hig2) && (depth_i >= min_depth) &&
                      (depth_i <= max_depth);
      const int doj = sipj_active && (r2 < hjg2) && (depth_j >= min_depth) &&
                      (depth_j <= max_depth);

#if defined(SWIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (sipi->ti_drift != e->ti_current)
        error("Particle sipi not drifted to current time");
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      /* Hit or miss? */
      if (doi && doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);

      } else if (doi) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_NONSYM_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);

      } else if (doj) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hj < h_min || hj >= h_max) error("Inappropriate h for this level!");
#endif

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM_SIDM(r2, dx, hj, hi, sipj, sipi, a, H);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF_SIDM);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param siparts_i The #sipart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param sicount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_SIDM_NAIVE(struct runner *r, const struct cell *restrict ci,
                               struct sipart *restrict siparts_i,
                               const int *ind, const int sicount,
                               const struct cell *restrict cj,
                               const double shift[3]) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  const int sicount_j = cj->sidm.count;
  struct sipart *restrict siparts_j = cj->sidm.parts;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the siparts_i. */
  for (int siid = 0; siid < sicount; siid++) {

    /* Get a hold of the ith part in ci. */
    struct sipart *restrict sipi = &siparts_i[ind[siid]];
    double sipix[3];
    for (int k = 0; k < 3; k++) sipix[k] = sipi->x[k] - shift[k];
    const float hi = sipi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!sipart_is_active(sipi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int sijd = 0; sijd < sicount_j; sijd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sijd];

      if (sipart_is_inhibited(sipj, e)) continue;

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = sipix[k] - sipj->x[k];
        r2 += dx[k] * dx[k];
      }

#if defined(SWIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (sipi->ti_drift != e->ti_current)
        error("Particle sipi not drifted to current time");
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM_SIDM(r2, dx, hi, sipj->h, sipi, sipj, a, H);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */
}

/**
 * @brief Determine which version of DOPAIR_SUBSET needs to be called depending
 * on the
 * orientation of the cells or whether DOPAIR_SUBSET needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param siparts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param sicount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR1_SUBSET_BRANCH_SIDM(struct runner *r,
                                const struct cell *restrict ci,
                                struct sipart *restrict siparts_i,
                                const int *ind, const int sicount,
                                struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->sidm.count == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  DOPAIR1_SUBSET_SIDM_NAIVE(r, ci, siparts_i, ind, sicount, cj, shift);
}

/**
 * @brief Compute the interactions between a cell, but only for the
 *      given indices in c.
 *
 * @param r The #runner.
 * @param i The #cell.
 * @param siparts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param sicount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_SIDM(struct runner *r, const struct cell *c,
                         struct sipart *restrict siparts, const int *ind,
                         const int sicount) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms and physical constants */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_cell = c->sidm.count;
  struct sipart *restrict siparts_j = c->sidm.parts;
  /* Loop over the parts in ci. */
  for (int siid = 0; siid < count_cell; siid++) {

    /* Get a hold of the ith part in ci. */
    struct sipart *sipi = &siparts[ind[siid]];
    const float sipix[3] = {(float)(sipi->x[0] - c->loc[0]),
                            (float)(sipi->x[1] - c->loc[1]),
                            (float)(sipi->x[2] - c->loc[2])};
    const float hi = sipi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    if (sipart_is_inhibited(sipi, e)) continue;

    /* Loop over the parts in cj. */
    for (int sijd = 0; sijd < count_cell; sijd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sijd];

      /* Skip oneself */
      if (sipi == sipj) continue;

      if (sipart_is_inhibited(sipj, e)) continue;

      const float hj = sipj->h;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - c->loc[0]),
                              (float)(sipj->x[1] - c->loc[1]),
                              (float)(sipj->x[2] - c->loc[2])};
      float dx[3] = {sipix[0] - sipjx[0], sipix[1] - sipjx[1],
                     sipix[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS)
      /* Check that particles have been drifted to the current time */
      if (sipi->ti_drift != e->ti_current)
        error("Particle sipi not drifted to current time");
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM_SIDM(r2, dx, hi, hj, sipi, sipj, a, H);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */
}

/**
 * @brief Determine which version of DOSELF_SUBSET needs to be called depending
 * on the optimisation level.

 * @param r The #runner.
 * @param ci The first #cell.
 * @param siparts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param sicount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BRANCH_SIDM(struct runner *r, const struct cell *ci,
                                struct sipart *restrict siparts, const int *ind,
                                const int sicount) {

  DOSELF1_SUBSET_SIDM(r, ci, siparts, ind, sicount);
}

/**
 * @brief Determine which version of DOSELF1 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOSELF1_BRANCH_SIDM(struct runner *r, const struct cell *c,
                         const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (c->sidm.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_sidm(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->sidm.h_max_old * kernel_gamma > c->dmin)
    if (!limit_max_h && c->sidm.h_max_active * kernel_gamma > c->dmin)
      error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");
#endif

  /* Check that cells are drifted. */
  if (!cell_are_sipart_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_SIDM_NAIVE(r, c, limit_min_h, limit_max_h);
#else
  DOSELF1_SIDM(r, c, limit_min_h, limit_max_h);
#endif
}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 */
void DOPAIR1_BRANCH_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                         const int limit_min_h, const int limit_max_h) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (ci->sidm.count == 0 || cj->sidm.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_sidm(ci, e) && !cell_is_active_sidm(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_sipart_drifted(ci, e) || !cell_are_sipart_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* No sorted intreactions here -> use the naive ones */
  DOPAIR1_SIDM_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
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
void DOSUB_PAIR1_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                      int recurse_below_h_max, const int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  if (!cell_is_active_sidm(ci, e) && !cell_is_active_sidm(cj, e)) return;
  if (ci->sidm.count == 0 || cj->sidm.count == 0) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->sidm.count < space_recurse_size_pair_sidm ||
      !cj->split || cj->sidm.count < space_recurse_size_pair_sidm) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH_SIDM(r, ci, cj, /*limit_h_min=*/0,
                        /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Both ci and cj are split */

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max && (!cell_can_recurse_in_subpair_sidm_task(ci) ||
                                 !cell_can_recurse_in_subpair_sidm_task(cj))) {
      DOPAIR1_BRANCH_SIDM(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_SIDM(r, ci->progeny[pid], cj->progeny[pjd],
                         recurse_below_h_max,
                         /*gettimer=*/0);
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR_SIDM);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param recurse_below_h_max Are we currently recursing at a level where we
 * violated the h < cell size condition.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SIDM(struct runner *r, struct cell *c, int recurse_below_h_max,
                      const int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  if (c->sidm.count == 0 || !cell_is_active_sidm(c, r->e)) return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->sidm.count < space_recurse_size_self_sidm) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH_SIDM(r, c, /*limit_h_min=*/0,
                        /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself_sidm_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->hydro.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH_SIDM(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1_SIDM(r, c->progeny[k], recurse_below_h_max, /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_SIDM(r, c->progeny[k], c->progeny[j],
                             recurse_below_h_max,
                             /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_SIDM);
}

/**
 * @brief Find which sub-cell of a cell contain the subset of particles given
 * by the list of indices.
 *
 * Will throw an error if the sub-cell can't be found.
 *
 * @param c The #cell
 * @param siparts An array of #sipart.
 * @param ind Index of the #sipart's in the particle array to find in the subs.
 */
struct cell *FIND_SUB_SIDM(const struct cell *const c,
                           const struct sipart *const siparts,
                           const int *const ind) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!c->split) error("Can't search for subs in a non-split cell");
#endif

  /* Find out in which sub-cell of ci the parts are.
   *
   * Note: We only need to check the first particle in the list */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL) {
      if (&siparts[ind[0]] >= &c->progeny[k]->sidm.parts[0] &&
          &siparts[ind[0]] <
              &c->progeny[k]->sidm.parts[c->progeny[k]->sidm.count]) {
        return c->progeny[k];
      }
    }
  }
  error("Invalid sub!");
  return NULL;
}

void DOSUB_PAIR1_SUBSET_SIDM(struct runner *r, struct cell *ci,
                             struct sipart *siparts, const int *ind,
                             const int sicount, struct cell *cj,
                             const int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (ci->sidm.count == 0 || cj->sidm.count == 0) return;
  if (!cell_is_active_sidm(ci, e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_pair_sidm_task(ci) &&
      cell_can_recurse_in_pair_sidm_task(cj)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_SIDM(ci, siparts, ind);

    /* Get the type of pair and flip ci/cj if needed. */
    double shift[3];
    const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_SUBSET_SIDM(r, ci->progeny[pid], siparts, ind, sicount,
                                cj->progeny[pjd],
                                /*gettimer=*/0);
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
        DOSUB_PAIR1_SUBSET_SIDM(r, cj->progeny[pjd], siparts, ind, sicount,
                                ci->progeny[pid],
                                /*gettimer=*/0);
    }

  }
  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_sidm(ci, e)) {

    /* Do any of the cells need to be drifted first? */
    if (!cell_are_sipart_drifted(cj, e)) error("Cell should be drifted!");

    DOPAIR1_SUBSET_BRANCH_SIDM(r, ci, siparts, ind, sicount, cj);
  }
}

void DOSUB_SELF1_SUBSET_SIDM(struct runner *r, struct cell *ci,
                             struct sipart *siparts, const int *ind,
                             const int sicount, const int gettimer) {

  const struct engine *e = r->e;

  /* Should we even bother? */
  if (ci->sidm.count == 0) return;
  if (!cell_is_active_sidm(ci, e)) return;

  /* Recurse? */
  if (ci->split && cell_can_recurse_in_self_sidm_task(ci)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_SIDM(ci, siparts, ind);

    /* Loop over all progeny. */
    DOSUB_SELF1_SUBSET_SIDM(r, sub, siparts, ind, sicount, /*gettimer=*/0);
    for (int j = 0; j < 8; j++)
      if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
        DOSUB_PAIR1_SUBSET_SIDM(r, sub, siparts, ind, sicount, ci->progeny[j],
                                /*gettimer=*/0);
  }

  /* Otherwise, compute self-interaction. */
  else
    DOSELF1_SUBSET_BRANCH_SIDM(r, ci, siparts, ind, sicount);
}