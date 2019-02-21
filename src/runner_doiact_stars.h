/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOSELF1_STARS(f) PASTE(runner_doself_stars, f)
#define DOSELF1_STARS _DOSELF1_STARS(FUNCTION)

#define _DO_NONSYM_PAIR1_STARS(f) PASTE(runner_do_nonsym_pair_stars, f)
#define DO_NONSYM_PAIR1_STARS _DO_NONSYM_PAIR1_STARS(FUNCTION)

#define _DOPAIR1_STARS(f) PASTE(runner_dopair_stars, f)
#define DOPAIR1_STARS _DOPAIR1_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_STARS(f) PASTE(runner_dopair_subset_stars, f)
#define DOPAIR1_SUBSET_STARS _DOPAIR1_SUBSET_STARS(FUNCTION)

#define _DOSELF1_SUBSET_STARS(f) PASTE(runner_doself_subset_stars, f)
#define DOSELF1_SUBSET_STARS _DOSELF1_SUBSET_STARS(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_doself_subset_branch_stars, f)
#define DOSELF1_SUBSET_BRANCH_STARS _DOSELF1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_STARS(f) \
  PASTE(runner_dopair_subset_branch_stars, f)
#define DOPAIR1_SUBSET_BRANCH_STARS _DOPAIR1_SUBSET_BRANCH_STARS(FUNCTION)

#define _DOSUB_SUBSET_STARS(f) PASTE(runner_dosub_subset_stars, f)
#define DOSUB_SUBSET_STARS _DOSUB_SUBSET_STARS(FUNCTION)

#define _DOSELF1_BRANCH_STARS(f) PASTE(runner_doself_branch_stars, f)
#define DOSELF1_BRANCH_STARS _DOSELF1_BRANCH_STARS(FUNCTION)

#define _DOPAIR1_BRANCH_STARS(f) PASTE(runner_dopair_branch_stars, f)
#define DOPAIR1_BRANCH_STARS _DOPAIR1_BRANCH_STARS(FUNCTION)

#define _DOSUB_PAIR1_STARS(f) PASTE(runner_dosub_pair_stars, f)
#define DOSUB_PAIR1_STARS _DOSUB_PAIR1_STARS(FUNCTION)

#define _DOSUB_SELF1_STARS(f) PASTE(runner_dosub_self_stars, f)
#define DOSUB_SELF1_STARS _DOSUB_SELF1_STARS(FUNCTION)

#define _IACT_STARS(f) PASTE(runner_iact_nonsym_stars, f)
#define IACT_STARS _IACT_STARS(FUNCTION)

/**
 * @brief Calculate the number density of #part around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_STARS(struct runner *r, struct cell *c, int timer) {
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;
  if (c->hydro.count == 0 && c->stars.count == 0) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int count = c->hydro.count;
  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts[sid];
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];
      const float hj = pj->h;

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

      if (r2 > 0.f && r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, si, pj, a, H);
      }
    } /* loop over the parts in ci. */
  }   /* loop over the sparts in ci. */
}

/**
 * @brief Calculate the number density of cj #part around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_STARS(struct runner *r, struct cell *restrict ci,
                           struct cell *restrict cj) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (!cell_is_active_stars(ci, e)) return;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

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

      if (r2 < hig2) IACT_STARS(r2, dx, hi, hj, si, pj, a, H);

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

void DOPAIR1_STARS(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj, int timer) {

  if (ci->stars.count != 0 && cj->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS(r, ci, cj);
  if (cj->stars.count != 0 && ci->hydro.count != 0)
    DO_NONSYM_PAIR1_STARS(r, cj, ci);
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
void DOPAIR1_SUBSET_STARS(struct runner *r, struct cell *restrict ci,
                          struct spart *restrict sparts_i, int *restrict ind,
                          int scount, struct cell *restrict cj,
                          const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < scount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *restrict spi = &sparts_i[ind[pid]];
    double spix[3];
    for (int k = 0; k < 3; k++) spix[k] = spi->x[k] - shift[k];
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

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = spix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS(r2, dx, hi, pj->h, spi, pj, a, H);
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

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;

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
      const float hj = pj->h;

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
      if (r2 > 0.f && r2 < hig2) {
        IACT_STARS(r2, dx, hi, hj, spi, pj, a, H);
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
                                 int *restrict ind, int scount) {

  DOSELF1_SUBSET_STARS(r, ci, sparts, ind, scount);
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
                                 int *restrict ind, int scount,
                                 struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  DOPAIR1_SUBSET_STARS(r, ci, sparts_i, ind, scount, cj, shift);
}

void DOSUB_SUBSET_STARS(struct runner *r, struct cell *ci, struct spart *sparts,
                        int *ind, int scount, struct cell *cj, int sid,
                        int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (!cell_is_active_stars(ci, e) &&
      (cj == NULL || !cell_is_active_stars(cj, e)))
    return;

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

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_stars_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET_STARS(r, sub, sparts, ind, scount, NULL, -1, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET_STARS(r, sub, sparts, ind, scount, ci->progeny[j], -1,
                             0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF1_SUBSET_BRANCH_STARS(r, ci, sparts, ind, scount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_stars_task(ci) &&
        cell_can_recurse_in_pair_stars_task(cj)) {

      /* Get the type of pair if not specified explicitly. */
      double shift[3] = {0.0, 0.0, 0.0};
      sid = space_getsid(s, &ci, &cj, shift);

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[4], sparts, ind, scount,
                               cj->progeny[3], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[3], sparts, ind, scount,
                               ci->progeny[4], -1, 0);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[2], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[2], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[1], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[6], sparts, ind, scount,
                               cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[5], sparts, ind, scount,
                               ci->progeny[6], -1, 0);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[1], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[1], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[1], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[1], -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[1], sparts, ind, scount,
                               cj->progeny[6], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[6], sparts, ind, scount,
                               ci->progeny[1], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[3], sparts, ind, scount,
                               cj->progeny[6], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[6], sparts, ind, scount,
                               ci->progeny[3], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[5], sparts, ind, scount,
                               cj->progeny[6], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[6], sparts, ind, scount,
                               ci->progeny[5], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[0], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[2], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[4], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET_STARS(r, ci->progeny[7], sparts, ind, scount,
                               cj->progeny[6], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET_STARS(r, cj->progeny[6], sparts, ind, scount,
                               ci->progeny[7], -1, 0);
          break;
      }

    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_stars(ci, e) || cell_is_active_stars(cj, e)) {

      /* Do any of the cells need to be drifted first? */
      if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");

      DOPAIR1_SUBSET_BRANCH_STARS(r, ci, sparts, ind, scount, cj);
    }

  } /* otherwise, pair interaction. */
}

/**
 * @brief Determine which version of DOSELF1_STARS needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_STARS(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->stars.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  DOSELF1_STARS(r, c, 1);
}

#define RUNNER_CHECK_SORT(TYPE, PART, cj, ci, sid)                          \
  ({                                                                        \
    const struct entry *restrict sort_j = cj->TYPE.sort[sid];               \
                                                                            \
    for (int pjd = 0; pjd < cj->TYPE.count; pjd++) {                        \
      const struct PART *p = &cj->TYPE.parts[sort_j[pjd].i];                \
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
            "cj->" #TYPE ".dx_max_sort_old=%e",                             \
            cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->TYPE.dx_max_sort, \
            cj->TYPE.dx_max_sort_old);                                      \
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
void DOPAIR1_BRANCH_STARS(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;
  const int ci_active = cell_is_active_stars(ci, e);
  const int cj_active = cell_is_active_stars(cj, e);
  const int do_ci = (ci->stars.count != 0 && cj->hydro.count != 0 && ci_active);
  const int do_cj = (cj->stars.count != 0 && ci->hydro.count != 0 && cj_active);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_spart_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_ci && (!(ci->stars.sorted & (1 << sid)) ||
                ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_ci && (!(cj->hydro.sorted & (1 << sid)) ||
                cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

  if (do_cj &&
      (!cell_are_part_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* Have the cells been sorted? */
  if (do_cj && (!(ci->hydro.sorted & (1 << sid)) ||
                ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin))
    error("Interacting unsorted cells.");

  if (do_cj && (!(cj->stars.sorted & (1 << sid)) ||
                cj->stars.dx_max_sort_old > space_maxreldx * cj->dmin))
    error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
  if (do_ci) {
    RUNNER_CHECK_SORT(hydro, part, cj, ci, sid);
    RUNNER_CHECK_SORT(stars, spart, ci, cj, sid);
  }

  if (do_cj) {
    RUNNER_CHECK_SORT(hydro, part, ci, cj, sid);
    RUNNER_CHECK_SORT(stars, spart, cj, ci, sid);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  DOPAIR1_STARS(r, ci, cj, 1);
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction linking the cells
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1_STARS(struct runner *r, struct cell *ci, struct cell *cj,
                       int sid, int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  int should_do = ci->stars.count != 0 && cj->hydro.count != 0 &&
                  cell_is_active_stars(ci, e);
  should_do |= cj->stars.count != 0 && ci->hydro.count != 0 &&
               cell_is_active_stars(cj, e);
  if (!should_do) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_stars_task(ci) &&
      cell_can_recurse_in_pair_stars_task(cj)) {

    /* Different types of flags. */
    switch (sid) {

      /* Regular sub-cell interactions of a single cell. */
      case 0: /* (  1 ,  1 ,  1 ) */
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        break;

      case 1: /* (  1 ,  1 ,  0 ) */
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[1], -1, 0);
        break;

      case 2: /* (  1 ,  1 , -1 ) */
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        break;

      case 3: /* (  1 ,  0 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[2], -1, 0);
        break;

      case 4: /* (  1 ,  0 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[0], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[1], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[2], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[3], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[3], -1, 0);
        break;

      case 5: /* (  1 ,  0 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[3], -1, 0);
        break;

      case 6: /* (  1 , -1 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[2], -1, 0);
        break;

      case 7: /* (  1 , -1 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[3], -1, 0);
        break;

      case 8: /* (  1 , -1 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[4], cj->progeny[3], -1, 0);
        break;

      case 9: /* (  0 ,  1 ,  1 ) */
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[4], -1, 0);
        break;

      case 10: /* (  0 ,  1 ,  0 ) */
        if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[0], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[4], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[1], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[4], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[5], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[5], -1, 0);
        break;

      case 11: /* (  0 ,  1 , -1 ) */
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[6], cj->progeny[5], -1, 0);
        break;

      case 12: /* (  0 ,  0 ,  1 ) */
        if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[1], cj->progeny[0], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[1], cj->progeny[2], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[1], cj->progeny[4], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[1], cj->progeny[6], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[2], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[3], cj->progeny[6], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[4], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[5], cj->progeny[6], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1_STARS(r, ci->progeny[7], cj->progeny[6], -1, 0);
        break;
    }

  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci = ci->stars.count != 0 && cj->hydro.count != 0 &&
                      cell_is_active_stars(ci, e);
    const int do_cj = cj->stars.count != 0 && ci->hydro.count != 0 &&
                      cell_is_active_stars(cj, e);

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (!cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");

      /* Do any of the cells need to be sorted first? */
      if (!(ci->stars.sorted & (1 << sid)) ||
          ci->stars.dx_max_sort_old > ci->dmin * space_maxreldx)
        error("Interacting unsorted cell (sparts).");

      if (!(cj->hydro.sorted & (1 << sid)) ||
          cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx)
        error("Interacting unsorted cell (parts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (parts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");

      /* Do any of the cells need to be sorted first? */
      if (!(ci->hydro.sorted & (1 << sid)) ||
          ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx)
        error("Interacting unsorted cell (parts).");

      if (!(cj->stars.sorted & (1 << sid)) ||
          cj->stars.dx_max_sort_old > cj->dmin * space_maxreldx)
        error("Interacting unsorted cell (sparts).");
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_STARS(r, ci, cj);
  }
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_STARS(struct runner *r, struct cell *ci, int gettimer) {

  /* Should we even bother? */
  if (ci->hydro.count == 0 || ci->stars.count == 0 ||
      !cell_is_active_stars(ci, r->e))
    return;

  /* Recurse? */
  if (cell_can_recurse_in_self_stars_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_STARS(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_STARS(r, ci->progeny[k], ci->progeny[j], -1, 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_spart_drifted(ci, r->e)) error("Interacting undrifted cell.");

    DOSELF1_BRANCH_STARS(r, ci);
  }
}
