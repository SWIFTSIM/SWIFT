
/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1_NOSORT(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  if (!cell_is_drifted(ci, e)) cell_drift_particles(ci, e);
  if (!cell_is_drifted(cj, e)) cell_drift_particles(cj, e);

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;

  if (cell_is_active(ci, e)) {

    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count_i; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[pid];
      if (!part_is_active(pi, e)) continue;
      const float hi = pi->h;

      double pix[3];
      for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      const float hig2 = hi * hi * kernel_gamma2;

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[pjd];

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {
          IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);
        }

      } /* loop over the parts in cj. */

    } /* loop over the parts in ci. */

  } /* Cell ci is active */

  if (cell_is_active(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pj = &parts_j[pjd];
      if (!part_is_active(pj, e)) continue;
      const float hj = pj->h;

      double pjx[3];
      for (int k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Loop over the parts in ci. */
      for (int pid = 0; pid < count_i; pid++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pi = &parts_i[pid];

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pjx[k] - pi->x[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {
          IACT_NONSYM(r2, dx, hj, pi->h, pj, pi);
        }

      } /* loop over the parts in ci. */

    } /* loop over the parts in cj. */

  } /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR2_NOSORT(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  if (!cell_is_drifted(ci, e)) cell_drift_particles(ci, e);
  if (!cell_is_drifted(cj, e)) cell_drift_particles(cj, e);

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  space_getsid_and_swap_cells(e->s, &ci, &cj, shift);

  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;

  if (cell_is_active(ci, e)) {

    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count_i; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[pid];
      if (!part_is_active(pi, e)) continue;
      const float hi = pi->h;

      double pix[3];
      for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      const float hig2 = hi * hi * kernel_gamma2;

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[pjd];
        const float hjg2 = pj->h * pj->h * kernel_gamma2;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hjg2) {
          IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);
        }

      } /* loop over the parts in cj. */

    } /* loop over the parts in ci. */

  } /* Cell ci is active */

  if (cell_is_active(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pj = &parts_j[pjd];
      if (!part_is_active(pj, e)) continue;
      const float hj = pj->h;

      double pjx[3];
      for (int k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Loop over the parts in ci. */
      for (int pid = 0; pid < count_i; pid++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pi = &parts_i[pid];
        const float hig2 = pi->h * pi->h * kernel_gamma2;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pjx[k] - pi->x[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2 || r2 < hig2) {
          IACT_NONSYM(r2, dx, hj, pi->h, pj, pi);
        }

      } /* loop over the parts in ci. */

    } /* loop over the parts in cj. */

  } /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR_SUBSET_NOSORT(struct runner *r, struct cell *restrict ci,
                          struct part *restrict parts_i, int *restrict ind,
                          int count, struct cell *restrict cj) {

  struct engine *e = r->e;

  TIMER_TIC;

  const int count_j = cj->count;
  struct part *restrict parts_j = cj->parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts_i. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[ind[pid]];
    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    if (!part_is_active(pi, e))
      error("Trying to correct smoothing length of inactive particle !");

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */

  TIMER_TOC(timer_dopair_subset);
}
