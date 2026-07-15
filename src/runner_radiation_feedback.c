/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner_radiation_feedback.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "runner.h"
#include "timers.h"

/**
 * @brief Top-level steering function for HII ionization feedback.
 *
 * This function recurses down the cell hierarchy until it reaches a level
 * where the cell size is comparable to the maximum HII radius of the stars
 * within it (the operational working level). If the cell is small enough,
 * it stops recursing and hands off control to the sub-task coordinator.
 *
 * @param r The #runner thread.
 * @param c The #cell to process.
 * @param timer If true, records the timing of this operation.
 */
void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer) {
#ifdef IONIZATION_FEEDBACK_LOOP
  struct engine *e = r->e;
  struct stars_props *star_props = e->stars_properties;

  /* Determine the search radius. Contrary to many implementations out there,
     we do not iterate if the search radius is too small. At the next step,
     the star will have a larger search radius since it's hii radius has
     increased.

     TODO: Implement retry if search radius is too small */
  const float r_hii_max =
      max(c->stars.h_hii_max_active, c->stars.h_max_active) * kernel_gamma;
  const float max_search_radius = star_props->max_HII_search_radius;
  const float interaction_limit =
      min(search_radius_factor * r_hii_max, max_search_radius);

  /* Anything to do here? */
  if (c->stars.count == 0 || c->hydro.count == 0 || !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->cellID != c->stars.radiation_level->cellID && timer == 1) {
    warning(
        "Not running the HII reionization task on the radiation level (c = "
        "%lld, "
        "c->stars.radiation_level = %lld)!",
        c->cellID, c->stars.radiation_level->cellID);
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  for (struct link *l = c->stars.radiation_in; l != NULL; l = l->next) {
    /* We have already handled the self case */
    if (l->t->type == task_type_self) continue;

    struct cell *cj = l->t->cj;
    struct cell *ci = l->t->ci;
    /* ci/cj may sit above their own hydro.super (radiation_level can
     * coarsen past it), in which case hydro.super is NULL -- print -1
     * rather than dereferencing it. */
    const long long ci_super_id =
        ci->hydro.super != NULL ? ci->hydro.super->cellID : -1;
    const long long cj_super_id =
        cj->hydro.super != NULL ? cj->hydro.super->cellID : -1;
    message(
        "[%lld, %lld] hydro super: %lld , %lld | radiation_level: %lld %lld",
        ci->cellID, cj->cellID, ci_super_id, cj_super_id,
        ci->stars.radiation_level->cellID, cj->stars.radiation_level->cellID);
  }
#endif

  TIMER_TIC;

  /* Is the cell split */
  if (c->split) {
    /* Keep recursing deeper into the super-cell hierarchy. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];
        runner_do_stars_hii_ionization_feedback(r, cp, 0);
        c->stars.h_hii_max = max(c->stars.h_hii_max, cp->stars.h_hii_max);
        c->stars.h_hii_max_active =
            max(c->stars.h_hii_max_active, cp->stars.h_hii_max_active);
        c->stars.h_max_active =
            max(c->stars.h_max_active, cp->stars.h_max_active);
      }
    }
  } else {
    /* We have reached the 'Working Level' */
    runner_dosub_stars_hii_ionization_feedback(r, c, interaction_limit);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Coordinator function to manage HII neighbor gathering and ionization
 * loops for active stars.
 *
 * This function acts at the working-level cell context. For each active star,
 * it allocates a temporary neighbor buffer, triggers a strict local
 * self-interaction search, and iterates over the scheduled engine task links
 * (`radiation_in`) to discover interacting neighbor branches. Once neighbors
 * are gathered, it *tags* them as ionized through atomic operations.
 *
 * A star that still has photons left after one pass retries in one of two
 * ways, chosen by why the pass fell short: if the neighbour buffer filled
 * up, more candidates may exist within the SAME search radius (bumped out
 * by the nearest-max_ngbs cut), so it retries there
 * (Stars:max_HII_retry_full_buffer); if the buffer did NOT fill up (every
 * reachable candidate was already found and ionized) the radius itself is
 * the bottleneck, so it is expanded
 * (Stars:max_HII_radius_expansion_tries, Stars:HII_radius_expansion_factor)
 * and the pass retried at the larger radius. Both retries are scoped to
 * the individual star inside the loop below: a star that exhausts its
 * budget on the first pass costs nothing extra, and only a star that
 * genuinely needs more reach pays for the additional passes.
 *
 * @param r The #runner thread.
 * @param c The working-level #cell containing the active stars to process.
 * @param interaction_limit The initial search range boundary for neighbor
 * mapping, before any retry-driven radius expansion.
 */
void runner_dosub_stars_hii_ionization_feedback(struct runner *r,
                                                struct cell *c,
                                                const float interaction_limit) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const struct feedback_props *feedback_props = e->feedback_props;

  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const double time = e->time;

  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  /* Anything to do here? */
  if (c->stars.count == 0 || c->hydro.count == 0 || !cell_is_active_stars(c, e))
    return;

  const struct stars_props *star_props = e->stars_properties;
  const int max_retry_full_buffer = star_props->max_HII_retry_full_buffer;
  const int max_radius_expansion_tries =
      star_props->max_HII_radius_expansion_tries;
  const float radius_expansion_factor = star_props->HII_radius_expansion_factor;
  const float max_search_radius = star_props->max_HII_search_radius;

  struct hii_neighbor ngb_buffer[max_ngbs];

  for (int i = 0; i < scount; i++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[i];

    /* Is this part within the timestep? */
    if (spart_is_inhibited(si, e)) continue;
    if (!spart_is_active(si, e)) continue;
    if (!feedback_is_HII_ionization_active(si, e)) continue;
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
    message("Star %lld can do ionization! r_hii = %e", si->id,
            si->h_hii * kernel_gamma);
#endif
#ifdef SWIFT_DEBUG_CHECKS
    /* Check that particles have been drifted to the current time */
    if (si->ti_drift != ti_current)
      error(
          "Particle si (%lld) not drifted to current time c = %lld, "
          "c->super = %lld",
          si->id, c->cellID, c->super->cellID);
#endif

    /* Logic: a pass can fall short of exhausting the star's photon budget
     * for two DIFFERENT reasons, which need two different retries:
     *
     * 1. The buffer filled up (count_found == max_ngbs): there may be MORE
     *    candidates within the SAME radius that got bumped out by the
     *    nearest-max_ngbs cut. Retry at the same radius (cheap: the
     *    volume searched doesn't grow) -- existing max_retry_full_buffer
     *    mechanism.
     * 2. The buffer did NOT fill up (every not-yet-ionized particle within
     *    the current radius was found and processed) but photons remain:
     *    the RADIUS itself is the bottleneck, not the buffer. Retrying at
     *    the same radius would just re-find nothing new. Only expanding
     *    the radius can make progress -- max_HII_radius_expansion_tries.
     *
     * Conflating these (retrying the same radius regardless of which
     * happened) leaves a star stalled whenever case 2 occurs: with no
     * further growth this pass, the only remaining way for the search
     * radius to grow is the NEXT rebuild's interaction_limit, which is
     * only pulled forward by h_hii itself (unchanged here, since nothing
     * new was found) or by the unrelated h_max term drifting up over
     * possibly many rebuild cycles. Expanding immediately keeps this
     * confined to the one star that actually needs it: neither retry
     * kind touches the outer loop over sparts, so a star that already
     * exhausts its budget on the first pass costs nothing extra. */
    float dynamic_search_radius = interaction_limit;
    int num_retry_full_buffer = 0;
    int num_radius_expansions = 0;
    while (1) {
      int count_found = 0;

      /* First loop over particles in the current cell */
      runner_doself_stars_hii_ionization_feedback(
          r, c, si, dynamic_search_radius, ngb_buffer, max_ngbs, &count_found);

      /* Now loop over particles in the neighboring cells via task links */
      for (struct link *l = c->stars.radiation_level->stars.radiation_in;
           l != NULL; l = l->next) {

        struct cell *cj = NULL;
        if (l->t->type == task_type_self) {
          cj = c->stars.radiation_level;
        } else {
          cj = (l->t->cj == c->stars.radiation_level) ? l->t->ci : l->t->cj;
        }

        /* Get the relative distance between the pairs, wrapping. */
        double shift[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < 3; k++) {
          if (cj->loc[k] - c->loc[k] < -e->s->dim[k] / 2)
            shift[k] = e->s->dim[k];
          else if (cj->loc[k] - c->loc[k] > e->s->dim[k] / 2)
            shift[k] = -e->s->dim[k];
        }

        /* Get the sorting index. */
        int sid = 0;
        for (int k = 0; k < 3; k++)
          sid = 3 * sid + ((cj->loc[k] - c->loc[k] + shift[k] < 0)   ? 0
                           : (cj->loc[k] - c->loc[k] + shift[k] > 0) ? 2
                                                                     : 1);

        /* Switch the cells around? */
        const int flipped = runner_flip[sid];
        sid = sortlistID[sid];

        runner_do_stars_hii_ionization_feedback_branch(
            r, c, cj, sid, flipped, shift, si, dynamic_search_radius,
            ngb_buffer, max_ngbs, &count_found);
      } /* Neighbour search */

      const int buffer_was_full = (count_found == max_ngbs);

      /***************************************************/
      /* It's time to sort the gas particles */
      if (count_found > 0) {
        /* Verify that the neighbors are properly sorted by distance */
        runner_do_stars_hii_ionization_feedback_check_sort(ngb_buffer,
                                                           count_found);

        const integertime_t ti_begin =
            get_integer_time_begin(e->ti_current - 1, si->time_bin);

        /* Now let's ionize the gas particles! */
        for (int k = 0; k < count_found; k++) {

          const int pixel = ngb_buffer[k].pixel;

          /* This particle's pixel has no photons left; skip it (not
             break -- other entries in the buffer may belong to a pixel
             that still has budget). */
          if (feedback_get_star_ionization_rate(si, pixel) <= 0.0) {
            continue;
          }

          struct part *pj = ngb_buffer[k].p;
          struct xpart *xpj = ngb_buffer[k].xp;
          const float r2 = ngb_buffer[k].r2;

          /* Do the ionization */
          feedback_iact_HII_ionization(si, pj, xpj, r2, pixel, phys_const,
                                       hydro_props, us, cosmo, cooling,
                                       feedback_props, ti_begin, time);

        } /* Loop over the sorted particles */
      }

      /* Fully exhausted: nothing more to do this pass. */
      if (feedback_get_star_ionization_rate_max(si) <= 0.0) break;

      if (buffer_was_full) {
        /* Case 1: more candidates may exist within the SAME radius. */
        if (num_retry_full_buffer >= max_retry_full_buffer) break;
        ++num_retry_full_buffer;
      } else {
        /* Case 2: genuinely nothing left within the current radius, but
         * photons remain -- expand the radius, bounded by the configured
         * ceiling. If we're already at that ceiling, expanding again
         * would just re-search the identical volume, so stop. */
        if (num_radius_expansions >= max_radius_expansion_tries) break;
        if (dynamic_search_radius >= max_search_radius) break;
        dynamic_search_radius = min(
            dynamic_search_radius * radius_expansion_factor, max_search_radius);
        ++num_radius_expansions;
      }
    }

    c->stars.h_hii_max = max(c->stars.h_hii_max, si->h_hii);
    c->stars.h_hii_max_active = max(c->stars.h_hii_max_active, si->h_hii);

    /*****************************************/
    /* Update the star after HII ionization */

    /* TODO: Move into a function */
    if (feedback_is_HII_ionization_active(si, e)) {
      /* Compute the times */
      double star_age_beg_step = 0;
      double dt_enrichment = 0;
      integertime_t ti_begin = 0;
      compute_time(si, with_cosmology, cosmo, &star_age_beg_step,
                   &dt_enrichment, &ti_begin, ti_current, time_base, time);

      /* Log when this HII region was (re)built */
      si->feedback_data.radiation.HII_region_last_rebuild = star_age_beg_step;
    }
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
    if (feedback_get_star_ionization_rate_max(si) <= 0.0) {
      message(
          "Star %lld has exhausted all its ionizing photons! r_hii = %e, "
          "num_retry_full_buffer = %d, num_radius_expansions = %d",
          si->id, si->h_hii * kernel_gamma, num_retry_full_buffer,
          num_radius_expansions);
    } else {
      message(
          "Star %lld has NOT exhausted all its ionizing photons! Remaining: "
          "%e, final_search_radius = %e, sp->h_hii = %e, "
          "num_retry_full_buffer = %d, num_radius_expansions = %d",
          si->id, feedback_get_star_ionization_rate_max(si),
          dynamic_search_radius, kernel_gamma * si->h_hii,
          num_retry_full_buffer, num_radius_expansions);
    }
#endif
  } /* Loop over sparts */
}

/**
 * @brief Downward recursion branch handler for pair cell interactions.
 *
 * This function evaluates a target cell pair configuration. It applies a
 * bounding-box optimization check along the sorting axis to determine if the
 * cell block is within reach. If the target cell is split, it continues down
 * the sub-tree hierarchy; otherwise, it resolves the leaf interaction using
 * either a sorted array approach or a naive fallback mechanism.
 *
 * @param r The #runner thread.
 * @param ci The local cell context containing the active source star particle.
 * @param cj The target neighbor cell containing gas particle candidates.
 * @param sid The sorting axis index between the pair.
 * @param flipped Flag indicating if the cell order is flipped relative to the
 * axis direction.
 * @param shift The periodic wrapping shift vector applied to the coordinates.
 * @param si The #spart (star) performing the feedback.
 * @param search_radius The maximum search distance around the star.
 * @param ngb_buffer The tracking array where valid gas neighbors are collected.
 * @param max_size The maximum capacity of the neighbor tracking array.
 * @param count_found (return) Total tracking count of validated gas neighbors
 * gathered.
 */
void runner_do_stars_hii_ionization_feedback_branch(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *ngb_buffer, int max_size,
    int *count_found) {

  /* OPTIMIZATION: Check if cj is completely out of reach before doing anything
   * else */
  if (!runner_hii_check_cell_can_be_reached(ci, cj, sid, flipped, shift, si,
                                            search_radius)) {
    return;
  }

  /* Don't reprocess cj or its child if it's ci. These particles were
     processed in the self interaction */
  if (ci == cj) return;

  if (cj->split) {

    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_branch(
            r, ci, cj->progeny[k], sid, flipped, shift, si, search_radius,
            ngb_buffer, max_size, count_found);
      }
    }
  } else {

    /* Let's first lock the cell */
    lock_lock(&cj->hydro.extra_sort_lock);

    /* The sorted-window prune is only a pre-filter: it projects gas and star
     * onto runner_shift[use_sid] (unit vector) and keeps particles within the
     * [di_min, di_max] slab; the exact r2 < search_radius^2 test inside the
     * dopair then decides ionization. Since the projection difference never
     * exceeds the 3D distance, the slab is conservative for ANY of the 13
     * sids, so the ionized set is invariant to which sid we use -- only speed
     * changes. We therefore prefer the carried coarse sid, but if cj is not
     * sorted along it (the common case for a progeny leaf, whose own pair
     * tasks use different sub-sids) we fall back to ANY sid cj is sorted
     * along rather than to the O(N) naive search. The freshness gate
     * (dx_max_sort_old) is per-cell, hence valid for any sid. */
    int use_sid = -1;
    int use_flipped = flipped;
    if (cj->hydro.dx_max_sort_old <= space_maxreldx * cj->dmin) {
      if (cj->hydro.sorted & (1 << sid)) {
        use_sid = sid; /* keep flipped: preserves near-end early-start */
      } else {
        for (int s = 0; s < 13; s++) {
          if (cj->hydro.sorted & (1 << s)) {
            use_sid = s;
            use_flipped = 0; /* forward scan is valid for any sorted sid */
            break;
          }
        }
      }
    }

    /* Unlock now that we have read the sort state. */
    if (lock_unlock(&cj->hydro.extra_sort_lock) != 0)
      error("Impossible to unlock cell!");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    use_sid = -1; /* force the naive path */
#endif
    if (use_sid < 0) {
      runner_dopair_naive_stars_hii_ionization_feedback(
          r, ci, cj, shift, si, search_radius, ngb_buffer, max_size,
          count_found);
    } else {
      runner_dopair_stars_hii_ionization_feedback(
          r, ci, cj, use_sid, use_flipped, shift, si, search_radius, ngb_buffer,
          max_size, count_found);
    }
  }
}

/**
 * @brief Gather gas particles within the HII radius from the star's own cell.
 *
 * Performs a brute-force search over all hydro particles in cell @p c that
 * fall within the HII smoothing length of the star @p si. Candidates are
 * added to the @p buffer if they are not already inhibited or ionized.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param search_radius The distance of potential gas candidates.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_doself_stars_hii_ionization_feedback(
    struct runner *r, struct cell *c, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found) {

  struct engine *e = r->e;
  const int count = c->hydro.count;

  /* Anything to do here?*/
  if (count == 0) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e))
    error("Interacting undrifted cell (hydro).");
  if (!cell_are_spart_drifted(c, e))
    error("Interacting undrifted cell (stars).");

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const float r2_max = search_radius * search_radius;
  const float six[3] = {si->x[0], si->x[1], si->x[2]};

  /* TODO: Use sorted cells to check that the particles/the cell are/is
     within the star's h_hii (or h if h_hii == 0) */

  /* Loop over the parts in c. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = &xparts[pjd];

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (!feedback_part_can_be_ionized(pj, xpj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that particles have been drifted to the current time */
    if (pj->ti_drift != e->ti_current)
      error(
          "Particle pj (%lld) not drifted to current time. c = %lld, "
          "c->super = %lld",
          pj->id, c->cellID, c->super->cellID);
#endif

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0]), (float)(pj->x[1]),
                          (float)(pj->x[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    /* Skipping the insert once a pixel is exhausted stops a dead pixel's
       candidates from crowding out other pixels' on retry. */
    const int pixel =
        runner_hii_get_pixel(dx, si->feedback_data.radiation.n_HII_pixels);
    if (r2 < r2_max && feedback_get_star_ionization_rate(si, pixel) > 0.0)
      runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, c,
                               pixel);
  } /* Loop in current cell */
}

/**
 * @brief Gather gas particles within the HII radius from a neighboring cell
 * using a naive O(N_gas) search.
 *
 * Similar to the 'self' version, but handles the distance calculations between
 * two different cells (@p ci and @p cj), including periodic boundary
 * conditions/wrapping.
 *
 * @param r The #runner thread.
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param search_radius The distance of potential gas candidates.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_dopair_naive_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const double shift[3],
    struct spart *si, const float search_radius, struct hii_neighbor *buffer,
    int max_size, int *count_found) {

  struct engine *e = r->e;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  if (cj != ci && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_dopair_naive_stars_hii_ionization_feedback(
            r, ci, cj->progeny[k], shift, si, search_radius, buffer, max_size,
            count_found);
      }
    }
  } else {

    /* Check that cells are drifted. */
    if (!cell_are_part_drifted(cj, e))
      error("Interacting undrifted cell (hydro).");
    if (!cell_are_spart_drifted(cj, e))
      error("Interacting undrifted cell (stars).");

    struct part *restrict parts_j = cj->hydro.parts;
    struct xpart *restrict xparts_j = cj->hydro.xparts;
    const float r2_max = search_radius * search_radius;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;
      if (!feedback_part_can_be_ionized(pj, xpj, e)) continue;

      /* message("[pair] Found %lld!", pj->id); */

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error(
            "Particle pj (%lld) not drifted to current time. c = %lld, "
            "c->super = %lld",
            pj->id, cj->cellID, cj->super->cellID);
#endif

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int pixel =
          runner_hii_get_pixel(dx, si->feedback_data.radiation.n_HII_pixels);
      if (r2 < r2_max && feedback_get_star_ionization_rate(si, pixel) > 0.0)
        runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj,
                                 pixel);
    } /* Loop in current cell */
  }
}

/**
 * @brief Gather gas particles within the HII radius from a neighboring cell
 * using the sorted cells to speed up the search.
 *
 * Projects the source star position along the target sorting axis vector @p
 * sid and analyzes the pre-sorted gas array layout of cell @p cj. It
 * implements directional parsing boundaries to read only the specific window
 * of gas particles capable of interacting, skipping distant arrays.
 *
 *
 * @param r The #runner thread.
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param search_radius The distance of potential gas candidates.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_dopair_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found) {

  struct engine *e = r->e;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(cj, e))
    error("Interacting undrifted cell (hydro).");
  if (!cell_are_spart_drifted(cj, e))
    error("Interacting undrifted cell (stars).");

  /* Get the hydro sorts for our gas cell */
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
  const float dx_max = cj->hydro.dx_max_sort + ci->stars.dx_max_sort;

  /* Project the star onto the sorting axis.
   * The gas particle projections (sort_j[].d) and the r2 distances below both
   * use ABSOLUTE gas positions (pj->x), so the star must be expressed in that
   * same absolute frame: the periodic shift brings cj adjacent to ci, hence a
   * gas particle at pj->x sits at effective position pj->x + shift, and the
   * correct minimum-image separation is (si->x - shift) - pj->x. The star is
   * therefore ALWAYS shifted by -shift, independent of `flipped` (which only
   * selects the sorted-scan direction below). */
  const double pix[3] = {si->x[0] - shift[0], si->x[1] - shift[1],
                         si->x[2] - shift[2]};

  double di_star = 0.0;
  for (int k = 0; k < 3; k++) {
    di_star += pix[k] * runner_shift[sid][k];
  }

  /* Compute bounds along the sorting axis */
  const double di_min = di_star - search_radius - dx_max;
  const double di_max = di_star + search_radius + dx_max;

  const double dj_min = sort_j[0].d;
  const double dj_max = sort_j[count_j - 1].d;

  /* Skip if the entire cell is completely out of reach along the sorted axis
   */
  if (di_max < dj_min || di_min > dj_max) {
    return;
  }

  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const float r2_max = search_radius * search_radius;
  const float six[3] = {(float)pix[0], (float)pix[1], (float)pix[2]};

  if (!flipped) {
    /* Star is on the 'left' relative to the axis direction.
     * We scan forward through the sorted list until particles are beyond
     * di_max. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d <= di_max; pjd++) {

      /* Skip particles that haven't reached the interaction zone yet */
      if (sort_j[pjd].d < di_min) continue;

      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj = &xparts_j[sort_j[pjd].i];

      if (part_is_inhibited(pj, e)) continue;
      if (!feedback_part_can_be_ionized(pj, xpj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != e->ti_current)
        error("Particle pj (%lld) not drifted to current time.", pj->id);
#endif

      const float pjx[3] = {(float)(pj->x[0]), (float)(pj->x[1]),
                            (float)(pj->x[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int pixel =
          runner_hii_get_pixel(dx, si->feedback_data.radiation.n_HII_pixels);
      if (r2 < r2_max && feedback_get_star_ionization_rate(si, pixel) > 0.0)
        runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj,
                                 pixel);
    }
  } else {
    /* Star is on the 'right' relative to the axis direction.
     * We scan backward through the sorted list until particles are below
     * di_min. */
    for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d >= di_min; pjd--) {

      /* Skip particles that are too far on the other side of the slab */
      if (sort_j[pjd].d > di_max) continue;

      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj = &xparts_j[sort_j[pjd].i];

      if (part_is_inhibited(pj, e)) continue;
      if (!feedback_part_can_be_ionized(pj, xpj, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (pj->ti_drift != e->ti_current)
        error("Particle pj (%lld) not drifted to current time.", pj->id);
#endif

      const float pjx[3] = {(float)(pj->x[0]), (float)(pj->x[1]),
                            (float)(pj->x[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int pixel =
          runner_hii_get_pixel(dx, si->feedback_data.radiation.n_HII_pixels);
      if (r2 < r2_max && feedback_get_star_ionization_rate(si, pixel) > 0.0)
        runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj,
                                 pixel);
    }
  }
}

/**
 * @brief Check if a star's interaction radius can reach any part of a
 * neighboring cell along a specific axis.
 *
 * @param ci The cell containing the star particle.
 * @param cj The neighboring gas cell to check.
 * @param sid The sorting axis index between the pair.
 * @param flipped Flag indicating if the cell order is flipped relative to the
 * axis direction.
 * @param shift The periodic wrapping shift vector.
 * @param si The star particle pointer.
 * @param search_radius The maximum HII ionization radius of the star.
 * @return 1 if the cell is within reach (or if we cannot determine it safely),
 * 0 if it can be completely skipped.
 */
__attribute__((always_inline)) INLINE int runner_hii_check_cell_can_be_reached(
    const struct cell *ci, const struct cell *cj, const int sid,
    const int flipped, const double shift[3], const struct spart *si,
    const float search_radius) {

  const float dx_max = cj->hydro.dx_max_sort + ci->stars.dx_max_sort;

  /* Project the star onto the sorting axis. The cell extent below (dj_min/
   * dj_max) is built from cj's ABSOLUTE geometry (cj->loc), so the star must
   * be expressed in that same absolute frame: shift brings cj adjacent to ci,
   * hence the star is ALWAYS shifted by -shift, independent of `flipped`.
   * (Must match the identical convention in
   * runner_dopair_stars_hii_ionization_feedback; a flipped-dependent sign here
   * offsets the reach window by 2*shift and wrongly skips wrapped neighbours.)
   */
  const double pix[3] = {si->x[0] - shift[0], si->x[1] - shift[1],
                         si->x[2] - shift[2]};

  double di_star = 0.0;
  for (int k = 0; k < 3; k++) {
    di_star += pix[k] * runner_shift[sid][k];
  }

  /* Star reach window along the axis direction */
  const double di_min = di_star - search_radius - dx_max;
  const double di_max = di_star + search_radius + dx_max;

  /* Calculate cell boundaries directly from geometry to remain safe if unsorted
   */
  double cj_center[3] = {cj->loc[0] + cj->width[0] * 0.5,
                         cj->loc[1] + cj->width[1] * 0.5,
                         cj->loc[2] + cj->width[2] * 0.5};

  double dj_center = 0.0;
  double dj_extent = 0.0;
  for (int k = 0; k < 3; k++) {
    dj_center += cj_center[k] * runner_shift[sid][k];
    dj_extent += cj->width[k] * 0.5 * fabs(runner_shift[sid][k]);
  }

  const double dj_min = dj_center - dj_extent;
  const double dj_max = dj_center + dj_extent;

  /* Skip if completely out of range */
  if (di_max < dj_min || di_min > dj_max) {
    return 0;
  }

  return 1;
}
