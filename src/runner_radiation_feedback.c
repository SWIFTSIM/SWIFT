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
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "hydro.h"
#include "stars.h"
#include "timers.h"

#define search_radius_factor 1.2f
#define max_ngbs 128

/**
 * @brief Top-level function for HII ionization feedback.
 *
 * This function recurses down the cell hierarchy until it reaches a level
 * where the cell size is comparable to the maximum HII radius of the stars
 * within it. At that point, it calls the 'branch' function to handle the
 * actual neighbor gathering.
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
  const float r_hii_max = max(c->stars.h_hii_max_old, c->stars.h_max_old) * kernel_gamma;
  const float max_search_radius = star_props->max_HII_search_radius;
  const float interaction_limit =
      min(search_radius_factor * r_hii_max, max_search_radius);
  const int can_recurse = c->split && (interaction_limit < 0.5f * c->dmin);

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

  /* for (struct link *l = c->stars.radiation_in; l != NULL; l = l->next) { */
  /*   /\* We have already handled the self case *\/ */
  /*   if (l->t->type == task_type_self) continue; */

  /*   struct cell *cj = l->t->cj; */
  /*   struct cell *ci = l->t->ci; */
  /*   message( */
  /*       "[%lld, %lld] hydro super: %lld , %lld | radiation_level: %lld %lld",
   */
  /*       ci->cellID, cj->cellID, ci->hydro.super->cellID, */
  /*       cj->hydro.super->cellID, ci->stars.radiation_level->cellID, */
  /*       cj->stars.radiation_level->cellID); */
  /* } */

  TIMER_TIC;

  /*
   * Decision: Should we process this cell or go deeper?
   * We stop recursing if:
   * 1. The cell is a leaf (!c->split).
   * 2. The cell is small enough that we can't efficiently split the
   *    interaction radius further (similar to Mosaics logic).
   * 3. The R_HII_max covers a significant fraction of the cell.
   */
  /* Is the cell split and not smaller than the smoothing length? */
  if (can_recurse) {
    /* Keep recursing deeper into the super-cell hierarchy. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];
        runner_do_stars_hii_ionization_feedback(r, cp, 0);
        c->stars.h_hii_max = max(c->stars.h_hii_max, cp->stars.h_hii_max);
        c->stars.h_max_active =
            max(c->stars.h_max_active, cp->stars.h_max_active);
      }
    }
  } else {
    /* We have reached the 'Working Level' */
    runner_do_stars_hii_ionization_feedback_branch(r, c, interaction_limit);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Branching point to coordinate HII ionization for stars in a cell.
 *
 * This function manages the neighbor gathering for all active stars in the
 * cell. It allocates a temporary buffer, gathers neighbors from the cell itself
 * (self) and its neighbors (pair) by climbing the cell hierarchy, and
 * prepares the gathered particles for sorting and ionization.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the stars to process.
 */
void runner_do_stars_hii_ionization_feedback_branch(
    struct runner *r, struct cell *c, const float interaction_limit) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  /* const struct feedback_props *fb_props = e->feedback_props; */
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;

  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  /* Anything to do here? */
  if (c->stars.count == 0 || c->hydro.count == 0 || !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* TODO: Reviser ces checks */
  /* Did we mess up the recursion? */
  if (interaction_limit > c->dmin && c != c->hydro.super)
    error("Cell (%lld) size (%e) smaller than HII interaction length (%e)",
          c->cellID, c->dmin, interaction_limit);
#endif

  /* TODO: Add multiple tries if the star has not exhausted its photons */
  struct hii_neighbor ngb_buffer[max_ngbs];

  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[sid];

    /* Is this part within the timestep? */
    if (spart_is_inhibited(si, e)) continue;
    if (!spart_is_active(si, e)) continue;
    if (!feedback_is_active(si, e)) continue; /* TODO: Do we want this? */
    if (!feedback_is_HII_ionization_active(si, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that particles have been drifted to the current time */
    if (si->ti_drift != e->ti_current)
      error(
          "Particle si (%lld) not drifted to current time c = %lld, "
          "c->super = %lld",
          si->id, c->cellID, c->super->cellID);
#endif

    /* Logic: If the buffer was full, there might be more neighbors just
     * outside the current R_max of the buffer. We retry until the buffer
     * is no longer maxed out OR the star runs out of photons. */
    int buffer_was_full = 0;
    do {
      int count_found = 0;
      buffer_was_full = 0;

      /***************************************************/
      /* First loop over particles in the current cell */
      runner_do_stars_hii_ionization_feedback_self(
          r, c->hydro.super, si, interaction_limit, ngb_buffer, max_ngbs,
          &count_found);

      /* Now loop over particles in the neighboring cells */
      for (struct link *l = c->hydro.super->stars.radiation_in; l != NULL;
           l = l->next) {
        /* We have already handled the self case */
        if (l->t->type == task_type_self) continue;
        struct cell *c_in = (l->t->cj == c->hydro.super) ? l->t->ci : l->t->cj;
        runner_do_stars_hii_ionization_feedback_pair(
            r, c, c_in, si, interaction_limit, ngb_buffer, max_ngbs,
            &count_found);
      } /* Neighbour search */

      /* Flag if we hit the limit */
      if (count_found == max_ngbs) buffer_was_full = 1;

      /***************************************************/
      /* It's time to sort the gas particles */
      if (count_found > 0) {
#ifdef SWIFT_DEBUG_CHECKS
        /* Verify that the neighbors are properly sorted by distance */
        for (int k = 0; k < count_found - 1; k++) {
          if (ngb_buffer[k].r2 > ngb_buffer[k + 1].r2) {
            error(
                "HII neighbor buffer not properly sorted! "
                "Index %d (r2=%e) is larger than index %d (r2=%e).",
                k, ngb_buffer[k].r2, k + 1, ngb_buffer[k + 1].r2);
          }
        }
#endif
        const integertime_t ti_begin =
            get_integer_time_begin(e->ti_current - 1, si->time_bin);

        /* Now let's ionize the gas particles! */
        for (int k = 0; k < count_found; k++) {

          /* No more photons to consume */
          if (feedback_get_star_ionization_rate(si) <= 0.0) {
            message("Star has exhausted all its ionizing photons!");
            break;
          }

          struct part *pj = ngb_buffer[k].p;
          struct xpart *xpj = ngb_buffer[k].xp;
          const float r2 = ngb_buffer[k].r2;

          /* Do the ionization */
          feedback_do_HII_ionization(si, pj, xpj, r2, phys_const, hydro_props,
                                     us, cosmo, cooling, ti_begin);

        } /* Loop over the sorted particles */
      }
      /* If the star still has photons and we previously filled the buffer,
       * we must go again to find the neighbors that were 'bumped out'. */
    } while (buffer_was_full && feedback_get_star_ionization_rate(si) > 0);

    c->stars.h_hii_max = max(c->stars.h_hii_max, si->h_hii);
    c->stars.h_max_active = max(c->stars.h_max_active, si->h_hii);
  } /* Loop over sparts */
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
void runner_do_stars_hii_ionization_feedback_self(
    struct runner *r, struct cell *c, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found) {

  struct engine *e = r->e;
  const int count = c->hydro.count;

  /* Anything to do here?*/
  if (count == 0) return;

  /* If cell is leaf, do the brute force search */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_self(
            r, c->progeny[k], si, search_radius, buffer, max_size, count_found);
      }
    }
  } else {

    /* Check that cells are drifted. */
    if (!cell_are_part_drifted(c, e))
      error("Interacting undrifted cell (hydro).");
    if (!cell_are_spart_drifted(c, e))
      error("Interacting undrifted cell (stars).");
    
    struct part *restrict parts = c->hydro.parts;
    struct xpart *restrict xparts = c->hydro.xparts;
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
      if (feedback_part_can_be_ionized(pj, xpj, e)) continue;

      /* message("[self] Found %lld!", pj->id); */

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

      runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, c);
    } /* Loop in current cell */
  }
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
void runner_do_stars_hii_ionization_feedback_pair_naive(
    struct runner *r, struct cell *ci, struct cell *cj, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_pair_naive(r, ci, cj->progeny[k], si,
                                                     search_radius, buffer,
                                                     max_size, count_found);
      }
    }
  } else {
    /* Get the relative distance between the pairs, wrapping. */
    double shift[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; k++) {
      if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
        shift[k] = e->s->dim[k];
      else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
        shift[k] = -e->s->dim[k];
    }

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

      runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj);
    } /* Loop in current cell */
  }
}

/**
 * @brief Gather gas particles within the HII radius from a neighboring cell.
 * 
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
void runner_do_stars_hii_ionization_feedback_pair(
    struct runner *r, struct cell *ci, struct cell *cj, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  if (cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_pair(r, ci, cj->progeny[k], si,
                                                     search_radius, buffer,
                                                     max_size, count_found);
      }
    }
  } else {
    /* Get the relative distance between the pairs, wrapping. */
    double shift[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; k++) {
      if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
        shift[k] = e->s->dim[k];
      else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
        shift[k] = -e->s->dim[k];
    }

    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;
      if (feedback_part_can_be_ionized(pj, xpj, e)) continue;

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

      runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj);
    } /* Loop in current cell */
  }
}

/**
 * @brief Maintain a sorted buffer by inserting a new neighbor at the correct
 * position.
 *
 * If the buffer is full, it replaces the furthest element if the new one is
 * closer.
 */
__attribute__((always_inline)) INLINE void runner_hii_buffer_insert(
    struct hii_neighbor *buffer, int max_size, int *count_found, float r2,
    struct part *p, struct xpart *xp, struct cell *c) {

  /* Case A: Buffer is not yet full */
  if (*count_found < max_size) {
    int i = *count_found - 1;
    /* Shift elements to make room (standard insertion) */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
    (*count_found)++;
  }
  /* Case B: Buffer is full, check if new particle is closer than the furthest
   */
  else if (r2 < buffer[max_size - 1].r2) {
    int i = max_size - 2;
    /* Shift elements to replace the furthest */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
  }
}
