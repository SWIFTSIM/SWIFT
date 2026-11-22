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

/* Standard headers. */
#include <stdlib.h>

/* This object's header. */
#include "runner_radiation_feedback.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "clocks.h"
#include "cycle.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "memuse.h"
#include "runner.h"
#include "timers.h"

#ifdef IONIZATION_FEEDBACK_LOOP
/**
 * @brief Grow a runner's HII maintenance scratch buffer to at least
 * @p needed entries.
 *
 * Doubling grow-on-demand, like gravity_cache_init's pattern
 * (runner_doiact_grav.c): never shrinks. Existing content is not
 * preserved; callers refill the buffer from scratch every maintenance pass.
 *
 * @param r The #runner whose buffer to grow.
 * @param needed The minimum required capacity, in entries.
 */
static INLINE void runner_hii_maintenance_buffer_ensure(struct runner *r,
                                                        int needed) {
  if (r->hii_maintenance_buffer_size >= needed) return;

  int new_size = 2 * r->hii_maintenance_buffer_size;
  if (new_size < needed) new_size = needed;
  if (new_size < max_ngbs) new_size = max_ngbs;

  struct hii_neighbor *new_buf = (struct hii_neighbor *)swift_realloc(
      "hii_maintenance_buffer", r->hii_maintenance_buffer,
      (size_t)new_size * sizeof(struct hii_neighbor));
  if (new_buf == NULL)
    error("Failed to grow the HII maintenance scratch buffer to %d entries.",
          new_size);

  r->hii_maintenance_buffer = new_buf;
  r->hii_maintenance_buffer_size = new_size;
}

#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
/**
 * @brief qsort() comparator for ascending (r2, id) order.
 *
 * id is the tiebreak: r2 alone is not a total order, since distinct
 * particles can sit at the identical squared distance. Explicit branches
 * rather than subtraction, since `a - b` on r2's magnitudes or on a
 * `long long` id either loses precision or overflows.
 *
 * Guarded like its only call site (the maintenance-buffer sort in
 * runner_dosub_stars_hii_ionization_feedback()): under NO_COOLING that call
 * is compiled out, and an unused static function trips -Werror
 * -Wunused-function.
 */
static int runner_hii_maintenance_compare(const void *a, const void *b) {
  const struct hii_neighbor *na = (const struct hii_neighbor *)a;
  const struct hii_neighbor *nb = (const struct hii_neighbor *)b;
  if (na->r2 < nb->r2) return -1;
  if (na->r2 > nb->r2) return 1;
  if (na->p->id < nb->p->id) return -1;
  if (na->p->id > nb->p->id) return 1;
  return 0;
}
#endif

/**
 * @brief The cone frame rotation of a star, or NULL when cones follow the box.
 *
 * @param si The #spart.
 * @param feedback_props The #feedback_props.
 * @param pass_index Number of passes the star has made before this one.
 * @param storage Memory for the matrix, used only when a rotation is returned.
 * @return @a storage filled with the rotation, or NULL.
 */
__attribute__((always_inline)) INLINE static const double *
runner_hii_star_cone_frame(const struct spart *si,
                           const struct feedback_props *feedback_props,
                           const int64_t pass_index, double storage[9]) {
  if (!feedback_get_HII_random_cone_frame(feedback_props)) return NULL;
  const int64_t epoch = runner_hii_cone_frame_epoch(
      pass_index, feedback_get_HII_cone_frame_redraw_period(feedback_props));
  runner_hii_cone_frame_from_id(si->id, epoch, storage);
  return storage;
}

/**
 * @brief Assign a gas candidate to its angular pixel and insert it into
 * the neighbour buffer if it is within range and that pixel still has
 * budget left.
 *
 * Shared by every gather call site (doself, dopair naive, dopair sorted
 * forward/backward), which otherwise each repeat this exact pixel +
 * range + budget + insert sequence.
 *
 * @param gather The #hii_gather_state of the traversal.
 * @param pj The candidate gas particle.
 * @param xpj The candidate's extended data.
 * @param c The #cell the candidate belongs to (for SWIFT_DEBUG_CHECKS).
 * @param dx Star-minus-particle position offset.
 * @param r2 Squared distance between star and candidate.
 * @param r2_max Squared search radius.
 */
__attribute__((always_inline)) INLINE static void
runner_hii_try_insert_candidate(struct hii_gather_state *gather,
                                struct part *pj, struct xpart *xpj,
                                struct cell *c, const float dx[3],
                                const float r2, const float r2_max) {

  if (r2 >= r2_max) return;

  struct spart *si = gather->si;
  struct hii_maintenance_context *ctx = gather->ctx;

#ifdef SWIFT_DEBUG_CHECKS
  /* Already-ionized gas is processed in place and passes no buffer; new
     candidates are buffered and pass no context. */
  if (gather->mode == hii_gather_already_ionized && ctx == NULL)
    error("Gathering already-ionized gas without a maintenance context.");
  if (gather->mode == hii_gather_new_candidates &&
      (gather->buffer == NULL || gather->max_size <= 0))
    error("Gathering new candidates without a usable neighbour buffer.");
#endif

  const int pixel = runner_hii_get_pixel(
      dx, feedback_get_star_HII_pixel_count(si), gather->frame);

  if (gather->mode == hii_gather_already_ionized) {
    /* Extent and count track what this star still HOLDS, not what it could
       afford this pass, and are updated live regardless of traversal order:
       gating them on the charge would shrink h_hii and drop the region's
       outskirts from the next pass's search radius. Only the charge itself
       (which shell lapses under a shortfall) is deferred to the caller's
       (r2, id)-ordered pass below. */
    ctx->r2_max_maintained = max(ctx->r2_max_maintained, r2);
    gather->count_found++;

    /* Buffer this candidate for (r2, id)-ordered charging once the whole
       traversal returns. The region has no reason to fit in max_ngbs, so
       the buffer grows to match rather than dropping outskirts. */
    runner_hii_maintenance_buffer_ensure(ctx->r, ctx->buffer_count + 1);
    struct hii_neighbor *slot =
        &ctx->r->hii_maintenance_buffer[ctx->buffer_count++];
    slot->r2 = r2;
    slot->p = pj;
    slot->xp = xpj;
    slot->pixel = pixel;
#ifdef SWIFT_DEBUG_CHECKS
    slot->c = c;
#endif
    return;
  }

  /* Skipping the insert once a pixel is exhausted stops a dead pixel's
     candidates from crowding out other pixels' on retry. */
  if (feedback_get_star_ionization_budget(si, pixel) > 0.0)
    runner_hii_buffer_insert(gather->buffer, gather->max_size,
                             &gather->count_found, r2, pj, xpj, c, pixel);
}

/**
 * @brief Does this gas particle belong to the set the current traversal is
 * gathering?
 *
 * @param pj The candidate gas particle.
 * @param xpj The candidate's extended data.
 * @param e The #engine.
 * @param star_id Id of the #spart performing the feedback.
 * @param mode Which set the traversal is after.
 */
__attribute__((always_inline)) INLINE static int
runner_hii_part_passes_gather_filter(const struct part *pj,
                                     const struct xpart *xpj,
                                     const struct engine *e,
                                     const long long star_id,
                                     enum hii_gather_mode mode) {

  if (mode == hii_gather_new_candidates)
    return feedback_part_can_be_ionized(pj, xpj, e);

  /* Only gas this star itself ionized: a particle inside two overlapping
     HII regions is owned, and paid for, by whichever star claimed it first. */
  return radiation_is_part_tagged_as_ionized(pj, xpj) &&
         radiation_get_part_ionized_star_id(pj, xpj) == star_id;
}

/**
 * @brief Window along a sorting axis that a star can reach.
 *
 * The star is shifted by -shift whatever the scan direction: shift brings the
 * neighbour cell next to the star's cell, so a gas particle at its absolute
 * position sits at pj->x + shift. The cell-level prune and the sorted scan
 * MUST share this window, and so share this function.
 *
 * @param si The #spart.
 * @param shift The periodic wrapping shift vector applied to the neighbour.
 * @param sid The sorting axis index between the pair.
 * @param search_radius The maximum search distance around the star.
 * @param dx_max Bound on how far the particles moved since the sort.
 * @param pix (return) The star position in the neighbour's absolute frame.
 * @param di_min (return) Lower end of the window along the axis.
 * @param di_max (return) Upper end of the window along the axis.
 */
__attribute__((always_inline)) INLINE static void runner_hii_axis_window(
    const struct spart *si, const double shift[3], const int sid,
    const float search_radius, const float dx_max, double pix[3],
    double *di_min, double *di_max) {

  for (int k = 0; k < 3; k++) pix[k] = si->x[k] - shift[k];

  double di_star = 0.0;
  for (int k = 0; k < 3; k++) {
    di_star += pix[k] * runner_shift[sid][k];
  }

  *di_min = di_star - search_radius - dx_max;
  *di_max = di_star + search_radius + dx_max;
}

/**
 * @brief Offer one gas particle of a cell to the gather.
 *
 * The particle is skipped if it is inhibited or outside the set the traversal
 * is after. Otherwise its distance to the star is measured and it goes to
 * #runner_hii_try_insert_candidate.
 *
 * @param gather The #hii_gather_state of the traversal.
 * @param pj The candidate gas particle.
 * @param xpj The candidate's extended data, or NULL for a foreign cell, which
 * has no local xparts.
 * @param c The #cell holding the candidate.
 * @param six The star position, in the frame of the candidate position
 * minus @p origin.
 * @param origin Position subtracted from the candidate's, or NULL for none.
 * @param r2_max Squared search radius.
 */
__attribute__((always_inline)) INLINE static void runner_hii_offer_candidate(
    struct hii_gather_state *gather, struct part *pj, struct xpart *xpj,
    struct cell *c, const float six[3], const double origin[3],
    const float r2_max) {

  const struct engine *e = gather->r->e;

  /* Early abort? */
  if (part_is_inhibited(pj, e)) return;
  if (!runner_hii_part_passes_gather_filter(pj, xpj, e, gather->si->id,
                                            gather->mode))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that particles have been drifted to the current time */
  if (pj->ti_drift != e->ti_current)
    error("Particle pj (%lld) not drifted to current time. c = %lld.", pj->id,
          c->cellID);
#endif

  /* Compute the pairwise distance. */
  const float pjx[3] = {
      (float)(origin == NULL ? pj->x[0] : pj->x[0] - origin[0]),
      (float)(origin == NULL ? pj->x[1] : pj->x[1] - origin[1]),
      (float)(origin == NULL ? pj->x[2] : pj->x[2] - origin[2])};
  const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  runner_hii_try_insert_candidate(gather, pj, xpj, c, dx, r2, r2_max);
}

/**
 * @brief Start-up of the subgrid radiation: print warnings about the HII
 * search-radius / box-size relationship.
 *
 * @param e The #engine.
 */
void feedback_radiation_startup_diagnostics(const struct engine *e) {

  /* Both warnings concern the HII search only */
  if (!(e->feedback_props->radiation_policy & radiation_policy_photoionization))
    return;

  /* Stars:HII_max_search_radius exceeding the periodic box's half-width is
     not itself fatal; the search radius just never gets to expand that
     far, which is indistinguishable from a genuinely converged result
     without an external comparison. Surface it at startup instead. */
  if (e->s->periodic && e->nodeID == 0) {
    const double half_width =
        0.5 * min3(e->s->dim[0], e->s->dim[1], e->s->dim[2]);
    if (e->stars_properties->HII_max_search_radius > half_width)
      warning(
          "Stars:HII_max_search_radius (%g) exceeds the periodic box's "
          "half-width (%g): the HII search radius will silently never "
          "reach the configured value, capping out at the box's own "
          "periodic wrap distance instead.",
          e->stars_properties->HII_max_search_radius, half_width);
  }

  /* Once a star's HII reach outgrows the top-level cell width there is no
     coarser grid to rebuild into, and only the stencil-completeness
     fallback keeps the rebuild criterion satisfiable. Configurations that
     can never reach it abort mid-run: on an unsatisfiable rebuild for a
     non-periodic box, or on space_regrid()'s "not enough cells" error for
     a periodic one. */
  if (e->nodeID == 0 && !space_radiation_top_stencil_can_cover_box(e->s))
    warning(
        "This box cannot coarsen to 3 top-level cells along every axis "
        "(non-cubic box, non-periodic, or Scheduler:min_top_level_cells "
        "> 3). Subgrid radiation then has no fallback if a star's HII "
        "reach outgrows the top-level cell width, and the run aborts: on "
        "an unsatisfiable rebuild if the box is non-periodic, or on "
        "space_regrid()'s 'not enough cells' error if it is periodic. "
        "Keep h_hii well below the top-level cell width, e.g. via "
        "Stars:HII_max_search_radius.");
}
#endif

/**
 * @brief Top-level steering function for HII ionization feedback.
 *
 * This function recurses into every split cell that holds an active star down
 * to a leaf, whatever the HII radius of the stars it holds. The leaf is the
 * working level: it hands off control to the sub-task coordinator, which
 * searches out to the radius computed here.
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

  /* Anything to do here? c->hydro.count == 0 is NOT included: a due star in
     a gas-free working-level cell must still reach
     runner_dosub_stars_hii_ionization_feedback() to have its
     attempt and its issued photons stamped (see there); only the recursion
     into progeny and the search itself may be skipped for lack of gas. */
  if (c->stars.count == 0 || !cell_is_active_stars(c, e)) return;

  const float r_hii_max =
      max(c->stars.h_hii_max_active, c->stars.h_max_active) * kernel_gamma;
  const float max_search_radius = star_props->HII_max_search_radius;
  const float interaction_limit =
      min(radiation_search_radius_factor * r_hii_max, max_search_radius);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->cellID != c->stars.radiation_level->cellID && timer == 1) {
    warning(
        "Not running the HII reionization task on the radiation level (c = "
        "%lld, "
        "c->stars.radiation_level = %lld)!",
        c->cellID, c->stars.radiation_level->cellID);
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Ordering assertion for the feedback_ghost activation fix in
   * cell_radiation_activate_supers_in(). Every hydro.super this gather
   * touches must have its stars.feedback_ghost already RUN this step
   * (ti_run == e->ti_current, not merely unskipped, since skip resets to 1
   * on completion). A trip means the wait edge wired at
   * engine_radiation_wire_super_deps() failed to hold. */
  for (struct link *l = c->stars.radiation_in; l != NULL; l = l->next) {
    struct cell *sides[2] = {l->t->ci, l->t->cj};
    for (int side = 0; side < 2; side++) {
      struct cell *cc = sides[side];
      if (cc == NULL) continue;
      if (cc->hydro.super == NULL) continue; /* radiation_level above super */
      struct cell *super = cc->hydro.super;
      if (super->stars.feedback_ghost == NULL) continue;

      /* Only meaningful if this super was covered by radiation activation
       * this step: cell_flag_do_hydro_drift marks that scope. An untouched
       * super legitimately carries a stale ti_run. */
      if (!cell_get_flag(super, cell_flag_do_hydro_drift)) continue;

      if (super->stars.feedback != NULL /* had a feedback task registered */
          && super->stars.feedback_ghost->ti_run != e->ti_current) {
        error(
            "hii_ionization_feedback about to touch hydro.super %lld whose "
            "stars.feedback_ghost has NOT run this step (ti_run=%lld, "
            "ti_current=%lld, skip=%d). Ordering edge severed.",
            super->cellID, super->stars.feedback_ghost->ti_run, e->ti_current,
            super->stars.feedback_ghost->skip);
      }
    }
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  for (struct link *l = c->stars.radiation_in; l != NULL; l = l->next) {
    /* We have already handled the self case */
    if (l->t->type == task_type_self) continue;

    struct cell *cj = l->t->cj;
    struct cell *ci = l->t->ci;
    /* ci/cj may sit above their own hydro.super (radiation_level can
     * coarsen past it), in which case hydro.super is NULL; print -1
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
 * @brief Run one neighbour traversal for a star: its own cell plus every
 * cell reachable through the radiation task links.
 *
 * Both phases of a rebuild pass use the same spatial search and differ only
 * in @p mode, so the traversal (including the periodic shift and sort-index
 * bookkeeping) lives here once instead of being duplicated per phase.
 *
 * @param c The working-level #cell holding the star.
 * @param gather The #hii_gather_state of the traversal.
 */
static void runner_hii_visit_neighbors(struct cell *c,
                                       struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;

  /* First loop over particles in the current cell */
  runner_doself_stars_hii_ionization_feedback(c, gather);

  /* Now loop over particles in the neighboring cells via task links */
  for (struct link *l = c->stars.radiation_level->stars.radiation_in; l != NULL;
       l = l->next) {

    struct cell *cj = NULL;
    if (l->t->type == task_type_self) {
      cj = c->stars.radiation_level;
    } else {
      cj = (l->t->cj == c->stars.radiation_level) ? l->t->ci : l->t->cj;
    }

#ifdef SWIFT_DEBUG_CHECKS
#ifdef WITH_MPI
    /* Coverage assert: every foreign radiation_in partner this gather
     * is about to walk must have been fully allocated and linked, with a
     * live rho recv keeping it refreshed. A trip here means the foreign
     * coverage extension (engine_addtasks_send/recv_hydro) did not reach
     * this cell. */
    if (cj->nodeID != engine_rank) {
      if (cj->hydro.count > 0 && cj->hydro.parts == NULL)
        error("Foreign radiation_in partner %lld has NULL parts.", cj->cellID);
      if (cell_get_recv(cj, task_subtype_rho) == NULL)
        error("Foreign radiation_in partner %lld has no live rho recv.",
              cj->cellID);
    }
#endif
#endif

    /* Get the relative distance between the pairs, wrapping. */
    double shift[3];
    runner_hii_get_pair_shift(e->s->dim, e->s->periodic, c->loc, cj->loc,
                              shift);

    /* Get the sorting index. */
    int sid = 0;
    for (int k = 0; k < 3; k++)
      sid = 3 * sid + ((cj->loc[k] - c->loc[k] + shift[k] < 0)   ? 0
                       : (cj->loc[k] - c->loc[k] + shift[k] > 0) ? 2
                                                                 : 1);

    /* Switch the cells around? */
    const int flipped = runner_flip[sid];
    sid = sortlistID[sid];

    runner_do_stars_hii_ionization_feedback_branch(c, cj, sid, flipped, shift,
                                                   gather);
  } /* Neighbour search */
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Counters of the claim walk of one star pass, printed at its end.
 */
struct hii_pass_diagnostics {

  /*! Newly claimed this pass, summed over every retry/expansion iteration.
      Read alongside the held count to tell growing, holding, or churning
      regions apart; also feeds the HII CHECKSUM line. */
  int count_claimed;

  /*! Order-independent hash of this pass's claimed particle ids (sum of
      id*constant), so it matches across runs that claim the same set via
      different traversal orders. Unsigned: overflow wraparound is defined,
      not UB. */
  unsigned long long claimed_id_hash;

  /*! Claims per angular pixel. */
  int pixel_claimed_count[HII_MAX_ANGULAR_PIXELS];

  /*! Last roll outcome per pixel: -1 never rolled, 0 lost, 1 won. Only the
      last outcome is kept, since the fixed per-pixel random draw makes every
      loss this pass identical. */
  int pixel_roll_state[HII_MAX_ANGULAR_PIXELS];

  /*! Photon budget of each pixel when the claim walk starts. */
  double pixel_budget_start[HII_MAX_ANGULAR_PIXELS];

  /*! Candidates offered to a cone. */
  int count_evaluated;

  /*! Fronts met. */
  int count_fronts;

  /*! Clock reading at the start of the claim walk. */
  ticks tic_claim_walk;
};

/**
 * @brief Reset the counters at the start of the claim walk of a star.
 *
 * @param diag The #hii_pass_diagnostics.
 * @param si The #spart.
 */
static void hii_pass_diagnostics_init(struct hii_pass_diagnostics *diag,
                                      const struct spart *si) {
  diag->count_claimed = 0;
  diag->claimed_id_hash = 0;
  for (int p = 0; p < HII_MAX_ANGULAR_PIXELS; p++) {
    diag->pixel_claimed_count[p] = 0;
  }
  for (int p = 0; p < feedback_get_star_HII_pixel_count(si); p++) {
    diag->pixel_roll_state[p] = -1;
    diag->pixel_budget_start[p] = feedback_get_star_ionization_budget(si, p);
  }
  diag->count_evaluated = 0;
  diag->count_fronts = 0;
  diag->tic_claim_walk = getticks();
}

/**
 * @brief Record one candidate offered to a cone.
 *
 * @param diag The #hii_pass_diagnostics.
 * @param pj The candidate.
 * @param pixel Angular pixel of the candidate.
 * @param outcome The #hii_candidate_outcome.
 * @param claimed_this_candidate Did this star claim the candidate?
 * @param was_tagged_before Was the candidate tagged before the offer?
 */
static void hii_pass_diagnostics_record_candidate(
    struct hii_pass_diagnostics *diag, const struct part *pj, const int pixel,
    const int outcome, const int claimed_this_candidate,
    const int was_tagged_before) {
  ++diag->count_evaluated;
  if (outcome == hii_candidate_front) ++diag->count_fronts;

  if (claimed_this_candidate) {
    ++diag->count_claimed;
    diag->claimed_id_hash += (unsigned long long)pj->id * 0x9E3779B97F4A7C15ULL;
    ++diag->pixel_claimed_count[pixel];
    diag->pixel_roll_state[pixel] = 1;
  } else if (!was_tagged_before) {
    diag->pixel_roll_state[pixel] = 0;
  }
}

/**
 * @brief Print the cost of the claim walk, the repeat-run comparator and the
 * per-pixel accounting of a star pass.
 *
 * @param diag The #hii_pass_diagnostics.
 * @param e The #engine.
 * @param si The #spart.
 * @param count_held Gas the star held at the start of the pass.
 */
static void hii_pass_diagnostics_print(const struct hii_pass_diagnostics *diag,
                                       const struct engine *e,
                                       const struct spart *si,
                                       const int count_held) {
  message("HII walk: star %lld evaluated %d fronts %d walk_ms %.4f", si->id,
          diag->count_evaluated, diag->count_fronts,
          clocks_from_ticks(getticks() - diag->tic_claim_walk));

  /* Repeat-run comparator: claimed_id_hash (a sum of id * constant) is
     order independent, but leftover_budget on the same line is a
     floating-point accumulation and is NOT; compare the hash, not the
     budget. */
  message(
      "HII CHECKSUM: star_id %lld ti_current %lld n_claimed %d "
      "claimed_id_hash %llu held %d leftover_budget %.17e",
      si->id, (long long)e->ti_current, diag->count_claimed,
      diag->claimed_id_hash, count_held,
      feedback_get_star_ionization_budget_total(si));

  /* Compact per-pixel claim-accounting line: pixel index, claims, photons
     spent/left, and the last probabilistic roll outcome. The per-star
     CHECKSUM line above aggregates across pixels and cannot show this.
     One line per star per pass regardless of pixel count. */
  {
    char pixel_summary[16384];
    size_t off = 0;
    const int n_pix = feedback_get_star_HII_pixel_count(si);
    for (int p = 0; p < n_pix && off < sizeof(pixel_summary); p++) {
      const double left = feedback_get_star_ionization_budget(si, p);
      const char roll = diag->pixel_roll_state[p] < 0
                            ? 'x'
                            : (diag->pixel_roll_state[p] ? 'W' : 'L');
      off +=
          (size_t)snprintf(pixel_summary + off, sizeof(pixel_summary) - off,
                           "%s[%d:c=%d,sp=%.17g,lf=%.17g,r=%c]",
                           p == 0 ? "" : " ", p, diag->pixel_claimed_count[p],
                           diag->pixel_budget_start[p] - left, left, roll);
    }
    message("HII pixels: star %lld %s", si->id, pixel_summary);
  }
}
#else
struct hii_pass_diagnostics {};
static void hii_pass_diagnostics_init(struct hii_pass_diagnostics *diag,
                                      const struct spart *si) {}
static void hii_pass_diagnostics_record_candidate(
    struct hii_pass_diagnostics *diag, const struct part *pj, const int pixel,
    const int outcome, const int claimed_this_candidate,
    const int was_tagged_before) {}
static void hii_pass_diagnostics_print(const struct hii_pass_diagnostics *diag,
                                       const struct engine *e,
                                       const struct spart *si,
                                       const int count_held) {}
#endif

/**
 * @brief Budget interval and stamps of one star pass.
 */
struct hii_star_pass {

  /*! Interval whose photons this pass issues, and whose recombinations they
      pay. */
  double dt_back;

  /*! Gap since the star's previous pass, which the tags of this pass
      outlive. */
  double dt_tag;

  /*! Age of the star at the start of its step, floored at 0. */
  double star_age_safe;

  /*! Number of passes the star made before this one. */
  int64_t pass_index;
};

/**
 * @brief Open the pass of a star: stamp the attempt and open the photon
 * budget of the interval the pass covers.
 *
 * @param e The #engine.
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param pass (return) The budget interval and stamps, set if this returns 1.
 * @return 0 if the pass has no photons to issue or the cell holds no gas, 1
 * otherwise.
 */
static int runner_hii_open_star_pass(const struct engine *e,
                                     const struct cell *c, struct spart *si,
                                     struct hii_star_pass *pass) {

  const struct cosmology *cosmo = e->cosmology;
  const struct feedback_props *feedback_props = e->feedback_props;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;
  const double time = e->time;

  /* Needed by the attempt/rebuild stamps below and the retry loop's
     expansion cap; computed once. */
  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin_star = 0;
  /* si->time_bin here, live: this task's ordering relative to the
     timestep task's overwrite (runner_do_timestep()) is not established
     by any dependency edge, so whether it is pre- or post-overwrite for
     a given star is untraced; passing it through unchanged keeps this
     call's behavior identical to before old_time_bin existed. */
  compute_time(si, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin_star, ti_current, time_base, time, si->time_bin);
  const double star_age_safe = max(star_age_beg_step, 0.0);

  const double dt_since_attempt =
      star_age_safe - feedback_get_star_HII_last_attempt(si);
  feedback_set_star_HII_last_attempt(si, star_age_safe);
  const int64_t pass_index = feedback_count_star_HII_pass(si);

  /* Each pass issues the photons of [issued age, star age], so every instant
     of the star's life is issued once. At age 0 nothing has been emitted yet:
     the first pass issues the floor interval ahead. Only a star of the ICs
     and every star in "rebuild every step" mode reach this: with a positive
     rebuild time a newborn star first passes at about one cadence of age. */
  const double dt_floor = (double)feedback_props->HII_rebuild_floor_Myr;
  const double issued_age = feedback_get_star_HII_photons_issued_age(si);
  const int is_emitting = feedback_get_star_ionization_rate(si, 0) > 0.0;
  double issue_until_age = max(issued_age, star_age_safe);
  if (is_emitting && issue_until_age == 0.0) issue_until_age = dt_floor;
  feedback_set_star_HII_photons_issued_age(si, issue_until_age);
  const double dt_issued = issue_until_age - issued_age;

  /* No photons to issue: the star has no rate yet (the start-up pass of a
     star in the ICs), or the pass lands inside the interval the first pass
     issued ahead. The region of the last pass stands, and the rate sample of
     the trapezoid stays at the start of the next interval. */
  if (!is_emitting || dt_issued <= 0.0) {
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
    message(
        "Star %lld: HII pass skipped, no photons to issue (age %e, "
        "photons issued up to age %e)",
        si->id, star_age_safe, issue_until_age);
#endif
    return 0;
  }

  /* Nothing to search for in a gas-free working-level cell: the photons of
     this interval escaped. See runner_do_stars_hii_ionization_feedback()
     for why this cell is reached at all. The rate sample moves to the end of
     the escaped interval. */
  if (c->hydro.count == 0) {
    feedback_resync_star_ionizing_photon_rate_cache(si);
    return 0;
  }

  /* This pass recounts the gas the star ionizes, so restart from zero.
     The gas-free skip above does not rebuild, and keeps the last mass. */
  feedback_reset_star_HII_mass(si);
  si->h_hii_request = 0.f;

  /* The ceiling is a backstop: a star pinned below its own cadence (e.g. by
     min_star_timestep in very dense gas) would otherwise get the photons of
     a long interval in one pass. */
  const double dt_nominal =
      feedback_get_star_HII_nominal_interval(feedback_props, dt_enrichment);
  const double dt_ceiling = HII_DT_BACK_MAX_INTERVALS * dt_nominal;
  const double dt_back = min(dt_issued, dt_ceiling);
  feedback_open_star_ionizing_photon_budget(si, dt_back);

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  if (dt_issued > dt_ceiling)
    message(
        "Star %lld: HII budget interval clamped, %e -> %e (the interval "
        "exceeds twice the nominal cadence despite gas being present; "
        "this star may be pinned below its own desired cadence)",
        si->id, dt_issued, dt_back);
#endif

  /* A tag must outlive the gap to the star's next pass, which the gap since
     its last attempt measures. For a first pass that gap is 0 and the floor
     stands in, so the tag lasts 2 x floor and the next pass comes one
     cadence later. The issued interval can be much shorter: the first pass
     past the floor interval in "rebuild every step" mode. */
  const double dt_gap = dt_since_attempt > 0.0 ? dt_since_attempt : dt_floor;

  pass->dt_back = dt_back;
  pass->dt_tag = min(dt_gap, dt_ceiling);
  pass->star_age_safe = star_age_safe;
  pass->pass_index = pass_index;
  return 1;
}

/**
 * @brief Charge the recombination losses of the gas the star already holds
 * ionized, before any photons are spent on new candidates.
 *
 * Collapses the star's extent to what it still holds. Does nothing when
 * cooling is disabled for tagged gas, since tags are then permanent.
 *
 * @param r The #runner thread.
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param frame Rotation of the cone frame (row-major 3x3), or NULL for the
 * box axes.
 * @param dt_back Interval whose photons this pass issues.
 * @param dt_tag Gap since the star's previous pass, which tags outlive.
 * @param interaction_limit The initial search range boundary.
 * @return Number of gas particles the star holds ionized.
 */
static int runner_hii_maintain_held_region(
    struct runner *r, struct cell *c, struct spart *si, const double *frame,
    const double dt_back, const double dt_tag, const float interaction_limit) {

  /* Held outside the #ifndef so it is in scope for the retry loop's
     expansion cap either way. Stays 0 when Phase 1 is compiled out, which
     correctly grants that build the full configured search reach. */
  int count_held = 0;

#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const double time = e->time;

  /* Phase 1: charge the recombination losses of the gas this star already
     holds ionized, before any photons are spent on new candidates. Skipped
     when cooling is disabled for tagged gas, since tags are then permanent
     by construction and there is nothing to renew. */
  struct hii_maintenance_context maintenance_ctx = {
      .phys_const = phys_const,
      .hydro_props = hydro_props,
      .us = us,
      .cosmo = cosmo,
      .cooling = cooling,
      .time = time,
      .dt_back = dt_back,
      .r2_max_maintained = 0.f,
      .photons_charged = 0.0,
      .r = r,
      .buffer_count = 0,
  };
  struct hii_gather_state held_gather = {.r = r,
                                         .si = si,
                                         .search_radius = interaction_limit,
                                         .buffer = NULL,
                                         .max_size = 0,
                                         .count_found = 0,
                                         .mode = hii_gather_already_ionized,
                                         .ctx = &maintenance_ctx,
                                         .frame = frame};
  runner_hii_visit_neighbors(c, &held_gather);
  count_held = held_gather.count_found;

#ifdef SWIFT_DEBUG_CHECKS
  /* hii_gather_already_ionized increments both counters together, one
     entry buffered per particle held. A mismatch means the buffer and
     the counters that drive it have silently gone out of sync. */
  if (count_held != maintenance_ctx.buffer_count)
    error("count_held (%d) != maintenance buffer_count (%d)", count_held,
          maintenance_ctx.buffer_count);
#endif

  /* Charge every buffered already-ionized candidate in ascending (r2, id)
     order: a budget shortfall then lapses the outermost shell first (the
     recession front), instead of a set determined by cell-traversal
     order. Every entry is charged regardless of pixel budget;
     feedback_iact_HII_maintain_ionized_part's own check handles that and
     unconditionally accrues mass_HII_region first, so skipping entries
     here would silently shrink that output field. */
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  const ticks tic_maintenance_sort = getticks();
#endif
  if (maintenance_ctx.buffer_count > 1) {
    qsort(r->hii_maintenance_buffer, maintenance_ctx.buffer_count,
          sizeof(struct hii_neighbor), runner_hii_maintenance_compare);
  }
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  const ticks maintenance_sort_ticks = getticks() - tic_maintenance_sort;
#endif
  const integertime_t ti_begin_maintenance =
      get_integer_time_begin(e->ti_current - 1, si->time_bin);
  for (int k = 0; k < maintenance_ctx.buffer_count; k++) {
    const struct hii_neighbor *nb = &r->hii_maintenance_buffer[k];
    const double charged = feedback_iact_HII_maintain_ionized_part(
        si, nb->p, nb->xp, nb->r2, nb->pixel, phys_const, hydro_props, us,
        cosmo, cooling, ti_begin_maintenance, time, dt_back, dt_tag);
    if (charged > 0.0) maintenance_ctx.photons_charged += charged;
  }
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message("HII maintenance sort: star %lld, n=%d, sort_ticks=%lld (%.4f ms)",
          si->id, maintenance_ctx.buffer_count,
          (long long)maintenance_sort_ticks,
          clocks_from_ticks(maintenance_sort_ticks));
#endif

  /* The region's extent is what this star still holds, not the high-water
     mark of everything it ever claimed: gas that drifted out of reach must
     stop holding h_hii (and with it the search radius and task depth) open.
     Phase 2 grows it back below for anything newly claimed.

     Holding nothing must collapse the extent; leaving a stale large
     h_hii keeps radiation_level coarse, so task_lock goes on taking a wide
     lock set for a star doing no work. But holding gas at exactly r2 == 0
     (a co-located particle) must NOT write 0, since cell_sanitize would
     then inflate h_hii to the whole cell. */
  const float h_hii_held =
      count_held == 0
          ? 0.f
          : (maintenance_ctx.r2_max_maintained > 0.f
                 ? sqrtf(maintenance_ctx.r2_max_maintained) * kernel_gamma_inv
                 : si->h_hii);

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message(
      "HII pass: star %lld dt_back = %e, held = %d, recomb photons = %.17g, "
      "h_hii %e -> %e",
      si->id, dt_back, count_held, maintenance_ctx.photons_charged, si->h_hii,
      h_hii_held);
#endif

  si->h_hii = h_hii_held;
#endif

  return count_held;
}

/**
 * @brief How the claim walk of a star pass ended.
 */
struct hii_claim_summary {

  /*! Same-radius retries made after the neighbour buffer filled up. */
  int num_retry_full_buffer;

  /*! Radius expansions made. */
  int num_radius_expansions;

  /*! Radius of the last search. */
  float final_search_radius;
};

/**
 * @brief Claim new candidates with the photons the star has left, retrying
 * while they last.
 *
 * @param r The #runner thread.
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param frame Rotation of the cone frame (row-major 3x3), or NULL for the
 * box axes.
 * @param dt_back Interval whose photons this pass issues.
 * @param dt_tag Gap since the star's previous pass, which tags outlive.
 * @param interaction_limit The initial search range boundary.
 * @param max_reachable_search_radius The neighbour cell margin of this pass.
 * @param count_held Gas the star already held at the start of the pass.
 * @param diag The #hii_pass_diagnostics.
 * @param summary (return) How the walk ended.
 */
static void runner_hii_claim_and_retry(
    struct runner *r, struct cell *c, struct spart *si, const double *frame,
    const double dt_back, const double dt_tag, const float interaction_limit,
    const float max_reachable_search_radius, const int count_held,
    struct hii_pass_diagnostics *diag, struct hii_claim_summary *summary) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const struct feedback_props *feedback_props = e->feedback_props;
  const double time = e->time;

  const struct stars_props *star_props = e->stars_properties;
  const int max_retry_full_buffer = star_props->HII_max_retry_full_buffer;
  const int max_radius_expansion_tries =
      star_props->HII_max_radius_expansion_tries;
  const float radius_expansion_factor = star_props->HII_radius_expansion_factor;
  const float max_search_radius = star_props->HII_max_search_radius;

  struct hii_neighbor ngb_buffer[max_ngbs];

  /* Logic: a pass can fall short of exhausting the star's photon budget
   * for two DIFFERENT reasons, which need two different retries:
   *
   * 1. The buffer filled up (count_found == max_ngbs): there may be MORE
   *    candidates within the SAME radius that got bumped out by the
   *    nearest-max_ngbs cut. Retry at the same radius (cheap: the
   *    volume searched doesn't grow); existing max_retry_full_buffer
   *    mechanism.
   * 2. The buffer did NOT fill up (every not-yet-ionized particle within
   *    the current radius was found and processed) but photons remain:
   *    the RADIUS itself is the bottleneck, not the buffer. Retrying at
   *    the same radius would just re-find nothing new. Only expanding
   *    the radius can make progress; HII_max_radius_expansion_tries.
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
  int num_empty_expansions = 0;

  hii_pass_diagnostics_init(diag, si);

  while (1) {
    /* Any candidate actually claimed during THIS attempt (any pixel). A
       same-radius retry (case 1 below) re-runs the identical search: if
       this stays 0, no pixel's budget or tag state moved, so the fixed
       (star, pixel, ti_begin) random draw guarantees the retry would
       reproduce this exact attempt; see the case-1 branch below. */
    int attempt_claimed_count = 0;
    /* Fronts met during this attempt: each one empties its cone, so a
       retry gathers past it. */
    int attempt_front_count = 0;

    struct hii_gather_state gather = {.r = r,
                                      .si = si,
                                      .search_radius = dynamic_search_radius,
                                      .buffer = ngb_buffer,
                                      .max_size = max_ngbs,
                                      .count_found = 0,
                                      .mode = hii_gather_new_candidates,
                                      .ctx = NULL,
                                      .frame = frame};
    runner_hii_visit_neighbors(c, &gather);
    const int count_found = gather.count_found;

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
           break; other entries in the buffer may belong to a pixel
           that still has budget). */
        if (feedback_get_star_ionization_budget(si, pixel) <= 0.0) {
          continue;
        }

        struct part *pj = ngb_buffer[k].p;
        struct xpart *xpj = ngb_buffer[k].xp;
        const float r2 = ngb_buffer[k].r2;
        const int was_tagged_before =
            radiation_is_part_tagged_as_ionized(pj, xpj);

        /* Do the ionization */
        const int outcome = feedback_iact_HII_ionization(
            si, pj, xpj, r2, pixel, phys_const, hydro_props, us, cosmo, cooling,
            feedback_props, ti_begin, time, dt_back, dt_tag);

        if (outcome == hii_candidate_front) ++attempt_front_count;

        /* A candidate starts this loop untagged, so tagged-with-si->id
           afterwards means this star just claimed it. A concurrent claim
           by another star's task tags a different id instead and is
           correctly excluded here. feedback_hii_claim_part always tags
           AND consumes budget together, so an untagged result here means
           the roll was lost or the candidate is a front. */
        const int claimed_this_candidate =
            !was_tagged_before &&
            radiation_is_part_tagged_as_ionized(pj, xpj) &&
            radiation_get_part_ionized_star_id(pj, xpj) == si->id;
        if (claimed_this_candidate) ++attempt_claimed_count;

        hii_pass_diagnostics_record_candidate(diag, pj, pixel, outcome,
                                              claimed_this_candidate,
                                              was_tagged_before);
      } /* Loop over the sorted particles */
    }

    /* Fully exhausted: nothing more to do this pass. */
    if (feedback_get_star_ionization_budget_max(si) <= 0.0) break;

    if (buffer_was_full) {
      /* Case 1: more candidates may exist within the SAME radius,
       * provided this attempt claimed at least one of them. If it
       * claimed none, no pixel's budget or tag state changed anywhere,
       * so a retry re-runs the identical search (same eligible set),
       * meets the identical per-pixel budgets, and every pixel's fixed
       * (star, pixel, ti_begin) random draw replays the identical
       * win/lose outcome; a provable no-op, not a heuristic. A front
       * also changes the state: it empties its cone, so a retry gathers
       * past it. */
      if (attempt_claimed_count == 0 && attempt_front_count == 0) break;
      if (num_retry_full_buffer >= max_retry_full_buffer) break;
      ++num_retry_full_buffer;
      num_empty_expansions = 0;
    } else {
      /* Case 2: genuinely nothing left within the current radius, but
       * photons remain: expand the radius unless one of the limits in
       * runner_hii_search_stop() binds. */
      const enum hii_search_stop stop = runner_hii_search_stop(
          dynamic_search_radius, max_search_radius, max_reachable_search_radius,
          num_radius_expansions, max_radius_expansion_tries, count_found,
          count_held, num_empty_expansions);

      if (stop == hii_search_stop_margin) {
        /* Only the cell layout stops this star, so ask the next passes,
           and the rebuild check once the radius outgrows the wired
           neighbours, for the radius it would have tried next. */
        si->h_hii_request = runner_hii_requested_radius(dynamic_search_radius,
                                                        radius_expansion_factor,
                                                        max_search_radius) /
                            (radiation_search_radius_factor * kernel_gamma);
#ifdef SWIFT_DEBUG_CHECKS
        atomic_inc(&e->radiation_reach_clamp_count);
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
        message(
            "Star %lld: radius expansion clamped at the reachable stencil "
            "bound (dynamic_search_radius = %e, max_reachable = %e), "
            "requesting %e.",
            si->id, dynamic_search_radius, max_reachable_search_radius,
            kernel_gamma * si->h_hii_request);
#endif
#endif
      }
      if (stop != hii_search_grow) break;

      if (count_found == 0 && count_held > 0)
        ++num_empty_expansions;
      else
        num_empty_expansions = 0;

      dynamic_search_radius =
          min3(dynamic_search_radius * radius_expansion_factor,
               max_search_radius, max_reachable_search_radius);
      ++num_radius_expansions;
    }
  }

  summary->num_retry_full_buffer = num_retry_full_buffer;
  summary->num_radius_expansions = num_radius_expansions;
  summary->final_search_radius = dynamic_search_radius;
}

/**
 * @brief Close the pass of a star: update the extent of the working-level
 * cell and stamp the rebuild.
 *
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param star_age_safe Age of the star at the start of its step, floored at 0.
 * @param summary How the claim walk ended.
 * @param diag The #hii_pass_diagnostics.
 * @param count_held Gas the star held at the start of the pass.
 */
static void runner_hii_close_star_pass(struct cell *c, struct spart *si,
                                       const double star_age_safe,
                                       const struct hii_claim_summary *summary,
                                       const struct hii_pass_diagnostics *diag,
                                       const int count_held) {

  const float h_search = stars_get_hii_search_length(si);
  c->stars.h_hii_max = max(c->stars.h_hii_max, h_search);
  c->stars.h_hii_max_active = max(c->stars.h_hii_max_active, h_search);

  /*****************************************/
  /* Update the star after HII ionization */

  /* Log when this HII region was (re)built. The clamped age is stamped, not
     the raw one: HII_region_last_rebuild is documented as never holding a
     negative value, and compute_time() can return one for a newborn. */
  feedback_set_star_HII_last_rebuild(si, star_age_safe);
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  if (feedback_get_star_ionization_budget_max(si) <= 0.0) {
    message(
        "Star %lld has exhausted all its ionizing photons! r_hii = %e, "
        "num_retry_full_buffer = %d, num_radius_expansions = %d",
        si->id, si->h_hii * kernel_gamma, summary->num_retry_full_buffer,
        summary->num_radius_expansions);
  } else {
    message(
        "Star %lld has NOT exhausted all its ionizing photons! Remaining: "
        "%e, final_search_radius = %e, sp->h_hii = %e, "
        "num_retry_full_buffer = %d, num_radius_expansions = %d",
        si->id, feedback_get_star_ionization_budget_max(si),
        summary->final_search_radius, kernel_gamma * si->h_hii,
        summary->num_retry_full_buffer, summary->num_radius_expansions);
  }
  message("Star %lld: claimed = %d this pass, held = %d from previous ones",
          si->id, diag->count_claimed, count_held);
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
 * (Stars:HII_max_retry_full_buffer); if the buffer did NOT fill up (every
 * reachable candidate was already found and ionized) the radius itself is
 * the bottleneck, so it is expanded
 * (Stars:HII_max_radius_expansion_tries, Stars:HII_radius_expansion_factor)
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
  const struct feedback_props *feedback_props = e->feedback_props;

  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  /* Anything to do here? c->hydro.count == 0 is deliberately NOT checked
     here; see the gas-free skip of runner_hii_open_star_pass(), and the
     matching comment in runner_do_stars_hii_ionization_feedback's own
     early-return. */
  if (c->stars.count == 0 || !cell_is_active_stars(c, e)) return;

  /* The gather only visits wired cells whatever the search radius says, so
   * the radius stops at this margin. A star stopped there asks for a wider
   * reach (h_hii_request) rather than waiting for an unrelated rebuild. */
  const float max_reachable_search_radius =
      interaction_limit + cell_radiation_neighbour_margin(c);

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
    if (si->ti_drift != e->ti_current)
      error(
          "Particle si (%lld) not drifted to current time c = %lld, "
          "c->super = %lld",
          si->id, c->cellID, c->super->cellID);
#endif

    struct hii_star_pass pass;
    if (!runner_hii_open_star_pass(e, c, si, &pass)) continue;

    double cone_frame_storage[9];
    const double *frame = runner_hii_star_cone_frame(
        si, feedback_props, pass.pass_index, cone_frame_storage);

    const int count_held = runner_hii_maintain_held_region(
        r, c, si, frame, pass.dt_back, pass.dt_tag, interaction_limit);

    struct hii_pass_diagnostics diag;
    struct hii_claim_summary summary;
    runner_hii_claim_and_retry(r, c, si, frame, pass.dt_back, pass.dt_tag,
                               interaction_limit, max_reachable_search_radius,
                               count_held, &diag, &summary);

    hii_pass_diagnostics_print(&diag, e, si, count_held);
    runner_hii_close_star_pass(c, si, pass.star_age_safe, &summary, &diag,
                               count_held);
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
 * @param ci The local cell context containing the active source star particle.
 * @param cj The target neighbor cell containing gas particle candidates.
 * @param sid The sorting axis index between the pair.
 * @param flipped Flag indicating if the cell order is flipped relative to the
 * axis direction.
 * @param shift The periodic wrapping shift vector applied to the coordinates.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_do_stars_hii_ionization_feedback_branch(
    struct cell *ci, struct cell *cj, const int sid, const int flipped,
    const double shift[3], struct hii_gather_state *gather) {

  /* OPTIMIZATION: Check if cj is completely out of reach before doing anything
   * else */
  if (!runner_hii_check_cell_can_be_reached(
          ci, cj, sid, flipped, shift, gather->si, gather->search_radius)) {
    return;
  }

  /* Don't reprocess cj or its child if it's ci. These particles were
     processed in the self interaction */
  if (ci == cj) return;

  if (cj->split) {

    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_branch(ci, cj->progeny[k], sid,
                                                       flipped, shift, gather);
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
     * sids, so the ionized set is invariant to which sid we use; only speed
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
      runner_dopair_naive_stars_hii_ionization_feedback(ci, cj, shift, gather);
    } else {
      runner_dopair_stars_hii_ionization_feedback(ci, cj, use_sid, use_flipped,
                                                  shift, gather);
    }
  }
}

/**
 * @brief Gather gas particles within the HII radius from the star's own cell.
 *
 * Performs a brute-force search over all hydro particles in cell @p c that
 * fall within the HII smoothing length of the star. Candidates are
 * added to the buffer if they are not already inhibited or ionized.
 *
 * @param c The #cell containing the gas particles.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_doself_stars_hii_ionization_feedback(
    struct cell *c, struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;
  struct spart *si = gather->si;
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
  const float r2_max = gather->search_radius * gather->search_radius;
  const float six[3] = {si->x[0], si->x[1], si->x[2]};

  /* TODO: Use sorted cells to check that the particles/the cell are/is
     within the star's h_hii (or h if h_hii == 0) */

  /* Loop over the parts in c. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. A foreign cell carries no
     * local xparts array; xpj must resolve to a genuine NULL there rather
     * than undefined NULL-based pointer arithmetic (feedback_get_
     * eligibility_temperature owns the xp == NULL contract downstream). */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = (xparts != NULL) ? &xparts[pjd] : NULL;

    runner_hii_offer_candidate(gather, pj, xpj, c, six, NULL, r2_max);
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
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param shift The periodic wrapping shift vector applied to @p cj.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_dopair_naive_stars_hii_ionization_feedback(
    struct cell *ci, struct cell *cj, const double shift[3],
    struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;
  struct spart *si = gather->si;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  if (cj != ci && cj->split) {
    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_dopair_naive_stars_hii_ionization_feedback(ci, cj->progeny[k],
                                                          shift, gather);
      }
    }
  } else {

    /* Gas comes from cj, the star from ci: check each for what it supplies. */
    if (!cell_are_part_drifted(cj, e))
      error("Interacting undrifted cell (hydro).");
    if (!cell_are_spart_drifted(ci, e))
      error("Interacting undrifted cell (stars).");

    struct part *restrict parts_j = cj->hydro.parts;
    struct xpart *restrict xparts_j = cj->hydro.xparts;
    const float r2_max = gather->search_radius * gather->search_radius;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. See the doself version above for
       * why xpj must be a genuine NULL when cj has no local xparts. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = (xparts_j != NULL) ? &xparts_j[pjd] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, cj->loc, r2_max);
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
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param sid The sorting axis index between the pair.
 * @param flipped Flag indicating if the cell order is flipped relative to the
 * axis direction.
 * @param shift The periodic wrapping shift vector applied to @p cj.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_dopair_stars_hii_ionization_feedback(
    struct cell *ci, struct cell *cj, const int sid, const int flipped,
    const double shift[3], struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;
  struct spart *si = gather->si;
  const float search_radius = gather->search_radius;
  const int count_j = cj->hydro.count;

  /* Anything to do here?*/
  if (count_j == 0) return;

  /* Gas comes from cj, the star from ci: check each for what it supplies. */
  if (!cell_are_part_drifted(cj, e))
    error("Interacting undrifted cell (hydro).");
  if (!cell_are_spart_drifted(ci, e))
    error("Interacting undrifted cell (stars).");

  /* Get the hydro sorts for our gas cell */
  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
  const float dx_max = cj->hydro.dx_max_sort + ci->stars.dx_max_sort;

  /* The gas particle projections (sort_j[].d) and the r2 distances below both
   * use ABSOLUTE gas positions (pj->x), so the star is expressed in that frame
   * by runner_hii_axis_window(). `flipped` only selects the scan direction. */
  double pix[3], di_min, di_max;
  runner_hii_axis_window(si, shift, sid, search_radius, dx_max, pix, &di_min,
                         &di_max);

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

      /* See the doself version (runner_doself_stars_hii_ionization_feedback)
       * for why xpj must be a genuine NULL when cj has no local xparts. */
      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj =
          (xparts_j != NULL) ? &xparts_j[sort_j[pjd].i] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, NULL, r2_max);
    }
  } else {
    /* Star is on the 'right' relative to the axis direction.
     * We scan backward through the sorted list until particles are below
     * di_min. */
    for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d >= di_min; pjd--) {

      /* Skip particles that are too far on the other side of the slab */
      if (sort_j[pjd].d > di_max) continue;

      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj =
          (xparts_j != NULL) ? &xparts_j[sort_j[pjd].i] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, NULL, r2_max);
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

  /* A gas particle can lie outside the box of cj by dx_outside_max, counted
   * from the last rebuild, which dx_max_sort (reset at each sort) can miss.
   * dx_outside_max is a maximum over the three axes: projected on the sorting
   * axis it grows by the sum of the absolute axis components. */
  const float axis_sum =
      (float)(fabs(runner_shift[sid][0]) + fabs(runner_shift[sid][1]) +
              fabs(runner_shift[sid][2]));
  const float dx_max =
      max(cj->hydro.dx_max_sort, cj->hydro.dx_outside_max * axis_sum) +
      ci->stars.dx_max_sort;

  /* The cell extent below (dj_min/dj_max) is built from cj's ABSOLUTE geometry
   * (cj->loc), which is the frame of runner_hii_axis_window(), also used by
   * the sorted scan. A flipped-dependent sign here would offset the reach
   * window by 2*shift and wrongly skip wrapped neighbours. */
  double pix[3], di_min, di_max;
  runner_hii_axis_window(si, shift, sid, search_radius, dx_max, pix, &di_min,
                         &di_max);

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
