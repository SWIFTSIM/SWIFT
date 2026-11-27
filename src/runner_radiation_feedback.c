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
 * The buffer doubles on demand and never shrinks. Its content is not kept:
 * every maintenance pass refills it.
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
 * The id breaks ties: distinct particles can sit at the same squared
 * distance. Branches replace a subtraction, which loses precision on r2 and
 * overflows on a `long long` id.
 *
 * Compiled out with its only caller, the maintenance-buffer sort, under
 * NO_COOLING, so that -Wunused-function stays quiet.
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
 * Every gather call site shares this sequence.
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
    /* Extent and count track what the star still holds, not what it can
       afford, so they do not depend on traversal order. Gating them on the
       charge would shrink h_hii and drop the region's outskirts from the next
       search radius. Only the charge waits for the caller's (r2, id) order. */
    ctx->r2_max_maintained = max(ctx->r2_max_maintained, r2);
    gather->count_found++;

    /* The held region need not fit in max_ngbs, so this buffer grows instead
       of dropping outskirts. */
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

  /* A pixel without photons must not crowd other pixels out of the buffer. */
  if (feedback_get_star_ionization_budget(si, pixel) > 0.0)
    runner_hii_buffer_insert(gather->buffer, gather->max_size,
                             &gather->count_found, r2, pj, xpj, c, pixel);
}

/**
 * @brief Count a gas particle held by another star if the search reached it
 * first and its pixel has photons left.
 *
 * @param gather The #hii_gather_state of the traversal.
 * @param dx Star-minus-particle position offset.
 * @param r2 Squared distance between star and particle.
 * @param r2_max Squared search radius.
 */
__attribute__((always_inline)) INLINE static void
runner_hii_count_held_by_other_star(struct hii_gather_state *gather,
                                    const float dx[3], const float r2,
                                    const float r2_max) {

  if (r2 >= r2_max || r2 <= gather->r2_others_min) return;

  struct spart *si = gather->si;
  const int pixel = runner_hii_get_pixel(
      dx, feedback_get_star_HII_pixel_count(si), gather->frame);
  if (feedback_get_star_ionization_budget(si, pixel) > 0.0)
    gather->count_held_by_others++;
}

/**
 * @brief Is this gas particle held by a star other than the one gathering new
 * candidates?
 *
 * @param pj The gas particle.
 * @param xpj The particle's extended data, or NULL for a foreign cell.
 * @param gather The #hii_gather_state of the traversal.
 */
__attribute__((always_inline)) INLINE static int
runner_hii_part_held_by_other_star(const struct part *pj,
                                   const struct xpart *xpj,
                                   const struct hii_gather_state *gather) {
  return gather->mode == hii_gather_new_candidates &&
         radiation_is_part_tagged_as_ionized(pj, xpj) &&
         radiation_get_part_ionized_star_id(pj, xpj) != gather->si->id;
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
 * #runner_hii_try_insert_candidate, or to
 * #runner_hii_count_held_by_other_star if another star holds it.
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

  /* Gas held by another star fails the filter but is still counted. */
  if (part_is_inhibited(pj, e)) return;
  const int held_by_other = runner_hii_part_held_by_other_star(pj, xpj, gather);
  if (!held_by_other && !runner_hii_part_passes_gather_filter(
                            pj, xpj, e, gather->si->id, gather->mode))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that particles have been drifted to the current time */
  if (pj->ti_drift != e->ti_current)
    error("Particle pj (%lld) not drifted to current time. c = %lld.", pj->id,
          c->cellID);
#endif

  const float pjx[3] = {
      (float)(origin == NULL ? pj->x[0] : pj->x[0] - origin[0]),
      (float)(origin == NULL ? pj->x[1] : pj->x[1] - origin[1]),
      (float)(origin == NULL ? pj->x[2] : pj->x[2] - origin[2])};
  const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

  if (held_by_other)
    runner_hii_count_held_by_other_star(gather, dx, r2, r2_max);
  else
    runner_hii_try_insert_candidate(gather, pj, xpj, c, dx, r2, r2_max);
}

/**
 * @brief Start-up of the subgrid radiation: print warnings about the HII
 * search-radius / box-size relationship.
 *
 * @param e The #engine.
 */
void feedback_radiation_startup_diagnostics(const struct engine *e) {

  if (!(e->feedback_props->radiation_policy & radiation_policy_photoionization))
    return;

  /* Not fatal, but a search radius capped by the box looks like a converged
     result, so say so at start-up. */
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
     coarser grid to rebuild into. Only the stencil-completeness fallback
     keeps the rebuild criterion satisfiable, and a box that cannot coarsen
     to 3 top-level cells has no such fallback. */
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

  /* A due star in a gas-free working-level cell must still reach
     runner_dosub_stars_hii_ionization_feedback() to get its attempt and its
     issued photons stamped, so the gas count is not tested here. */
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
  /* Every hydro.super this gather touches must have run its
   * stars.feedback_ghost this step (ti_run == e->ti_current, as skip resets
   * to 1 on completion). A trip means the wait edge wired in
   * engine_radiation_wire_super_deps() failed to hold. */
  for (struct link *l = c->stars.radiation_in; l != NULL; l = l->next) {
    struct cell *sides[2] = {l->t->ci, l->t->cj};
    for (int side = 0; side < 2; side++) {
      struct cell *cc = sides[side];
      if (cc == NULL) continue;
      if (cc->hydro.super == NULL) continue; /* radiation_level above super */
      struct cell *super = cc->hydro.super;
      if (super->stars.feedback_ghost == NULL) continue;

      /* Only a super covered by the radiation activation this step
       * (cell_flag_do_hydro_drift) is checked: others keep a stale ti_run. */
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
    /* The self case is handled elsewhere. */
    if (l->t->type == task_type_self) continue;

    struct cell *cj = l->t->cj;
    struct cell *ci = l->t->ci;
    /* hydro.super is NULL when radiation_level is coarser than it. */
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

  if (c->split) {
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
    /* Working level reached. */
    runner_dosub_stars_hii_ionization_feedback(r, c, interaction_limit);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Run one neighbour traversal for a star: its own cell plus every
 * cell reachable through the radiation task links.
 *
 * Both phases of a rebuild pass share this search and differ only in the
 * gather mode.
 *
 * @param c The working-level #cell holding the star.
 * @param gather The #hii_gather_state of the traversal.
 */
static void runner_hii_visit_neighbors(struct cell *c,
                                       struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;

  /* The star's own cell. */
  runner_doself_stars_hii_ionization_feedback(c, gather);

  /* The neighbouring cells, through the task links. */
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
    /* Every foreign radiation_in partner must be allocated and linked, with a
     * live rho recv keeping it fresh. A trip means the foreign coverage
     * extension (engine_addtasks_send/recv_hydro) missed this cell. */
    if (cj->nodeID != engine_rank) {
      if (cj->hydro.count > 0 && cj->hydro.parts == NULL)
        error("Foreign radiation_in partner %lld has NULL parts.", cj->cellID);
      if (cell_get_recv(cj, task_subtype_rho) == NULL)
        error("Foreign radiation_in partner %lld has no live rho recv.",
              cj->cellID);
    }
#endif
#endif

    /* Periodic shift between the pair. */
    double shift[3];
    runner_hii_get_pair_shift(e->s->dim, e->s->periodic, c->loc, cj->loc,
                              shift);

    /* Sorting index. */
    int sid = 0;
    for (int k = 0; k < 3; k++)
      sid = 3 * sid + ((cj->loc[k] - c->loc[k] + shift[k] < 0)   ? 0
                       : (cj->loc[k] - c->loc[k] + shift[k] > 0) ? 2
                                                                 : 1);

    /* The sorted axis may point from cj to c. */
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

  /* Compare runs by claimed_id_hash, which is order independent. The
     leftover_budget on the same line is a floating-point sum and is not. */
  message(
      "HII CHECKSUM: star_id %lld ti_current %lld n_claimed %d "
      "claimed_id_hash %llu held %d leftover_budget %.17e",
      si->id, (long long)e->ti_current, diag->count_claimed,
      diag->claimed_id_hash, count_held,
      feedback_get_star_ionization_budget_total(si));

  /* One line per star and pass: for each pixel, the claims, the photons spent
     and left, and the last roll (W won, L lost, x none). The CHECKSUM line
     sums over pixels and cannot show this. */
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

  /*! feedback_get_star_HII_carry_tree_reach() before the budget opened:
      the tree reach of the last coverage stop, signed. */
  float last_tree_reach;
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

  double star_age_beg_step = 0;
  double dt_enrichment = 0;
  integertime_t ti_begin_star = 0;
  /* The live si->time_bin: no dependency orders this task against the
     overwrite in runner_do_timestep(), so it may be either value. */
  compute_time(si, with_cosmology, cosmo, &star_age_beg_step, &dt_enrichment,
               &ti_begin_star, ti_current, time_base, time, si->time_bin);
  const double star_age_safe = max(star_age_beg_step, 0.0);

  const double dt_since_attempt =
      star_age_safe - feedback_get_star_HII_last_attempt(si);
  feedback_set_star_HII_last_attempt(si, star_age_safe);
  const int64_t pass_index = feedback_count_star_HII_pass(si);

  /* Only a star of the ICs and every star in "rebuild every step" mode reach
     the floor interval: with a positive rebuild time a newborn star first
     passes at about one cadence of age. */
  const double dt_floor = (double)feedback_props->HII_rebuild_floor_Myr;
  const double issued_age = feedback_get_star_HII_photons_issued_age(si);
  const int is_emitting = feedback_get_star_ionization_rate(si, 0) > 0.0;
  const double issue_until_age = feedback_get_star_HII_issue_until_age(
      issued_age, star_age_safe, dt_floor, is_emitting);
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
  const float last_tree_reach = feedback_get_star_HII_carry_tree_reach(si);
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
  pass->last_tree_reach = last_tree_reach;
  return 1;
}

#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
/**
 * @brief Extent h_hii of the region a star still holds ionized.
 *
 * Compiled out with its only caller under NO_COOLING, so that
 * -Wunused-function stays quiet.
 *
 * The extent is what the star holds now, not the largest it ever had, so gas
 * that drifted away stops holding h_hii open. Holding nothing collapses it: a
 * stale h_hii keeps radiation_level coarse and the lock set wide. Gas at
 * exactly r2 == 0 keeps the current extent, since writing 0 would make
 * cell_sanitize inflate h_hii to the whole cell.
 *
 * @param count_held Number of gas particles the star holds ionized.
 * @param r2_max_held Squared distance to the star of the farthest of them.
 * @param h_hii_now The current h_hii of the star.
 * @return The h_hii the star holds after this pass.
 */
__attribute__((always_inline)) INLINE static float runner_hii_held_extent(
    const int count_held, const float r2_max_held, const float h_hii_now) {
  if (count_held == 0) return 0.f;
  if (r2_max_held > 0.f) return sqrtf(r2_max_held) * kernel_gamma_inv;
  return h_hii_now;
}
#endif

/**
 * @brief Charge the recombination losses of the gas the star already holds
 * ionized, before any photons are spent on new candidates.
 *
 * The held gas is gathered out to runner_hii_held_search_radius(), not to the
 * start radius. Collapses the star's extent to what it still holds. Does
 * nothing when cooling is disabled for tagged gas, since tags are then
 * permanent.
 *
 * @param r The #runner thread.
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param frame Rotation of the cone frame (row-major 3x3), or NULL for the
 * box axes.
 * @param dt_back Interval whose photons this pass issues.
 * @param dt_tag Gap since the star's previous pass, which tags outlive.
 * @param covered_radius Radius inside which the pass visits every cell that
 * held gas at the last rebuild.
 * @param held_radius (return) Distance to the star of the farthest gas it
 * holds, 0 if it holds none.
 * @return Number of gas particles the star holds ionized.
 */
static int runner_hii_maintain_held_region(
    struct runner *r, struct cell *c, struct spart *si, const double *frame,
    const double dt_back, const double dt_tag, const float covered_radius,
    float *held_radius) {

  /* Stays 0 when the maintenance below is compiled out, which grants the full
     configured search reach. */
  int count_held = 0;
  *held_radius = 0.f;

#ifndef IONIZATION_FEEDBACK_DEBUG_NO_COOLING
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const double time = e->time;

  /* Skipped when cooling is disabled for tagged gas: tags are then permanent
     and there is nothing to renew. */
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

  const float held_search_radius = runner_hii_held_search_radius(
      covered_radius, e->stars_properties->HII_max_search_radius);

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  const ticks tic_held_gather = getticks();
#endif
  struct hii_gather_state held_gather = {.r = r,
                                         .si = si,
                                         .search_radius = held_search_radius,
                                         .buffer = NULL,
                                         .max_size = 0,
                                         .count_found = 0,
                                         .mode = hii_gather_already_ionized,
                                         .ctx = &maintenance_ctx,
                                         .frame = frame};
  runner_hii_visit_neighbors(c, &held_gather);
  count_held = held_gather.count_found;
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  const ticks held_gather_ticks = getticks() - tic_held_gather;
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* The gather buffers one entry per particle held, so both counts agree. */
  if (count_held != maintenance_ctx.buffer_count)
    error("count_held (%d) != maintenance buffer_count (%d)", count_held,
          maintenance_ctx.buffer_count);
#endif

  /* Charge the held gas in ascending (r2, id) order: a budget shortfall then
     lapses the outermost shell first, not a set fixed by the cell traversal.
     Every entry is charged whatever its pixel budget, because
     feedback_iact_HII_maintain_ionized_part() checks that itself and first
     accrues mass_HII_region, so skipping one would shrink that field. */
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
  message(
      "HII maintenance sort: star %lld, n=%d, sort_ticks=%lld (%.4f ms), "
      "gather_ms %.4f",
      si->id, maintenance_ctx.buffer_count, (long long)maintenance_sort_ticks,
      clocks_from_ticks(maintenance_sort_ticks),
      clocks_from_ticks(held_gather_ticks));
#endif

  /* The claim walk grows the extent again for new gas. */
  const float h_hii_held = runner_hii_held_extent(
      count_held, maintenance_ctx.r2_max_maintained, si->h_hii);

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message(
      "HII pass: star %lld dt_back = %e, held = %d, recomb photons = %.17g, "
      "h_hii %e -> %e",
      si->id, dt_back, count_held, maintenance_ctx.photons_charged, si->h_hii,
      h_hii_held);
#endif

  si->h_hii = h_hii_held;
  if (count_held > 0) *held_radius = sqrtf(maintenance_ctx.r2_max_maintained);
#endif

  return count_held;
}

/*! HII coverage stops of this rank, the same summed over the ranks at the
    end of each step, and that sum at the last report. All cumulative since
    the run or its last restart: they are not in the restart file. */
static struct radiation_hii_coverage_stats radiation_hii_coverage_local;
static struct radiation_hii_coverage_stats radiation_hii_coverage;
static struct radiation_hii_coverage_stats radiation_hii_coverage_last_report;

/**
 * @brief Sum the HII coverage stops of all the ranks. Collective.
 */
void feedback_radiation_hii_coverage_reduce(void) {
#ifdef WITH_MPI
  if (MPI_Allreduce(
          &radiation_hii_coverage_local, &radiation_hii_coverage,
          sizeof(struct radiation_hii_coverage_stats) / sizeof(double),
          MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to aggregate the HII coverage stops.");
#else
  radiation_hii_coverage = radiation_hii_coverage_local;
#endif
}

/**
 * @brief Print the HII coverage stops of the run when a step added some.
 * Rank 0 only.
 */
void feedback_radiation_hii_coverage_report(void) {

  /* A pass stopped at the edge of the covered gas asks for a wider tree, and
   * released photons mean the rebuild could not widen it. Cut start radii
   * alone print nothing. Photons left at one-probe, tries and ceiling stops
   * escape, and their share of the issued photons measures the search stall. */
  const struct radiation_hii_coverage_stats *cov = &radiation_hii_coverage;
  const struct radiation_hii_coverage_stats *cov_last =
      &radiation_hii_coverage_last_report;
  const int new_coverage_stops =
      cov->carried_passes + cov->released_passes >
      cov_last->carried_passes + cov_last->released_passes;
  const int new_one_probe_stops =
      cov->one_probe_stops > cov_last->one_probe_stops;
  const int new_growth_limit_stops =
      cov->growth_limit_stops > cov_last->growth_limit_stops;

  if (new_coverage_stops)
    message(
        "HII coverage stops: %.0f carried (+%.0f), %.0f released (+%.0f) "
        "with %.3e photons (+%.3e), %.0f start radii cut (+%.0f).",
        cov->carried_passes, cov->carried_passes - cov_last->carried_passes,
        cov->released_passes, cov->released_passes - cov_last->released_passes,
        cov->released_photons,
        cov->released_photons - cov_last->released_photons,
        cov->start_radius_cuts,
        cov->start_radius_cuts - cov_last->start_radius_cuts);
  if (new_one_probe_stops)
    message(
        "HII one-probe stops with photons left: %.0f (+%.0f) with %.3e "
        "photons (+%.3e).",
        cov->one_probe_stops, cov->one_probe_stops - cov_last->one_probe_stops,
        cov->one_probe_photons,
        cov->one_probe_photons - cov_last->one_probe_photons);
  if (new_growth_limit_stops)
    message(
        "HII tries or ceiling stops with photons left: %.0f (+%.0f) with "
        "%.3e photons (+%.3e).",
        cov->growth_limit_stops,
        cov->growth_limit_stops - cov_last->growth_limit_stops,
        cov->growth_limit_photons,
        cov->growth_limit_photons - cov_last->growth_limit_photons);
  if (new_coverage_stops || new_one_probe_stops || new_growth_limit_stops)
    radiation_hii_coverage_last_report = *cov;
}

/**
 * @brief Keep or release the unspent photons of a pass stopped at the covered
 * radius, and count the outcome.
 *
 * The decision compares tree reaches (cell_radiation_tree_reach()), not stop
 * radii: a stop radius moves with the star, the tree reach only at a rebuild.
 *
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param feedback_props The #feedback_props.
 * @param last_tree_reach feedback_get_star_HII_carry_tree_reach() before the
 * budget of this pass opened.
 * @param stats (return) The coverage-stop counters of this rank.
 * @return The photons released (escaped), 0 if they are kept.
 */
double runner_hii_carry_at_coverage_stop(
    const struct cell *c, struct spart *si,
    const struct feedback_props *feedback_props, const float last_tree_reach,
    struct radiation_hii_coverage_stats *stats) {

  const double released = feedback_carry_star_HII_unspent_photons(
      si, feedback_props, cell_radiation_tree_reach(c), last_tree_reach);

  const int carried = feedback_get_star_HII_carry_tree_reach(si) > 0.f;
  if (carried) {
    atomic_add_d(&stats->carried_passes, 1.);
  } else {
    atomic_add_d(&stats->released_passes, 1.);
    atomic_add_d(&stats->released_photons, released);
  }

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  message(
      "Star %lld: tree reach %e (last coverage stop %e), carried %.17g, "
      "released %.17g.",
      si->id, cell_radiation_tree_reach(c), last_tree_reach,
      carried ? feedback_get_star_ionization_budget_total(si) : 0., released);
#endif
  return released;
}

/**
 * @brief Count a pass of a star holding gas that its one empty probe stopped
 * with photons left. Those photons escape.
 *
 * Under SWIFT_DEBUG_CHECKS_VERBOSE, also search once more out to the reach of
 * the walk and report a stall (runner_hii_is_stall()): a candidate the
 * stopped walk could have reached but did not.
 *
 * @param r The #runner thread.
 * @param c The working-level #cell holding the star.
 * @param si The #spart.
 * @param frame Rotation of the cone frame (row-major 3x3), or NULL for the
 * box axes.
 * @param radius The radius the walk stopped at.
 * @param covered_radius Radius inside which the pass visits every cell that
 * held gas at the last rebuild.
 * @param buffer Scratch neighbour buffer of #max_ngbs entries.
 */
static void runner_hii_count_one_probe_stop(struct runner *r, struct cell *c,
                                            struct spart *si,
                                            const double *frame,
                                            const float radius,
                                            const float covered_radius,
                                            struct hii_neighbor *buffer) {

  struct radiation_hii_coverage_stats *stats = &radiation_hii_coverage_local;
  const double photons_left = feedback_get_star_ionization_budget_total(si);
  atomic_add_d(&stats->one_probe_stops, 1.);
  atomic_add_d(&stats->one_probe_photons, photons_left);

#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
  const float max_search_radius = r->e->stars_properties->HII_max_search_radius;
  const float reach =
      runner_hii_held_search_radius(covered_radius, max_search_radius);
  if (radius >= reach) return;

  /* The walk found no candidate inside radius, so any found here lies
     beyond it. The gather only reads the gas. */
  struct hii_gather_state gather = {.r = r,
                                    .si = si,
                                    .search_radius = reach,
                                    .buffer = buffer,
                                    .max_size = max_ngbs,
                                    .count_found = 0,
                                    .mode = hii_gather_new_candidates,
                                    .ctx = NULL,
                                    .frame = frame};
  runner_hii_visit_neighbors(c, &gather);
  if (gather.count_found > 0 &&
      runner_hii_is_stall(radius, sqrtf(buffer[0].r2), covered_radius,
                          max_search_radius))
    message(
        "HII stall signature: star %lld stopped by the one-probe rule at %e "
        "below the walk reach %e; %d candidates inside it, nearest at %e; "
        "photons left %.17g",
        si->id, radius, reach, gather.count_found, sqrtf(buffer[0].r2),
        photons_left);
#endif
}

/**
 * @brief Count a pass that the tries or the ceiling stopped with photons left.
 * Those photons escape.
 *
 * @param si The #spart.
 */
static void runner_hii_count_growth_limit_stop(const struct spart *si) {
  struct radiation_hii_coverage_stats *stats = &radiation_hii_coverage_local;
  atomic_add_d(&stats->growth_limit_stops, 1.);
  atomic_add_d(&stats->growth_limit_photons,
               feedback_get_star_ionization_budget_total(si));
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
 * @brief Offer the gathered candidates to a star in distance order, each
 * pixel spending its own photon budget.
 *
 * @param r The #runner thread.
 * @param si The #spart.
 * @param buffer The gathered candidates, checked for sort order here.
 * @param count_found Number of entries in @p buffer.
 * @param dt_back Interval whose photons this pass issues.
 * @param dt_tag Gap since the star's previous pass, which tags outlive.
 * @param diag The #hii_pass_diagnostics.
 * @param claimed_count (in, out) Candidates claimed by the star, over all
 * pixels.
 * @param front_count (in, out) Fronts met: each empties its cone.
 */
static void runner_hii_claim_candidates(
    struct runner *r, struct spart *si, struct hii_neighbor *buffer,
    const int count_found, const double dt_back, const double dt_tag,
    struct hii_pass_diagnostics *diag, int *claimed_count, int *front_count) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;
  const struct feedback_props *feedback_props = e->feedback_props;
  const double time = e->time;

  runner_do_stars_hii_ionization_feedback_check_sort(buffer, count_found);

  const integertime_t ti_begin =
      get_integer_time_begin(e->ti_current - 1, si->time_bin);

  for (int k = 0; k < count_found; k++) {

    const int pixel = buffer[k].pixel;

    /* Skip, do not break: later entries may sit in a pixel with photons. */
    if (feedback_get_star_ionization_budget(si, pixel) <= 0.0) continue;

    struct part *pj = buffer[k].p;
    struct xpart *xpj = buffer[k].xp;
    const float r2 = buffer[k].r2;
    const int was_tagged_before = radiation_is_part_tagged_as_ionized(pj, xpj);

    const int outcome = feedback_iact_HII_ionization(
        si, pj, xpj, r2, pixel, phys_const, hydro_props, us, cosmo, cooling,
        feedback_props, ti_begin, time, dt_back, dt_tag);

    if (outcome == hii_candidate_front) ++(*front_count);

    /* The candidate starts untagged, so a tag with si->id means this star
       claimed it. A claim by another star tags a different id. An
       untagged result means a lost roll or a front, because
       feedback_hii_claim_part tags and spends budget together. */
    const int claimed_this_candidate =
        !was_tagged_before && radiation_is_part_tagged_as_ionized(pj, xpj) &&
        radiation_get_part_ionized_star_id(pj, xpj) == si->id;
    if (claimed_this_candidate) ++(*claimed_count);

    hii_pass_diagnostics_record_candidate(
        diag, pj, pixel, outcome, claimed_this_candidate, was_tagged_before);
  }
}

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
 * @param start_radius Radius of the first search
 * (runner_hii_claim_start_radius()), before the cut to @p covered_radius.
 * @param covered_radius Radius inside which the pass visits every cell that
 * held gas at the last rebuild; the search never goes beyond it.
 * @param last_tree_reach feedback_get_star_HII_carry_tree_reach() before the
 * budget of this pass opened.
 * @param count_held Gas the star already held at the start of the pass.
 * @param held_radius Distance to the star of the farthest gas it holds, 0 if
 * it holds none.
 * @param diag The #hii_pass_diagnostics.
 * @param summary (return) How the walk ended.
 */
static void runner_hii_claim_and_retry(
    struct runner *r, struct cell *c, struct spart *si, const double *frame,
    const double dt_back, const double dt_tag, const float start_radius,
    const float covered_radius, const float last_tree_reach,
    const int count_held, const float held_radius,
    struct hii_pass_diagnostics *diag, struct hii_claim_summary *summary) {

  struct engine *e = r->e;
  const struct feedback_props *feedback_props = e->feedback_props;

  const struct stars_props *star_props = e->stars_properties;
  const int max_retry_full_buffer = star_props->HII_max_retry_full_buffer;
  const int max_radius_expansion_tries =
      star_props->HII_max_radius_expansion_tries;
  const float radius_expansion_factor = star_props->HII_radius_expansion_factor;
  const float max_search_radius = star_props->HII_max_search_radius;

  struct hii_neighbor ngb_buffer[max_ngbs];

  /* A region whose held gas has flowed out since the last rebuild can start
     beyond the covered radius; the walk then starts at the covered radius. */
  if (start_radius > covered_radius)
    atomic_add_d(&radiation_hii_coverage_local.start_radius_cuts, 1.);

  /* A pass can leave photons unspent for two reasons, and each needs its own
   * retry:
   *
   * 1. The buffer filled up: the nearest-max_ngbs cut may have bumped out
   *    more candidates within the same radius. Retry at that radius
   *    (Stars:HII_max_retry_full_buffer).
   * 2. The buffer did not fill up: every candidate within the radius was
   *    found, so the radius is the limit and a retry would find nothing
   *    new. Expand it (Stars:HII_max_radius_expansion_tries).
   *
   * Without case 2 the radius would only grow at a later rebuild, through
   * h_hii or the drift of h_max. Both retries stay inside the loop of the
   * one star that needs them. */
  float dynamic_search_radius = min(start_radius, covered_radius);
  int num_retry_full_buffer = 0;
  int num_radius_expansions = 0;
  int num_empty_expansions = 0;

  /* Gas of other stars counts only where no earlier search of this walk
     reached: inside the held region it says nothing about the front. */
  float others_min_radius = held_radius;

  hii_pass_diagnostics_init(diag, si);

  while (1) {
    /* Candidates claimed in this attempt, over all pixels. See case 1 below. */
    int attempt_claimed_count = 0;
    /* Fronts met in this attempt: each empties its cone, so a retry gathers
       past it. */
    int attempt_front_count = 0;

    struct hii_gather_state gather = {
        .r = r,
        .si = si,
        .search_radius = dynamic_search_radius,
        .buffer = ngb_buffer,
        .max_size = max_ngbs,
        .count_found = 0,
        .mode = hii_gather_new_candidates,
        .ctx = NULL,
        .frame = frame,
        .r2_others_min = others_min_radius * others_min_radius,
        .count_held_by_others = 0};
    runner_hii_visit_neighbors(c, &gather);
    const int count_found = gather.count_found;
    const int count_met = runner_hii_probe_count_met(&gather);

    const int buffer_was_full = (count_found == max_ngbs);

    if (count_found > 0)
      runner_hii_claim_candidates(r, si, ngb_buffer, count_found, dt_back,
                                  dt_tag, diag, &attempt_claimed_count,
                                  &attempt_front_count);

    /* Fully exhausted: nothing more to do this pass. */
    if (feedback_get_star_ionization_budget_max(si) <= 0.0) break;

    if (buffer_was_full) {
      /* Case 1: retry only if this attempt claimed something or met a front.
       * Otherwise no budget or tag changed, so the retry would repeat the same
       * search, and the fixed (star, pixel, ti_begin) random draw of each
       * pixel would repeat the same outcome. */
      if (attempt_claimed_count == 0 && attempt_front_count == 0) break;
      if (num_retry_full_buffer >= max_retry_full_buffer) break;
      ++num_retry_full_buffer;
      num_empty_expansions = 0;
    } else {
      /* Case 2: expand the radius unless a limit of runner_hii_search_stop()
       * binds. */
      const enum hii_search_stop stop = runner_hii_search_stop(
          dynamic_search_radius, max_search_radius, covered_radius,
          num_radius_expansions, max_radius_expansion_tries, count_met,
          count_held, num_empty_expansions);

      if (stop == hii_search_stop_covered_radius) {
        /* Only the cell layout stops this star, so ask the next passes,
           and the rebuild check once the radius outgrows the wired
           neighbours, for the radius it would have tried next. The
           unspent photons are kept for the next pass if the tree reach grew
           since the last coverage stop, else released as escaped. That pass
           may run before any rebuild. */
        si->h_hii_request = runner_hii_requested_radius(dynamic_search_radius,
                                                        radius_expansion_factor,
                                                        max_search_radius) /
                            (radiation_search_radius_factor * kernel_gamma);
#ifdef SWIFT_DEBUG_CHECKS_VERBOSE
        message(
            "Star %lld: search stopped at the covered radius %e, "
            "requesting %e.",
            si->id, dynamic_search_radius, kernel_gamma * si->h_hii_request);
#endif
        runner_hii_carry_at_coverage_stop(c, si, feedback_props,
                                          last_tree_reach,
                                          &radiation_hii_coverage_local);
      }
      if (stop == hii_search_stop_one_probe)
        runner_hii_count_one_probe_stop(r, c, si, frame, dynamic_search_radius,
                                        covered_radius, ngb_buffer);
      if (stop == hii_search_stop_tries || stop == hii_search_stop_ceiling)
        runner_hii_count_growth_limit_stop(si);
      if (stop != hii_search_grow) break;

      num_empty_expansions = runner_hii_next_empty_probes(count_met, count_held,
                                                          num_empty_expansions);

      others_min_radius = dynamic_search_radius;
      dynamic_search_radius =
          min3(dynamic_search_radius * radius_expansion_factor,
               max_search_radius, covered_radius);
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

  /* Stamp the clamped age: HII_region_last_rebuild never holds a negative
     value, and compute_time() can return one for a newborn. */
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
 * @brief Run the HII pass of every active star of a working-level cell.
 *
 * Each star opens its photon budget, charges the gas it already holds, then
 * claims new candidates, tagging them as ionized through atomic operations.
 * A star with photons left retries at the same radius if the neighbour buffer
 * filled up, and at a larger radius otherwise (see
 * runner_hii_claim_and_retry()).
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

  /* The gas count is not tested here: see runner_hii_open_star_pass(). */
  if (c->stars.count == 0 || !cell_is_active_stars(c, e)) return;

  for (int i = 0; i < scount; i++) {

    struct spart *si = &sparts[i];

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

    /* The gather only visits linked cells whatever the search radius says,
       so no search of this pass goes beyond this radius. A star stopped
       there asks for a wider reach (h_hii_request). */
    const float covered_radius =
        cell_radiation_covered_radius(c, si->x, e->s->dim, e->s->periodic);

    float held_radius = 0.f;
    const int count_held = runner_hii_maintain_held_region(
        r, c, si, frame, pass.dt_back, pass.dt_tag, covered_radius,
        &held_radius);

    const float max_search_radius = e->stars_properties->HII_max_search_radius;
    const float start_radius = runner_hii_claim_start_radius(
        interaction_limit, held_radius, max_search_radius);

#ifdef SWIFT_DEBUG_CHECKS
    /* The held gas lies inside the reach of its gather, and the walk starts
       beyond it unless the ceiling binds. */
    const float held_reach =
        runner_hii_held_search_radius(covered_radius, max_search_radius);
    if (held_radius > held_reach * (1.f + 1e-6f))
      error("Star %lld: held gas at %e beyond the held-gas reach %e.", si->id,
            held_radius, held_reach);
    if (start_radius < interaction_limit ||
        start_radius < min(radiation_search_radius_factor * held_radius,
                           max_search_radius))
      error("Star %lld: claim walk starts at %e, inside the held gas at %e.",
            si->id, start_radius, held_radius);
#endif

    struct hii_pass_diagnostics diag;
    struct hii_claim_summary summary;
    runner_hii_claim_and_retry(r, c, si, frame, pass.dt_back, pass.dt_tag,
                               start_radius, covered_radius,
                               pass.last_tree_reach, count_held, held_radius,
                               &diag, &summary);

    hii_pass_diagnostics_print(&diag, e, si, count_held);
    runner_hii_close_star_pass(c, si, pass.star_age_safe, &summary, &diag,
                               count_held);
  } /* Loop over sparts */
}

/**
 * @brief Gather from a neighbour cell, recursing into its progeny.
 *
 * A cell out of reach along the sorting axis is skipped. A leaf is scanned
 * through its sort if it has a usable one, naively otherwise.
 *
 * @param ci The cell holding the star.
 * @param cj The neighbour cell holding gas candidates.
 * @param sid The sorting axis index between the pair.
 * @param flipped Is the cell order flipped relative to the axis direction?
 * @param shift The periodic wrapping shift vector applied to @p cj.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_do_stars_hii_ionization_feedback_branch(
    struct cell *ci, struct cell *cj, const int sid, const int flipped,
    const double shift[3], struct hii_gather_state *gather) {

  /* Skip a cell that is entirely out of reach. */
  if (!runner_hii_check_cell_can_be_reached(
          ci, cj, sid, flipped, shift, gather->si, gather->search_radius)) {
    return;
  }

  /* The self interaction already processed ci. */
  if (ci == cj) return;

  if (cj->split) {

    for (int k = 0; k < 8; k++) {
      if (cj->progeny[k] != NULL) {
        runner_do_stars_hii_ionization_feedback_branch(ci, cj->progeny[k], sid,
                                                       flipped, shift, gather);
      }
    }
  } else {

    lock_lock(&cj->hydro.extra_sort_lock);

    /* The sorted window only pre-filters: the r2 test decides. A projection
     * never exceeds the 3D distance, so the window is conservative along any
     * of the 13 sids and the claimed set does not depend on the sid, only the
     * speed does. Prefer the carried sid. If cj is not sorted along it, as for
     * most progeny leaves, use any sid cj is sorted along rather than the O(N)
     * naive search. The freshness gate (dx_max_sort_old) is per cell, so it
     * holds for any sid. */
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

    /* The sort state is read. */
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
 * Every gas particle of the cell is offered to the gather.
 *
 * @param c The #cell containing the gas particles.
 * @param gather The #hii_gather_state of the traversal.
 */
void runner_doself_stars_hii_ionization_feedback(
    struct cell *c, struct hii_gather_state *gather) {

  struct engine *e = gather->r->e;
  struct spart *si = gather->si;
  const int count = c->hydro.count;

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

    /* A foreign cell has no xparts array: xpj must be NULL there, not NULL
     * plus an offset (feedback_get_eligibility_temperature() owns the xp ==
     * NULL contract). */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = (xparts != NULL) ? &xparts[pjd] : NULL;

    runner_hii_offer_candidate(gather, pj, xpj, c, six, NULL, r2_max);
  } /* Loop in current cell */
}

/**
 * @brief Gather gas particles within the HII radius from a neighbouring cell,
 * offering every particle of it.
 *
 * Distances account for the periodic shift between @p ci and @p cj.
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

      /* xpj is NULL when cj has no xparts, as in the self loop. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = (xparts_j != NULL) ? &xparts_j[pjd] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, cj->loc, r2_max);
    } /* Loop in current cell */
  }
}

/**
 * @brief Gather gas particles within the HII radius from a neighbouring cell,
 * scanning only the window of its sort along @p sid that the star can reach.
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

  if (count_j == 0) return;

  /* Gas comes from cj, the star from ci: check each for what it supplies. */
  if (!cell_are_part_drifted(cj, e))
    error("Interacting undrifted cell (hydro).");
  if (!cell_are_spart_drifted(ci, e))
    error("Interacting undrifted cell (stars).");

  const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
  const float dx_max = cj->hydro.dx_max_sort + ci->stars.dx_max_sort;

  /* Projections and distances both use absolute gas positions, so
   * runner_hii_axis_window() puts the star in that frame. `flipped` only
   * selects the scan direction. */
  double pix[3], di_min, di_max;
  runner_hii_axis_window(si, shift, sid, search_radius, dx_max, pix, &di_min,
                         &di_max);

  const double dj_min = sort_j[0].d;
  const double dj_max = sort_j[count_j - 1].d;

  /* Skip a cell out of reach along the sorted axis. */
  if (di_max < dj_min || di_min > dj_max) {
    return;
  }

  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const float r2_max = search_radius * search_radius;
  const float six[3] = {(float)pix[0], (float)pix[1], (float)pix[2]};

  if (!flipped) {
    /* Scan forward until the particles pass di_max. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d <= di_max; pjd++) {

      if (sort_j[pjd].d < di_min) continue;

      /* xpj is NULL when cj has no xparts, as in the self loop. */
      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj =
          (xparts_j != NULL) ? &xparts_j[sort_j[pjd].i] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, NULL, r2_max);
    }
  } else {
    /* Scan backward until the particles fall below di_min. */
    for (int pjd = count_j - 1; pjd >= 0 && sort_j[pjd].d >= di_min; pjd--) {

      if (sort_j[pjd].d > di_max) continue;

      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      struct xpart *restrict xpj =
          (xparts_j != NULL) ? &xparts_j[sort_j[pjd].i] : NULL;

      runner_hii_offer_candidate(gather, pj, xpj, cj, six, NULL, r2_max);
    }
  }
}

/**
 * @brief Can a star's search radius reach any part of a neighbouring cell
 * along a sorting axis?
 *
 * @param ci The cell containing the star particle.
 * @param cj The neighboring gas cell to check.
 * @param sid The sorting axis index between the pair.
 * @param flipped Flag indicating if the cell order is flipped relative to the
 * axis direction.
 * @param shift The periodic wrapping shift vector.
 * @param si The star particle pointer.
 * @param search_radius The search radius around the star.
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

  /* The cell extent below uses the absolute geometry of cj, the frame of
   * runner_hii_axis_window() and of the sorted scan. A sign that depended on
   * flipped would offset the window by 2*shift and skip wrapped neighbours. */
  double pix[3], di_min, di_max;
  runner_hii_axis_window(si, shift, sid, search_radius, dx_max, pix, &di_min,
                         &di_max);

  /* The geometry, not the sort, bounds the cell: cj may be unsorted. */
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

  if (di_max < dj_min || di_min > dj_max) {
    return 0;
  }

  return 1;
}
