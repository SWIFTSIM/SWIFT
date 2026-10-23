/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_RADIATION_FEEDBACK_H
#define SWIFT_RUNNER_RADIATION_FEEDBACK_H

/* Config parameters. */
#include <config.h>
#include <math.h>

#if defined(HAVE_CHEALPIX)
#include <chealpix.h>
#endif

/* Local headers. */
#include "error.h"
#include "feedback.h"
#include "hydro.h"
#include "minmax.h"
#include "random.h"
#include "stars.h"

/* Avoid cyclic inclusions. */
struct runner;
struct cell;
struct engine;

/* radiation_search_radius_factor lives in cell_stars.h: cell.h's rebuild
 * and split predicates need the same value and cannot depend on this
 * feedback-scheme header. */

/**
 * @brief Does this engine policy need the radiation subgrid task layer
 * (the per-region hii_ionization_feedback task and its super pointer)?
 *
 * True only when GEAR ionization feedback is compiled in and both the
 * stars and feedback policies are active.
 *
 * @param with_stars Is engine_policy_stars set?
 * @param with_feedback Is engine_policy_feedback set?
 */
#ifdef IONIZATION_FEEDBACK_LOOP
__attribute__((always_inline)) INLINE static int
feedback_radiation_subgrid_needed(int with_stars, int with_feedback) {
  return with_stars && with_feedback;
}
#else
__attribute__((always_inline)) INLINE static int
feedback_radiation_subgrid_needed(int with_stars, int with_feedback) {
  return 0;
}
#endif

/**
 * @brief Does this engine policy need the radiation gather task layer
 * (radiation_in/radiation_out self and pair tasks)?
 *
 * True only when GEAR ionization feedback is compiled in and both the
 * hydro and stars policies are active. Unlike
 * #feedback_radiation_subgrid_needed, the feedback policy is not required:
 * the geometric gather task graph is built from hydro + stars alone, with
 * the ionization task itself gated separately.
 *
 * @param with_hydro Is engine_policy_hydro set?
 * @param with_stars Is engine_policy_stars set?
 */
#ifdef IONIZATION_FEEDBACK_LOOP
__attribute__((always_inline)) INLINE static int
feedback_radiation_gather_tasks_needed(int with_hydro, int with_stars) {
  return with_hydro && with_stars;
}
#else
__attribute__((always_inline)) INLINE static int
feedback_radiation_gather_tasks_needed(int with_hydro, int with_stars) {
  return 0;
}
#endif

/**
 * @brief Result of offering one untagged candidate to a cone.
 */
enum hii_candidate_outcome {

  /*! The candidate is claimed and its cost charged. */
  hii_candidate_claimed = 0,

  /*! Already ionized: nothing charged, the cone goes on. */
  hii_candidate_passed = 1,

  /*! Not claimed: the candidate is an ionization front that absorbs photons. */
  hii_candidate_front = 2,
};

/**
 * @brief Decide whether a cone claims an untagged candidate.
 *
 * The first candidate a cone cannot afford ends the cone: it is claimed on a
 * won roll, and it is a front that absorbs the rest of the budget otherwise
 * (Hopkins et al. 2018, Appendix E). The win probability is budget / cost, so
 * the expected spend equals the budget. The roll is the one-per-cone draw of
 * the caller, so a cone rolls once per pass.
 *
 * @param deterministic 1 if the boundary candidate is always claimed.
 * @param budget Photons left in the cone.
 * @param cost Photons the candidate costs: its neutral atoms plus one
 * interval of recombinations.
 * @param random_number The cone's draw in [0, 1), read only if the candidate
 * is unaffordable.
 * @return The #hii_candidate_outcome.
 */
__attribute__((always_inline)) INLINE static enum hii_candidate_outcome
runner_hii_candidate_outcome(const int deterministic, const double budget,
                             const double cost, const float random_number) {

  if (cost <= budget) return hii_candidate_claimed;
  if (deterministic) return hii_candidate_claimed;

  const float proba = budget / cost;
  if (random_number <= proba) return hii_candidate_claimed;

  return hii_candidate_front;
}

/**
 * @brief Decide whether a held particle whose upkeep was paid in part keeps
 * its ionization tag.
 *
 * The tag is kept with probability paid / owed, so the expected photons spent
 * on the particle equal what it costs. A fully paid upkeep always renews.
 *
 * @param paid Photons charged for the upkeep.
 * @param owed Photons the upkeep costs.
 * @param random_number A draw in [0, 1), read only if @p paid < @p owed.
 * @return 1 if the tag is renewed.
 */
__attribute__((always_inline)) INLINE static int runner_hii_upkeep_renewed(
    const double paid, const double owed, const double random_number) {
  if (paid >= owed) return 1;
  return random_number < paid / owed;
}

/**
 * @brief Limit that ends the radius growth of an HII pass.
 */
enum hii_search_stop {

  /*! No limit binds: the radius grows again. */
  hii_search_grow = 0,

  /*! Stars:HII_max_radius_expansion_tries is used up. */
  hii_search_stop_tries = 1,

  /*! The radius is at Stars:HII_max_search_radius. */
  hii_search_stop_ceiling = 2,

  /*! The radius is at the neighbour cell margin, and nothing else binds. */
  hii_search_stop_margin = 3,

  /*! A star that holds gas has had its one empty probe. */
  hii_search_stop_one_probe = 4,
};

/**
 * @brief Decide whether an HII pass whose search found no more candidates
 * grows its radius, and if not, which limit stops it.
 *
 * A star already holding an equilibrium region gets one empty probe past it
 * and no more: a pixel left with a sliver of budget never ends the pass, so
 * it would otherwise spend every expansion on ever larger empty volumes. A
 * star holding nothing keeps its full configured reach, the case the
 * expansion exists for (a fresh star, or one in a cavity whose nearest gas is
 * several expansions away). The neighbour cell margin is reported only when
 * no other limit binds, so it marks a star that the cell layout alone holds
 * back.
 *
 * @param radius The radius just searched.
 * @param max_search_radius Stars:HII_max_search_radius.
 * @param max_reachable_radius The neighbour cell margin of this pass.
 * @param num_expansions Expansions done so far in this pass.
 * @param max_expansions Stars:HII_max_radius_expansion_tries.
 * @param count_found Candidates found at @p radius.
 * @param count_held Gas the star already held at the start of the pass.
 * @param num_empty_expansions Consecutive empty probes of a star holding gas.
 */
__attribute__((always_inline)) INLINE static enum hii_search_stop
runner_hii_search_stop(const float radius, const float max_search_radius,
                       const float max_reachable_radius,
                       const int num_expansions, const int max_expansions,
                       const int count_found, const int count_held,
                       const int num_empty_expansions) {

  if (num_expansions >= max_expansions) return hii_search_stop_tries;
  if (radius >= max_search_radius) return hii_search_stop_ceiling;
  const int probe_used =
      count_found == 0 && count_held > 0 && num_empty_expansions >= 1;
  if (probe_used) return hii_search_stop_one_probe;
  if (radius >= max_reachable_radius) return hii_search_stop_margin;
  return hii_search_grow;
}

/**
 * @brief Search radius a pass stopped by the neighbour cell margin asks for:
 * the next radius it would have tried.
 *
 * @param radius The radius just searched.
 * @param expansion_factor Stars:HII_radius_expansion_factor.
 * @param max_search_radius Stars:HII_max_search_radius.
 */
__attribute__((always_inline)) INLINE static float runner_hii_requested_radius(
    const float radius, const float expansion_factor,
    const float max_search_radius) {
  const float next = radius * expansion_factor;
  return next < max_search_radius ? next : max_search_radius;
}

/* Per-pass capacity of the HII neighbour search buffer (stack-allocated,
 * compile-time constant by design -- see runner_dosub_stars_hii_ionization_
 * feedback()). Sized to comfortably cover a single rebuild's typical new
 * shell of gas at gas_mass=0.01-resolution test problems without needing
 * many retries (empirically ~a few thousand particles per shell in
 * StromgrenSphere); denser regions (e.g. clumps in a cosmological zoom-in)
 * rely on Stars:HII_max_retry_full_buffer (struct stars_props, runtime-
 * configurable) to cover the rest across several passes at the same
 * search radius. */
#define max_ngbs 1024

/**
 * @brief Temporary structure to store gas particles found within the HII
 * interaction radius of a star.
 *
 * This is used to gather and sort potential ionization candidates before
 * performing the feedback.
 */
struct hii_neighbor {

  /*! Squared distance between the star and the gas particle. */
  float r2;

  /*! Pointer to the gas particle data. */
  struct part *p;

  /*! Pointer to the gas particle extra data. */
  struct xpart *xp;

  /*! Angular pixel this particle was assigned to (always 0 while angular
      splitting is disabled). */
  int pixel;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Pointer to the cell this particle belongs to */
  struct cell *c;
#endif
};

/**
 * @brief Which gas particles a neighbour traversal is looking for.
 *
 * The two phases of a rebuild pass share one spatial search but want
 * disjoint particle sets, so the traversal is parameterised rather than
 * duplicated.
 */
enum hii_gather_mode {

  /*! Not-yet-ionized gas, buffered and sorted by distance so the nearest
      candidates get first claim on the star's photons. */
  hii_gather_new_candidates,

  /*! Gas this star already ionized, processed as it is found: an unbounded
      set that must be visited in full, so it cannot go through the
      fixed-capacity buffer. */
  hii_gather_already_ionized,
};

/**
 * @brief Everything feedback_iact_HII_maintain_ionized_part() needs, carried
 * down the traversal to the point where already-ionized gas is found.
 *
 * Only populated for #hii_gather_already_ionized; NULL otherwise.
 */
struct hii_maintenance_context {
  const struct phys_const *phys_const;
  const struct hydro_props *hydro_props;
  const struct unit_system *us;
  const struct cosmology *cosmo;
  const struct cooling_function_data *cooling;

  /*! Current simulation time. */
  double time;

  /*! Time elapsed since this star's previous HII rebuild pass. */
  double dt_back;

  /*! (return) Largest squared distance at which this star still holds gas
      tagged and owned, i.e. the region's extent this pass. Independent of
      what the budget could actually pay for, so it does not depend on
      traversal order. */
  float r2_max_maintained;

  /*! (return) Photons charged to balance recombination this pass. Filled in
      by the caller after the traversal returns and the buffered candidates
      below have been sorted and charged in (r2, id) order; 0 during the
      traversal itself. */
  double photons_charged;

  /*! Runner performing this pass. Already-ionized candidates found during
      the traversal are appended to r->hii_maintenance_buffer instead of
      being charged immediately, so the whole set can be sorted into
      ascending (r2, id) order first. */
  struct runner *r;

  /*! Number of candidates appended to r->hii_maintenance_buffer so far this
      pass. */
  int buffer_count;
};

int runner_hii_check_cell_can_be_reached(const struct cell *ci,
                                         const struct cell *cj, const int sid,
                                         const int flipped,
                                         const double shift[3],
                                         const struct spart *si,
                                         const float search_radius);

/**
 * @brief Top-level steering function for HII ionization feedback, called
 * unconditionally from the runner's task dispatch (runner_main.c).
 *
 * A no-op when GEAR ionization feedback is compiled out: the task type
 * still exists in the enum, but with_radiation_tasks
 * (#feedback_radiation_gather_tasks_needed) never creates the task, so
 * this is dead code reachable only by construction, not at runtime.
 */
#ifdef IONIZATION_FEEDBACK_LOOP
void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer);
#else
__attribute__((always_inline)) INLINE static void
runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                        int timer) {}
#endif

/**
 * @brief Print startup warnings about the HII search-radius / box-size
 * relationship, called unconditionally from engine_config().
 *
 * A no-op when GEAR ionization feedback is compiled out.
 *
 * @param e The #engine.
 */
#ifdef IONIZATION_FEEDBACK_LOOP
void feedback_radiation_startup_diagnostics(const struct engine *e);
#else
__attribute__((always_inline)) INLINE static void
feedback_radiation_startup_diagnostics(const struct engine *e) {}
#endif

void runner_dosub_stars_hii_ionization_feedback(struct runner *r,
                                                struct cell *c,
                                                const float interaction_limit);
void runner_do_stars_hii_ionization_feedback_branch(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *ngb_buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx, const double *frame);

void runner_doself_stars_hii_ionization_feedback(
    struct runner *r, struct cell *c, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx, const double *frame);
void runner_dopair_naive_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const double shift[3],
    struct spart *si, const float search_radius, struct hii_neighbor *buffer,
    int max_size, int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx, const double *frame);
void runner_dopair_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx, const double *frame);

/**
 * @brief Shift that brings the neighbour cell next to the star's cell.
 *
 * @param dim The size of the simulation box.
 * @param periodic Is the box periodic?
 * @param loc_c Lower corner of the star's cell.
 * @param loc_cj Lower corner of the neighbour cell.
 * @param shift (return) Vector to add to the neighbour's coordinates.
 */
__attribute__((always_inline)) INLINE static void runner_hii_get_pair_shift(
    const double dim[3], const int periodic, const double loc_c[3],
    const double loc_cj[3], double shift[3]) {

  for (int k = 0; k < 3; k++) {
    const double dx = loc_cj[k] - loc_c[k];
    if (periodic && dx < -dim[k] / 2)
      shift[k] = dim[k];
    else if (periodic && dx > dim[k] / 2)
      shift[k] = -dim[k];
    else
      shift[k] = 0.0;
  }
}

/**
 * @brief (r2, id) < (r2, id): is @p a_r2/@p a_id strictly before
 * @p b_r2/@p b_id in the buffer's sort order?
 *
 * id is the tiebreak: r2 alone is not a total order (distinct particles can
 * sit at the identical squared distance, e.g. a symmetric IC), and the
 * insertion order must be reproducible across runs and threads.
 */
__attribute__((always_inline)) INLINE static int runner_hii_r2_id_before(
    float a_r2, long long a_id, float b_r2, long long b_id) {
  return (a_r2 < b_r2) || (a_r2 == b_r2 && a_id < b_id);
}

/**
 * @brief Maintain a buffer sorted in ascending (r2, id) order by inserting a
 * new neighbor at the correct position.
 *
 * If the buffer is full, it replaces the furthest element if the new one
 * sorts before it.
 */
__attribute__((always_inline)) INLINE static void runner_hii_buffer_insert(
    struct hii_neighbor *buffer, int max_size, int *count_found, float r2,
    struct part *p, struct xpart *xp, struct cell *c, int pixel) {

  /* Case A: Buffer is not yet full */
  if (*count_found < max_size) {
    int i = *count_found - 1;
    /* Shift elements to make room (standard insertion) */
    while (i >= 0 &&
           runner_hii_r2_id_before(r2, p->id, buffer[i].r2, buffer[i].p->id)) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
    buffer[i + 1].pixel = pixel;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
    (*count_found)++;
  }
  /* Case B: Buffer is full, check if new particle sorts before the furthest
   */
  else if (runner_hii_r2_id_before(r2, p->id, buffer[max_size - 1].r2,
                                   buffer[max_size - 1].p->id)) {
    int i = max_size - 2;
    /* Shift elements to replace the furthest */
    while (i >= 0 &&
           runner_hii_r2_id_before(r2, p->id, buffer[i].r2, buffer[i].p->id)) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
    buffer[i + 1].pixel = pixel;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
  }
}

/**
 * @brief Index of the cone frame a star uses at one of its passes.
 *
 * @param pass_index Number of passes the star has made before this one.
 * @param redraw_period Passes per frame, or 0 for a frame kept for the star's
 * life.
 * @return pass_index / redraw_period, or 0 when @a redraw_period <= 0.
 */
__attribute__((always_inline)) INLINE static int64_t
runner_hii_cone_frame_epoch(const int64_t pass_index, const int redraw_period) {
  if (redraw_period <= 0) return 0;
  return pass_index / redraw_period;
}

/**
 * @brief Build the rotation of a star's HEALPix cone frame from its ID and
 * the frame epoch.
 *
 * Shoemake's uniform random unit quaternion, built from three numbers drawn
 * with random_unit_interval(), gives a rotation uniform on SO(3). The result
 * depends only on the star ID, the epoch and the compile-time random seed, so
 * it is the same on every rank, after a restart and in any processing order.
 *
 * @param id The star ID.
 * @param epoch The frame epoch (see #runner_hii_cone_frame_epoch).
 * @param R (return) The rotation matrix, row-major.
 */
__attribute__((always_inline)) INLINE static void runner_hii_cone_frame_from_id(
    const int64_t id, const int64_t epoch, double R[9]) {

  const integertime_t ti = 3 * (integertime_t)epoch;
  const double u1 =
      random_unit_interval(id, ti, random_number_HII_frame_rotation);
  const double u2 =
      random_unit_interval(id, ti + 1, random_number_HII_frame_rotation);
  const double u3 =
      random_unit_interval(id, ti + 2, random_number_HII_frame_rotation);

  const double a = sqrt(1.0 - u1);
  const double b = sqrt(u1);
  const double qx = a * sin(2.0 * M_PI * u2);
  const double qy = a * cos(2.0 * M_PI * u2);
  const double qz = b * sin(2.0 * M_PI * u3);
  const double qw = b * cos(2.0 * M_PI * u3);

  R[0] = 1.0 - 2.0 * (qy * qy + qz * qz);
  R[1] = 2.0 * (qx * qy - qz * qw);
  R[2] = 2.0 * (qx * qz + qy * qw);
  R[3] = 2.0 * (qx * qy + qz * qw);
  R[4] = 1.0 - 2.0 * (qx * qx + qz * qz);
  R[5] = 2.0 * (qy * qz - qx * qw);
  R[6] = 2.0 * (qx * qz - qy * qw);
  R[7] = 2.0 * (qy * qz + qx * qw);
  R[8] = 1.0 - 2.0 * (qx * qx + qy * qy);
}

/**
 * @brief Compute the angular (HEALPix) pixel a gas candidate falls into,
 * as seen from the star.
 *
 * Hard single-pixel assignment (no fractional/multi-pixel overlap): a
 * candidate belongs to exactly one pixel, never split across several.
 *
 * @param dx Star-minus-particle position offset (dx = x_star - x_particle).
 * @param n_HII_pixels Number of active pixels for this star (1 = spherical,
 * HEALPix disabled).
 * @param frame Rotation of the cone frame (row-major 3x3, see
 * #runner_hii_cone_frame_from_id), or NULL for cones fixed to the box axes.
 * @return Pixel index in [0, n_HII_pixels).
 */
__attribute__((always_inline)) INLINE static int runner_hii_get_pixel(
    const float dx[3], int n_HII_pixels, const double *frame) {

  if (n_HII_pixels <= 1) return 0;

#if defined(HAVE_CHEALPIX)
  /* n_HII_pixels = 12*nside^2 for every value radiation_init() accepts
     (src/feedback/GEAR/radiation.c, validated against the build's
     HII_MAX_ANGULAR_PIXELS); invert it to recover nside for
     vec2pix_ring64, which takes nside, not a pixel count. Only exact
     values reach here, so this is exact integer arithmetic, not a
     rounded approximation, for any nside the build supports. */
  const int nside = (int)lround(sqrt(n_HII_pixels / 12.0));
#ifdef SWIFT_DEBUG_CHECKS
  if (n_HII_pixels != 12 * nside * nside)
    error(
        "n_HII_pixels=%d does not correspond to any nside supported by "
        "radiation_init() -- got non-integer or unrecognised nside=%d.",
        n_HII_pixels, nside);
#endif
  /* A particle exactly on the star gives a zero-length direction, which is
     undefined for vec2pix_ring64 -- fall back to pixel 0 (measure-zero in
     continuous positions, but must not read/write out of bounds). */
  if (dx[0] == 0.0f && dx[1] == 0.0f && dx[2] == 0.0f) return 0;

  /* Pixels are defined around the star, so the direction of interest is
     from the star to the particle -- the negation of dx. vec2pix_ring64
     does not require a normalized vector. */
  double dir[3] = {-(double)dx[0], -(double)dx[1], -(double)dx[2]};
  if (frame != NULL) {
    const double d0 = dir[0], d1 = dir[1], d2 = dir[2];
    dir[0] = frame[0] * d0 + frame[1] * d1 + frame[2] * d2;
    dir[1] = frame[3] * d0 + frame[4] * d1 + frame[5] * d2;
    dir[2] = frame[6] * d0 + frame[7] * d1 + frame[8] * d2;
  }
  int64_t ipix = 0;
  vec2pix_ring64(nside, dir, &ipix);
  return (int)ipix;
#else
  error(
      "n_HII_pixels > 1 but HAVE_CHEALPIX is not defined -- this should be "
      "unreachable (radiation_init() validates this at startup).");
  return 0;
#endif
}

/**
 * @brief Verify that the gathered HII neighbor buffer is properly sorted in
 * ascending (r2, id) order, matching #runner_hii_buffer_insert's tiebreak.
 *
 * @param ngb_buffer The array of gathered #hii_neighbor structures.
 * @param count_found The number of valid neighbors stored in the buffer.
 */
__attribute__((always_inline)) INLINE static void
runner_do_stars_hii_ionization_feedback_check_sort(
    const struct hii_neighbor *ngb_buffer, const int count_found) {
#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < count_found - 1; k++) {
    /* Strict: (r2, id) is a total order over distinct particles (unique
       ids), so two adjacent entries failing this means a genuine ordering
       bug or a duplicate particle. */
    if (!runner_hii_r2_id_before(ngb_buffer[k].r2, ngb_buffer[k].p->id,
                                 ngb_buffer[k + 1].r2,
                                 ngb_buffer[k + 1].p->id)) {
      error(
          "HII neighbor buffer not properly (r2, id)-sorted! "
          "Index %d (r2=%e, id=%lld) does not sort before index %d "
          "(r2=%e, id=%lld).",
          k, ngb_buffer[k].r2, ngb_buffer[k].p->id, k + 1, ngb_buffer[k + 1].r2,
          ngb_buffer[k + 1].p->id);
    }
  }
#endif
}
#endif /* SWIFT_RUNNER_RADIATION_FEEDBACK_H */
