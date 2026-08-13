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
#ifndef SWIFT_RUNNER_RADIATION_FEEBDACK_H
#define SWIFT_RUNNER_RADIATION_FEEBDACK_H

/* Config parameters. */
#include <config.h>
#include <math.h>

#if defined(HAVE_CHEALPIX)
#include <chealpix.h>
#endif

/* Local headers. */
#include "error.h"
#include "hydro.h"
#include "stars.h"

/* Avoid cyclic inclusions. */
struct runner;
struct cell;

/* radiation_search_radius_factor lives in cell_stars.h: cell.h's rebuild
 * and split predicates need the same value and cannot depend on this
 * feedback-scheme header. */

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

void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer);
void runner_dosub_stars_hii_ionization_feedback(struct runner *r,
                                                struct cell *c,
                                                const float interaction_limit);
void runner_do_stars_hii_ionization_feedback_branch(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *ngb_buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx);

void runner_doself_stars_hii_ionization_feedback(
    struct runner *r, struct cell *c, struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx);
void runner_dopair_naive_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const double shift[3],
    struct spart *si, const float search_radius, struct hii_neighbor *buffer,
    int max_size, int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx);
void runner_dopair_stars_hii_ionization_feedback(
    struct runner *r, struct cell *ci, struct cell *cj, const int sid,
    const int flipped, const double shift[3], struct spart *si,
    const float search_radius, struct hii_neighbor *buffer, int max_size,
    int *count_found, enum hii_gather_mode mode,
    struct hii_maintenance_context *ctx);

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
 * @brief Compute the angular (HEALPix) pixel a gas candidate falls into,
 * as seen from the star.
 *
 * Hard single-pixel assignment (no fractional/multi-pixel overlap): a
 * candidate belongs to exactly one pixel, never split across several.
 *
 * @param dx Star-minus-particle position offset (dx = x_star - x_particle).
 * @param n_HII_pixels Number of active pixels for this star (1 = spherical,
 * HEALPix disabled).
 * @return Pixel index in [0, n_HII_pixels).
 */
__attribute__((always_inline)) INLINE static int runner_hii_get_pixel(
    const float dx[3], int n_HII_pixels) {

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
  const double dir[3] = {-(double)dx[0], -(double)dx[1], -(double)dx[2]};
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
#endif /* SWIFT_RUNNER_RADIATION_FEEBDACK_H */
