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
#ifndef SWIFT_RADIATION_HII_GEAR_H
#define SWIFT_RADIATION_HII_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_hii.h
 * @brief Physics of the HII regions that does not depend on the runner: the
 * claim of a candidate by a cone, the renewal of a held particle's tag, and
 * the cone frame and angular pixel of a candidate.
 */

/* Config parameters. */
#include <config.h>
#include <math.h>

#if defined(HAVE_CHEALPIX)
#include <chealpix.h>
#endif

/* Local headers. */
#include "error.h"
#include "inline.h"
#include "random.h"

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
        "radiation_init(): got non-integer or unrecognised nside=%d.",
        n_HII_pixels, nside);
#endif
  /* A particle exactly on the star gives a zero-length direction, which is
     undefined for vec2pix_ring64; fall back to pixel 0 (measure-zero in
     continuous positions, but must not read/write out of bounds). */
  if (dx[0] == 0.0f && dx[1] == 0.0f && dx[2] == 0.0f) return 0;

  /* Pixels are defined around the star, so the direction of interest is
     from the star to the particle; the negation of dx. vec2pix_ring64
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
      "n_HII_pixels > 1 but HAVE_CHEALPIX is not defined: this should be "
      "unreachable (radiation_init() validates this at startup).");
  return 0;
#endif
}

#endif /* SWIFT_RADIATION_HII_GEAR_H */
