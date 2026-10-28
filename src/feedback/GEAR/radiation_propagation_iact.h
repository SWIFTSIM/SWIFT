/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H
#define SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_propagation_iact.h
 * @brief Gas-gas gradient- and force-loop hooks of the packet transport of the
 * ISRF. See theory/GEAR/Radiation/04_isrf_packet_transport.tex.
 */

#include "radiation.h"
#include "radiation_packet.h"

#ifdef GEAR_SUBGRID_RADIATION_ISRF

/**
 * @brief Gradient-loop interaction between two particles (symmetric): the
 * gradient accumulators of the packet transport of both are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor (unused, kept for interface conformance).
 * @param H Current Hubble parameter (unused, same reason).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  if (!pi->feedback_data.ISRF_packet_transport) return;

  const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
  radiation_packet_gradient_accumulate(dx_ij, hi, pi, pj);
  radiation_packet_gradient_accumulate(dx, hj, pj, pi);
}

/**
 * @brief Gradient-loop interaction between two particles (non-symmetric):
 * only particle i's accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor (unused, kept for interface conformance).
 * @param H Current Hubble parameter (unused, same reason).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H) {

  if (!pi->feedback_data.ISRF_packet_transport) return;

  const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
  radiation_packet_gradient_accumulate(dx_ij, hi, pi, pj);
}

/**
 * @brief Force-loop packet hop between two particles both updated by this
 * visit: each books its share into the other and gathers the other's field.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param book_i_to_j Book the share of i into j (j is local).
 * @param book_j_to_i Book the share of j into i (i is local).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_packet_hop(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int book_i_to_j, const int book_j_to_i) {

  if (!pi->feedback_data.ISRF_packet_transport) return;

  const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
  if (book_i_to_j) radiation_packet_book_share(dx_ij, hi, pi, pj, a, H);
  if (book_j_to_i) radiation_packet_book_share(dx, hj, pj, pi, a, H);
  radiation_packet_gather_field(r2, hi, pi, pj);
  radiation_packet_gather_field(r2, hj, pj, pi);
}

/**
 * @brief Force-loop packet hop from the updated particle i into j; j's own
 * field joins the kernel mean of i.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Sender, active.
 * @param pj Receiver, active or not.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param book Book the share (j is local and this is its booking visit).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_packet_hop(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const float a,
                                   const float H, const int book) {

  if (!pi->feedback_data.ISRF_packet_transport) return;

  radiation_packet_gather_field(r2, hi, pi, pj);
  if (!book) return;
  const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
  radiation_packet_book_share(dx_ij, hi, pi, pj, a, H);
}

#else /* GEAR_SUBGRID_RADIATION_ISRF */

/* Without the part: the hydro loops keep their calls, which do nothing. */

__attribute__((always_inline)) INLINE static void runner_iact_isrf_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H) {}

__attribute__((always_inline)) INLINE static void runner_iact_isrf_packet_hop(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int book_i_to_j, const int book_j_to_i) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_packet_hop(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const float a,
                                   const float H, const int book) {}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
