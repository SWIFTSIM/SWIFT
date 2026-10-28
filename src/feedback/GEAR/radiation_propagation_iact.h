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

/* The reach of a sender is its own kernel support (inner part) and, for an
   active receiver, the rest of the kernel support of the receiver (outer
   part). With a one-sided gradient loop, the gradient loop adds the outer
   part of a sender in the visit of the receiver. With a two-sided gradient
   loop the visit of the sender reaches every receiver, so the activity of the
   receiver plays no role. */
#ifdef EXTRA_HYDRO_LOOP_TYPE2
#define RADIATION_PACKET_RECEIVER_ACTIVE(active) 1
#define RADIATION_PACKET_GRADIENT_WANT RADIATION_PACKET_REACH_ALL
#else
#define RADIATION_PACKET_RECEIVER_ACTIVE(active) (active)
#define RADIATION_PACKET_GRADIENT_WANT RADIATION_PACKET_REACH_INNER
#endif

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
  radiation_packet_gradient_accumulate(dx_ij, hi, hj, 1,
                                       RADIATION_PACKET_REACH_ALL, pi, pj);
  radiation_packet_gradient_accumulate(dx, hj, hi, 1,
                                       RADIATION_PACKET_REACH_ALL, pj, pi);
}

/**
 * @brief Gradient-loop interaction between two particles (non-symmetric):
 * particle i is updated and gathers its inner part of the reach. If j is
 * active, its outer part of the reach, which this visit of i is the only one to
 * see, is added to the accumulators of j.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (only its outer part is added, if active).
 * @param a Current scale factor (unused, kept for interface conformance).
 * @param H Current Hubble parameter (unused, same reason).
 * @param pj_active Is j active in this step?
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H, const int pj_active) {

  if (pi->feedback_data.ISRF_packet_transport) {
    const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
    radiation_packet_gradient_accumulate(
        dx_ij, hi, hj, RADIATION_PACKET_RECEIVER_ACTIVE(pj_active),
        RADIATION_PACKET_GRADIENT_WANT, pi, pj);
  }

#ifndef EXTRA_HYDRO_LOOP_TYPE2
  if (pj_active && pj->feedback_data.ISRF_packet_transport)
    radiation_packet_gradient_accumulate(dx, hj, hi, 1,
                                         RADIATION_PACKET_REACH_OUTER, pj, pi);
#endif
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
  if (book_i_to_j)
    radiation_packet_book_share(dx_ij, hi, hj, 1, RADIATION_PACKET_REACH_ALL,
                                pi, pj, a, H);
  if (book_j_to_i)
    radiation_packet_book_share(dx, hj, hi, 1, RADIATION_PACKET_REACH_ALL, pj,
                                pi, a, H);
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
 * @param pj_active Is j active in this step? It decides whether the outer part
 * of the reach of i is booked.
 * @param reverse Book the outer part of the share of j into i. For a sender j
 * that is active and foreign: the rank of i is the only one that books it.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_packet_hop(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const float a,
                                   const float H, const int book,
                                   const int pj_active, const int reverse) {

  if (pi->feedback_data.ISRF_packet_transport) {
    radiation_packet_gather_field(r2, hi, pi, pj);
    if (book) {
      const float dx_ij[3] = {-dx[0], -dx[1], -dx[2]};
      radiation_packet_book_share(dx_ij, hi, hj,
                                  RADIATION_PACKET_RECEIVER_ACTIVE(pj_active),
                                  RADIATION_PACKET_REACH_ALL, pi, pj, a, H);
    }
  }

#ifndef EXTRA_HYDRO_LOOP_TYPE2
  if (reverse && pj->feedback_data.ISRF_packet_transport)
    radiation_packet_book_share(dx, hj, hi, 1, RADIATION_PACKET_REACH_OUTER, pj,
                                pi, a, H);
#endif
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
                                 const float H, const int pj_active) {}

__attribute__((always_inline)) INLINE static void runner_iact_isrf_packet_hop(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const int book_i_to_j, const int book_j_to_i) {}

__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_packet_hop(const float r2, const float dx[3],
                                   const float hi, const float hj,
                                   struct part *restrict pi,
                                   struct part *restrict pj, const float a,
                                   const float H, const int book,
                                   const int pj_active, const int reverse) {}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
