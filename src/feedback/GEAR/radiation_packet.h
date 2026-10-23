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
#ifndef SWIFT_RADIATION_PACKET_GEAR_H
#define SWIFT_RADIATION_PACKET_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_packet.h
 * @brief Packet transport of the ISRF bands: each active gas particle hands
 * its packet to its downstream neighbours, once per own step.
 *
 * The gradient loop sums the sender's transmission weights `S`, the force
 * loop books each share `w/S` into the receiver. Both loops call the one
 * out-of-line #radiation_packet_transmission_weight, so the shares of a
 * sender sum to 1 up to rounding. A share is booked by the rank that owns the
 * receiver, once.
 */

#include "error.h"
#include "inline.h"
#include "part.h"

#include <math.h>

#ifdef GEAR_SUBGRID_RADIATION_ISRF

/*! Cosine below which the transmission weight is exactly zero. */
#define RADIATION_PACKET_MU_MIN 1e-3f

/*! Optical depth above which a share is fully absorbed. */
#define RADIATION_PACKET_MAX_OPTICAL_DEPTH 80.f

float radiation_packet_transmission_weight(const float dx[3], const float h_s,
                                           const float n[3], const float power);
void radiation_packet_end_gradient(struct part *p);
void radiation_packet_end_force(struct part *p);

/**
 * @brief Does the particle hold no energy in any band?
 *
 * @param pk The packet.
 */
__attribute__((always_inline)) INLINE static int radiation_packet_is_empty(
    const struct feedback_isrf_packet_data *pk) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    if (pk->band[b].energy != 0.f) return 0;
  return 1;
}

/**
 * @brief Gradient loop: add one downstream neighbour to the sender's weight
 * sum and forward advance.
 *
 * @param dx Comoving vector from the sender to the neighbour.
 * @param h_s Comoving smoothing length of the sender.
 * @param s The sender.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_gradient_accumulate(const float dx[3], const float h_s,
                                     struct part *restrict s) {

  struct feedback_isrf_packet_data *pk = &s->feedback_data.isrf_packet;
  if (radiation_packet_is_empty(pk)) return;

  const float w = radiation_packet_transmission_weight(dx, h_s, pk->direction,
                                                       pk->cone_power);
  if (w == 0.f) return;

  pk->weight_sum += w;
  pk->forward_advance += w * (dx[0] * pk->direction[0] +
                              dx[1] * pk->direction[1] +
                              dx[2] * pk->direction[2]);
}

/**
 * @brief Force loop: book the sender's share for one receiver.
 *
 * The share is attenuated by the receiver's dust rate over the increase of
 * the distance to the sender's origin. Energy and photon number lose the
 * same fraction, so the mean photon energy of a band is unchanged.
 *
 * @param dx Comoving vector from the sender to the receiver.
 * @param h_s Comoving smoothing length of the sender.
 * @param s The sender.
 * @param r The receiver.
 * @param a Current scale factor.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_book_share(
    const float dx[3], const float h_s, const struct part *restrict s,
    struct part *restrict r, const float a) {

  const struct feedback_isrf_packet_data *pk_s = &s->feedback_data.isrf_packet;
  if (radiation_packet_is_empty(pk_s)) return;

  const float w = radiation_packet_transmission_weight(
      dx, h_s, pk_s->direction, pk_s->cone_power);
  if (w == 0.f) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (!(pk_s->weight_sum > 0.f))
    error("Particle %lld books a share with weight sum %e.", s->id,
          pk_s->weight_sum);
#endif
  const float fraction = w / pk_s->weight_sum;

  /* Offset from the receiver to the sender's origin. */
  const float d_s = pk_s->origin_distance;
  const float offset[3] = {-dx[0] - d_s * pk_s->direction[0],
                           -dx[1] - d_s * pk_s->direction[1],
                           -dx[2] - d_s * pk_s->direction[2]};
  const float d_r = sqrtf(offset[0] * offset[0] + offset[1] * offset[1] +
                          offset[2] * offset[2]);
  float path = d_r - d_s;
  if (path < 0.f) path = 0.f;
  const float path_phys = a * path;

  struct feedback_isrf_packet_data *pk_r = &r->feedback_data.isrf_packet;
  float delivered_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const struct feedback_isrf_packet_band_data *src = &pk_s->band[b];
    struct feedback_isrf_packet_band_data *dst = &pk_r->band[b];

    const float tau = r->feedback_data.isrf_operator[b].kappa * path_phys;
    float absorbed_fraction = 1.f;
    if (tau < RADIATION_PACKET_MAX_OPTICAL_DEPTH)
      absorbed_fraction = -expm1f(-tau);

    const float share_energy = src->energy * fraction;
    const float share_photons = src->photons * fraction;
    const float absorbed_energy = share_energy * absorbed_fraction;
    const float absorbed_photons = share_photons * absorbed_fraction;
    const float transmitted_energy = share_energy - absorbed_energy;

    dst->inbox_energy += transmitted_energy;
    dst->inbox_photons += share_photons - absorbed_photons;
    dst->absorbed_energy += (double)absorbed_energy;
    dst->absorbed_photons += (double)absorbed_photons;
    delivered_energy += transmitted_energy;
  }

  for (int k = 0; k < 3; k++)
    pk_r->inbox_origin_offset[k] += delivered_energy * offset[k];

#ifdef SWIFT_DEBUG_CHECKS
  pk_r->shares_received++;
#endif
}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_PACKET_GEAR_H */
