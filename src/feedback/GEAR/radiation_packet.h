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
#include <string.h>

#ifdef GEAR_SUBGRID_RADIATION_ISRF

/*! Cosine below which the transmission weight is exactly zero. */
#define RADIATION_PACKET_MU_MIN 1e-3f

/*! Optical depth above which a share is fully absorbed. */
#define RADIATION_PACKET_MAX_OPTICAL_DEPTH 80.f

/*! Photon number (as `N E_lo`) below which a share is booked as absorbed:
    near the float minimum 1.2e-38 the number would flush to zero before the
    energy, and the mean photon energy would leave the band. */
#define RADIATION_PACKET_MIN_PHOTONS 1e-35f

float radiation_packet_transmission_weight(const float dx[3], const float h_s,
                                           const float n[3], const float power);
int radiation_packet_emission_cone(const long long star_id,
                                   const integertime_t ti_begin,
                                   const float dx[3]);
void radiation_packet_merge(struct part *p, const integertime_t ti_current,
                            const integertime_t ti_step,
                            const double time_base);
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
  pk->forward_advance +=
      w * (dx[0] * pk->direction[0] + dx[1] * pk->direction[1] +
           dx[2] * pk->direction[2]);
}

/**
 * @brief Deliver one share to a receiver: dust absorption over the path, then
 * the inbox, or the reservoir when the source has the coarser step.
 *
 * Energy and photon number lose the same fraction, so the mean photon energy
 * of a band is unchanged.
 *
 * @param r The receiver.
 * @param energy Energy of the share per band, before absorption.
 * @param photons Photon number of the share per band, as `N E_lo`.
 * @param offset Comoving offset from the receiver to the share's origin.
 * @param path_phys Physical path length over which the share is absorbed.
 * @param to_reservoir Is the source's step coarser than the receiver's?
 * @param source_step_end_ti End of the source's step.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_deliver(
    struct part *restrict r, const float energy[ISRF_OPERATOR_COUNT],
    const float photons[ISRF_OPERATOR_COUNT], const float offset[3],
    const float path_phys, const int to_reservoir,
    const integertime_t source_step_end_ti) {

  struct feedback_isrf_packet_data *pk_r = &r->feedback_data.isrf_packet;
  float delivered_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *dst = &pk_r->band[b];

    const float tau = r->feedback_data.isrf_operator[b].kappa * path_phys;
    float absorbed_fraction = 1.f;
    if (tau < RADIATION_PACKET_MAX_OPTICAL_DEPTH)
      absorbed_fraction = -expm1f(-tau);

    float absorbed_energy = energy[b] * absorbed_fraction;
    float absorbed_photons = photons[b] * absorbed_fraction;
    float transmitted_energy = energy[b] - absorbed_energy;
    float transmitted_photons = photons[b] - absorbed_photons;
    if (transmitted_photons < RADIATION_PACKET_MIN_PHOTONS) {
      absorbed_energy = energy[b];
      absorbed_photons = photons[b];
      transmitted_energy = 0.f;
      transmitted_photons = 0.f;
    }

    if (to_reservoir) {
      dst->reservoir_energy += transmitted_energy;
      dst->reservoir_photons += transmitted_photons;
    } else {
      dst->inbox_energy += transmitted_energy;
      dst->inbox_photons += transmitted_photons;
    }
    dst->absorbed_energy += (double)absorbed_energy;
    dst->absorbed_photons += (double)absorbed_photons;
    delivered_energy += transmitted_energy;
  }

  if (to_reservoir) {
    for (int k = 0; k < 3; k++)
      pk_r->reservoir_origin_offset[k] += delivered_energy * offset[k];
    if (source_step_end_ti > pk_r->reservoir_end_ti)
      pk_r->reservoir_end_ti = source_step_end_ti;
  } else {
    for (int k = 0; k < 3; k++)
      pk_r->inbox_origin_offset[k] += delivered_energy * offset[k];
  }

#ifdef SWIFT_DEBUG_CHECKS
  pk_r->shares_received++;
#endif
}

/**
 * @brief Force loop: book the sender's share for one receiver.
 *
 * The share is absorbed over the increase of the distance to the sender's
 * origin, so path lengths telescope along a chain of hops.
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

  const float w = radiation_packet_transmission_weight(dx, h_s, pk_s->direction,
                                                       pk_s->cone_power);
  if (w == 0.f) return;

#ifdef SWIFT_DEBUG_CHECKS
  if (!(pk_s->weight_sum > 0.f))
    error("Particle %lld books a share with weight sum %e.", s->id,
          pk_s->weight_sum);
#endif
  const float fraction = w / pk_s->weight_sum;

  const float d_s = pk_s->origin_distance;
  const float offset[3] = {-dx[0] - d_s * pk_s->direction[0],
                           -dx[1] - d_s * pk_s->direction[1],
                           -dx[2] - d_s * pk_s->direction[2]};
  const float d_r = sqrtf(offset[0] * offset[0] + offset[1] * offset[1] +
                          offset[2] * offset[2]);
  float path = d_r - d_s;
  if (path < 0.f) path = 0.f;

  float energy[ISRF_OPERATOR_COUNT], photons[ISRF_OPERATOR_COUNT];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    energy[b] = pk_s->band[b].energy * fraction;
    photons[b] = pk_s->band[b].photons * fraction;
  }

  radiation_packet_deliver(r, energy, photons, offset, a * path,
                           s->time_bin > r->time_bin, pk_s->step_end_ti);
}

/**
 * @brief Zero a star's packet emission state and escape ledger, once at the
 * star's creation.
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_first_init_star(struct spart *sp) {
  memset(&sp->feedback_data.radiation.packet, 0,
         sizeof(sp->feedback_data.radiation.packet));
}

/**
 * @brief Zero a star's emission-cone weights before its density loop.
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_init_star_cones(struct spart *sp) {
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++)
    sp->feedback_data.radiation.packet.cone_weight[c] = 0.f;
}

/**
 * @brief Star density loop: add a gas neighbour to its emission cone.
 *
 * @param si The star.
 * @param dx Comoving vector from the gas to the star (si - pj).
 * @param weight Kernel weight times gas mass, as in the feedback loop.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_star_density(
    struct spart *si, const float dx[3], const float weight,
    const integertime_t ti_current) {
  const integertime_t ti_begin =
      get_integer_time_begin(ti_current, si->time_bin);
  const int c = radiation_packet_emission_cone(si->id, ti_begin, dx);
  si->feedback_data.radiation.packet.cone_weight[c] += weight;
}

/**
 * @brief Feedback loop: give a gas neighbour its share of the star's cone.
 *
 * The star's energy of this step is split equally between the cones, and
 * within a cone by #radiation_packet_star_density's weight. The share is
 * absorbed over the distance to the star.
 *
 * @param si The star.
 * @param pj The gas neighbour.
 * @param dx Comoving vector from the gas to the star (si - pj).
 * @param r Comoving distance.
 * @param weight Kernel weight times gas mass, as in the density loop.
 * @param Delta_t The star's feedback step.
 * @param a Current scale factor.
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_star_emit(
    const struct spart *si, struct part *pj, const float dx[3], const float r,
    const float weight, const float Delta_t, const float a,
    const integertime_t ti_current) {

  if (weight == 0.f) return;
  const integertime_t ti_step = get_integer_timestep(si->time_bin);
  const integertime_t ti_begin =
      get_integer_time_begin(ti_current, si->time_bin);
  const int c = radiation_packet_emission_cone(si->id, ti_begin, dx);
  const float cone_weight = si->feedback_data.radiation.packet.cone_weight[c];
  /* Possible only at the kernel edge, where the weight underflows in one loop
     and not the other: that cone's light was booked as escaped. */
  if (!(cone_weight > 0.f)) return;
  const double fraction =
      (double)weight / ((double)cone_weight * RADIATION_PACKET_EMISSION_CONES);

  float energy[ISRF_OPERATOR_COUNT], photons[ISRF_OPERATOR_COUNT];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const enum radiation_isrf_moment m = radiation_isrf_operator_owner[b];
    const double E =
        si->feedback_data.radiation.L_band[m] * (double)Delta_t * fraction;
    energy[b] = (float)E;
    photons[b] =
        (float)(E * si->feedback_data.radiation.packet.photon_fraction[b]);
  }

  radiation_packet_deliver(pj, energy, photons, dx, a * r,
                           si->time_bin > pj->time_bin, ti_current + ti_step);
}

void radiation_packet_star_book_escape(struct spart *sp, const float Delta_t);

#else /* GEAR_SUBGRID_RADIATION_ISRF */

__attribute__((always_inline)) INLINE static void
radiation_packet_first_init_star(struct spart *sp) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_init_star_cones(struct spart *sp) {}

__attribute__((always_inline)) INLINE static void radiation_packet_star_density(
    struct spart *si, const float dx[3], const float weight,
    const integertime_t ti_current) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_star_book_escape(struct spart *sp, const float Delta_t) {}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_PACKET_GEAR_H */
