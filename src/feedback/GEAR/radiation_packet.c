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
#include <config.h>

/* This object's header. */
#include "radiation_packet.h"

/* Local includes. */
#include "dimension.h"
#include "kernel_hydro.h"

#include <math.h>

/**
 * @brief Transmission weight `W(r, h_s) max(0, mu)^p` of a sender towards
 * one neighbour, with `mu` the cosine between the hop and the packet
 * direction.
 *
 * Exactly 0 for `mu <= RADIATION_PACKET_MU_MIN`, for `r = 0` and beyond the
 * sender's kernel support. Kept out of line so that the gradient loop, which
 * sums the weights, and the force loop, which divides by that sum, evaluate
 * the same compiled arithmetic. The distance is computed here too: a caller's
 * `r2` may be contracted differently at each call site.
 *
 * @param dx Comoving vector from the sender to the neighbour.
 * @param h_s Comoving smoothing length of the sender.
 * @param n Unit direction of the packet.
 * @param power Power `p` of the cone.
 * @return The weight.
 */
__attribute__((noinline)) float radiation_packet_transmission_weight(
    const float dx[3], const float h_s, const float n[3], const float power) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float H2 = h_s * h_s * kernel_gamma2;
  if (!(r2 > 0.f) || r2 >= H2) return 0.f;

  const float r = sqrtf(r2);
  const float mu = (dx[0] * n[0] + dx[1] * n[1] + dx[2] * n[2]) / r;
  if (mu <= RADIATION_PACKET_MU_MIN) return 0.f;

  const float h_inv = 1.f / h_s;
  float w;
  kernel_eval(r * h_inv, &w);
  return w * pow_dimension(h_inv) * powf(mu, power);
}

/**
 * @brief Start of an active step: the inbox and the due part of the
 * reservoir become the held packet.
 *
 * The origin is the energy-weighted mean of the arriving origins and the
 * direction points from it to the particle. The reservoir releases the part
 * due in the step that starts now, `dt / t_rem`, all of it once `t_rem <= dt`.
 *
 * @param p The particle, active.
 * @param ti_current Start of the particle's new step.
 * @param ti_step Length of the particle's step.
 * @param time_base Physical length of one integer time step.
 */
void radiation_packet_merge(struct part *p, const integertime_t ti_current,
                            const integertime_t ti_step,
                            const double time_base) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
#ifdef SWIFT_DEBUG_CHECKS
  if (!radiation_packet_is_empty(pk))
    error("Particle %lld still holds a packet at its merge.", p->id);
#endif
  pk->step_end_ti = ti_current + ti_step;

  float release = 1.f;
  if (pk->reservoir_end_ti > ti_current) {
    const double t_rem =
        (double)(pk->reservoir_end_ti - ti_current) * time_base;
    const double dt = (double)ti_step * time_base;
    if (t_rem > dt) release = (float)(dt / t_rem);
  }

  float total_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    float res_energy = band->reservoir_energy;
    float res_photons = band->reservoir_photons;
    if (release < 1.f &&
        res_photons * release >= RADIATION_PACKET_MIN_PHOTONS) {
      res_energy *= release;
      res_photons *= release;
      band->reservoir_energy -= res_energy;
      band->reservoir_photons -= res_photons;
    } else {
      band->reservoir_energy = 0.f;
      band->reservoir_photons = 0.f;
    }
    band->energy = band->inbox_energy + res_energy;
    band->photons = band->inbox_photons + res_photons;
    band->inbox_energy = 0.f;
    band->inbox_photons = 0.f;
    total_energy += band->energy;
  }

  float offset[3];
  for (int k = 0; k < 3; k++) {
    float res_offset = pk->reservoir_origin_offset[k];
    if (release < 1.f) {
      res_offset *= release;
      pk->reservoir_origin_offset[k] -= res_offset;
    } else {
      pk->reservoir_origin_offset[k] = 0.f;
    }
    offset[k] = pk->inbox_origin_offset[k] + res_offset;
    pk->inbox_origin_offset[k] = 0.f;
  }

  if (!(total_energy > 0.f)) return;

  for (int k = 0; k < 3; k++) offset[k] /= total_energy;
  const float d = sqrtf(offset[0] * offset[0] + offset[1] * offset[1] +
                        offset[2] * offset[2]);
  pk->origin_distance = d;
  if (d > 0.f) {
    for (int k = 0; k < 3; k++) pk->direction[k] = -offset[k] / d;
  } else {
    /* Source on the particle: any direction is as good. */
    pk->direction[0] = 0.f;
    pk->direction[1] = 0.f;
    pk->direction[2] = 1.f;
  }
}

/**
 * @brief End of the gradient loop: normalise the forward advance, or let a
 * packet with no downstream neighbour escape into the ledger.
 *
 * @param p The particle.
 */
void radiation_packet_end_gradient(struct part *p) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  if (pk->weight_sum > 0.f) {
    pk->forward_advance /= pk->weight_sum;
    return;
  }

  pk->forward_advance = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    band->escaped_energy += (double)band->energy;
    band->escaped_photons += (double)band->photons;
    band->energy = 0.f;
    band->photons = 0.f;
  }
}

/**
 * @brief End of the force loop: the held packet has left in shares.
 *
 * @param p The particle.
 */
void radiation_packet_end_force(struct part *p) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    pk->band[b].energy = 0.f;
    pk->band[b].photons = 0.f;
  }
}
