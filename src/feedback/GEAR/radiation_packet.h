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
 * loop books each share `w/S` into the receiver. The weight carries the
 * receiver's volume `mass/rho`, divided inside it. Both loops call the one
 * out-of-line #radiation_packet_transmission_weight, so the shares of a
 * sender sum to 1 up to rounding. A share is booked by the rank that owns the
 * receiver, once.
 */

#include "cosmology.h"
#include "error.h"
#include "hydro.h"
#include "inline.h"
#include "kernel_hydro.h"
#include "part.h"
#include "star_gas_weights.h"

#include <math.h>
#include <string.h>

#ifdef SWIFT_DEBUG_CHECKS
#include "atomic.h"
#endif

#ifdef GEAR_SUBGRID_RADIATION_ISRF

/*! Cosine below which the transmission weight is exactly zero. */
#define RADIATION_PACKET_MU_MIN 1e-3f

/*! Optical depth above which a share is fully absorbed. */
#define RADIATION_PACKET_MAX_OPTICAL_DEPTH 80.f

/*! Photon number (as `N E_lo`) below which a share is booked as absorbed:
    near the float minimum 1.2e-38 the number would flush to zero before the
    energy, and the mean photon energy would leave the band. */
#define RADIATION_PACKET_MIN_PHOTONS 1e-35f

/*! Distance over the star's kernel support below which a gas particle
    receives the star's light. */
#define RADIATION_PACKET_EMISSION_EDGE STAR_GAS_WEIGHT_EDGE

/*! Is the speed of the light of a share limited? Set by
    #radiation_packet_set_light_speed. */
extern int radiation_packet_light_speed_limited;

float radiation_packet_transmission_weight(const float dx[3], const float h_s,
                                           const float n[3], const float power,
                                           const float mass, const float rho);
int radiation_packet_emission_cone(const long long star_id,
                                   const integertime_t ti_begin,
                                   const double dx[3]);
int radiation_packet_star_gas_weight(const struct spart *si,
                                     const struct part *pj, const float dx[3],
                                     const integertime_t ti_begin, double d[3],
                                     double *q, double *weight);
void radiation_packet_set_box(const double dim[3], const int periodic);
void radiation_packet_set_redshift(const double speed_of_light,
                                   const double weight_pe,
                                   const double weight_lw,
                                   const double photon_weight_lw);
void radiation_packet_set_light_speed(const double limit_speed,
                                      const double time_base,
                                      const int cosmological);
double radiation_packet_arrival_delay(const double path_phys, const double H);
double radiation_packet_arrival_delay_receding(const double path_phys,
                                               const double H,
                                               const double closing_speed);
integertime_t radiation_packet_bin_end(
    const struct feedback_isrf_packet_data *pk, const int j);
integertime_t radiation_packet_far_end(
    const struct feedback_isrf_packet_data *pk, const int c);
void radiation_packet_step_diagnostics(const struct engine *e);
void radiation_packet_bin_counters(long long *widened, long long *beyond,
                                   long long *clamped);
void radiation_packet_deposit(struct feedback_isrf_packet_data *pk,
                              struct feedback_isrf_packet_slot_data *sl,
                              const float energy[ISRF_OPERATOR_COUNT],
                              const float photons[ISRF_OPERATOR_COUNT],
                              const double first, const double last);
void radiation_packet_redshift_share(const double H, const double path_phys,
                                     double energy[ISRF_OPERATOR_COUNT],
                                     double photons[ISRF_OPERATOR_COUNT],
                                     double lost_energy[ISRF_OPERATOR_COUNT],
                                     double lost_photons[ISRF_OPERATOR_COUNT]);
void radiation_packet_merge(struct part *p, const integertime_t ti_current,
                            const integertime_t ti_step, const double time_base,
                            const struct cosmology *cosmo);
void radiation_packet_end_gradient(struct part *p, const float a, const float c,
                                   const float dt, const float tau);
void radiation_packet_end_force(struct part *p);
int radiation_packet_octant(const float frame[9], const float direction[3]);
int radiation_packet_slot(const float frame[9], const float direction[3]);
double radiation_packet_slot_border_distance(const float frame[9],
                                             const float direction[3]);
void radiation_packet_couple_field(struct part *p);
void radiation_packet_inject_test_beams(
    struct part *p, const float beams[RADIATION_PACKET_TEST_BEAMS][11],
    const integertime_t ti_current);

/**
 * @brief Exponential running mean of the own field over the time constant
 * `tau`: `m + (x - m) dt / max(tau, dt)`.
 *
 * Only the field is averaged, never the packet, so the ledger is untouched.
 * Explicit branches, not a clamp: under -ffast-math a NaN would not be
 * absorbed, and `m + (x - m) * 1` is not bit-exact.
 *
 * @param m Mean at the particle's previous step.
 * @param x Value of this step.
 * @param dt Length of this step.
 * @param tau Time constant. 0 or less is off: the mean is `x`.
 * @return The new mean. A step of length 0 or less keeps `m`.
 */
__attribute__((always_inline)) INLINE static float radiation_packet_read_mean(
    const float m, const float x, const float dt, const float tau) {
  if (!(tau > 0.f) || dt >= tau) return x;
  if (!(dt > 0.f)) return m;
  return m + (x - m) * (dt / tau);
}

/**
 * @brief Does the slot hold no energy in any band?
 *
 * @param sl The slot.
 */
__attribute__((always_inline)) INLINE static int radiation_packet_slot_is_empty(
    const struct feedback_isrf_packet_slot_data *sl) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    if (sl->band[b].energy != 0.f) return 0;
  return 1;
}

/**
 * @brief Does the particle hold no energy in any band of any slot?
 *
 * @param pk The packet.
 */
__attribute__((always_inline)) INLINE static int radiation_packet_is_empty(
    const struct feedback_isrf_packet_data *pk) {
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    if (!radiation_packet_slot_is_empty(&pk->slot[s])) return 0;
  return 1;
}

/**
 * @brief Gradient loop: add one downstream neighbour to the weight sum and
 * forward advance of each non-empty slot of the sender.
 *
 * @param dx Comoving vector from the sender to the neighbour.
 * @param h_s Comoving smoothing length of the sender.
 * @param s The sender.
 * @param r The neighbour.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_gradient_accumulate(const float dx[3], const float h_s,
                                     struct part *restrict s,
                                     const struct part *restrict r) {

  struct feedback_isrf_packet_data *pk = &s->feedback_data.isrf_packet;
  for (int k = 0; k < RADIATION_PACKET_SLOTS; k++) {
    struct feedback_isrf_packet_slot_data *sl = &pk->slot[k];
    if (radiation_packet_slot_is_empty(sl)) continue;

    const float w = radiation_packet_transmission_weight(
        dx, h_s, sl->direction, pk->cone_power, hydro_get_mass(r),
        r->feedback_data.rho_prev);
    if (w == 0.f) continue;

    sl->weight_sum += w;
    sl->forward_advance +=
        w * (dx[0] * sl->direction[0] + dx[1] * sl->direction[1] +
             dx[2] * sl->direction[2]);
  }
}

/**
 * @brief Deliver one share to a receiver: dust absorption over the path, then
 * the inbox, or the reservoir when the source has the coarser step, of the
 * slot that holds the direction from the share's origin to the
 * receiver.
 *
 * Energy and photon number lose the same fraction to dust, so the mean photon
 * energy of a band is unchanged by it. In an expanding box the share first
 * loses energy and photons to the expansion over the light travel time of its
 * path, #radiation_packet_redshift_share. The share arrives and is split in
 * double: a float product near 1e-38 flushes to zero under -ffast-math, which
 * would leave an energy without its photon number.
 *
 * The receiver may drift while a share waits in its inbox. The offset to the
 * origin and the arrival time are fixed at the delivery, so for a receiver
 * that moves away from the origin at the speed `u` the light front moves at
 * `v + u`. The limit speed `v` is therefore taken as a fraction of `c` below
 * `1 - u_max / c`, with `u_max` the largest speed of the gas.
 *
 * @param r The receiver.
 * @param frame Rotation matrix of the octant frame, as
 * #feedback_props.ISRF_packet_octant_frame.
 * @param energy Energy of the share per band, before absorption.
 * @param photons Photon number of the share per band, as `N E_lo`.
 * @param offset Comoving offset from the receiver to the share's origin.
 * @param path_phys Physical path length over which the share is absorbed.
 * @param H Current Hubble rate. 0 or less for a static box.
 * @param to_reservoir Is the source's step coarser than the receiver's?
 * @param source_step_end_ti End of the source's step. Read only without a
 * speed limit: the reservoir is released over it.
 * @param arrival Integer time at which the first light of the share reaches
 * the receiver. Read only when the speed of the light is limited: the share
 * then waits in its arrival bin (#radiation_packet_deposit) until the
 * receiver's merge at or after the end of the bin.
 * @param arrival_last Integer time at which the last light of the share
 * reaches the receiver. The light arrives evenly in between. The same as
 * @p arrival for a pulse.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_deliver(
    struct part *restrict r, const float frame[9],
    const double energy[ISRF_OPERATOR_COUNT],
    const double photons[ISRF_OPERATOR_COUNT], const float offset[3],
    const float path_phys, const float H, const int to_reservoir,
    const integertime_t source_step_end_ti, const double arrival,
    const double arrival_last) {

#ifdef SWIFT_DEBUG_CHECKS
  if (r->time_bin == time_bin_inhibited)
    error("Packet share booked into the removed gas particle %lld.", r->id);
  float frame_norm = 0.f;
  for (int k = 0; k < 9; k++) frame_norm += frame[k] * frame[k];
  if (!(fabsf(frame_norm - 3.f) < 1e-3f))
    error("Packet share booked with an unset octant frame (norm %e).",
          frame_norm);
#endif

  const float from_origin[3] = {-offset[0], -offset[1], -offset[2]};
  struct feedback_isrf_packet_data *pk_r = &r->feedback_data.isrf_packet;
  struct feedback_isrf_packet_slot_data *sl =
      &pk_r->slot[radiation_packet_slot(frame, from_origin)];

  /* In a static box the arrays are used as given. */
  const double *en = energy;
  const double *ph = photons;
  double expanded_energy[ISRF_OPERATOR_COUNT];
  double expanded_photons[ISRF_OPERATOR_COUNT];
  double lost_energy[ISRF_OPERATOR_COUNT];
  double lost_photons[ISRF_OPERATOR_COUNT];
  if (H > 0.f) {
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      expanded_energy[b] = energy[b];
      expanded_photons[b] = photons[b];
    }
    radiation_packet_redshift_share(H, path_phys, expanded_energy,
                                    expanded_photons, lost_energy,
                                    lost_photons);
    en = expanded_energy;
    ph = expanded_photons;
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      pk_r->band[b].redshifted_energy += lost_energy[b];
      pk_r->band[b].redshifted_photons += lost_photons[b];
    }
  }

  float transmitted_energy_b[ISRF_OPERATOR_COUNT];
  float transmitted_photons_b[ISRF_OPERATOR_COUNT];
  float delivered_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *ledger = &pk_r->band[b];

    const float tau = r->feedback_data.isrf_operator[b].kappa * path_phys;
    double absorbed_fraction = 1.0;
    if (tau < RADIATION_PACKET_MAX_OPTICAL_DEPTH)
      absorbed_fraction = -expm1((double)-tau);

    double absorbed_energy = en[b] * absorbed_fraction;
    double absorbed_photons = ph[b] * absorbed_fraction;
    float transmitted_energy = (float)(en[b] - absorbed_energy);
    float transmitted_photons = (float)(ph[b] - absorbed_photons);
    if (transmitted_photons < RADIATION_PACKET_MIN_PHOTONS) {
      absorbed_energy = en[b];
      absorbed_photons = ph[b];
      transmitted_energy = 0.f;
      transmitted_photons = 0.f;
    }
    transmitted_energy_b[b] = transmitted_energy;
    transmitted_photons_b[b] = transmitted_photons;
    ledger->absorbed_energy += absorbed_energy;
    ledger->absorbed_photons += absorbed_photons;
    delivered_energy += transmitted_energy;
  }

  if (radiation_packet_light_speed_limited) {
    /* The light waits in its arrival bin until the receiver's merge at or
       after the end of the bin. */
    if (delivered_energy > 0.f)
      radiation_packet_deposit(pk_r, sl, transmitted_energy_b,
                               transmitted_photons_b, arrival, arrival_last);
    for (int k = 0; k < 3; k++)
      sl->inbox_origin_offset[k] += (double)delivered_energy * offset[k];
  } else if (to_reservoir) {
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      sl->band[b].reservoir_energy += transmitted_energy_b[b];
      sl->band[b].reservoir_photons += transmitted_photons_b[b];
    }
    for (int k = 0; k < 3; k++)
      sl->reservoir_origin_offset[k] += (double)delivered_energy * offset[k];
    if (source_step_end_ti > sl->reservoir_end_ti)
      sl->reservoir_end_ti = source_step_end_ti;
  } else {
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      sl->band[b].bin_energy[0] += transmitted_energy_b[b];
      sl->band[b].bin_photons[0] += transmitted_photons_b[b];
    }
    for (int k = 0; k < 3; k++)
      sl->inbox_origin_offset[k] += (double)delivered_energy * offset[k];
  }

#ifdef SWIFT_DEBUG_CHECKS
  pk_r->shares_received++;
#endif
}

/**
 * @brief Force loop: book the shares of each non-empty slot of the sender for
 * one receiver.
 *
 * Each share is absorbed over the increase of the distance to the slot's
 * origin, so path lengths telescope along a chain of hops. The octant frame
 * is the sender's copy.
 *
 * @param dx Comoving vector from the sender to the receiver.
 * @param h_s Comoving smoothing length of the sender.
 * @param s The sender.
 * @param r The receiver.
 * @param a Current scale factor.
 * @param H Current Hubble rate. 0 in a static box.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_book_share(
    const float dx[3], const float h_s, const struct part *restrict s,
    struct part *restrict r, const float a, const float H) {

  const struct feedback_isrf_packet_data *pk_s = &s->feedback_data.isrf_packet;
  for (int k = 0; k < RADIATION_PACKET_SLOTS; k++) {
    const struct feedback_isrf_packet_slot_data *sl = &pk_s->slot[k];
    if (radiation_packet_slot_is_empty(sl)) continue;

    const float w = radiation_packet_transmission_weight(
        dx, h_s, sl->direction, pk_s->cone_power, hydro_get_mass(r),
        r->feedback_data.rho_prev);
    if (w == 0.f) continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (!(sl->weight_sum > 0.f))
      error("Particle %lld books a share with weight sum %e in slot %d.", s->id,
            sl->weight_sum, k);
#endif
    const float fraction = w / sl->weight_sum;

    const float d_s = sl->origin_distance;
    const float offset[3] = {-dx[0] - d_s * sl->direction[0],
                             -dx[1] - d_s * sl->direction[1],
                             -dx[2] - d_s * sl->direction[2]};
    const float d_r = sqrtf(offset[0] * offset[0] + offset[1] * offset[1] +
                            offset[2] * offset[2]);
    float path = d_r - d_s;
    if (path < 0.f) path = 0.f;

    double energy[ISRF_OPERATOR_COUNT], photons[ISRF_OPERATOR_COUNT];
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      energy[b] = (double)sl->band[b].energy * (double)fraction;
      photons[b] = (double)sl->band[b].photons * (double)fraction;
    }

    /* With a limit, the light of a coarser sender is spread over one step of
       the sender from its first arrival time. The delay is the light time of
       the separation of the two particles: by the triangle inequality the
       light cannot reach the receiver sooner, whatever the mean origin of
       the packet. */
    double arrival = 0., arrival_last = 0.;
    const integertime_t source_end_ti = pk_s->step_end_ti;
    if (radiation_packet_light_speed_limited) {
      const double separation =
          sqrt((double)dx[0] * dx[0] + (double)dx[1] * dx[1] +
               (double)dx[2] * dx[2]);
      double closing = 0.;
      if (separation > 0.)
        closing = ((double)r->v[0] * dx[0] + (double)r->v[1] * dx[1] +
                   (double)r->v[2] * dx[2]) /
                  (separation * (double)a);
      arrival = sl->arrival + radiation_packet_arrival_delay_receding(
                                  a * separation, H, closing);
      arrival_last = s->time_bin > r->time_bin
                         ? arrival + (double)get_integer_timestep(s->time_bin)
                         : arrival;
    }
    radiation_packet_deliver(r, pk_s->octant_frame, energy, photons, offset,
                             a * path, H, s->time_bin > r->time_bin,
                             source_end_ti, arrival, arrival_last);
  }
}

/**
 * @brief Zero the kernel-mean accumulators of an active particle.
 *
 * @param p The particle.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_reset_field_mean(struct part *p) {
  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  pk->ngb_field_weight_sum = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    pk->band[b].ngb_mean_field_energy = 0.f;
    pk->band[b].ngb_mean_field_photons = 0.f;
  }
}

/**
 * @brief Force loop: add a neighbour's own field to the kernel mean of the
 * particle updated in this visit.
 *
 * @param r2 Comoving square distance.
 * @param h Comoving smoothing length of the updated particle.
 * @param pi The updated particle.
 * @param pj The neighbour.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_gather_field(
    const float r2, const float h, struct part *restrict pi,
    const struct part *restrict pj) {

  const float u = sqrtf(r2) / h;
  if (!(u < kernel_gamma)) return;
  float w;
  kernel_eval(u, &w);
  w *= hydro_get_mass(pj);
  struct feedback_isrf_packet_data *pk = &pi->feedback_data.isrf_packet;
  pk->ngb_field_weight_sum += w;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    pk->band[b].ngb_mean_field_energy +=
        w * pj->feedback_data.isrf_packet.band[b].field_energy;
    pk->band[b].ngb_mean_field_photons +=
        w * pj->feedback_data.isrf_packet.band[b].field_photons;
  }
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
 * @brief Zero a star's emission weights before its density loop.
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_init_star_cones(struct spart *sp) {
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++)
    sp->feedback_data.radiation.packet.cone_weight[c] = 0.0;
#ifdef SWIFT_DEBUG_CHECKS
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++) {
    sp->feedback_data.radiation.packet.cone_count[c] = 0;
    sp->feedback_data.radiation.packet.cone_edge_count[c] = 0;
    sp->feedback_data.radiation.packet.cone_id_sum[c] = 0;
    sp->feedback_data.radiation.packet.cone_count_delivered[c] = 0;
    sp->feedback_data.radiation.packet.cone_id_sum_delivered[c] = 0;
    sp->feedback_data.radiation.packet.cone_weight_delivered[c] = 0.0;
  }
  sp->feedback_data.radiation.packet.h_density = sp->h;
  sp->feedback_data.radiation.packet.h_feedback = 0.f;
  sp->feedback_data.radiation.packet.cones_booked = 0;
#endif
}

/**
 * @brief Weight sum of one emission cone at which its light is fully
 * delivered, for the star's current support.
 *
 * @param sp The #spart.
 */
__attribute__((always_inline)) INLINE static double radiation_packet_cone_floor(
    const struct spart *sp) {
  return star_gas_cone_floor((double)kernel_gamma * (double)sp->h,
                             RADIATION_PACKET_EMISSION_CONES);
}

/**
 * @brief Star density loop: add a gas neighbour's weight to the weight sum of
 * its emission cone, if it is a receiver.
 *
 * @param si The star.
 * @param pj The gas neighbour.
 * @param dx Comoving vector from the gas to the star (si - pj).
 * @param ti_current Current integer time.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_star_density(
    struct spart *si, const struct part *pj, const float dx[3],
    const integertime_t ti_current) {
  const integertime_t ti_begin =
      get_integer_time_begin(ti_current, si->time_bin);
  double d[3], q, weight;
  const int c =
      radiation_packet_star_gas_weight(si, pj, dx, ti_begin, d, &q, &weight);
  if (weight > 0.0) si->feedback_data.radiation.packet.cone_weight[c] += weight;
#ifdef SWIFT_DEBUG_CHECKS
  if (weight > 0.0) {
    si->feedback_data.radiation.packet.cone_count[c]++;
    si->feedback_data.radiation.packet.cone_id_sum[c] += pj->id;
  } else if (q < 1.0) {
    si->feedback_data.radiation.packet.cone_edge_count[c]++;
  }
#endif
}

/**
 * @brief Feedback loop: give a gas neighbour its share of the star's light.
 *
 * Each of the 12 cones holds `L dt / 12`. A cone whose weight sum is below the
 * cone floor lets the rest of its twelfth escape. The light that does not
 * escape is shared by all receivers of the star in proportion to their
 * weights. The share is absorbed over the distance to the star.
 *
 * @param si The star.
 * @param pj The gas neighbour.
 * @param dx Comoving vector from the gas to the star (si - pj).
 * @param Delta_t The star's feedback step. It is the step that has just
 * ended: the feedback of a step is applied at its end, and #compute_time gives
 * the begin of the step that finished at @p ti_current
 * (#feedback_prepare_radiation_feedback caches that step). The light was
 * emitted over `[ti_current - ti_step, ti_current]`, so it reaches a receiver
 * at the distance `r` evenly between `ti_current - ti_step + delay` and
 * `ti_current + delay`, with `delay` the light time of `r`.
 * @param a Current scale factor.
 * @param H Current Hubble rate. 0 in a static box.
 * @param ti_current Current integer time, the end of the star's step.
 * @param frame Rotation matrix of the octant frame, as
 * #feedback_props.ISRF_packet_octant_frame.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_star_emit(
    const struct spart *si, struct part *pj, const float dx[3],
    const float Delta_t, const float a, const float H,
    const integertime_t ti_current, const float frame[9]) {

  const integertime_t ti_step = get_integer_timestep(si->time_bin);
  const integertime_t ti_begin =
      get_integer_time_begin(ti_current, si->time_bin);
  double d[3], q, weight;
  const int c =
      radiation_packet_star_gas_weight(si, pj, dx, ti_begin, d, &q, &weight);
  if (!(weight > 0.0)) return;
  const double S = si->feedback_data.radiation.packet.cone_weight[c];
#ifdef SWIFT_DEBUG_CHECKS
  /* Several feedback tasks of one star may run at once. */
  struct spart *s_dbg = (struct spart *)si;
  atomic_inc(&s_dbg->feedback_data.radiation.packet.cone_count_delivered[c]);
  atomic_add(&s_dbg->feedback_data.radiation.packet.cone_id_sum_delivered[c],
             pj->id);
  atomic_add_d(&s_dbg->feedback_data.radiation.packet.cone_weight_delivered[c],
               weight);
  atomic_max_f(&s_dbg->feedback_data.radiation.packet.h_feedback, si->h);
#endif
  /* Possible only if this gas was not a receiver in the density loop. */
  if (!(S > 0.0)) return;
  const double fraction = star_gas_global_share(
      weight, si->feedback_data.radiation.packet.cone_weight,
      radiation_packet_cone_floor(si), RADIATION_PACKET_EMISSION_CONES);

  double energy[ISRF_OPERATOR_COUNT], photons[ISRF_OPERATOR_COUNT];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const enum radiation_isrf_moment m = radiation_isrf_operator_owner[b];
    energy[b] =
        si->feedback_data.radiation.L_band[m] * (double)Delta_t * fraction;
    photons[b] =
        energy[b] * si->feedback_data.radiation.packet.photon_fraction[b];
  }

  const float offset[3] = {(float)d[0], (float)d[1], (float)d[2]};
  const float r = (float)sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double delay = radiation_packet_arrival_delay(a * r, H);
  /* Without a limit the reservoir spreads the light over the star's next
     step. */
  const integertime_t window_end =
      ti_current + (radiation_packet_light_speed_limited ? 0 : ti_step) +
      (integertime_t)ceil(delay);
  radiation_packet_deliver(pj, frame, energy, photons, offset, a * r, H,
                           si->time_bin > pj->time_bin, window_end,
                           (double)(ti_current - ti_step) + delay,
                           (double)ti_current + delay);
}

void radiation_packet_star_book_emission(struct spart *sp, const float Delta_t);
#ifdef SWIFT_DEBUG_CHECKS
void radiation_packet_check_star_cones(const struct spart *sp,
                                       const integertime_t ti_current);
#endif

/**
 * @brief Book the packet of a gas particle that leaves the gas, converted
 * into or swallowed by another particle, into that particle's ledger.
 *
 * The light in flight (held, inbox and reservoir, in every slot) is lost with
 * the gas and is booked as removed. The gas's absorbed, escaped and redshifted
 * ledgers are carried over.
 *
 * @param p The gas particle.
 * @param l The ledger of the star or sink.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_book_removed_gas(const struct part *p,
                                  struct feedback_isrf_packet_removed_data *l) {
  const struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const struct feedback_isrf_packet_band_data *d = &pk->band[b];
    for (int k = 0; k < RADIATION_PACKET_SLOTS; k++) {
      const struct feedback_isrf_packet_slot_band_data *t =
          &pk->slot[k].band[b];
      double waiting_energy = 0., waiting_photons = 0.;
      for (int i = 0; i < RADIATION_PACKET_BINS; i++) {
        waiting_energy += (double)t->bin_energy[i];
        waiting_photons += (double)t->bin_photons[i];
      }
      if (radiation_packet_light_speed_limited) {
        for (int i = 0; i < RADIATION_PACKET_FAR_CELLS; i++) {
          waiting_energy += (double)t->far_energy[i];
          waiting_photons += (double)t->far_photons[i];
        }
      } else {
        waiting_energy += (double)t->reservoir_energy;
        waiting_photons += (double)t->reservoir_photons;
      }
      l->removed_energy[b] += (double)t->energy + waiting_energy;
      l->removed_photons[b] += (double)t->photons + waiting_photons;
    }
    l->inherited_energy[b] +=
        d->absorbed_energy + d->escaped_energy + d->redshifted_energy;
    l->inherited_photons[b] +=
        d->absorbed_photons + d->escaped_photons + d->redshifted_photons;
  }
}

/**
 * @brief Add the ledger of a swallowed sink to the sink that swallows it.
 *
 * @param l The ledger of the surviving sink.
 * @param l_swallowed The ledger of the swallowed sink.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_add_removed(
    struct feedback_isrf_packet_removed_data *l,
    const struct feedback_isrf_packet_removed_data *l_swallowed) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    l->removed_energy[b] += l_swallowed->removed_energy[b];
    l->removed_photons[b] += l_swallowed->removed_photons[b];
    l->inherited_energy[b] += l_swallowed->inherited_energy[b];
    l->inherited_photons[b] += l_swallowed->inherited_photons[b];
  }
}

/**
 * @brief Gas converted into a star: book its packet to the star, after the
 * star's #radiation_packet_first_init_star.
 *
 * @param p The gas particle.
 * @param sp The new star.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_book_removed_gas_to_star(const struct part *p,
                                          struct spart *sp) {
  radiation_packet_book_removed_gas(
      p, &sp->feedback_data.radiation.packet.removed);
}

/**
 * @brief Gas converted into or swallowed by a sink: book its packet to the
 * sink.
 *
 * @param p The gas particle.
 * @param s The sink.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_book_removed_gas_to_sink(const struct part *p,
                                          struct sink *s) {
  radiation_packet_book_removed_gas(p, &s->feedback_data.isrf_packet);
}

/**
 * @brief A sink swallows another: the survivor takes its ledger.
 *
 * @param si The surviving sink.
 * @param sj The swallowed sink.
 */
__attribute__((always_inline)) INLINE static void radiation_packet_sink_merger(
    struct sink *si, const struct sink *sj) {
  radiation_packet_add_removed(&si->feedback_data.isrf_packet,
                               &sj->feedback_data.isrf_packet);
}

/**
 * @brief Zero a sink's ledger, at its creation or first initialisation.
 *
 * @param s The sink.
 */
__attribute__((always_inline)) INLINE static void
radiation_packet_first_init_sink(struct sink *s) {
  memset(&s->feedback_data.isrf_packet, 0,
         sizeof(s->feedback_data.isrf_packet));
}

#else /* GEAR_SUBGRID_RADIATION_ISRF */

__attribute__((always_inline)) INLINE static void
radiation_packet_book_removed_gas_to_star(const struct part *p,
                                          struct spart *sp) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_book_removed_gas_to_sink(const struct part *p,
                                          struct sink *s) {}

__attribute__((always_inline)) INLINE static void radiation_packet_sink_merger(
    struct sink *si, const struct sink *sj) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_first_init_sink(struct sink *s) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_first_init_star(struct spart *sp) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_init_star_cones(struct spart *sp) {}

__attribute__((always_inline)) INLINE static void radiation_packet_star_density(
    struct spart *si, const struct part *pj, const float dx[3],
    const integertime_t ti_current) {}

__attribute__((always_inline)) INLINE static void
radiation_packet_star_book_emission(struct spart *sp, const float Delta_t) {}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */

#endif /* SWIFT_RADIATION_PACKET_GEAR_H */
