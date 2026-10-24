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
#include "engine.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "random.h"
#include "star_gas_weights.h"

#if defined(HAVE_CHEALPIX)
#include <chealpix.h>
#endif

#include <float.h>
#include <math.h>

/* The periodic box, set once at start-up by radiation_packet_set_box(). */
static double radiation_packet_box_dim[3] = {0.0, 0.0, 0.0};
static int radiation_packet_box_periodic = 0;
static int radiation_packet_box_is_set = 0;

/**
 * @brief Give the star emission the box it needs to wrap star-gas vectors.
 *
 * @param dim The box size.
 * @param periodic Is the box periodic?
 */
void radiation_packet_set_box(const double dim[3], const int periodic) {
  for (int k = 0; k < 3; k++) radiation_packet_box_dim[k] = dim[k];
  radiation_packet_box_periodic = periodic;
  radiation_packet_box_is_set = 1;
}

/**
 * @brief Start-up of the packet emission, from engine_config().
 *
 * @param e The #engine.
 */
void radiation_packet_startup(const struct engine *e) {
  radiation_packet_set_box(e->s->dim, e->s->periodic);
}

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
 * @brief Emission cone of a gas neighbour of a star: the HEALPix nside 1
 * pixel of the direction from the star to the gas, in a random frame drawn
 * from the star's ID and the start of its step.
 *
 * Kept out of line, so that the star density loop that sums the cones and
 * the feedback loop that fills them, on any rank, use the same arithmetic.
 *
 * @param star_id ID of the star.
 * @param ti_begin Start of the star's step.
 * @param dx Comoving vector from the gas to the star (star - gas).
 * @return The cone, in [0, RADIATION_PACKET_EMISSION_CONES).
 */
__attribute__((noinline)) int radiation_packet_emission_cone(
    const long long star_id, const integertime_t ti_begin, const double dx[3]) {

#if defined(HAVE_CHEALPIX)
  if (dx[0] == 0.0 && dx[1] == 0.0 && dx[2] == 0.0) return 0;

  const integertime_t ti = 3 * ti_begin;
  const double u1 =
      random_unit_interval(star_id, ti, random_number_ISRF_emission_frame);
  const double u2 =
      random_unit_interval(star_id, ti + 1, random_number_ISRF_emission_frame);
  const double u3 =
      random_unit_interval(star_id, ti + 2, random_number_ISRF_emission_frame);

  /* Uniform random rotation from a unit quaternion (Shoemake 1992). */
  const double a = sqrt(1.0 - u1);
  const double b = sqrt(u1);
  const double qx = a * sin(2.0 * M_PI * u2);
  const double qy = a * cos(2.0 * M_PI * u2);
  const double qz = b * sin(2.0 * M_PI * u3);
  const double qw = b * cos(2.0 * M_PI * u3);

  const double d0 = -dx[0], d1 = -dx[1], d2 = -dx[2];
  double dir[3];
  dir[0] = (1.0 - 2.0 * (qy * qy + qz * qz)) * d0 +
           2.0 * (qx * qy - qz * qw) * d1 + 2.0 * (qx * qz + qy * qw) * d2;
  dir[1] = 2.0 * (qx * qy + qz * qw) * d0 +
           (1.0 - 2.0 * (qx * qx + qz * qz)) * d1 +
           2.0 * (qy * qz - qx * qw) * d2;
  dir[2] = 2.0 * (qx * qz - qy * qw) * d0 + 2.0 * (qy * qz + qx * qw) * d1 +
           (1.0 - 2.0 * (qx * qx + qy * qy)) * d2;

  int64_t ipix = 0;
  vec2pix_ring64(1, dir, &ipix);
  return (int)ipix;
#else
  error("The ISRF packet emission needs the HEALPix C library.");
  return 0;
#endif
}

/**
 * @brief Emission cone and volume weight of a gas neighbour of a star.
 *
 * The star-gas vector is taken in double from the absolute positions; the
 * runner's float @p dx only picks the periodic image, from which it differs
 * by about 0 or about one box length. Both star loops call this one
 * out-of-line function on the same positions, so they agree bit for bit on
 * the cone, on the weight and on which gas is a receiver, whatever path the
 * runner took to each pair.
 *
 * @param si The star.
 * @param pj The gas neighbour.
 * @param dx The runner's comoving vector from the gas to the star.
 * @param ti_begin Start of the star's step.
 * @param d (return) Comoving vector from the gas to the star, in double.
 * @param q (return) Distance over the star's kernel support.
 * @param weight (return) Volume weight, 0 when the gas is not a receiver:
 * beyond #RADIATION_PACKET_EMISSION_EDGE or with no density.
 * @return The emission cone.
 */
__attribute__((noinline)) int radiation_packet_star_gas_weight(
    const struct spart *si, const struct part *pj, const float dx[3],
    const integertime_t ti_begin, double d[3], double *q, double *weight) {

  if (!radiation_packet_box_is_set)
    error("The packet emission was not given the box size.");

  for (int k = 0; k < 3; k++) {
    d[k] = si->x[k] - pj->x[k];
    if (radiation_packet_box_periodic) {
      const double wrap = d[k] - (double)dx[k];
      if (wrap > 0.5 * radiation_packet_box_dim[k])
        d[k] -= radiation_packet_box_dim[k];
      else if (wrap < -0.5 * radiation_packet_box_dim[k])
        d[k] += radiation_packet_box_dim[k];
    }
  }
  const double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  const double u = r / (double)si->h;
  *q = u * (double)kernel_gamma_inv;
  *weight = 0.0;
  if (*q < RADIATION_PACKET_EMISSION_EDGE)
    *weight = star_gas_volume_weight(u, (double)hydro_get_mass(pj),
                                     (double)hydro_get_comoving_density(pj));
  return radiation_packet_emission_cone(si->id, ti_begin, d);
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

  /* One decision for every band and the origin sum: a part below the photon
     cutoff releases the whole reservoir, and its origin sum with it. */
  int partial = release < 1.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT && partial; b++) {
    const float res_photons = pk->band[b].reservoir_photons;
    if (res_photons > 0.f &&
        res_photons * release < RADIATION_PACKET_MIN_PHOTONS)
      partial = 0;
  }

  float total_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    float res_energy = band->reservoir_energy;
    float res_photons = band->reservoir_photons;
    if (partial) {
      res_energy *= release;
      res_photons *= release;
      band->reservoir_energy -= res_energy;
      band->reservoir_photons -= res_photons;
    } else {
      band->reservoir_energy = 0.f;
      band->reservoir_photons = 0.f;
    }
    /* Added, not assigned: a packet still held would otherwise be lost. */
    const float arriving = band->inbox_energy + res_energy;
    band->energy += arriving;
    band->photons += band->inbox_photons + res_photons;
    band->inbox_energy = 0.f;
    band->inbox_photons = 0.f;
    total_energy += arriving;
  }

  /* Double: the product of a share near the transmission cutoff with a small
     offset component is below FLT_MIN and would flush to zero. */
  double offset[3];
  for (int k = 0; k < 3; k++) {
    double res_offset = pk->reservoir_origin_offset[k];
    if (partial) {
      res_offset *= (double)release;
      pk->reservoir_origin_offset[k] -= res_offset;
    } else {
      pk->reservoir_origin_offset[k] = 0.0;
    }
    offset[k] = pk->inbox_origin_offset[k] + res_offset;
    pk->inbox_origin_offset[k] = 0.0;
  }

  if (!(total_energy > 0.f)) return;

  for (int k = 0; k < 3; k++) offset[k] /= (double)total_energy;
  const double d = sqrt(offset[0] * offset[0] + offset[1] * offset[1] +
                        offset[2] * offset[2]);
  pk->origin_distance = (float)d;
  if (d > 0.0) {
    for (int k = 0; k < 3; k++) pk->direction[k] = (float)(-offset[k] / d);
  } else {
    /* Source on the particle: any direction is as good. */
    pk->direction[0] = 0.f;
    pk->direction[1] = 0.f;
    pk->direction[2] = 1.f;
  }
}

/**
 * @brief End of the gradient loop: normalise the forward advance and read the
 * own field off the held packet, or let a packet with no downstream neighbour
 * escape into the ledger.
 *
 * The held packet is the light that crosses the particle in one step, over
 * the forward advance, so its specific energy density is
 * `(E / m) l_par / (c dt)`.
 *
 * @param p The particle.
 * @param a Current scale factor.
 * @param c Speed of light.
 * @param dt Physical length of the particle's step.
 */
void radiation_packet_end_gradient(struct part *p, const float a, const float c,
                                   const float dt) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  if (pk->weight_sum > 0.f) {
    pk->forward_advance /= pk->weight_sum;
    float factor = 0.f;
    const float m = hydro_get_mass(p);
    if (dt > 0.f && m > 0.f) factor = a * pk->forward_advance / (m * c * dt);
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      pk->band[b].field_energy = pk->band[b].energy * factor;
      pk->band[b].field_photons = pk->band[b].photons * factor;
    }
    return;
  }

  pk->forward_advance = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    pk->band[b].field_energy = 0.f;
    pk->band[b].field_photons = 0.f;
  }
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    band->escaped_energy += (double)band->energy;
    band->escaped_photons += (double)band->photons;
    band->energy = 0.f;
    band->photons = 0.f;
  }
}

/**
 * @brief End of the force loop: the held packet has left in shares, and the
 * kernel mean of the own fields is completed with the particle's own term.
 *
 * @param p The particle.
 */
void radiation_packet_end_force(struct part *p) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  float w_self;
  kernel_eval(0.f, &w_self);
  const float mw = hydro_get_mass(p) * w_self;
  const float norm = 1.f / (pk->ngb_field_weight_sum + mw);
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    band->energy = 0.f;
    band->photons = 0.f;
    band->ngb_mean_field_energy =
        (band->ngb_mean_field_energy + mw * band->field_energy) * norm;
    band->ngb_mean_field_photons =
        (band->ngb_mean_field_photons + mw * band->field_photons) * norm;
  }
}

/**
 * @brief Book a star's emission of this feedback step, and the twelfths of
 * the sky in its empty cones as escaped, once per step on the star's own
 * rank.
 *
 * A star with no receiver has 12 empty cones: all its light escapes.
 *
 * @param sp The star, after its density loop.
 * @param Delta_t The star's feedback step.
 */
void radiation_packet_star_book_emission(struct spart *sp,
                                         const float Delta_t) {
#ifdef SWIFT_DEBUG_CHECKS
  sp->feedback_data.radiation.packet.cones_booked = 1;
#endif

  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const enum radiation_isrf_moment m = radiation_isrf_operator_owner[b];
    const double E = sp->feedback_data.radiation.L_band[m] * (double)Delta_t;
    sp->feedback_data.radiation.packet.emitted_energy[b] += E;
    sp->feedback_data.radiation.packet.emitted_photons[b] +=
        E * sp->feedback_data.radiation.packet.photon_fraction[b];
  }

  const int n_empty =
      RADIATION_PACKET_EMISSION_CONES - radiation_packet_occupied_cones(sp);
  if (n_empty == 0) return;

  const double fraction = (double)n_empty / RADIATION_PACKET_EMISSION_CONES;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const enum radiation_isrf_moment m = radiation_isrf_operator_owner[b];
    const double E =
        sp->feedback_data.radiation.L_band[m] * (double)Delta_t * fraction;
    sp->feedback_data.radiation.packet.escaped_energy[b] += E;
    sp->feedback_data.radiation.packet.escaped_photons[b] +=
        E * sp->feedback_data.radiation.packet.photon_fraction[b];
  }
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Debugging checks: report whether the feedback loop delivered to the
 * same receivers as the density loop that summed the weights.
 *
 * Both loops take each weight from radiation_packet_star_gas_weight() on the
 * same positions, so the weight sums can differ only by summation order. One
 * line per star step reports the receivers, the empty cones and those of them
 * that hold gas beyond #RADIATION_PACKET_EMISSION_EDGE only, the relative
 * difference of the weight sums, the worst membership difference (density
 * minus feedback loop), the step's PE energy and the delivered fraction of the
 * star's light minus one.
 *
 * @param sp The star, after its feedback loop.
 * @param ti_current Current integer time, for the message.
 */
void radiation_packet_check_star_cones(const struct spart *sp,
                                       const integertime_t ti_current) {
  if (!sp->feedback_data.radiation.packet.cones_booked) return;
  int n_receivers = 0, n_empty = 0, n_edge_only = 0, n_diff = 0;
  long long id_diff = 0;
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++) {
    const int n = sp->feedback_data.radiation.packet.cone_count[c];
    const int n_d = sp->feedback_data.radiation.packet.cone_count_delivered[c];
    const long long id = sp->feedback_data.radiation.packet.cone_id_sum[c];
    const long long id_d =
        sp->feedback_data.radiation.packet.cone_id_sum_delivered[c];
    n_receivers += n;
    if (!sp->feedback_data.radiation.packet.cone_occupied[c]) {
      n_empty++;
      if (sp->feedback_data.radiation.packet.cone_edge_count[c] > 0)
        n_edge_only++;
    }
    /* The worst membership difference: the largest count difference, else
       a swap at the same count. */
    if (abs(n - n_d) > abs(n_diff) ||
        (n_diff == 0 && id_diff == 0 && id != id_d)) {
      n_diff = n - n_d;
      id_diff = id - id_d;
    }
  }
  const double S = sp->feedback_data.radiation.packet.weight_sum;
  const double S_d = sp->feedback_data.radiation.packet.weight_sum_delivered;
  const double dS = S > 0.0 ? (S_d - S) / S : 0.0;
  const double covered = (double)(RADIATION_PACKET_EMISSION_CONES - n_empty) /
                         RADIATION_PACKET_EMISSION_CONES;
  const double E_pe = sp->feedback_data.radiation.L_band[ISRF_MOMENT_PE] *
                      (double)sp->feedback_data.radiation.Delta_t;
  message(
      "Packet emission of star %lld at ti_current %lld: %d receivers, %d "
      "empty cones (%d with edge gas only), weight sum relative difference "
      "%e, h %e and %e, membership difference %d gas (ID sum %lld, density "
      "minus feedback loop), PE energy %e, delivered fraction minus one %e.",
      sp->id, ti_current, n_receivers, n_empty, n_edge_only, dS,
      sp->feedback_data.radiation.packet.h_density,
      sp->feedback_data.radiation.packet.h_feedback, n_diff, id_diff, E_pe,
      covered * dS);
}
#endif
