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
#include "parser.h"
#include "random.h"
#include "star_gas_weights.h"

#if defined(HAVE_CHEALPIX)
#include <chealpix.h>
#endif

#include <float.h>
#include <math.h>
#include <stdint.h>

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
 * @brief Octant of the fixed frame that holds a direction: the reception slot
 * of the light that travels along it.
 *
 * Bit `k` of the result is set when the component `k` of the direction in the
 * frame is zero or positive, so a direction on a border goes to the positive
 * side whatever the sign of its zero. Directions more than 90 degrees apart
 * are never in one octant; exactly 90 degrees can be. The octant of the
 * opposite direction is `7 - octant`, except on a border. Kept out of line, so
 * that every caller bins with the same compiled arithmetic. The components are
 * formed in double.
 *
 * @param frame Rotation matrix (row-major) from the box to the frame.
 * @param direction Direction in the box frame. Need not be a unit vector.
 * @return The octant, in [0, RADIATION_PACKET_SLOTS).
 */
__attribute__((noinline)) int radiation_packet_octant(
    const float frame[9], const float direction[3]) {

  int octant = 0;
  for (int k = 0; k < 3; k++) {
    const double q = (double)frame[3 * k + 0] * (double)direction[0] +
                     (double)frame[3 * k + 1] * (double)direction[1] +
                     (double)frame[3 * k + 2] * (double)direction[2];
    if (q >= 0.0) octant |= 1 << k;
  }
  return octant;
}

/**
 * @brief Build the rotation matrix of the octant frame.
 *
 * The matrix is `Rz(a_z) Ry(a_y) Rx(a_x)` of the three angles. A positive
 * seed offset composes it with a rotation drawn uniformly from that seed
 * (Shoemake 1992), so that repeats of a run use independent frames. Offset 0
 * is the frame of the angles alone.
 *
 * @param angles_deg Rotation angles about x, y and z, in degrees.
 * @param seed_offset Seed of the extra rotation. 0 for none.
 * @param frame (return) The matrix, row-major.
 */
void radiation_packet_build_octant_frame(const float angles_deg[3],
                                         const int seed_offset,
                                         float frame[9]) {

  const double to_rad = M_PI / 180.0;
  const double cx = cos(angles_deg[0] * to_rad),
               sx = sin(angles_deg[0] * to_rad);
  const double cy = cos(angles_deg[1] * to_rad),
               sy = sin(angles_deg[1] * to_rad);
  const double cz = cos(angles_deg[2] * to_rad),
               sz = sin(angles_deg[2] * to_rad);

  const double Rx[3][3] = {{1., 0., 0.}, {0., cx, -sx}, {0., sx, cx}};
  const double Ry[3][3] = {{cy, 0., sy}, {0., 1., 0.}, {-sy, 0., cy}};
  const double Rz[3][3] = {{cz, -sz, 0.}, {sz, cz, 0.}, {0., 0., 1.}};

  double A[3][3], R[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      A[i][j] = 0.;
      for (int k = 0; k < 3; k++) A[i][j] += Ry[i][k] * Rx[k][j];
    }
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      R[i][j] = 0.;
      for (int k = 0; k < 3; k++) R[i][j] += Rz[i][k] * A[k][j];
    }

  if (seed_offset > 0) {
    const double u1 =
        random_unit_interval(seed_offset, 0, random_number_ISRF_octant_frame);
    const double u2 =
        random_unit_interval(seed_offset, 1, random_number_ISRF_octant_frame);
    const double u3 =
        random_unit_interval(seed_offset, 2, random_number_ISRF_octant_frame);
    const double a = sqrt(1.0 - u1), b = sqrt(u1);
    const double qx = a * sin(2.0 * M_PI * u2), qy = a * cos(2.0 * M_PI * u2);
    const double qz = b * sin(2.0 * M_PI * u3), qw = b * cos(2.0 * M_PI * u3);
    const double Q[3][3] = {
        {1.0 - 2.0 * (qy * qy + qz * qz), 2.0 * (qx * qy - qz * qw),
         2.0 * (qx * qz + qy * qw)},
        {2.0 * (qx * qy + qz * qw), 1.0 - 2.0 * (qx * qx + qz * qz),
         2.0 * (qy * qz - qx * qw)},
        {2.0 * (qx * qz - qy * qw), 2.0 * (qy * qz + qx * qw),
         1.0 - 2.0 * (qx * qx + qy * qy)}};
    double C[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        C[i][j] = 0.;
        for (int k = 0; k < 3; k++) C[i][j] += Q[i][k] * R[k][j];
      }
    memcpy(R, C, sizeof(R));
  }

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) frame[3 * i + j] = (float)R[i][j];
}

/**
 * @brief Read the test beams of the packet transport.
 *
 * A beam is 11 numbers: the lower and the upper corner of a box (3 + 3), its
 * direction (3) and the energy per particle in the PE and in the LW band
 * (2), all in internal units. Absent or with zero energy, a beam is off.
 * The direction is normalised. An active beam needs a box with upper corner
 * above the lower one on every axis. Finiteness is tested on the bit pattern,
 * which -ffinite-math-only cannot fold away.
 *
 * @param params The parameter file.
 * @param beams (return) The beams.
 * @return The number of active beams.
 */
int radiation_packet_read_test_beams(
    struct swift_params *params, float beams[RADIATION_PACKET_TEST_BEAMS][11]) {

  int n_active = 0;
  memset(beams, 0, sizeof(float) * 11 * RADIATION_PACKET_TEST_BEAMS);
  parser_get_opt_param_float_array(
      params, "GEARFeedback:ISRF_packet_test_beam_1", 11, beams[0]);
  parser_get_opt_param_float_array(
      params, "GEARFeedback:ISRF_packet_test_beam_2", 11, beams[1]);

  for (int i = 0; i < RADIATION_PACKET_TEST_BEAMS; i++) {
    float *bm = beams[i];
    for (int k = 0; k < 11; k++) {
      uint32_t bits;
      memcpy(&bits, &bm[k], sizeof(bits));
      volatile uint32_t opaque = bits;
      if (((opaque >> 23) & 0xFFu) == 0xFFu)
        error("GEARFeedback:ISRF_packet_test_beam_%d: non-finite entry %d.",
              i + 1, k);
    }
    if (bm[9] < 0.f || bm[10] < 0.f)
      error("GEARFeedback:ISRF_packet_test_beam_%d: negative energy.", i + 1);
    if (bm[9] == 0.f && bm[10] == 0.f) continue;
    const float norm = sqrtf(bm[6] * bm[6] + bm[7] * bm[7] + bm[8] * bm[8]);
    if (!(norm > 0.f))
      error("GEARFeedback:ISRF_packet_test_beam_%d: zero direction.", i + 1);
    for (int k = 6; k < 9; k++) bm[k] /= norm;
    for (int k = 0; k < 3; k++)
      if (!(bm[3 + k] > bm[k]))
        error(
            "GEARFeedback:ISRF_packet_test_beam_%d: the box has upper corner "
            "%g not above lower corner %g on axis %d.",
            i + 1, bm[3 + k], bm[k], k);
    n_active++;
    message(
        "ISRF packet test beam %d is ACTIVE: energy %g (PE) and %g (LW) per "
        "particle in the box, injected once at each particle's first merge, "
        "in no ledger.",
        i + 1, bm[9], bm[10]);
  }
  return n_active;
}

/**
 * @brief Give the particles inside a test beam box their beam pulse, once.
 *
 * The pulse goes into the inbox of the slot of the beam direction, as a share
 * whose origin lies one smoothing length behind the particle along the
 * direction. Its photon number is its energy, which is the mean photon energy
 * of the band at the lower edge. The energy is not in any ledger of the
 * particle: a test adds it up from the parameters. Called before the merge of
 * the first step of the particle.
 *
 * @param p The particle.
 * @param beams The beams, see #radiation_packet_read_test_beams.
 */
void radiation_packet_inject_test_beams(
    struct part *p, const float beams[RADIATION_PACKET_TEST_BEAMS][11]) {

  for (int i = 0; i < RADIATION_PACKET_TEST_BEAMS; i++) {
    const float *bm = beams[i];
    if (bm[9] == 0.f && bm[10] == 0.f) continue;
    int inside = 1;
    for (int k = 0; k < 3; k++)
      if (!(p->x[k] >= (double)bm[k] && p->x[k] < (double)bm[3 + k]))
        inside = 0;
    if (!inside) continue;

    const double energy[ISRF_OPERATOR_COUNT] = {bm[9], bm[10]};
    const float offset[3] = {-p->h * bm[6], -p->h * bm[7], -p->h * bm[8]};
    radiation_packet_deliver(p, p->feedback_data.isrf_packet.octant_frame,
                             energy, energy, offset, 0.f, 0, 0);
  }
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
 * @brief Emission cone and weight of a gas neighbour of a star.
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
 * @param weight (return) Tent weight, 0 when the gas is not a receiver:
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
    *weight = star_gas_tent_weight(*q, (double)hydro_get_mass(pj),
                                   (double)hydro_get_comoving_density(pj));
  return radiation_packet_emission_cone(si->id, ti_begin, d);
}

/**
 * @brief Start of an active step, for one slot: the inbox and the due part of
 * the reservoir become the held packet.
 *
 * The origin is the energy-weighted mean of the arriving origins and the
 * direction points from it to the particle. The reservoir releases the part
 * due in the step that starts now, `dt / t_rem`, all of it once `t_rem <= dt`.
 *
 * @param sl The slot.
 * @param ti_current Start of the particle's new step.
 * @param ti_step Length of the particle's step.
 * @param time_base Physical length of one integer time step.
 * @param slot_index Index of the slot, for the debugging check.
 * @param frame Octant frame, for the debugging check.
 */
static void radiation_packet_merge_slot(
    struct feedback_isrf_packet_slot_data *sl, const integertime_t ti_current,
    const integertime_t ti_step, const double time_base, const int slot_index,
    const float frame[9]) {

  float release = 1.f;
  if (sl->reservoir_end_ti > ti_current) {
    const double t_rem =
        (double)(sl->reservoir_end_ti - ti_current) * time_base;
    const double dt = (double)ti_step * time_base;
    if (t_rem > dt) release = (float)(dt / t_rem);
  }

  /* One decision for every band and the origin sum: a part below the photon
     cutoff releases the whole reservoir, and its origin sum with it. */
  int partial = release < 1.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT && partial; b++) {
    const float res_photons = sl->band[b].reservoir_photons;
    if (res_photons > 0.f &&
        res_photons * release < RADIATION_PACKET_MIN_PHOTONS)
      partial = 0;
  }

  float total_energy = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
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
    double res_offset = sl->reservoir_origin_offset[k];
    if (partial) {
      res_offset *= (double)release;
      sl->reservoir_origin_offset[k] -= res_offset;
    } else {
      sl->reservoir_origin_offset[k] = 0.0;
    }
    offset[k] = sl->inbox_origin_offset[k] + res_offset;
    sl->inbox_origin_offset[k] = 0.0;
  }

  if (!(total_energy > 0.f)) return;

  for (int k = 0; k < 3; k++) offset[k] /= (double)total_energy;
  const double d = sqrt(offset[0] * offset[0] + offset[1] * offset[1] +
                        offset[2] * offset[2]);
  sl->origin_distance = (float)d;
  if (d > 0.0) {
    for (int k = 0; k < 3; k++) sl->direction[k] = (float)(-offset[k] / d);
#ifdef SWIFT_DEBUG_CHECKS
    /* Every share of the slot came from its octant, an octant is convex, so
       the mean origin lies in it too. A direction too close to a border may
       flip by rounding to float. */
    double q_min = 1.;
    for (int k = 0; k < 3; k++) {
      double q = 0.;
      for (int j = 0; j < 3; j++)
        q += (double)frame[3 * k + j] * (double)sl->direction[j];
      if (fabs(q) < q_min) q_min = fabs(q);
    }
    if (q_min > 1e-5 &&
        radiation_packet_octant(frame, sl->direction) != slot_index)
      error("The direction of slot %d is in octant %d.", slot_index,
            radiation_packet_octant(frame, sl->direction));
#endif
  } else {
    /* Mean origin on the particle: any direction is as good. It can lie
       outside the slot's octant, which is harmless: it is not a mean of
       directions, and downstream shares are binned by their own geometry. */
    sl->direction[0] = 0.f;
    sl->direction[1] = 0.f;
    sl->direction[2] = 1.f;
  }
}

/**
 * @brief Start of an active step: in every slot, the inbox and the due part of
 * the reservoir become the held packet.
 *
 * Each share was binned into its slot when it arrived, so the slots stay
 * apart here and merge nothing across octants.
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

  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    radiation_packet_merge_slot(&pk->slot[s], ti_current, ti_step, time_base, s,
                                pk->octant_frame);
}

/**
 * @brief End of the gradient loop: normalise the forward advance of each slot
 * and read the own field off the held packets, or let a slot with no
 * downstream neighbour escape into the ledger.
 *
 * The held packet of a slot is the light that crosses the particle in one
 * step, over the slot's forward advance, so its specific energy density is
 * `(E / m) l_par / (c dt)`. The own field is the sum of that value over the
 * slots, then its running mean over `tau` (#radiation_packet_read_mean), so
 * the kernel mean gathers the mean.
 *
 * @param p The particle.
 * @param a Current scale factor.
 * @param c Speed of light.
 * @param dt Physical length of the particle's step.
 * @param tau Time constant of the running mean. 0 or less is off.
 */
void radiation_packet_end_gradient(struct part *p, const float a, const float c,
                                   const float dt, const float tau) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  const float m = hydro_get_mass(p);

  /* A zero frame copy would send every share to one slot. Once per step. */
  float frame_norm = 0.f;
  for (int k = 0; k < 9; k++)
    frame_norm += pk->octant_frame[k] * pk->octant_frame[k];
  if (!(fabsf(frame_norm - 3.f) < 1e-3f))
    error("Particle %lld has an unset octant frame (norm %e) at its step.",
          p->id, frame_norm);

  double field_energy[ISRF_OPERATOR_COUNT], field_photons[ISRF_OPERATOR_COUNT];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    field_energy[b] = 0.;
    field_photons[b] = 0.;
  }

  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
    struct feedback_isrf_packet_slot_data *sl = &pk->slot[s];
    if (sl->weight_sum > 0.f) {
      sl->forward_advance /= sl->weight_sum;
      float factor = 0.f;
      if (dt > 0.f && m > 0.f) factor = a * sl->forward_advance / (m * c * dt);
      for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
        field_energy[b] += (double)(sl->band[b].energy * factor);
        field_photons[b] += (double)(sl->band[b].photons * factor);
      }
      continue;
    }

    sl->forward_advance = 0.f;
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
      pk->band[b].escaped_energy += (double)band->energy;
      pk->band[b].escaped_photons += (double)band->photons;
      band->energy = 0.f;
      band->photons = 0.f;
    }
  }

  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    pk->band[b].field_energy = radiation_packet_read_mean(
        pk->band[b].field_energy, (float)field_energy[b], dt, tau);
    pk->band[b].field_photons = radiation_packet_read_mean(
        pk->band[b].field_photons, (float)field_photons[b], dt, tau);
  }
}

/**
 * @brief End of the force loop: the held packets have left in shares, and the
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
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      pk->slot[s].band[b].energy = 0.f;
      pk->slot[s].band[b].photons = 0.f;
    }
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_band_data *band = &pk->band[b];
    band->ngb_mean_field_energy =
        (band->ngb_mean_field_energy + mw * band->field_energy) * norm;
    band->ngb_mean_field_photons =
        (band->ngb_mean_field_photons + mw * band->field_photons) * norm;
  }
}

/**
 * @brief Write the kernel mean of the own field, summed over the slots, into
 * the band moments that Grackle reads: PE to PE, LW to LW and the LW photon
 * number (`N E_lo`) to LW_PHOTON.
 *
 * Overwrites the moments, including any value read from the initial
 * conditions. Call after #radiation_packet_end_force, which completes the
 * kernel mean.
 *
 * @param p The particle.
 */
void radiation_packet_couple_field(struct part *p) {

  const struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    p->feedback_data.isrf_moment[radiation_isrf_operator_owner[b]].u =
        (double)pk->band[b].ngb_mean_field_energy;
  p->feedback_data.isrf_moment[ISRF_MOMENT_LW_PHOTON].u =
      (double)pk->band[ISRF_OPERATOR_LW].ngb_mean_field_photons;
}

/**
 * @brief Book a star's emission of this feedback step, and the light of each
 * cone that its receivers do not take as escaped, once per step on the star's
 * own rank.
 *
 * An empty cone loses its whole twelfth, so a star with no receiver loses all
 * its light.
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

  const double floor = radiation_packet_cone_floor(sp);
  double fraction = 0.0;
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++)
    fraction +=
        star_gas_cone_escape(sp->feedback_data.radiation.packet.cone_weight[c],
                             floor, RADIATION_PACKET_EMISSION_CONES);
  if (fraction == 0.0) return;

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
 * same receivers as the density loop that summed the cone weights.
 *
 * Both loops take each weight from radiation_packet_star_gas_weight() on the
 * same positions, so the cone weight sums can differ only by summation order.
 * One line per star step reports the receivers, the empty cones and those of
 * them that hold gas beyond #RADIATION_PACKET_EMISSION_EDGE only, the largest
 * relative difference of a cone weight sum, the worst membership difference
 * (density minus feedback loop), the step's PE energy and the delivered plus
 * escaped fraction of the star's light, minus one.
 *
 * @param sp The star, after its feedback loop.
 * @param ti_current Current integer time, for the message.
 */
void radiation_packet_check_star_cones(const struct spart *sp,
                                       const integertime_t ti_current) {
  if (!sp->feedback_data.radiation.packet.cones_booked) return;
  const double floor = radiation_packet_cone_floor(sp);
  int n_receivers = 0, n_empty = 0, n_edge_only = 0, n_diff = 0;
  long long id_diff = 0;
  double dS = 0.0, accounted = 0.0;
  for (int c = 0; c < RADIATION_PACKET_EMISSION_CONES; c++) {
    const int n = sp->feedback_data.radiation.packet.cone_count[c];
    const int n_d = sp->feedback_data.radiation.packet.cone_count_delivered[c];
    const long long id = sp->feedback_data.radiation.packet.cone_id_sum[c];
    const long long id_d =
        sp->feedback_data.radiation.packet.cone_id_sum_delivered[c];
    const double S = sp->feedback_data.radiation.packet.cone_weight[c];
    const double S_d =
        sp->feedback_data.radiation.packet.cone_weight_delivered[c];
    n_receivers += n;
    if (!(S > 0.0)) {
      n_empty++;
      if (sp->feedback_data.radiation.packet.cone_edge_count[c] > 0)
        n_edge_only++;
    } else {
      dS = fmax(dS, fabs(S_d - S) / S);
    }
    /* The delivered shares are those of the feedback loop's weight sum. */
    accounted +=
        star_gas_cone_share(S_d, S, floor, RADIATION_PACKET_EMISSION_CONES) +
        star_gas_cone_escape(S, floor, RADIATION_PACKET_EMISSION_CONES);
    /* The worst membership difference: the largest count difference, else
       a swap at the same count. */
    if (abs(n - n_d) > abs(n_diff) ||
        (n_diff == 0 && id_diff == 0 && id != id_d)) {
      n_diff = n - n_d;
      id_diff = id - id_d;
    }
  }
  const double E_pe = sp->feedback_data.radiation.L_band[ISRF_MOMENT_PE] *
                      (double)sp->feedback_data.radiation.Delta_t;
  message(
      "Packet emission of star %lld at ti_current %lld: %d receivers, %d "
      "empty cones (%d with edge gas only), cone weight sum largest relative "
      "difference %e, h %e and %e, membership difference %d gas (ID sum %lld, "
      "density minus feedback loop), PE energy %e, delivered fraction minus "
      "one %e.",
      sp->id, ti_current, n_receivers, n_empty, n_edge_only, dS,
      sp->feedback_data.radiation.packet.h_density,
      sp->feedback_data.radiation.packet.h_feedback, n_diff, id_diff, E_pe,
      accounted - 1.0);
}
#endif
