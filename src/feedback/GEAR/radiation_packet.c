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

/* The periodic box, set at start-up by radiation_packet_set_box(). */
static double radiation_packet_box_dim[3] = {0.0, 0.0, 0.0};
static int radiation_packet_box_periodic = 0;
static int radiation_packet_box_is_set = 0;

/**
 * @brief Give the star emission the box, to wrap star-gas vectors.
 *
 * @param dim The box size.
 * @param periodic Is the box periodic?
 */
void radiation_packet_set_box(const double dim[3], const int periodic) {
  for (int k = 0; k < 3; k++) radiation_packet_box_dim[k] = dim[k];
  radiation_packet_box_periodic = periodic;
  radiation_packet_box_is_set = 1;
}

/* The redshift constants, set at start-up by radiation_packet_set_redshift().
 */
static double radiation_packet_speed_of_light = 0.0;
static double radiation_packet_weight_energy[ISRF_OPERATOR_COUNT] = {1.0, 1.0};
static double radiation_packet_weight_photon_lw = 0.0;
static int radiation_packet_redshift_is_set = 0;

/**
 * @brief Give the redshift of a share the speed of light and the band-edge
 * weights.
 *
 * @param speed_of_light The speed of light, internal units.
 * @param weight_pe #feedback_props.band_edge_weight_pe.
 * @param weight_lw #feedback_props.band_edge_weight_lw.
 * @param photon_weight_lw #feedback_props.band_edge_photon_weight_lw.
 */
void radiation_packet_set_redshift(const double speed_of_light,
                                   const double weight_pe,
                                   const double weight_lw,
                                   const double photon_weight_lw) {
  radiation_packet_speed_of_light = speed_of_light;
  radiation_packet_weight_energy[ISRF_OPERATOR_PE] = weight_pe;
  radiation_packet_weight_energy[ISRF_OPERATOR_LW] = weight_lw;
  radiation_packet_weight_photon_lw = photon_weight_lw;
  radiation_packet_redshift_is_set = 1;
}

/* The speed limit of the light of a share, set once at start-up by
   radiation_packet_set_light_speed(). */
int radiation_packet_light_speed_limited = 0;
static double radiation_packet_arrival_per_length = 0.0;
static double radiation_packet_limit_speed = 0.0;
static int radiation_packet_ticks_are_ln_a = 0;

#ifdef SWIFT_DEBUG_CHECKS_ISRF
/*! Number of times the far cells of a particle became wider, of pieces of
    light that fell beyond the widest far cells, and of deliveries to a
    receiver that moves away faster than the clamp allows, since the last
    report. */
static long long radiation_packet_coarsening_count = 0;
static long long radiation_packet_overflow_count = 0;
static long long radiation_packet_clamp_count = 0;

/*! Count one widening of the far cells of a particle. */
static void radiation_packet_count_coarsening(void) {
  __sync_fetch_and_add(&radiation_packet_coarsening_count, 1);
}

/*! Count one piece of light that fell beyond the widest far cells. */
static void radiation_packet_count_overflow(void) {
  __sync_fetch_and_add(&radiation_packet_overflow_count, 1);
}

/*! Count one delivery to a receiver that moves away faster than the clamp
    allows. */
static void radiation_packet_count_clamp(void) {
  __sync_fetch_and_add(&radiation_packet_clamp_count, 1);
}
#else
static void radiation_packet_count_coarsening(void) {}
static void radiation_packet_count_overflow(void) {}
static void radiation_packet_count_clamp(void) {}
#endif

/**
 * @brief Limit the speed of the light of a share, or lift the limit.
 *
 * A share booked over the physical path `L` arrives `L / limit_speed` after
 * its sender's arrival time. That time is counted in integer times: a tick is
 * `time_base` in a static box and a step of `ln a` in a cosmological one,
 * where the delay is `H L / limit_speed` ticks of `ln a`, the rate of the
 * start of the hop.
 *
 * @param limit_speed The limit speed in internal units. 0 or less lifts the
 * limit.
 * @param time_base Length of one integer time, in time units (static) or in
 * `ln a` (cosmological).
 * @param cosmological Are the integer times steps of `ln a`?
 */
void radiation_packet_set_light_speed(const double limit_speed,
                                      const double time_base,
                                      const int cosmological) {
  radiation_packet_light_speed_limited = limit_speed > 0.0;
  radiation_packet_limit_speed = limit_speed > 0.0 ? limit_speed : 0.0;
  radiation_packet_arrival_per_length =
      limit_speed > 0.0 ? 1.0 / (limit_speed * time_base) : 0.0;
  radiation_packet_ticks_are_ln_a = cosmological;
}

/**
 * @brief Time, in integer times, that the light needs for a path at the limit
 * speed.
 *
 * In a cosmological run the delay uses the Hubble rate at the start of the
 * hop. With a falling rate (`q > 0`) the true delay in `ln a` is a little
 * shorter and the share arrives late by a second-order term. In an
 * accelerating phase (`q < 0`) the rate grows, the true delay is longer, and
 * the share arrives early by the term `zeta^2 q / 2` of `ln a`, with
 * `zeta = H path / c`. It is negligible for a hop of a few smoothing lengths
 * (`zeta` of order 1e-4 or less).
 *
 * @param path_phys Physical length of the path.
 * @param H Current Hubble rate. Read only in a cosmological run.
 * @return The delay, 0 if the speed is not limited.
 */
double radiation_packet_arrival_delay(const double path_phys, const double H) {
  if (!radiation_packet_light_speed_limited) return 0.0;
  return path_phys * radiation_packet_arrival_per_length *
         (radiation_packet_ticks_are_ln_a ? H : 1.0);
}

/**
 * @brief Light time of a path to a receiver that moves away from the sender.
 *
 * The offset to the origin and the arrival time are fixed at the delivery.
 * A receiver that moves away at the speed `u` along the path is reached by
 * the front at the speed `v - u`, so the delay is the one of
 * #radiation_packet_arrival_delay times `v / (v - u)`. A receiver that moves
 * towards the sender is treated as at rest. The factor is limited to
 * #RADIATION_PACKET_RECEDING_MAX_FACTOR: when the receiver moves away faster
 * than the fraction 1 - 1 / factor of `v`, the light is delivered at that
 * delay, which is earlier than the front reaches the receiver in the frame of
 * the origin. The limit applies in the gas frame. Each such delivery is
 * counted and reported by #radiation_packet_step_diagnostics.
 *
 * @param path_phys Physical length of the path.
 * @param H Current Hubble rate. Read only in a cosmological run.
 * @param closing_speed Physical speed of the receiver away from the sender,
 * along the path. Negative values are taken as 0.
 * @return The delay, 0 if the speed is not limited.
 */
double radiation_packet_arrival_delay_receding(const double path_phys,
                                               const double H,
                                               const double closing_speed) {
  const double base = radiation_packet_arrival_delay(path_phys, H);
  if (!(closing_speed > 0.0)) return base;
  const double v = radiation_packet_limit_speed;
  if (v - closing_speed < (1. / RADIATION_PACKET_RECEDING_MAX_FACTOR) * v) {
    radiation_packet_count_clamp();
    return base * RADIATION_PACKET_RECEDING_MAX_FACTOR;
  }
  return base * v / (v - closing_speed);
}

/**
 * @brief Start-up of the packet emission, from engine_config().
 *
 * @param e The #engine.
 */
void radiation_packet_startup(const struct engine *e) {
  if (e->feedback_props == NULL) return;
  radiation_packet_set_box(e->s->dim, e->s->periodic);
  const struct feedback_props *fp = e->feedback_props;
  if (fp->ISRF_packet_light_speed == isrf_packet_light_speed_unlimited)
    radiation_packet_set_light_speed(0.0, e->time_base, 0);
  else
    radiation_packet_set_light_speed(
        (double)fp->ISRF_packet_light_speed_fraction_of_c *
            e->physical_constants->const_speed_light_c,
        e->time_base, e->policy & engine_policy_cosmology);
  radiation_packet_set_redshift(e->physical_constants->const_speed_light_c,
                                e->feedback_props->band_edge_weight_pe,
                                e->feedback_props->band_edge_weight_lw,
                                e->feedback_props->band_edge_photon_weight_lw);
}

/**
 * @brief Redshift of a share over its path in an expanding box.
 *
 * With the redshift depth `zeta = H path_phys / c`, the energy of band `b`
 * falls as `exp(-w_b zeta)`, `w_b` the band-edge weight. The LW energy loses `w
 * - 1` of its `w` parts through the 11.2 eV edge into the PE band, as the fluid
 * transport does. The LW photon number falls as `exp(-w_N zeta)` and the
 * photons it loses enter the PE band in units of the PE lower edge. The PE
 * photon number loses photons at the rate `(w_PE - 1) E / (N E_lo)`, from its
 * own mean photon energy, which has no table value. Photons are `N E_lo`.
 *
 * @param H Current Hubble rate. 0 or less changes nothing.
 * @param path_phys Physical path length of the share.
 * @param energy (in/out) Energy of the share per band.
 * @param photons (in/out) Photon number of the share per band, `N E_lo`.
 * @param lost_energy (return) Energy removed from each band, minus the energy
 * it gained from the other band.
 * @param lost_photons (return) Photon number removed, as #lost_energy.
 */
void radiation_packet_redshift_share(const double H, const double path_phys,
                                     double energy[ISRF_OPERATOR_COUNT],
                                     double photons[ISRF_OPERATOR_COUNT],
                                     double lost_energy[ISRF_OPERATOR_COUNT],
                                     double lost_photons[ISRF_OPERATOR_COUNT]) {

  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    lost_energy[b] = 0.0;
    lost_photons[b] = 0.0;
  }
  if (!radiation_packet_redshift_is_set)
    error("The packet redshift was not given the speed of light.");
  const double zeta = H * path_phys / radiation_packet_speed_of_light;
  if (!(zeta > 0.0)) return;

  const double w_pe = radiation_packet_weight_energy[ISRF_OPERATOR_PE];
  const double w_lw = radiation_packet_weight_energy[ISRF_OPERATOR_LW];

  /* The PE photon rate is fixed by the energy and photons before the step. */
  double loss_pe_photons = 0.0;
  if (energy[ISRF_OPERATOR_PE] > 0.0 && photons[ISRF_OPERATOR_PE] > 0.0) {
    const double w_n_pe =
        (w_pe - 1.0) * energy[ISRF_OPERATOR_PE] / photons[ISRF_OPERATOR_PE];
    loss_pe_photons = -photons[ISRF_OPERATOR_PE] * expm1(-w_n_pe * zeta);
  }
  const double loss_pe = -energy[ISRF_OPERATOR_PE] * expm1(-w_pe * zeta);
  const double loss_lw = -energy[ISRF_OPERATOR_LW] * expm1(-w_lw * zeta);
  const double loss_lw_photons =
      -photons[ISRF_OPERATOR_LW] *
      expm1(-radiation_packet_weight_photon_lw * zeta);

  /* Gain of the PE band: what the LW band loses through its lower edge. */
  const double gain_pe = loss_lw * (w_lw - 1.0) / w_lw;
  const double gain_pe_photons =
      loss_lw_photons *
      (RADIATION_PE_BAND_LOWER_EDGE_CGS / RADIATION_LW_BAND_LOWER_EDGE_CGS);

  energy[ISRF_OPERATOR_PE] += gain_pe - loss_pe;
  energy[ISRF_OPERATOR_LW] -= loss_lw;
  photons[ISRF_OPERATOR_PE] += gain_pe_photons - loss_pe_photons;
  photons[ISRF_OPERATOR_LW] -= loss_lw_photons;

  lost_energy[ISRF_OPERATOR_PE] = loss_pe - gain_pe;
  lost_energy[ISRF_OPERATOR_LW] = loss_lw;
  lost_photons[ISRF_OPERATOR_PE] = loss_pe_photons - gain_pe_photons;
  lost_photons[ISRF_OPERATOR_LW] = loss_lw_photons;
}

/**
 * @brief Transmission weight `V W(r, h_s) max(0, mu)^p` of a sender towards
 * one neighbour, with `mu` the cosine between the hop and the packet
 * direction and `V = mass / rho` the comoving volume of the neighbour. The
 * expansion factor of the comoving density is common to all neighbours and
 * cancels in the shares.
 *
 * The sum over the neighbours is the SPH quadrature of the integral of
 * `W mu^p`, so the share to a side follows the solid angle and not the number
 * of neighbours there.
 *
 * Exactly 0 for `mu <= RADIATION_PACKET_MU_MIN`, for `r = 0` and beyond the
 * kernel support of the sender. Kept out of line so that the gradient loop and
 * the force loop evaluate the same compiled arithmetic, `mass / rho` and the
 * distance included: a caller's `r2` may be contracted differently at each
 * call site.
 *
 * @param dx Comoving vector from the sender to the neighbour.
 * @param h_s Comoving smoothing length of the sender.
 * @param n Unit direction of the packet.
 * @param power Power `p` of the cone.
 * @param mass Mass of the neighbour.
 * @param rho Comoving density snapshot `rho_prev` of the neighbour.
 * @return The weight.
 */
__attribute__((noinline)) float radiation_packet_transmission_weight(
    const float dx[3], const float h_s, const float n[3], const float power,
    const float mass, const float rho) {

  const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  const float H2 = h_s * h_s * kernel_gamma2;
  if (!(r2 > 0.f) || r2 >= H2) return 0.f;

  const float r = sqrtf(r2);
  const float mu = (dx[0] * n[0] + dx[1] * n[1] + dx[2] * n[2]) / r;
  if (mu <= RADIATION_PACKET_MU_MIN) return 0.f;

  const float h_inv = 1.f / h_s;
  float w;
  kernel_eval(r * h_inv, &w);
  return (mass / rho) * w * pow_dimension(h_inv) * powf(mu, power);
}

/**
 * @brief Octant of the fixed frame that holds a direction: the reception slot
 * of the light that travels along it.
 *
 * Bit `k` of the result is set when the component `k` of the direction in the
 * frame is zero or positive, so a direction on a border goes to the positive
 * side whatever the sign of its zero. Directions more than 90 degrees apart are
 * never in one octant; exactly 90 degrees can be. The octant of the opposite
 * direction is `7 - octant`, except on a border. Kept out of line, so that
 * every caller bins with the same compiled arithmetic. The components are
 * formed in double.
 *
 * @param frame Rotation matrix (row-major) from the box to the frame.
 * @param direction Direction in the box frame. Need not be a unit vector.
 * @return The octant, in [0, 8).
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
 * @brief Reception slot that holds a direction: the cell of the partition of
 * the sphere, fixed in the octant frame, that contains it.
 *
 * The cells have equal solid angle and are convex cones, so the mean of
 * directions in one cell is in it. With 8 slots the cells are the octants
 * (#radiation_packet_octant). With 12 they are the faces of a rhombic
 * dodecahedron: slot `4 m + 2 a + b`, `m` the frame axis of smallest absolute
 * coordinate (the lowest on a tie), `a` and `b` 1 for a coordinate that is zero
 * or positive on the other two axes in axis order. With 16 slot `o + 8 h`, `o`
 * the octant and `h` 1 when the absolute second frame coordinate exceeds the
 * absolute first. A tie goes to the same side whatever the sign of a zero.
 * Kept out of line, as #radiation_packet_octant.
 *
 * @param frame Rotation matrix (row-major) from the box to the frame.
 * @param direction Direction in the box frame. Need not be a unit vector.
 * @return The slot, in [0, RADIATION_PACKET_SLOTS).
 */
__attribute__((noinline)) int radiation_packet_slot(const float frame[9],
                                                    const float direction[3]) {

#if RADIATION_PACKET_SLOTS == 8
  return radiation_packet_octant(frame, direction);
#else
  double q[3];
  for (int k = 0; k < 3; k++)
    q[k] = (double)frame[3 * k + 0] * (double)direction[0] +
           (double)frame[3 * k + 1] * (double)direction[1] +
           (double)frame[3 * k + 2] * (double)direction[2];
#if RADIATION_PACKET_SLOTS == 16
  int slot = 0;
  for (int k = 0; k < 3; k++)
    if (q[k] >= 0.0) slot |= 1 << k;
  if (fabs(q[1]) > fabs(q[0])) slot |= 8;
  return slot;
#else
  int m = 0;
  if (fabs(q[1]) < fabs(q[m])) m = 1;
  if (fabs(q[2]) < fabs(q[m])) m = 2;
  const int i = (m == 0) ? 1 : 0;
  const int j = (m == 2) ? 1 : 2;
  return 4 * m + 2 * (q[i] >= 0.0) + (q[j] >= 0.0);
#endif
#endif
}

/**
 * @brief Distance of a direction from the nearest border of its slot, as the
 * sine of the angle to the nearest border plane.
 *
 * Used by the debug check of #radiation_packet_set_origin to skip directions
 * that rounding to float can move across a border. For 12 slots it is
 * conservative: it takes the border between every pair of absolute coordinates.
 *
 * @param frame Rotation matrix (row-major) from the box to the frame.
 * @param direction Direction in the box frame.
 * @return The distance, 0 on a border, at most 1.
 */
double radiation_packet_slot_border_distance(const float frame[9],
                                             const float direction[3]) {

  double a[3], norm2 = 0.0;
  for (int k = 0; k < 3; k++) {
    const double q = (double)frame[3 * k + 0] * (double)direction[0] +
                     (double)frame[3 * k + 1] * (double)direction[1] +
                     (double)frame[3 * k + 2] * (double)direction[2];
    a[k] = fabs(q);
    norm2 += q * q;
  }
  const double norm = sqrt(norm2);
  if (!(norm > 0.0)) return 0.0;
  double d = 1.0;
#if RADIATION_PACKET_SLOTS == 8 || RADIATION_PACKET_SLOTS == 16
  for (int k = 0; k < 3; k++) d = fmin(d, a[k]);
#endif
#if RADIATION_PACKET_SLOTS == 16
  d = fmin(d, fabs(a[0] - a[1]) * M_SQRT1_2);
#elif RADIATION_PACKET_SLOTS == 12
  for (int k = 0; k < 3; k++)
    d = fmin(d, fabs(a[k] - a[(k + 1) % 3]) * M_SQRT1_2);
#endif
  return d / norm;
}

/**
 * @brief Build the rotation matrix of the octant frame.
 *
 * The matrix is `Rz(a_z) Ry(a_y) Rx(a_x)`. A positive seed offset composes it
 * with a rotation drawn uniformly from that seed (Shoemake 1992), so that
 * repeats of a run use independent frames.
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
 * A beam is #RADIATION_PACKET_TEST_BEAM_FIELDS numbers in internal units, in
 * the order of the RADIATION_PACKET_TEST_BEAM_* offsets. Absent or with zero
 * energy, a beam is off. The direction is normalised. Finiteness is tested on
 * the bit pattern, which -ffinite-math-only cannot fold away.
 *
 * @param params The parameter file.
 * @param beams (return) The beams.
 * @return The number of active beams.
 */
int radiation_packet_read_test_beams(
    struct swift_params *params,
    float beams[RADIATION_PACKET_TEST_BEAMS]
               [RADIATION_PACKET_TEST_BEAM_FIELDS]) {

  int n_active = 0;
  memset(beams, 0,
         sizeof(float) * RADIATION_PACKET_TEST_BEAM_FIELDS *
             RADIATION_PACKET_TEST_BEAMS);
  parser_get_opt_param_float_array(params,
                                   "GEARFeedback:ISRF_packet_test_beam_1",
                                   RADIATION_PACKET_TEST_BEAM_FIELDS, beams[0]);
  parser_get_opt_param_float_array(params,
                                   "GEARFeedback:ISRF_packet_test_beam_2",
                                   RADIATION_PACKET_TEST_BEAM_FIELDS, beams[1]);

  for (int i = 0; i < RADIATION_PACKET_TEST_BEAMS; i++) {
    float *bm = beams[i];
    float *lower = &bm[RADIATION_PACKET_TEST_BEAM_LOWER];
    float *upper = &bm[RADIATION_PACKET_TEST_BEAM_UPPER];
    float *direction = &bm[RADIATION_PACKET_TEST_BEAM_DIRECTION];
    const float energy_PE = bm[RADIATION_PACKET_TEST_BEAM_ENERGY_PE];
    const float energy_LW = bm[RADIATION_PACKET_TEST_BEAM_ENERGY_LW];

    for (int k = 0; k < RADIATION_PACKET_TEST_BEAM_FIELDS; k++) {
      uint32_t bits;
      memcpy(&bits, &bm[k], sizeof(bits));
      volatile uint32_t opaque = bits;
      if (((opaque >> 23) & 0xFFu) == 0xFFu)
        error("GEARFeedback:ISRF_packet_test_beam_%d: non-finite entry %d.",
              i + 1, k);
    }
    if (energy_PE < 0.f || energy_LW < 0.f)
      error("GEARFeedback:ISRF_packet_test_beam_%d: negative energy.", i + 1);
    if (energy_PE == 0.f && energy_LW == 0.f) continue;
    const float norm =
        sqrtf(direction[0] * direction[0] + direction[1] * direction[1] +
              direction[2] * direction[2]);
    if (!(norm > 0.f))
      error("GEARFeedback:ISRF_packet_test_beam_%d: zero direction.", i + 1);
    for (int k = 0; k < 3; k++) direction[k] /= norm;
    for (int k = 0; k < 3; k++)
      if (!(upper[k] > lower[k]))
        error(
            "GEARFeedback:ISRF_packet_test_beam_%d: the box has upper corner "
            "%g not above lower corner %g on axis %d.",
            i + 1, upper[k], lower[k], k);
    n_active++;
    message(
        "ISRF packet test beam %d is ACTIVE: energy %g (PE) and %g (LW) per "
        "particle in the box, injected once at each particle's first merge, "
        "in no ledger.",
        i + 1, energy_PE, energy_LW);
  }
  return n_active;
}

/**
 * @brief Give the particles inside a test beam box their beam pulse, once.
 *
 * The pulse goes into the inbox of the slot of the beam direction, as a share
 * whose origin lies one smoothing length behind the particle. Its photon number
 * is its energy (mean photon energy at the band lower edge). The energy is in
 * no ledger: a test adds it up from the parameters. Called before the first
 * merge of the particle.
 *
 * @param p The particle.
 * @param beams The beams, see #radiation_packet_read_test_beams.
 * @param ti_current Current integer time, the arrival time of the pulse.
 */
void radiation_packet_inject_test_beams(
    struct part *p,
    const float beams[RADIATION_PACKET_TEST_BEAMS]
                     [RADIATION_PACKET_TEST_BEAM_FIELDS],
    const integertime_t ti_current) {

  for (int i = 0; i < RADIATION_PACKET_TEST_BEAMS; i++) {
    const float *bm = beams[i];
    const float *lower = &bm[RADIATION_PACKET_TEST_BEAM_LOWER];
    const float *upper = &bm[RADIATION_PACKET_TEST_BEAM_UPPER];
    const float *direction = &bm[RADIATION_PACKET_TEST_BEAM_DIRECTION];
    const double energy[ISRF_OPERATOR_COUNT] = {
        bm[RADIATION_PACKET_TEST_BEAM_ENERGY_PE],
        bm[RADIATION_PACKET_TEST_BEAM_ENERGY_LW]};
    if (energy[0] == 0. && energy[1] == 0.) continue;

    int inside = 1;
    for (int k = 0; k < 3; k++)
      if (!(p->x[k] >= (double)lower[k] && p->x[k] < (double)upper[k]))
        inside = 0;
    if (!inside) continue;

    const float offset[3] = {-p->h * direction[0], -p->h * direction[1],
                             -p->h * direction[2]};
    radiation_packet_deliver(p, p->feedback_data.isrf_packet.octant_frame,
                             energy, energy, offset, 0.f, 0.f, 0, 0,
                             (double)ti_current, (double)ti_current);
  }
}

/**
 * @brief Emission cone of a gas neighbour of a star: the HEALPix nside 1
 * pixel of the direction from the star to the gas, in a random frame drawn
 * from the star's ID and the start of its step.
 *
 * Kept out of line, so that the density loop and the feedback loop of the star,
 * on any rank, use the same arithmetic.
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
 * float @p dx of the runner only picks the periodic image. Both star loops call
 * this out-of-line function on the same positions, so they agree bit for bit on
 * the cone, the weight and the receivers, whatever path the runner took.
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

#ifdef SWIFT_DEBUG_CHECKS_ISRF
/**
 * @brief Counters of the arrival bins since the last report, for tests.
 *
 * @param widened (return) Times the far cells of a particle became wider.
 * @param beyond (return) Pieces of light beyond the widest far cells.
 * @param clamped (return) Deliveries where the clamp of the delay of a
 * receding receiver acted (#radiation_packet_arrival_delay_receding).
 */
void radiation_packet_bin_counters(long long *widened, long long *beyond,
                                   long long *clamped) {
  *widened = radiation_packet_coarsening_count;
  *beyond = radiation_packet_overflow_count;
  *clamped = radiation_packet_clamp_count;
}

/**
 * @brief Print, once per step, how many times the far arrival cells had to
 * widen, how many pieces of light were beyond the widest far cells and how
 * many deliveries went to a receiver that moved away faster than the clamp of
 * the delay allows, summed over the ranks, then reset the counters of every
 * rank. Prints nothing if all are zero. Collective: called by every rank from
 * engine_step().
 *
 * @param e The #engine.
 */
void radiation_packet_step_diagnostics(const struct engine *e) {
  if (!radiation_packet_light_speed_limited) return;
  long long counts[3] = {radiation_packet_coarsening_count,
                         radiation_packet_overflow_count,
                         radiation_packet_clamp_count};
  radiation_packet_coarsening_count = 0;
  radiation_packet_overflow_count = 0;
  radiation_packet_clamp_count = 0;
#ifdef WITH_MPI
  long long total[3];
  MPI_Reduce(counts, total, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  if (e->nodeID == 0)
    for (int k = 0; k < 3; k++) counts[k] = total[k];
#endif
  if (e->nodeID == 0 && (counts[0] > 0 || counts[1] > 0 || counts[2] > 0))
    message(
        "ISRF packet arrival, step %d: far cells widened %lld times, %lld "
        "pieces of light beyond the widest far cells, %lld deliveries to a "
        "receiver that moves away faster than %g of the limit speed (the "
        "light is delivered early).",
        e->step, counts[0], counts[1], counts[2],
        1. - 1. / RADIATION_PACKET_RECEDING_MAX_FACTOR);
}
#endif /* SWIFT_DEBUG_CHECKS_ISRF */

/**
 * @brief Largest float not above a distance, for a stored distance below the
 * end of a bin or a cell: a smaller distance is a later time, which is the
 * safe side. The result is at most @p x and within one float step of it.
 *
 * @param x The distance.
 */
static float radiation_packet_round_down(const double x) {
  /* Through memory: -ffast-math folds (double)(float)x to x, so that the float
     is rounded to the nearest and never corrected downwards. */
  volatile float rounded = (float)x;
  float f = rounded;
  if ((double)f > x) f = nextafterf(f, -FLT_MAX);
  return f;
}

/**
 * @brief Add a piece of a share to a near arrival bin.
 *
 * @param sl The slot.
 * @param j The bin.
 * @param energy Energy of the piece per band.
 * @param photons Photon number of the piece per band.
 * @param below_end Distance from the end of the bin to the latest arrival time
 * of the piece, rounded down.
 */
static void radiation_packet_add_to_bin(
    struct feedback_isrf_packet_slot_data *sl, const int j,
    const float energy[ISRF_OPERATOR_COUNT],
    const float photons[ISRF_OPERATOR_COUNT], const float below_end) {
  float held = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    held += sl->band[b].bin_energy[j];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    sl->band[b].bin_energy[j] += energy[b];
    sl->band[b].bin_photons[j] += photons[b];
  }
  sl->bin_stamp_below_end[j] =
      held > 0.f ? fminf(sl->bin_stamp_below_end[j], below_end) : below_end;
}

/**
 * @brief Add a piece of a share to a far arrival cell.
 *
 * The first piece in a cell keeps its interval: its light arrives evenly
 * between @p lo and @p hi. A second piece makes the cell hold several shares,
 * of which the cell keeps one box (the later first time, the later last time)
 * if both are boxes inside the cell, and otherwise only the latest arrival
 * time. A piece that carries light of
 * earlier pieces (below the photon cutoff) keeps the interval of its own
 * light: the carried light arrived before it, so it leaves no earlier than it
 * has arrived.
 *
 * @param sl The slot.
 * @param c The cell.
 * @param energy Energy of the piece per band.
 * @param photons Photon number of the piece per band.
 * @param lo First arrival time of the piece.
 * @param hi Last arrival time of the piece.
 * @param end End of the cell.
 */
static void radiation_packet_add_to_far(
    struct feedback_isrf_packet_slot_data *sl, const int c,
    const float energy[ISRF_OPERATOR_COUNT],
    const float photons[ISRF_OPERATOR_COUNT], const double lo, const double hi,
    const double end) {
  float held = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    held += sl->band[b].far_energy[c];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    sl->band[b].far_energy[c] += energy[b];
    sl->band[b].far_photons[c] += photons[b];
  }
  const float last_below = radiation_packet_round_down(end - hi);
  if (held > 0.f) {
    /* Two boxes make one box from the later of the first times to the later
       of the last times: the light of each share arrives, along its line,
       at or above that box line, so the cell never releases more than has
       arrived. A cell that holds light beyond its end keeps only the latest
       arrival time. */
    const int box = sl->far_first_below[c] >= 0.f && hi <= end;
    sl->far_first_below[c] = box ? fminf(sl->far_first_below[c],
                                         radiation_packet_round_down(end - lo))
                                 : -1.f;
    sl->far_last_below[c] = fminf(sl->far_last_below[c], last_below);
  } else {
    /* Light that arrives beyond the end of the cell cannot be kept as a box
       inside it. */
    sl->far_first_below[c] =
        hi <= end ? radiation_packet_round_down(end - lo) : -1.f;
    sl->far_last_below[c] = last_below;
  }
}

/**
 * @brief Width of the near arrival bins: the length of the step.
 *
 * @param pk The packet of the receiver.
 */
static integertime_t radiation_packet_step_width(
    const struct feedback_isrf_packet_data *pk) {
  return pk->step_ti > 0 ? pk->step_ti : 1;
}

/**
 * @brief End of a near arrival bin on the receiver's step lattice.
 *
 * The bins are one step wide and lie on the multiples of the step. Bin 0 ends
 * at the first multiple at or after the end of the current step, and bin `j`
 * at `j` steps later. A lattice that does not move with the step keeps what
 * waits in a bin in place, so that its light comes nearer the merge time as
 * time passes. The bins never become wider: light that arrives beyond the
 * last bin goes to the far cells.
 *
 * @param pk The packet of the receiver.
 * @param j The bin.
 * @return The end of the bin, an integer time.
 */
integertime_t radiation_packet_bin_end(
    const struct feedback_isrf_packet_data *pk, const int j) {
  const integertime_t w = radiation_packet_step_width(pk);
  return ((pk->step_end_ti + w - 1) / w) * w + (integertime_t)j * w;
}

/*! Largest width of a far cell, in integer times. */
#define RADIATION_PACKET_MAX_FAR_WIDTH ((integertime_t)1 << 60)

/*! Largest far scale: the shift `s + 1` stays below the width of an integer
    time. */
#define RADIATION_PACKET_MAX_FAR_SCALE (8 * (int)sizeof(integertime_t) - 2)

/**
 * @brief Largest far scale for a step: the width stays below
 * #RADIATION_PACKET_MAX_FAR_WIDTH.
 *
 * @param w The length of the step.
 */
static int radiation_packet_far_scale_limit(const integertime_t w) {
  int s = 0;
  while (s < RADIATION_PACKET_MAX_FAR_SCALE &&
         w <= (RADIATION_PACKET_MAX_FAR_WIDTH >> (s + 1)))
    s++;
  return s;
}

/**
 * @brief End of a far arrival cell for a given scale.
 *
 * The cells are `2^scale` steps wide and lie on the multiples of that width.
 * Cell 0 ends at the first multiple after the end of the last near bin, so
 * that it holds the light that arrives beyond the near bins, and cell `c` at
 * `c` widths later.
 *
 * @param pk The packet of the receiver.
 * @param scale The far scale.
 * @param c The cell.
 * @return The end of the cell, an integer time.
 */
static integertime_t radiation_packet_far_end_at(
    const struct feedback_isrf_packet_data *pk, const int scale, const int c) {
  const integertime_t width = radiation_packet_step_width(pk) << scale;
  const integertime_t reach =
      radiation_packet_bin_end(pk, RADIATION_PACKET_BINS - 1);
  return (reach / width + 1 + (integertime_t)c) * width;
}

/**
 * @brief End of a far arrival cell.
 *
 * @param pk The packet of the receiver.
 * @param c The cell.
 * @return The end of the cell, an integer time.
 */
integertime_t radiation_packet_far_end(
    const struct feedback_isrf_packet_data *pk, const int c) {
  return radiation_packet_far_end_at(pk, pk->far_scale, c);
}

/*! Number of arrival containers of a slot: the near bins, then the far cells.
 */
#define RADIATION_PACKET_PIECES \
  (RADIATION_PACKET_BINS + RADIATION_PACKET_FAR_CELLS)

/**
 * @brief Ends of the near bins and of the far cells of a receiver, as doubles.
 *
 * Gives the same values as #radiation_packet_bin_end and
 * #radiation_packet_far_end with two divisions in all.
 *
 * @param pk The packet of the receiver.
 * @param ub (return) The ends: the near bins first, then the far cells.
 */
static void radiation_packet_piece_ends(
    const struct feedback_isrf_packet_data *pk,
    double ub[RADIATION_PACKET_PIECES]) {
  const integertime_t w = radiation_packet_step_width(pk);
  const integertime_t e0 = ((pk->step_end_ti + w - 1) / w) * w;
  for (int j = 0; j < RADIATION_PACKET_BINS; j++)
    ub[j] = (double)(e0 + (integertime_t)j * w);
  const integertime_t reach =
      e0 + (integertime_t)(RADIATION_PACKET_BINS - 1) * w;
  const integertime_t width = w << pk->far_scale;
  const integertime_t f0 = (reach / width + 1) * width;
  for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
    ub[RADIATION_PACKET_BINS + c] = (double)(f0 + (integertime_t)c * width);
}

/*! Light of one container (a near bin or a far cell) that waits to be
    placed again: it arrives evenly between `first` and `last`. */
struct radiation_packet_record {
  float energy[ISRF_OPERATOR_COUNT];
  float photons[ISRF_OPERATOR_COUNT];
  double first;
  double last;
};

/**
 * @brief Does the slot wait for light in any near bin or far cell?
 *
 * @param sl The slot.
 */
static int radiation_packet_slot_is_pending(
    const struct feedback_isrf_packet_slot_data *sl) {
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    for (int j = 0; j < RADIATION_PACKET_BINS; j++)
      if (sl->band[b].bin_energy[j] != 0.f || sl->band[b].bin_photons[j] != 0.f)
        return 1;
    for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
      if (sl->band[b].far_energy[c] != 0.f || sl->band[b].far_photons[c] != 0.f)
        return 1;
  }
  return 0;
}

/**
 * @brief Take the light of a near bin out of a slot as a pulse at its latest
 * arrival time.
 *
 * @param sl The slot.
 * @param j The bin.
 * @param end The end of the bin, on the lattice that the light was stored in.
 * @param rec (return) The light.
 * @return Did the bin hold light?
 */
static int radiation_packet_take_bin(struct feedback_isrf_packet_slot_data *sl,
                                     const int j, const integertime_t end,
                                     struct radiation_packet_record *rec) {
  int occupied = 0;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
    rec->energy[b] = band->bin_energy[j];
    rec->photons[b] = band->bin_photons[j];
    if (rec->energy[b] != 0.f || rec->photons[b] != 0.f) occupied = 1;
    band->bin_energy[j] = 0.f;
    band->bin_photons[j] = 0.f;
  }
  rec->last = radiation_packet_time_later((double)end -
                                          (double)sl->bin_stamp_below_end[j]);
  rec->first = rec->last;
  sl->bin_stamp_below_end[j] = 0.f;
  return occupied;
}

/**
 * @brief Take the light of a far cell out of a slot: a box if the cell holds
 * one share, a pulse at its latest arrival time if it holds several.
 *
 * @param sl The slot.
 * @param c The cell.
 * @param end The end of the cell, on the lattice that the light was stored in.
 * @param rec (return) The light.
 * @return Did the cell hold light?
 */
static int radiation_packet_take_far(struct feedback_isrf_packet_slot_data *sl,
                                     const int c, const integertime_t end,
                                     struct radiation_packet_record *rec) {
  int occupied = 0;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
    rec->energy[b] = band->far_energy[c];
    rec->photons[b] = band->far_photons[c];
    if (rec->energy[b] != 0.f || rec->photons[b] != 0.f) occupied = 1;
    band->far_energy[c] = 0.f;
    band->far_photons[c] = 0.f;
  }
  rec->last =
      radiation_packet_time_later((double)end - (double)sl->far_last_below[c]);
  rec->first = sl->far_first_below[c] >= 0.f
                   ? radiation_packet_time_later((double)end -
                                                 (double)sl->far_first_below[c])
                   : rec->last;
  sl->far_first_below[c] = 0.f;
  sl->far_last_below[c] = 0.f;
  return occupied;
}

/**
 * @brief Put a share into the arrival bins and cells of a slot, on the
 * current lattices of the receiver.
 *
 * The near bins hold the light that arrives by the end of the bin: bin 0 the
 * light that arrives by its end `E`, bin `j` the light in `(E + (j - 1) w, E +
 * j w]`, `w` the step. The receiver releases a bin whole at its first merge
 * at or after the end of the bin, so nothing leaves before it has arrived. The
 * far cells hold the light beyond the last bin. A share arrives evenly between
 * @p first and @p last. A pulse goes whole into its bin or cell. A longer
 * share is split by its overlap with the bins and cells, in double: each piece
 * is the difference of two cumulative values, so no piece is negative, and the
 * sum of the pieces is the share. A piece below the photon cutoff in a band
 * joins the next later piece, never an earlier one. Light beyond the last far
 * cell goes to it and counts as overflow.
 *
 * @param sl The slot.
 * @param ub The ends of the near bins and the far cells of the receiver
 * (#radiation_packet_piece_ends).
 * @param energy Energy of the share per band.
 * @param photons Photon number of the share per band.
 * @param first First arrival time of the share.
 * @param last Last arrival time of the share.
 */
static void radiation_packet_place(struct feedback_isrf_packet_slot_data *sl,
                                   const double ub[RADIATION_PACKET_PIECES],
                                   const float energy[ISRF_OPERATOR_COUNT],
                                   const float photons[ISRF_OPERATOR_COUNT],
                                   const double first, const double last) {

  const int n = RADIATION_PACKET_PIECES;

  /* Cumulative fraction of the share in the pieces 0..j, and beyond (j = n). */
  double cum[RADIATION_PACKET_PIECES + 1];
  if (!(last > first)) {
    int j = 0;
    while (j < n && ub[j] < first) j++;
    for (int i = 0; i <= n; i++) cum[i] = i >= j ? 1. : 0.;
  } else {
    const double length = last - first;
    for (int i = 0; i < n; i++)
      cum[i] = fmin(1., fmax(0., (fmin(last, ub[i]) - first) / length));
    cum[n] = 1.;
  }
  int last_piece = n;
  for (int i = 0; i <= n; i++)
    if (cum[i] >= 1.) {
      last_piece = i;
      break;
    }

  /* A piece that is not the last is rounded down, so that no part of the
     share is placed in an earlier piece than its arrival allows. The last
     piece takes the rest. */
  float carry_e[ISRF_OPERATOR_COUNT], carry_n[ISRF_OPERATOR_COUNT];
  double given_e[ISRF_OPERATOR_COUNT], given_n[ISRF_OPERATOR_COUNT];
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    carry_e[b] = 0.f;
    carry_n[b] = 0.f;
    given_e[b] = 0.;
    given_n[b] = 0.;
  }

  for (int j = 0; j <= last_piece; j++) {
    float pe[ISRF_OPERATOR_COUNT], pn[ISRF_OPERATOR_COUNT];
    int any = 0;
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      const double share = cum[j] - (j > 0 ? cum[j - 1] : 0.);
      if (j == last_piece) {
        pe[b] = (float)((double)energy[b] - given_e[b]);
        pn[b] = (float)((double)photons[b] - given_n[b]);
      } else {
        pe[b] = radiation_packet_round_down((double)energy[b] * share);
        pn[b] = radiation_packet_round_down((double)photons[b] * share);
      }
      given_e[b] += (double)pe[b];
      given_n[b] += (double)pn[b];
      pe[b] += carry_e[b];
      pn[b] += carry_n[b];
      carry_e[b] = 0.f;
      carry_n[b] = 0.f;
      /* A piece below the cutoff joins the next later piece. */
      if (j != last_piece && pn[b] > 0.f &&
          pn[b] < RADIATION_PACKET_MIN_PHOTONS) {
        carry_e[b] = pe[b];
        carry_n[b] = pn[b];
        pe[b] = 0.f;
        pn[b] = 0.f;
      }
      if (pe[b] != 0.f || pn[b] != 0.f) any = 1;
    }
    if (!any) continue;
    if (j < RADIATION_PACKET_BINS) {
      const double stamp = fmin(last, ub[j]);
      radiation_packet_add_to_bin(sl, j, pe, pn,
                                  radiation_packet_round_down(ub[j] - stamp));
    } else if (j < n) {
      radiation_packet_add_to_far(sl, j - RADIATION_PACKET_BINS, pe, pn,
                                  fmax(first, ub[j - 1]), fmin(last, ub[j]),
                                  ub[j]);
    } else {
      radiation_packet_add_to_far(sl, RADIATION_PACKET_FAR_CELLS - 1, pe, pn,
                                  fmax(first, ub[n - 1]), last, ub[n - 1]);
      radiation_packet_count_overflow();
    }
  }
}

/**
 * @brief Move the waiting light of a particle to the lattices of its new step
 * and, when needed, to another far scale.
 *
 * Call it after the step of the packet changed (or, with @p extra_last, before
 * a share is placed whose light arrives beyond the far cells). The far scale
 * becomes the smallest at which all waiting light, and the share, fits the
 * far cells. Light that waits in a near bin keeps its latest arrival time,
 * light in a far cell that holds one share keeps its interval, and light in a
 * far cell that holds several shares becomes a pulse at its latest arrival
 * time: it is never moved earlier. What the new lattices cannot hold in the
 * near bins goes to the far cells.
 *
 * @param pk The packet of the particle.
 * @param old_first End of bin 0 on the previous lattice.
 * @param old_width Width of the previous near bins.
 * @param old_far_end Ends of the far cells on the previous lattice.
 * @param old_scale The previous far scale.
 * @param extra_last Last arrival time of a share about to be placed, or a
 * negative value for none.
 */
static void radiation_packet_regrid(
    struct feedback_isrf_packet_data *pk, const integertime_t old_first,
    const integertime_t old_width,
    const integertime_t old_far_end[RADIATION_PACKET_FAR_CELLS],
    const int old_scale, const double extra_last) {

  const int n = RADIATION_PACKET_BINS;
  const integertime_t w = radiation_packet_step_width(pk);
  const integertime_t new_first = radiation_packet_bin_end(pk, 0);

  /* The latest arrival time of everything that waits. */
  double need = extra_last;
  int pending = extra_last >= 0.;
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
    const struct feedback_isrf_packet_slot_data *sl = &pk->slot[s];
    if (!radiation_packet_slot_is_pending(sl)) continue;
    pending = 1;
    for (int j = 0; j < n; j++) {
      float e = 0.f;
      for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
        e += fabsf(sl->band[b].bin_energy[j]) +
             fabsf(sl->band[b].bin_photons[j]);
      if (e != 0.f)
        need = fmax(need, (double)(old_first + (integertime_t)j * old_width) -
                              (double)sl->bin_stamp_below_end[j]);
    }
    for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++) {
      float e = 0.f;
      for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
        e += fabsf(sl->band[b].far_energy[c]) +
             fabsf(sl->band[b].far_photons[c]);
      if (e != 0.f)
        need =
            fmax(need, (double)old_far_end[c] - (double)sl->far_last_below[c]);
    }
  }
  if (!pending) {
    pk->far_scale = 0;
    return;
  }

  /* The smallest far scale at which it fits. */
  const int limit = radiation_packet_far_scale_limit(w);
  int scale = 0;
  if (need > (double)radiation_packet_bin_end(pk, n - 1))
    while (scale < limit &&
           (double)radiation_packet_far_end_at(
               pk, scale, RADIATION_PACKET_FAR_CELLS - 1) < need)
      scale++;

  const int near_changed = new_first != old_first || w != old_width;
  const int far_changed =
      scale != old_scale || w != old_width ||
      radiation_packet_far_end_at(pk, scale, 0) != old_far_end[0];
  pk->far_scale = (unsigned char)scale;
  if (!near_changed && !far_changed) return;

  double ub[RADIATION_PACKET_PIECES];
  radiation_packet_piece_ends(pk, ub);
  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++) {
    struct feedback_isrf_packet_slot_data *sl = &pk->slot[s];
    if (!radiation_packet_slot_is_pending(sl)) continue;
    struct radiation_packet_record rec[RADIATION_PACKET_PIECES];
    int count = 0;
    if (near_changed)
      for (int j = 0; j < n; j++)
        if (radiation_packet_take_bin(
                sl, j, old_first + (integertime_t)j * old_width, &rec[count]))
          count++;
    if (far_changed)
      for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
        if (radiation_packet_take_far(sl, c, old_far_end[c], &rec[count]))
          count++;
    for (int i = 0; i < count; i++)
      radiation_packet_place(sl, ub, rec[i].energy, rec[i].photons,
                             rec[i].first, rec[i].last);
  }
}

/**
 * @brief Put a share into the arrival bins and cells of a slot.
 *
 * If the last light of the share is beyond the last far cell, the far cells
 * become wider until it fits (#radiation_packet_regrid), up to a limit. The
 * near bins never widen, so light that arrives soon is released as it would
 * be without the far light. See #radiation_packet_place for the placing.
 *
 * Limit of the far light. Two far cells cannot keep the separate arrival times
 * of many shares, so several shares in a cell leave together when the cell is
 * due (they are never early). A steady source is steady only up to a delay.
 * The largest |held - mean| / mean of a receiver fed by one steady sender of
 * its own step is 0 for delays up to 6 receiver steps, 1 at 7, 3 at 8, 7 at
 * 10, 31 at 25 and 63 at 40; for a sender 2 or 4 times coarser it is 0 up to 4
 * steps, about 1 at 4.5 to 6 and up to 3 at 7 to 8. It does not occur at the
 * default limit of 0.99 c, nor in the fixed mode with margin 0.5 except for the
 * outermost neighbours.
 *
 * @param pk The packet of the receiver.
 * @param sl The slot.
 * @param energy Energy of the share per band.
 * @param photons Photon number of the share per band.
 * @param first First arrival time of the share.
 * @param last Last arrival time of the share.
 */
void radiation_packet_deposit(struct feedback_isrf_packet_data *pk,
                              struct feedback_isrf_packet_slot_data *sl,
                              const float energy[ISRF_OPERATOR_COUNT],
                              const float photons[ISRF_OPERATOR_COUNT],
                              const double first, const double last) {

  double ub[RADIATION_PACKET_PIECES];
  radiation_packet_piece_ends(pk, ub);
  if (last < (double)max_nr_timesteps &&
      last > ub[RADIATION_PACKET_PIECES - 1] &&
      pk->far_scale <
          radiation_packet_far_scale_limit(radiation_packet_step_width(pk))) {
    integertime_t old_far_end[RADIATION_PACKET_FAR_CELLS];
    for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
      old_far_end[c] = radiation_packet_far_end(pk, c);
    radiation_packet_regrid(pk, radiation_packet_bin_end(pk, 0),
                            radiation_packet_step_width(pk), old_far_end,
                            pk->far_scale, last);
    radiation_packet_count_coarsening();
    radiation_packet_piece_ends(pk, ub);
  }
  radiation_packet_place(sl, ub, energy, photons, first, last);
}

/**
 * @brief Book the bands of a bucket whose whole photon number is below the
 * cutoff as absorbed, as the delivery does.
 *
 * @param energy (in/out) Energy of the bucket in a band.
 * @param photons (in/out) Photon number of the bucket in a band.
 * @param ledger The ledger of the band.
 */
static void radiation_packet_absorb_below_cutoff(
    float *energy, float *photons,
    struct feedback_isrf_packet_band_data *ledger) {
  /* A negative value cannot be released either. */
  if ((*photons != 0.f || *energy != 0.f) &&
      *photons < RADIATION_PACKET_MIN_PHOTONS) {
    radiation_packet_ledger_absorbed(ledger, (double)*energy, (double)*photons);
    *energy = 0.f;
    *photons = 0.f;
  }
}

/**
 * @brief Set the origin and the direction of the held packet of a slot from
 * the energy-weighted sum of the origin offsets.
 *
 * @param sl The slot.
 * @param offset Sum of energy times offset of what was released.
 * @param total_energy Energy released.
 * @param slot_index Index of the slot, for the debugging check.
 * @param frame Octant frame, for the debugging check.
 */
static void radiation_packet_set_origin(
    struct feedback_isrf_packet_slot_data *sl, double offset[3],
    const float total_energy, const int slot_index, const float frame[9]) {
  for (int k = 0; k < 3; k++) offset[k] /= (double)total_energy;
  const double d = sqrt(offset[0] * offset[0] + offset[1] * offset[1] +
                        offset[2] * offset[2]);
  sl->origin_distance = (float)d;
  if (d > 0.0) {
    for (int k = 0; k < 3; k++) sl->direction[k] = (float)(-offset[k] / d);
#ifdef SWIFT_DEBUG_CHECKS
    /* A cell is convex, so the mean origin of its shares lies in it too. A
       direction near a border may flip by rounding to float. */
    if (radiation_packet_slot_border_distance(frame, sl->direction) > 1e-5 &&
        radiation_packet_slot(frame, sl->direction) != slot_index)
      error("The direction of slot %d is in slot %d.", slot_index,
            radiation_packet_slot(frame, sl->direction));
#else
    (void)slot_index;
    (void)frame;
#endif
  } else {
    /* Mean origin on the particle: any direction is as good. It can lie
       outside the cell of the slot, which is harmless: shares are binned by
       their own geometry. */
    sl->direction[0] = 0.f;
    sl->direction[1] = 0.f;
    sl->direction[2] = 1.f;
  }
}

/**
 * @brief Start of an active step without a speed limit, for one slot: the
 * inbox and the due part of the reservoir become the held packet.
 *
 * The origin is the energy-weighted mean of the arriving origins and the
 * direction points from it to the particle. The reservoir releases the part due
 * in the step that ends now, `dt / t_rem`, with `dt` and `t_rem` measured from
 * the start of that step, and all of it once `t_rem <= dt`. A share with `n dt`
 * to run is met by `n` merges that release `1/n`, `1/(n-1)` of the rest, ...,
 * so each takes `1/n` of it. Same with proper-time steps.
 *
 * Without a speed limit every share arrives in bin 0, the inbox.
 *
 * @param sl The slot.
 * @param ti_current Time of the merge, the end of the particle's last step.
 * @param ti_step Length of that last step.
 * @param time_base Physical length of one integer time step. Under cosmology
 * it is a step in `ln a`, so @p cosmo is used instead.
 * @param cosmo The cosmological model, or NULL in a static box.
 * @param slot_index Index of the slot, for the debugging check.
 * @param frame Octant frame, for the debugging check.
 */
static void radiation_packet_merge_slot_unlimited(
    struct feedback_isrf_packet_slot_data *sl, const integertime_t ti_current,
    const integertime_t ti_step, const double time_base,
    const struct cosmology *cosmo, const int slot_index, const float frame[9]) {

  float release = 1.f;
  if (sl->reservoir_end_ti > ti_current) {
    const integertime_t ti_begin = ti_current - ti_step;
    double t_rem, dt;
    if (cosmo != NULL) {
      t_rem = cosmology_get_delta_time(cosmo, ti_begin, sl->reservoir_end_ti);
      dt = cosmology_get_delta_time(cosmo, ti_begin, ti_current);
    } else {
      t_rem = (double)(sl->reservoir_end_ti - ti_begin) * time_base;
      dt = (double)ti_step * time_base;
    }
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
    const float arriving = band->bin_energy[0] + res_energy;
    band->energy += arriving;
    band->photons += band->bin_photons[0] + res_photons;
    band->bin_energy[0] = 0.f;
    band->bin_photons[0] = 0.f;
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
  radiation_packet_set_origin(sl, offset, total_energy, slot_index, frame);
}

/**
 * @brief Far cell of a slot at the start of an active step with a speed
 * limit: release the part of its light that has arrived to the held packet.
 *
 * A far cell that holds one share releases the fraction of it that has
 * arrived, as the light arrives evenly. A far cell that holds several shares
 * releases all of it when its latest light has arrived, which is never early
 * for any of them.
 *
 * @param sl The slot.
 * @param ledger The ledger of the particle, per band.
 * @param c The far cell.
 * @param far_end End of the far cell on the previous lattice.
 * @param now Time of the merge.
 * @param total_energy (in, out) Energy released to the held packet, summed.
 * @param stamp (in, out) Latest arrival time of what is released.
 */
static void radiation_packet_release_far_cell(
    struct feedback_isrf_packet_slot_data *sl,
    struct feedback_isrf_packet_band_data ledger[ISRF_OPERATOR_COUNT],
    const int c, const integertime_t far_end, const double now,
    float *total_energy, double *stamp) {

  float cell_total = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    cell_total += sl->band[b].far_energy[c];
  if (!(cell_total > 0.f)) return;

  const double end = (double)far_end;
  const double last =
      radiation_packet_time_later(end - (double)sl->far_last_below[c]);
  const int single = sl->far_first_below[c] >= 0.f;
  const double first =
      single ? radiation_packet_time_later(end - (double)sl->far_first_below[c])
             : last;
  double f = 0.;
  if (now >= last)
    f = 1.;
  else if (single && last > first && now > first)
    f = (now - first) / (last - first);

  double f_b[ISRF_OPERATOR_COUNT];
  int all = 1;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
    double fb = f;
    const double ph = (double)band->far_photons[c];
    /* A part below the cutoff waits. A remainder below the cutoff is
       absorbed with the release, as it cannot be kept. The part that
       leaves is rounded down, so that it is never above the part due. */
    if (ph > 0. && fb < 1.) {
      if (ph * fb < (double)RADIATION_PACKET_MIN_PHOTONS)
        fb = 0.;
      else if (ph * (1. - fb) < (double)RADIATION_PACKET_MIN_PHOTONS) {
        radiation_packet_ledger_absorbed(
            &ledger[b], (double)band->far_energy[c] * (1. - fb),
            ph * (1. - fb));
        fb = 1.;
        band->far_energy[c] =
            radiation_packet_round_down((double)band->far_energy[c] * f);
        band->far_photons[c] = radiation_packet_round_down(ph * f);
      }
    }
    f_b[b] = fb;
    if (fb < 1.) all = 0;
  }
  float e = 0.f;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
    float re = band->far_energy[c], rp = band->far_photons[c];
    if (f_b[b] < 1.) {
      re = radiation_packet_round_down((double)re * f_b[b]);
      rp = radiation_packet_round_down((double)rp * f_b[b]);
    }
    band->far_energy[c] -= re;
    band->far_photons[c] -= rp;
    band->energy += re;
    band->photons += rp;
    e += re;
  }
  if (e > 0.f) {
    *total_energy += e;
    *stamp = fmax(*stamp, fmin(last, now));
  }
  if (all) {
    sl->far_first_below[c] = 0.f;
    sl->far_last_below[c] = 0.f;
  } else if (single && f > 0.) {
    /* What is left arrives evenly from now to the same end. */
    sl->far_first_below[c] = radiation_packet_round_down(end - now);
  }
}

/**
 * @brief Start of an active step with a speed limit, for one slot: the near
 * bins that have ended and the part of the far cells that has arrived become
 * the held packet.
 *
 * A bin is released whole at the first merge at or after its end
 * (#radiation_packet_place). Times of the bins are integers and are compared
 * with the merge time as integers. A far cell that holds one share releases
 * the fraction of it that has arrived, as the light arrives evenly. A far cell
 * that holds several shares releases all of it when its latest light has
 * arrived, which is never early for any of them. The caller then moves what
 * waits to the lattices of the new step (#radiation_packet_regrid).
 *
 * The held packet is stamped with the latest arrival time of what is released,
 * so that no share booked from it leaves before its light. The origin of the
 * packet is the mean origin of everything that waited in the slot.
 *
 * @param sl The slot.
 * @param ledger The ledger of the particle, per band.
 * @param old_first End of bin 0 on the lattice of the previous step.
 * @param old_width Width of the near bins of the previous step.
 * @param old_far_end Ends of the far cells on the previous lattice.
 * @param ti_current Time of the merge, the end of the particle's last step.
 * @param slot_index Index of the slot, for the debugging check.
 * @param frame Octant frame, for the debugging check.
 */
static void radiation_packet_merge_slot_limited(
    struct feedback_isrf_packet_slot_data *sl,
    struct feedback_isrf_packet_band_data ledger[ISRF_OPERATOR_COUNT],
    const integertime_t old_first, const integertime_t old_width,
    const integertime_t old_far_end[RADIATION_PACKET_FAR_CELLS],
    const integertime_t ti_current, const int slot_index,
    const float frame[9]) {

  const int n = RADIATION_PACKET_BINS;
  const double now = (double)ti_current;

  /* A band whose whole bin or cell is below the cutoff is absorbed before
     anything else. */
  double pending = 0.;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
    for (int j = 0; j < n; j++) {
      radiation_packet_absorb_below_cutoff(&band->bin_energy[j],
                                           &band->bin_photons[j], &ledger[b]);
      pending += band->bin_energy[j];
    }
    for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++) {
      radiation_packet_absorb_below_cutoff(&band->far_energy[c],
                                           &band->far_photons[c], &ledger[b]);
      pending += band->far_energy[c];
    }
  }
  double held_before = 0.;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++)
    held_before += sl->band[b].energy;

  float total_energy = 0.f;
  double stamp = -1e300;

  /* The near bins that have ended. */
  for (int j = 0; j < n; j++) {
    const integertime_t end = old_first + (integertime_t)j * old_width;
    if (end > ti_current) continue;
    float e = 0.f;
    for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
      struct feedback_isrf_packet_slot_band_data *band = &sl->band[b];
      band->energy += band->bin_energy[j];
      band->photons += band->bin_photons[j];
      e += band->bin_energy[j];
      band->bin_energy[j] = 0.f;
      band->bin_photons[j] = 0.f;
    }
    if (e > 0.f) {
      total_energy += e;
      stamp =
          fmax(stamp, radiation_packet_time_later(
                          (double)end - (double)sl->bin_stamp_below_end[j]));
    }
  }

  /* The far cells. */
  for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
    radiation_packet_release_far_cell(sl, ledger, c, old_far_end[c], now,
                                      &total_energy, &stamp);

  /* Origin: the mean origin of everything that waited. */
  double offset[3] = {0., 0., 0.};
  if (total_energy > 0.f && pending > 0.) {
    const double fraction = fmin(1., (double)total_energy / pending);
    for (int k = 0; k < 3; k++) {
      offset[k] = sl->inbox_origin_offset[k] * fraction;
      sl->inbox_origin_offset[k] -= offset[k];
    }
  }
  double left = 0.;
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    for (int j = 0; j < n; j++) left += sl->band[b].bin_energy[j];
    for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
      left += sl->band[b].far_energy[c];
  }
  if (!(left > 0.))
    for (int k = 0; k < 3; k++) sl->inbox_origin_offset[k] = 0.0;

  if (!(total_energy > 0.f)) return;
  sl->arrival = held_before > 0. ? fmax(sl->arrival, stamp) : stamp;
  radiation_packet_set_origin(sl, offset, total_energy, slot_index, frame);
}

/**
 * @brief Start of an active step: in every slot, the arrival bins that have
 * ended and the due part of the far cells or the reservoir become the held
 * packet.
 *
 * Each share was binned into its slot on arrival, so nothing merges across
 * slots. With a speed limit what still waits moves to the lattices of the new
 * step.
 *
 * @param p The particle, active.
 * @param ti_current Time of the merge, the end of the particle's last step.
 * @param ti_step Length of the particle's last step. It sets the release of
 * the reservoir and stands in for the next step in `step_end_ti`.
 * @param time_base Physical length of one integer time step. Under cosmology
 * it is a step in `ln a`, so @p cosmo is used instead.
 * @param cosmo The cosmological model, or NULL in a static box.
 */
void radiation_packet_merge(struct part *p, const integertime_t ti_current,
                            const integertime_t ti_step, const double time_base,
                            const struct cosmology *cosmo) {

  struct feedback_isrf_packet_data *pk = &p->feedback_data.isrf_packet;
#ifdef SWIFT_DEBUG_CHECKS
  if (!radiation_packet_is_empty(pk))
    error("Particle %lld still holds a packet at its merge.", p->id);
#endif

  if (!radiation_packet_light_speed_limited) {
    pk->step_end_ti = ti_current + ti_step;
    pk->step_ti = ti_step;
    for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
      radiation_packet_merge_slot_unlimited(&pk->slot[s], ti_current, ti_step,
                                            time_base, cosmo, s,
                                            pk->octant_frame);
    return;
  }

  const integertime_t old_first = radiation_packet_bin_end(pk, 0);
  const integertime_t old_width = radiation_packet_step_width(pk);
  const int old_scale = pk->far_scale;
  integertime_t old_far_end[RADIATION_PACKET_FAR_CELLS];
  for (int c = 0; c < RADIATION_PACKET_FAR_CELLS; c++)
    old_far_end[c] = radiation_packet_far_end(pk, c);

  for (int s = 0; s < RADIATION_PACKET_SLOTS; s++)
    radiation_packet_merge_slot_limited(&pk->slot[s], pk->band, old_first,
                                        old_width, old_far_end, ti_current, s,
                                        pk->octant_frame);
  pk->step_end_ti = ti_current + ti_step;
  pk->step_ti = ti_step;
  radiation_packet_regrid(pk, old_first, old_width, old_far_end, old_scale,
                          -1.);
}

/**
 * @brief End of the gradient loop: normalise the forward advance of each slot
 * and read the own field off the held packets, or let a slot with no
 * downstream neighbour escape into the ledger.
 *
 * The held packet of a slot is the light that crosses the particle in one
 * step, over the forward advance of the slot, so its specific energy density is
 * `(E / m) l_par / (c dt)`. The own field is the sum over the slots, then its
 * running mean over `tau` (#radiation_packet_read_mean).
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
  if (!(fabsf(frame_norm - 3.f) < RADIATION_PACKET_FRAME_NORM_TOLERANCE))
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
      radiation_packet_ledger_escaped(&pk->band[b], (double)band->energy,
                                      (double)band->photons);
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
 * @brief Write the kernel mean of the own field into the band moments that
 * Grackle reads: PE to PE, LW to LW and the LW photon number (`N E_lo`) to
 * LW_PHOTON.
 *
 * Overwrites the moments, including any value from the initial conditions. Call
 * after #radiation_packet_end_force.
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
 * @brief Book the emission of a star for this feedback step, and as escaped the
 * light that its receivers do not take, once per step on the rank of the star.
 * Only the ledger is booked, so this does nothing without
 * SWIFT_DEBUG_CHECKS_ISRF.
 *
 * An empty cone loses its whole twelfth.
 *
 * @param sp The star, after its density loop.
 * @param Delta_t The star's feedback step.
 */
void radiation_packet_star_book_emission(struct spart *sp,
                                         const float Delta_t) {
#ifdef SWIFT_DEBUG_CHECKS
  sp->feedback_data.radiation.packet.cones_booked = 1;
#endif

#ifdef SWIFT_DEBUG_CHECKS_ISRF
  for (int b = 0; b < ISRF_OPERATOR_COUNT; b++) {
    const enum radiation_isrf_moment m = radiation_isrf_operator_owner[b];
    const double E = sp->feedback_data.radiation.L_band[m] * (double)Delta_t;
    sp->feedback_data.radiation.packet.emitted_energy[b] += E;
    sp->feedback_data.radiation.packet.emitted_photons[b] +=
        E * sp->feedback_data.radiation.packet.photon_fraction[b];
  }

  const double floor = radiation_packet_cone_floor(sp);
  const double fraction =
      star_gas_escape_fraction(sp->feedback_data.radiation.packet.cone_weight,
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
#endif
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Debugging checks: report whether the feedback loop delivered to the
 * same receivers as the density loop that summed the cone weights.
 *
 * Both loops take each weight from radiation_packet_star_gas_weight() on the
 * same positions, so the cone weight sums can differ only by summation order.
 * One line per star step reports the receivers, the empty cones (and those
 * with gas beyond #RADIATION_PACKET_EMISSION_EDGE only), the largest relative
 * difference of a cone weight sum, the worst membership difference (density
 * minus feedback loop), the PE energy of the step and the delivered plus
 * escaped fraction of the light, minus one.
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
  double dS = 0.0, S_d_sum = 0.0;
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
    S_d_sum += S_d;
    /* The worst membership difference: the largest count difference, else
       a swap at the same count. */
    if (abs(n - n_d) > abs(n_diff) ||
        (n_diff == 0 && id_diff == 0 && id != id_d)) {
      n_diff = n - n_d;
      id_diff = id - id_d;
    }
  }
  /* The delivered shares are those of the feedback loop's weight sum. */
  const double accounted =
      star_gas_global_share(S_d_sum,
                            sp->feedback_data.radiation.packet.cone_weight,
                            floor, RADIATION_PACKET_EMISSION_CONES) +
      star_gas_escape_fraction(sp->feedback_data.radiation.packet.cone_weight,
                               floor, RADIATION_PACKET_EMISSION_CONES);
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
