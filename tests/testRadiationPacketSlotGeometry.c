/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2026.
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

/**
 * @file testRadiationPacketSlotGeometry.c
 * @brief Geometry of the reception slots of the GEAR ISRF packet transport,
 * for the slot count of the build (8, 12, 16, 24 or 48).
 *
 * Checks the compiled radiation_packet_slot() against an independent
 * reference: the nearest of K generator directions (8, 12, 24, 48) or the
 * azimuth sector (16). The checks are: every direction maps to one slot in
 * range; each slot is one reference cell; the cells hold equal solid angle;
 * positive combinations of directions of one slot stay in it; the opposite
 * direction is in another slot; the result depends only on the bits of the
 * frame and of the direction; directions on a border get a slot whose closed
 * cell holds them, whatever the sign of a zero; the border distance never
 * exceeds the true one. A planted wrong slot function must fail the same
 * gates.
 */
#include <config.h>

/* Some standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

#if defined(GEAR_SUBGRID_RADIATION_ISRF)

#include "feedback/GEAR/radiation_isrf.h"
#include "feedback/GEAR/radiation_packet.h"

#define K RADIATION_PACKET_SLOTS

/*! Directions drawn per frame in the partition check. */
#define N_SAMPLES 4800000

/*! Bar of a slot population, in binomial standard deviations. */
#define POPULATION_SIGMA 5.0

/*! Directions closer than this to a border (sine of the angle) are not
    compared with the reference: float rounding of the frame product can move
    them across. */
#define BORDER_SKIP 1e-7

/** @brief Signature of a slot function under test. */
typedef int (*slot_function)(const float frame[9], const float direction[3]);

/**
 * @brief Next value of a splitmix64 stream, mapped to [0, 1).
 *
 * @param state The stream state.
 * @return A number in [0, 1).
 */
static double test_uniform(uint64_t *state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z ^= z >> 31;
  return (z >> 11) * 0x1.0p-53;
}

/**
 * @brief Uniform unit vector, rounded to float.
 *
 * @param state The stream state.
 * @param v (return) The vector.
 */
static void test_direction(uint64_t *state, float v[3]) {
  const double z = 2. * test_uniform(state) - 1.;
  const double phi = 2. * M_PI * test_uniform(state);
  const double s = sqrt(fmax(0., 1. - z * z));
  v[0] = (float)(s * cos(phi));
  v[1] = (float)(s * sin(phi));
  v[2] = (float)z;
}

/*! Generator directions of the reference cells: the cell of a generator is
    the set of directions closer to it than to any other generator. */
static double generators[48][3];
static int n_generators = 0;

/**
 * @brief Fill #generators for the slot count of the build: the cube vertices
 * (8), the cube edge midpoints (12), the permutations of (0, +-1, +-2) (24)
 * or of (+-1, +-2, +-3) (48). 16 uses the azimuth rule instead.
 */
static void build_generators(void) {
  n_generators = 0;
  const int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                          {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  double base[3];
  if (K == 8) {
    for (int s = 0; s < 8; s++)
      for (int k = 0; k < 3; k++) generators[s][k] = ((s >> k) & 1) ? 1. : -1.;
    n_generators = 8;
    return;
  }
  if (K == 12) {
    base[0] = 0.;
    base[1] = 1.;
    base[2] = 1.;
  } else if (K == 24) {
    base[0] = 0.;
    base[1] = 1.;
    base[2] = 2.;
  } else if (K == 48) {
    base[0] = 1.;
    base[1] = 2.;
    base[2] = 3.;
  } else {
    return;
  }
  for (int p = 0; p < 6; p++)
    for (int s = 0; s < 8; s++) {
      double g[3];
      for (int k = 0; k < 3; k++) {
        const double v = base[perm[p][k]];
        g[k] = ((s >> k) & 1) ? -v : v;
      }
      int dup = 0;
      for (int n = 0; n < n_generators && !dup; n++)
        dup = g[0] == generators[n][0] && g[1] == generators[n][1] &&
              g[2] == generators[n][2];
      if (!dup) {
        if (n_generators >= 48) error("Too many generators.");
        memcpy(generators[n_generators++], g, sizeof(g));
      }
    }
}

/**
 * @brief Coordinates of a float direction in the frame, in long double.
 *
 * @param F The frame, row-major.
 * @param d The direction.
 * @param q (return) The frame coordinates.
 */
static void frame_coordinates(const float F[9], const float d[3],
                              long double q[3]) {
  for (int k = 0; k < 3; k++)
    q[k] = (long double)F[3 * k] * d[0] + (long double)F[3 * k + 1] * d[1] +
           (long double)F[3 * k + 2] * d[2];
}

/**
 * @brief Reference cell of a direction and its distance from the cell border.
 *
 * @param F The frame.
 * @param d The direction.
 * @param dist (return) Sine of the angle to the nearest border plane.
 * @return The reference cell.
 */
static int ref_cell(const float F[9], const float d[3], long double *dist) {
  long double q[3];
  frame_coordinates(F, d, q);
  const long double norm = sqrtl(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
  if (K == 16) {
    const long double phi = atan2l(q[1], q[0]);
    const long double t = phi < 0.L ? phi + 2.L * M_PI : phi;
    const int sector = (int)(t / (M_PI / 4.L)) % 8;
    long double a0 = fabsl(q[0]), a1 = fabsl(q[1]);
    long double m = fabsl(q[2]);
    if (a0 < m) m = a0;
    if (a1 < m) m = a1;
    const long double h = fabsl(a0 - a1) * 0.7071067811865475244L;
    if (h < m) m = h;
    *dist = norm > 0.L ? m / norm : 0.L;
    return sector + 8 * (q[2] >= 0.L);
  }
  int best = 0;
  long double best_dot = -1e30L;
  for (int n = 0; n < n_generators; n++) {
    const long double v = generators[n][0] * q[0] + generators[n][1] * q[1] +
                          generators[n][2] * q[2];
    if (v > best_dot) {
      best_dot = v;
      best = n;
    }
  }
  long double m = 1e30L;
  for (int n = 0; n < n_generators; n++) {
    if (n == best) continue;
    long double w[3], w2 = 0.L;
    for (int k = 0; k < 3; k++) {
      w[k] = (long double)generators[best][k] - generators[n][k];
      w2 += w[k] * w[k];
    }
    const long double s = (w[0] * q[0] + w[1] * q[1] + w[2] * q[2]) / sqrtl(w2);
    if (s < m) m = s;
  }
  *dist = norm > 0.L ? m / norm : 0.L;
  return best;
}

/**
 * @brief Whether the closed reference cell of a candidate slot holds a
 * direction: its generator is within rounding of the best one.
 *
 * @param F The frame.
 * @param d The direction.
 * @param cell The reference cell of the slot.
 * @return 1 if the closed cell holds the direction.
 */
static int closed_cell_holds(const float F[9], const float d[3],
                             const int cell) {
  long double q[3];
  frame_coordinates(F, d, q);
  if (K == 16) {
    long double dist;
    const int c = ref_cell(F, d, &dist);
    if (c == cell) return 1;
    return dist < 1e-12L;
  }
  long double best = -1e30L;
  for (int n = 0; n < n_generators; n++) {
    const long double v = generators[n][0] * q[0] + generators[n][1] * q[1] +
                          generators[n][2] * q[2];
    if (v > best) best = v;
  }
  const long double v = generators[cell][0] * q[0] +
                        generators[cell][1] * q[1] + generators[cell][2] * q[2];
  return v >= best - 1e-12L;
}

/** @brief Production slot function. */
static int production_slot(const float frame[9], const float direction[3]) {
  return radiation_packet_slot(frame, direction);
}

/**
 * @brief Planted wrong slot function: the production one with the larger and
 * the smaller of two coordinates swapped in the sign bit, as a wrong index
 * formula would. Must fail the partition gates.
 */
static int planted_slot(const float frame[9], const float direction[3]) {
  const int s = radiation_packet_slot(frame, direction);
  const float flipped[3] = {direction[0], -direction[1], direction[2]};
  if (radiation_packet_slot(frame, flipped) % 2 == 0) return s ^ 1;
  return s;
}

/**
 * @brief The partition gates for one frame and one slot function.
 *
 * @param name Name of the frame, for the messages.
 * @param F The frame.
 * @param slot The slot function under test.
 * @param seed Seed of the direction stream.
 * @param verbose Print the results.
 * @return The number of failed gates.
 */
static int check_partition(const char *name, const float F[9],
                           const slot_function slot, const uint64_t seed,
                           const int verbose) {
  uint64_t st = seed;
  long count[K] = {0};
  int cell_of_slot[K];
  for (int s = 0; s < K; s++) cell_of_slot[s] = -1;
  long out_of_range = 0, cell_mismatch = 0, skipped = 0, antipode = 0;
  long border_overstated = 0;
  for (long i = 0; i < N_SAMPLES; i++) {
    float d[3];
    test_direction(&st, d);
    const int s = slot(F, d);
    if (s < 0 || s >= K) {
      out_of_range++;
      continue;
    }
    count[s]++;
    long double dist;
    const int c = ref_cell(F, d, &dist);
    const double bd = radiation_packet_slot_border_distance(F, d);
    if (!(bd >= 0.) || (long double)bd > dist + 1e-9L) border_overstated++;
    if (dist <= BORDER_SKIP) {
      skipped++;
      continue;
    }
    if (cell_of_slot[s] < 0) cell_of_slot[s] = c;
    if (cell_of_slot[s] != c) cell_mismatch++;
    const float e[3] = {-d[0], -d[1], -d[2]};
    if (slot(F, e) == s) antipode++;
  }

  int failed = 0;
  /* Each slot is one reference cell, and two slots never share one. */
  int duplicate_cells = 0, empty = 0;
  for (int s = 0; s < K; s++) {
    if (cell_of_slot[s] < 0) empty++;
    for (int r = 0; r < s; r++)
      if (cell_of_slot[s] >= 0 && cell_of_slot[s] == cell_of_slot[r])
        duplicate_cells++;
  }
  if (out_of_range || cell_mismatch || duplicate_cells || empty) failed++;
  if (antipode) failed++;
  if (border_overstated) failed++;

  /* Equal solid angle: the bar is POPULATION_SIGMA binomial deviations. */
  const double p = 1. / K;
  const double sigma = sqrt(N_SAMPLES * p * (1. - p));
  double worst = 0.;
  for (int s = 0; s < K; s++) {
    const double dev = fabs((double)count[s] - N_SAMPLES * p) / sigma;
    if (!(dev <= worst)) worst = dev;
  }
  if (!(worst <= POPULATION_SIGMA)) failed++;
  if (skipped > N_SAMPLES / 1000) failed++;

  if (verbose)
    message(
        "Frame %s, K = %d: %d directions, %ld out of range, %ld off their "
        "slot's cell, %d slots sharing a cell, %d slots never hit, %ld "
        "antipodes in the same slot, %ld border distances above the true "
        "one, %ld near a border; worst slot population %.2f sigma (bar %.1f, "
        "%.2f %% of the mean).",
        name, K, N_SAMPLES, out_of_range, cell_mismatch, duplicate_cells, empty,
        antipode, border_overstated, skipped, worst, POPULATION_SIGMA,
        100. * POPULATION_SIGMA * sigma / (N_SAMPLES * p));
  return failed;
}

/**
 * @brief Convexity: a positive combination of two or three directions of one
 * slot is in that slot, unless it is within rounding of a border.
 *
 * @param name Name of the frame.
 * @param F The frame.
 * @param slot The slot function under test.
 * @param verbose Print the result.
 * @return The number of failed gates.
 */
static int check_convexity(const char *name, const float F[9],
                           const slot_function slot, const int verbose) {
  uint64_t st = 0x1234567ULL;
  long tested = 0, violations = 0;
  const long n = 2000000;
  for (long i = 0; i < n; i++) {
    float a[3], b[3], c[3];
    test_direction(&st, a);
    test_direction(&st, b);
    test_direction(&st, c);
    /* Pull b and c towards a so that many triples share a slot. */
    const double pull = test_uniform(&st);
    for (int k = 0; k < 3; k++) {
      b[k] = (float)(a[k] + pull * (b[k] - a[k]));
      c[k] = (float)(a[k] + pull * (c[k] - a[k]));
    }
    const int s = slot(F, a);
    if (slot(F, b) != s || slot(F, c) != s) continue;
    const double wa = test_uniform(&st), wb = test_uniform(&st),
                 wc = (i % 2) ? test_uniform(&st) : 0.;
    float m[3];
    for (int k = 0; k < 3; k++)
      m[k] = (float)(wa * a[k] + wb * b[k] + wc * c[k]);
    long double dist;
    ref_cell(F, m, &dist);
    if (dist <= BORDER_SKIP) continue;
    tested++;
    if (slot(F, m) != s) violations++;
  }
  if (verbose)
    message(
        "Frame %s, K = %d: %ld positive combinations of one slot's "
        "directions, %ld outside that slot.",
        name, K, tested, violations);
  return (violations != 0) + (tested < 100000);
}

/**
 * @brief Purity: the slot depends only on the bits of the frame and of the
 * direction. Copies at other addresses, repeated calls and scaling by powers
 * of two (exact in float and in the double products) give the same slot.
 *
 * @param F The frame.
 */
static void check_purity(const float F[9]) {
  float *Fc = malloc(9 * sizeof(float));
  if (Fc == NULL) error("Allocation failed.");
  memcpy(Fc, F, 9 * sizeof(float));
  uint64_t st = 99;
  long differ = 0, n = 0;
  for (int i = 0; i < 200000; i++) {
    float d[3];
    test_direction(&st, d);
    /* Every hundredth direction is put exactly on a border of the identity
       frame, where a tie rule decides. */
    if (i % 100 == 0) d[1] = (i % 200 == 0) ? d[0] : -d[0];
    const int s = radiation_packet_slot(F, d);
    float *dc = malloc(3 * sizeof(float));
    if (dc == NULL) error("Allocation failed.");
    for (int e = -20; e <= 20; e += 5) {
      for (int k = 0; k < 3; k++) dc[k] = ldexpf(d[k], e);
      if (radiation_packet_slot(Fc, dc) != s) differ++;
      if (radiation_packet_slot(F, dc) != s) differ++;
      n += 2;
    }
    free(dc);
  }
  free(Fc);
  if (differ != 0)
    error("The slot changed for %ld of %ld copies or exact rescalings.", differ,
          n);
  message("Purity: %ld copies and exact rescalings, all in the same slot.", n);
}

/**
 * @brief Directions on borders in the identity frame: every vector with
 * coordinates in {-2, -1, -0, +0, 1, 2}. The slot must be a closed cell that
 * holds the direction, must not depend on the sign of a zero, and the zero
 * vector gets a slot in range.
 */
static void check_ties(void) {
  const float identity[9] = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
  const float vals[6] = {-2.f, -1.f, -0.f, 0.f, 1.f, 2.f};
  int cell_of_slot[K];
  for (int s = 0; s < K; s++) cell_of_slot[s] = -1;
  /* Map the slots to reference cells with generic directions first. */
  uint64_t st = 7;
  for (int i = 0; i < 1000000; i++) {
    float d[3];
    test_direction(&st, d);
    long double dist;
    const int c = ref_cell(identity, d, &dist);
    if (dist > BORDER_SKIP)
      cell_of_slot[radiation_packet_slot(identity, d)] = c;
  }
  for (int s = 0; s < K; s++)
    if (cell_of_slot[s] < 0)
      error("Slot %d never hit in the identity frame.", s);

  int checked = 0, on_border = 0;
  for (int a = 0; a < 6; a++)
    for (int b = 0; b < 6; b++)
      for (int c = 0; c < 6; c++) {
        const float d[3] = {vals[a], vals[b], vals[c]};
        const int s = radiation_packet_slot(identity, d);
        if (s < 0 || s >= K) error("Slot %d out of range.", s);
        /* The same vector with every zero of the other sign. */
        const float z[3] = {d[0] == 0.f ? -d[0] : d[0],
                            d[1] == 0.f ? -d[1] : d[1],
                            d[2] == 0.f ? -d[2] : d[2]};
        if (radiation_packet_slot(identity, z) != s)
          error("The sign of a zero moved (%g, %g, %g) to another slot.", d[0],
                d[1], d[2]);
        if (d[0] == 0.f && d[1] == 0.f && d[2] == 0.f) continue;
        if (!closed_cell_holds(identity, d, cell_of_slot[s]))
          error(
              "Border direction (%g, %g, %g) is in slot %d, whose cell does "
              "not hold it.",
              d[0], d[1], d[2], s);
        long double dist;
        ref_cell(identity, d, &dist);
        if (dist < 1e-12L) on_border++;
        checked++;
      }
  if (on_border == 0) error("No direction of the tie table is on a border.");
  message(
      "Ties: %d directions (%d on a border) in a closed cell, the sign of "
      "a zero never matters.",
      checked, on_border);
}

/**
 * @brief Hand table of the slot index formula in the identity frame, ties
 * included, for 24 and 48 slots.
 */
static void check_hand_table(void) {
#if K == 24 || K == 48
  const float identity[9] = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
  const struct {
    float d[3];
    int s24, s48;
  } table[] = {
      {{-3.f, 2.f, 1.f}, 3, 7},    {{1.f, -2.f, 3.f}, 20, 41},
      {{0.5f, -4.f, -1.f}, 8, 17}, {{1.f, 1.f, 0.f}, 7, 15},
      {{1.f, 1.f, 1.f}, 7, 15},    {{2.f, -2.f, 1.f}, 6, 13},
      {{-0.f, -0.f, 1.f}, 23, 47}, {{0.f, 0.f, 1.f}, 23, 47},
      {{0.f, 0.f, 0.f}, 7, 15},    {{-0.f, 0.f, -0.f}, 7, 15},
      {{0.f, 0.f, -1.f}, 19, 39},  {{-1.f, -1.f, -1.f}, 2, 4},
  };
  for (size_t i = 0; i < sizeof(table) / sizeof(table[0]); i++) {
    const int got = radiation_packet_slot(identity, table[i].d);
    const int want = (K == 24) ? table[i].s24 : table[i].s48;
    if (got != want)
      error("Slot of (%g, %g, %g): got %d, expected %d.", table[i].d[0],
            table[i].d[1], table[i].d[2], got, want);
  }
  message("Hand table, %zu directions: ok.", sizeof(table) / sizeof(table[0]));
#endif
}

/**
 * @brief Exact solid angle of one cell of 24 and 48 slots from its vertices
 * (Van Oosterom and Strackee 1983): 4 pi / K.
 */
static void check_exact_cell_area(void) {
#if K == 24 || K == 48
  const double r3 = 1. / sqrt(3.), r2 = 1. / sqrt(2.);
  /* 24: face centre and two cube vertices; 48: face centre, edge midpoint and
     cube vertex. */
  const double a[3] = {0., 0., 1.};
  const double b24[3] = {r3, r3, r3}, c24[3] = {-r3, r3, r3};
  const double b48[3] = {0., r2, r2}, c48[3] = {r3, r3, r3};
  const double *b = (K == 24) ? b24 : b48, *c = (K == 24) ? c24 : c48;
  const double triple = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                        a[1] * (b[0] * c[2] - b[2] * c[0]) +
                        a[2] * (b[0] * c[1] - b[1] * c[0]);
  const double ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const double bc = b[0] * c[0] + b[1] * c[1] + b[2] * c[2];
  const double ca = c[0] * a[0] + c[1] * a[1] + c[2] * a[2];
  const double omega = 2. * atan2(fabs(triple), 1. + ab + bc + ca);
  if (!(fabs(omega - 4. * M_PI / K) < 1e-12))
    error("Cell solid angle %.15f, expected 4 pi / %d = %.15f.", omega, K,
          4. * M_PI / K);
  message("Cell solid angle %.12f = 4 pi / %d.", omega, K);
#endif
}

/**
 * @brief Smallest angle between a box axis and a slot border in a frame.
 *
 * @param F The frame.
 * @return The angle, in degrees.
 */
static double box_axis_margin_deg(const float F[9]) {
  double worst = 90.;
  for (int k = 0; k < 6; k++) {
    float d[3] = {0.f, 0.f, 0.f};
    d[k % 3] = (k < 3) ? 1.f : -1.f;
    long double dist;
    ref_cell(F, d, &dist);
    const double deg = asin(fmin(1., (double)dist)) * 180. / M_PI;
    if (deg < worst) worst = deg;
  }
  return worst;
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

  build_generators();
  if (K != 16 && n_generators != K)
    error("%d reference generators for %d slots.", n_generators, K);

  const float identity[9] = {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
  float generic[9], seeded[9], default_frame[9];
  const float angles[3] = {37.f, 23.f, 61.f};
  const float default_angles[3] = {24.f, 36.f, 66.f};
  radiation_packet_build_octant_frame(angles, 0, generic);
  radiation_packet_build_octant_frame(angles, 5, seeded);
  radiation_packet_build_octant_frame(default_angles, 0, default_frame);

  /* The gates must fail on a planted wrong slot function. */
  const int planted = check_partition("generic", generic, planted_slot, 11, 0) +
                      check_convexity("generic", generic, planted_slot, 0);
  if (planted == 0)
    error("The partition gates passed a planted wrong slot function.");
  message("A planted wrong slot function fails %d gates: ok.", planted);

  const struct {
    const char *name;
    const float *F;
  } frames[] = {{"identity", identity},
                {"generic", generic},
                {"seeded", seeded},
                {"default", default_frame}};
  int failed = 0;
  for (size_t f = 0; f < sizeof(frames) / sizeof(frames[0]); f++) {
    failed += check_partition(frames[f].name, frames[f].F, production_slot,
                              1000 + f, 1);
    failed += check_convexity(frames[f].name, frames[f].F, production_slot, 1);
  }
  if (failed) error("%d partition or convexity gates failed.", failed);

  check_purity(generic);
  check_purity(identity);
  check_ties();
  check_hand_table();
  check_exact_cell_area();

  message(
      "Default frame (24, 36, 66): box axes at least %.2f degrees from a "
      "slot border.",
      box_axis_margin_deg(default_frame));

  message("testRadiationPacketSlotGeometry, K = %d: all gates passed.", K);
  return 0;
}

#else

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  return 0;
}

#endif /* GEAR_SUBGRID_RADIATION_ISRF */
