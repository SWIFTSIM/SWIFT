/*******************************************************************************
 * This file is part of SWIFT.
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

/* Some standard headers. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "runner_radiation_feedback.h"
#include "swift.h"

/*! Number of star IDs in the matrix checks. */
#define N_IDS_MATRIX 100000

/*! Number of directions and of star IDs in the occupancy checks. */
#define N_SAMPLES 1000000

/*! Largest deviation of R R^T from the identity and of det R from 1. */
#define ORTHONORMAL_BAR 1e-12

/*! One-sided normal quantile for p = 0.001. */
#define Z_P001 3.0902323

#if defined(HAVE_CHEALPIX)

/**
 * @brief Next value of a splitmix64 stream, mapped to [0, 1).
 *
 * Independent of the generator under test.
 *
 * @param state The stream state.
 * @return A number in [0, 1).
 */
static double test_uniform(uint64_t *state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15ULL);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  z = z ^ (z >> 31);
  return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

/**
 * @brief Draw a direction uniform on the sphere.
 *
 * @param state The stream state.
 * @param dx (return) The direction, single precision, as the code passes it.
 */
static void test_direction(uint64_t *state, float dx[3]) {
  const double z = 2.0 * test_uniform(state) - 1.0;
  const double phi = 2.0 * M_PI * test_uniform(state);
  const double s = sqrt(1.0 - z * z);
  dx[0] = (float)(s * cos(phi));
  dx[1] = (float)(s * sin(phi));
  dx[2] = (float)z;
}

/**
 * @brief Critical value of chi-square at p = 0.001 (Wilson-Hilferty).
 *
 * @param dof The degrees of freedom.
 * @return The value a chi-square with @a dof degrees of freedom exceeds with
 * probability 0.001.
 */
static double chi2_critical(const int dof) {
  const double a = 2.0 / (9.0 * dof);
  const double t = 1.0 - a + Z_P001 * sqrt(a);
  return dof * t * t * t;
}

/**
 * @brief Chi-square of pixel counts against equal expected counts.
 *
 * @param counts The counts per pixel.
 * @param n_pix The number of pixels.
 * @param total The total count.
 * @return The chi-square statistic.
 */
static double chi2_flat(const long *counts, const int n_pix, const long total) {
  const double expected = (double)total / n_pix;
  double chi2 = 0.0;
  for (int p = 0; p < n_pix; p++) {
    const double d = counts[p] - expected;
    chi2 += d * d / expected;
  }
  return chi2;
}

#endif /* HAVE_CHEALPIX */

/**
 * @brief Print and score one check that must pass or must fail.
 *
 * @param name What is checked.
 * @param value The measured statistic.
 * @param bar The bar.
 * @param expect_fail Is this a control that must exceed the bar?
 * @return 1 if the check did not behave as required, 0 otherwise.
 */
static int report(const char *name, const double value, const double bar,
                  const int expect_fail) {
  const int exceeds = !(value <= bar);
  const int bad = (exceeds != expect_fail);
  printf("%s %s: %g (bar %g, %s)\n", bad ? "FAIL" : "ok  ", name, value, bar,
         expect_fail ? "control, must exceed" : "must not exceed");
  return bad;
}

/**
 * @brief Largest departure of a 3x3 matrix from a proper rotation.
 *
 * @param R The matrix, row-major.
 * @return The larger of max|R R^T - I| and |det R - 1|.
 */
static double departure_from_rotation(const double R[9]) {
  double worst = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int k = 0; k < 3; k++) s += R[3 * i + k] * R[3 * j + k];
      worst = fmax(worst, fabs(s - (i == j ? 1.0 : 0.0)));
    }
  const double det = R[0] * (R[4] * R[8] - R[5] * R[7]) -
                     R[1] * (R[3] * R[8] - R[5] * R[6]) +
                     R[2] * (R[3] * R[7] - R[4] * R[6]);
  return fmax(worst, fabs(det - 1.0));
}

#if defined(HAVE_CHEALPIX)

/**
 * @brief A rotation about z only, by an angle fixed by the ID.
 *
 * Deliberately biased: it is the positive control of the uniformity check.
 *
 * @param id The star ID.
 * @param R (return) The matrix.
 */
static void rotation_about_z_only(const int64_t id, double R[9]) {
  const double angle =
      2.0 * M_PI *
      random_unit_interval(id, 0, random_number_HII_frame_rotation);
  const double c = cos(angle), s = sin(angle);
  const double M[9] = {c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0};
  memcpy(R, M, sizeof(M));
}

/**
 * @brief The identity matrix.
 *
 * @param id Unused.
 * @param R (return) The matrix.
 */
static void rotation_identity(const int64_t id, double R[9]) {
  (void)id;
  const double M[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  memcpy(R, M, sizeof(M));
}

#endif /* HAVE_CHEALPIX */

/**
 * @brief Matrix checks: orthonormality, determinant, determinism.
 *
 * @return The number of failed checks.
 */
static int test_matrices(void) {
  int fail = 0;
  double worst = 0.0;
  double R[9], R2[9];
  int deterministic = 1;
  int distinct = 0;
  for (int64_t id = 1; id <= N_IDS_MATRIX; id++) {
    runner_hii_cone_frame_from_id(id, R);
    worst = fmax(worst, departure_from_rotation(R));
    runner_hii_cone_frame_from_id(id, R2);
    if (memcmp(R, R2, sizeof(R)) != 0) deterministic = 0;
    runner_hii_cone_frame_from_id(id + 1, R2);
    if (memcmp(R, R2, sizeof(R)) != 0) distinct++;
  }
  fail += report("rotation matrix departure from orthonormal, det 1", worst,
                 ORTHONORMAL_BAR, 0);

  /* Control: a matrix scaled by 1.001 is not a rotation. */
  runner_hii_cone_frame_from_id(12345, R);
  for (int k = 0; k < 9; k++) R[k] *= 1.001;
  fail += report("control: scaled matrix departure", departure_from_rotation(R),
                 ORTHONORMAL_BAR, 1);

  fail += report("same ID gives two different matrices (1 = yes)",
                 deterministic ? 0.0 : 1.0, 0.0, 0);
  fail += report("neighbouring IDs give different matrices (missing)",
                 (double)(N_IDS_MATRIX - distinct), 0.0, 0);
  return fail;
}

#if defined(HAVE_CHEALPIX)

/**
 * @brief The pixel of a direction in the box frame, without the helper.
 *
 * @param dx The offset star minus particle.
 * @param nside The HEALPix nside.
 * @return The pixel.
 */
static int pixel_reference(const float dx[3], const int nside) {
  const double dir[3] = {-(double)dx[0], -(double)dx[1], -(double)dx[2]};
  int64_t ipix = 0;
  vec2pix_ring64(nside, dir, &ipix);
  return (int)ipix;
}

/**
 * @brief Occupancy of pixels by uniform directions rotated by one star's R.
 *
 * @param rot The rotation builder.
 * @param id The star ID.
 * @param nside The HEALPix nside.
 * @param shear Squeeze z by 0.5 after the rotation (control of a non-rotation).
 * @return The chi-square of the occupancy against a flat one.
 */
static double occupancy_chi2(void (*rot)(int64_t, double *), const int64_t id,
                             const int nside, const int shear) {
  const int n_pix = 12 * nside * nside;
  long *counts = calloc(n_pix, sizeof(long));
  double R[9];
  rot(id, R);
  if (shear)
    for (int k = 6; k < 9; k++) R[k] *= 0.5;
  uint64_t state = 0x1234567ULL + (uint64_t)id;
  float dx[3];
  for (int i = 0; i < N_SAMPLES; i++) {
    test_direction(&state, dx);
    counts[runner_hii_get_pixel(dx, n_pix, R)]++;
  }
  const double chi2 = chi2_flat(counts, n_pix, N_SAMPLES);
  free(counts);
  return chi2;
}

/**
 * @brief Occupancy of the image of one fixed axis under R(id) over many IDs.
 *
 * A rotation uniform on SO(3) sends a fixed axis to a point uniform on the
 * sphere, so the pixel counts over many IDs are flat.
 *
 * @param rot The rotation builder.
 * @param nside The HEALPix nside.
 * @return The chi-square of the occupancy against a flat one.
 */
static double stacked_chi2(void (*rot)(int64_t, double *), const int nside) {
  const int n_pix = 12 * nside * nside;
  long *counts = calloc(n_pix, sizeof(long));
  /* A generic axis, not along a box axis or a diagonal. The sign convention
     of runner_hii_get_pixel negates dx, so -dx is the axis itself. */
  const double n = sqrt(0.3 * 0.3 + 0.5 * 0.5 + 0.81 * 0.81);
  const float dx[3] = {(float)(-0.3 / n), (float)(-0.5 / n),
                       (float)(-0.81 / n)};
  double R[9];
  for (int64_t id = 1; id <= N_SAMPLES; id++) {
    rot(id, R);
    counts[runner_hii_get_pixel(dx, n_pix, R)]++;
  }
  const double chi2 = chi2_flat(counts, n_pix, N_SAMPLES);
  free(counts);
  return chi2;
}

/**
 * @brief The pixel of a direction after an independent rotation in extended
 * precision.
 *
 * @param dx The offset star minus particle.
 * @param nside The HEALPix nside.
 * @param R The rotation, row-major.
 * @return The pixel.
 */
static int pixel_rotated_reference(const float dx[3], const int nside,
                                   const double R[9]) {
  long double v[3];
  for (int i = 0; i < 3; i++)
    v[i] = -(long double)dx[0] * R[3 * i] - (long double)dx[1] * R[3 * i + 1] -
           (long double)dx[2] * R[3 * i + 2];
  const double dir[3] = {(double)v[0], (double)v[1], (double)v[2]};
  int64_t ipix = 0;
  vec2pix_ring64(nside, dir, &ipix);
  return (int)ipix;
}

/**
 * @brief Pixel checks: uniformity, flag off, identity, one frame for any dx.
 *
 * @return The number of failed checks.
 */
static int test_pixels(void) {
  int fail = 0;
  char name[160];
  const int nsides[3] = {1, 2, 4};

  for (int s = 0; s < 3; s++) {
    const int nside = nsides[s];
    const int n_pix = 12 * nside * nside;
    const double bar = chi2_critical(n_pix - 1);

    snprintf(name, sizeof(name), "nside %d, uniform directions, one star's R",
             nside);
    fail +=
        report(name, occupancy_chi2(runner_hii_cone_frame_from_id, 7, nside, 0),
               bar, 0);
    snprintf(name, sizeof(name),
             "nside %d, control: squeezed (non-rotation) matrix", nside);
    fail +=
        report(name, occupancy_chi2(runner_hii_cone_frame_from_id, 7, nside, 1),
               bar, 1);

    snprintf(name, sizeof(name), "nside %d, stacked image of an axis over IDs",
             nside);
    fail += report(name, stacked_chi2(runner_hii_cone_frame_from_id, nside),
                   bar, 0);
    snprintf(name, sizeof(name),
             "nside %d, control: rotation about z only, stacked", nside);
    fail += report(name, stacked_chi2(rotation_about_z_only, nside), bar, 1);
    snprintf(name, sizeof(name), "nside %d, control: identity, stacked", nside);
    fail += report(name, stacked_chi2(rotation_identity, nside), bar, 1);
  }

  /* Flag off: a NULL frame gives the box-frame pixel, bit for bit. The
     identity matrix, the negative control of the rotation, leaves every pixel
     unchanged. The rotated pixel is also compared with one from an independent
     rotation in long double, for dx and for -dx (so both signs of the offset
     go through the same frame). */
  for (int s = 0; s < 3; s++) {
    const int nside = nsides[s];
    const int n_pix = 12 * nside * nside;
    uint64_t state = 99;
    long diff_null = 0, diff_identity = 0, bad_rotated = 0, bad_range = 0;
    double R[9], I[9];
    rotation_identity(0, I);
    float dx[3], mdx[3];
    for (int i = 0; i < 100000; i++) {
      test_direction(&state, dx);
      runner_hii_cone_frame_from_id(1 + i % 977, R);
      const int ref = pixel_reference(dx, nside);
      if (runner_hii_get_pixel(dx, n_pix, NULL) != ref) diff_null++;
      if (runner_hii_get_pixel(dx, n_pix, I) != ref) diff_identity++;
      const int p = runner_hii_get_pixel(dx, n_pix, R);
      for (int k = 0; k < 3; k++) mdx[k] = -dx[k];
      const int q = runner_hii_get_pixel(mdx, n_pix, R);
      if (p < 0 || p >= n_pix) bad_range++;
      if (p != pixel_rotated_reference(dx, nside, R)) bad_rotated++;
      if (q != pixel_rotated_reference(mdx, nside, R)) bad_rotated++;
    }
    snprintf(name, sizeof(name),
             "nside %d, flag off: pixels differing from box", nside);
    fail += report(name, (double)diff_null, 0.0, 0);
    snprintf(name, sizeof(name),
             "nside %d, identity frame: pixels changed (negative control)",
             nside);
    fail += report(name, (double)diff_identity, 0.0, 0);
    snprintf(name, sizeof(name),
             "nside %d, rotated: pixels differing from an independent rotation",
             nside);
    fail += report(name, (double)bad_rotated, 0.0, 0);
    snprintf(name, sizeof(name), "nside %d, rotated: pixel out of range",
             nside);
    fail += report(name, (double)bad_range, 0.0, 0);
  }

  /* Control: a pixel must change when the frame is on (the check above only
     proves that nothing changes when it is off). */
  {
    uint64_t state = 5;
    long changed = 0;
    double R[9];
    float dx[3];
    runner_hii_cone_frame_from_id(31, R);
    for (int i = 0; i < 10000; i++) {
      test_direction(&state, dx);
      if (runner_hii_get_pixel(dx, 12, R) != pixel_reference(dx, 1)) changed++;
    }
    fail += report(
        "control: share of pixels left unmoved by a random frame (nside 1)",
        1.0 - (double)changed / 10000, 0.5, 0);
  }
  return fail;
}

#endif /* HAVE_CHEALPIX */

/**
 * @brief 64-bit digest of the rotation matrices and pixels of a list of IDs.
 *
 * @param id_shift Added to every ID before the helper. 0 for the real code;
 * non-zero emulates a seed that depends on the rank (control).
 * @return An FNV-1a digest of the bits of every matrix and pixel.
 */
static uint64_t digest_of_ids(const int64_t id_shift) {
  uint64_t h = 1469598103934665603ULL;
  double R[9];
  for (int64_t id = 1; id <= 20000; id++) {
    runner_hii_cone_frame_from_id(id + id_shift, R);
    unsigned char bytes[sizeof(R)];
    memcpy(bytes, R, sizeof(R));
    for (size_t k = 0; k < sizeof(R); k++)
      h = (h ^ bytes[k]) * 1099511628211ULL;
#if defined(HAVE_CHEALPIX)
    const float dx[3] = {0.3f * (float)(id % 7) - 1.f, 0.2f, 0.1f + (float)id};
    h = (h ^ (uint64_t)runner_hii_get_pixel(dx, 48, R)) * 1099511628211ULL;
#endif
  }
  return h;
}

/**
 * @brief Rank independence: the same IDs give bit-identical results in two
 * rank contexts, and across real MPI ranks when run under mpirun.
 *
 * @param rank The MPI rank (0 without MPI).
 * @param n_ranks The number of ranks (1 without MPI).
 * @return The number of failed checks.
 */
static int test_rank_independence(const int rank, const int n_ranks) {
  int fail = 0;
  const uint64_t a = digest_of_ids(0);
  const uint64_t b = digest_of_ids(0);
  fail += report("rank context A against B: digests differ (1 = yes)",
                 a == b ? 0.0 : 1.0, 0.0, 0);
  /* Control: a seed that depends on the rank must give another digest. */
  const uint64_t c = digest_of_ids(1000003LL * (rank + 1));
  fail += report("control: rank-dependent seed gives another digest (1 = yes)",
                 a == c ? 0.0 : 1.0, 0.0, 1);

#ifdef WITH_MPI
  uint64_t lo = a, hi = a;
  MPI_Allreduce(&a, &lo, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&a, &hi, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  fail += report("MPI ranks: digest min differs from max (1 = yes)",
                 lo == hi ? 0.0 : 1.0, 0.0, 0);
  if (rank == 0) printf("      MPI ranks in this run: %d\n", n_ranks);
#else
  (void)n_ranks;
#endif
  return fail;
}

#if defined(HAVE_CHEALPIX)

/**
 * @brief Write the pixels of 1e6 directions to a file, for a comparison of two
 * builds of this test (with and without FMA contraction).
 *
 * Also counts, per nside, the directions whose pixel changes under a one-ULP
 * perturbation of any input component, i.e. those next to a pixel boundary.
 *
 * @param path The output file.
 * @return 0.
 */
static int dump_pixels(const char *path) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) error("Cannot open %s", path);
  const int nsides[3] = {1, 2, 4};
  double R[9];
  for (int s = 0; s < 3; s++) {
    const int n_pix = 12 * nsides[s] * nsides[s];
    uint64_t state = 4242 + s;
    float dx[3];
    for (int i = 0; i < N_SAMPLES; i++) {
      test_direction(&state, dx);
      runner_hii_cone_frame_from_id(1 + i % 4093, R);
      const int p = runner_hii_get_pixel(dx, n_pix, R);
      uint8_t near_boundary = 0;
      for (int k = 0; k < 3 && !near_boundary; k++) {
        float up[3] = {dx[0], dx[1], dx[2]};
        up[k] = nextafterf(dx[k], 1e30f);
        float dn[3] = {dx[0], dx[1], dx[2]};
        dn[k] = nextafterf(dx[k], -1e30f);
        if (runner_hii_get_pixel(up, n_pix, R) != p ||
            runner_hii_get_pixel(dn, n_pix, R) != p)
          near_boundary = 1;
      }
      const uint8_t rec[2] = {(uint8_t)p, near_boundary};
      fwrite(rec, 1, 2, f);
    }
  }
  fclose(f);
  return 0;
}

#endif

int main(int argc, char *argv[]) {

  int rank = 0, n_ranks = 1;
#ifdef WITH_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
#endif

  int fail = 0;

#if defined(HAVE_CHEALPIX)
  if (argc == 3 && strcmp(argv[1], "--dump") == 0) {
    dump_pixels(argv[2]);
    goto done;
  }
#else
  (void)argc;
  (void)argv;
#endif

  fail += test_matrices();
#if defined(HAVE_CHEALPIX)
  fail += test_pixels();
#else
  printf("skip pixel checks: HAVE_CHEALPIX is not defined\n");
#endif
  fail += test_rank_independence(rank, n_ranks);

#if defined(HAVE_CHEALPIX)
done:
#endif
#ifdef WITH_MPI
  MPI_Finalize();
#endif
  if (fail) printf("%d check(s) failed\n", fail);
  return fail ? 1 : 0;
}
