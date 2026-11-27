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

/**
 * @file testCompression.c
 * @brief Round-trip test of all the lossy compression filters.
 *
 * For every lossy compression scheme, and every field type the scheme can be
 * applied to, we write an array to an HDF5 file using the same set-up as the
 * snapshot writers (lossy filter + shuffle + deflate + fletcher32), read it
 * back, and check that:
 *  - the values agree with the originals to the accuracy promised by the
 *    filter (exactly for the Nbit filters),
 *  - the on-disk type has the expected class (and, for Nbit, the same
 *    signedness and the expected precision),
 *  - the filter name returned by set_hdf5_lossy_compression() is correct.
 */

#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Includes. */
#include "swift.h"

#ifdef HAVE_HDF5

#include <hdf5.h>

#define test_filename "testCompression.hdf5"
#define test_N 10000
#define test_chunk 1024

/* Amplitude of the values used to test the DScale filters. Keep the floats
 * close to unity so that float round-off stays below 10^-6. */
#define dscale_amplitude_float 1.
#define dscale_amplitude_double 100.

/* Simple deterministic 64-bit random number generator (splitmix64). */
static unsigned long long rng_state = 0x123456789abcdefULL;
static unsigned long long rng_next(void) {
  unsigned long long z = (rng_state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

/* Uniform random double in [0, 1). */
static double rng_uniform(void) {
  return (rng_next() >> 11) * (1.0 / 9007199254740992.0);
}

/* Classification of the filters. */
static int is_dscale(const enum lossy_compression_schemes c) {
  return c >= compression_write_d_scale_1 && c <= compression_write_d_scale_6;
}
static int dscale_digits(const enum lossy_compression_schemes c) {
  return 1 + (int)(c - compression_write_d_scale_1);
}
static int is_nbit(const enum lossy_compression_schemes c) {
  return c >= compression_write_Nbit_32 && c <= compression_write_Nbit_56;
}
static int nbit_bits(const enum lossy_compression_schemes c) {
  switch (c) {
    case compression_write_Nbit_32:
      return 32;
    case compression_write_Nbit_36:
      return 36;
    case compression_write_Nbit_40:
      return 40;
    case compression_write_Nbit_44:
      return 44;
    case compression_write_Nbit_48:
      return 48;
    case compression_write_Nbit_56:
      return 56;
    default:
      error("Not an Nbit filter");
      return 0;
  }
}

/* Number of mantissa bits kept by the mantissa-type filters (0 otherwise). */
static int mantissa_bits(const enum lossy_compression_schemes c) {
  switch (c) {
    case compression_write_d_mantissa_9:
    case compression_write_f_mantissa_9:
      return 9;
    case compression_write_d_mantissa_13:
    case compression_write_f_mantissa_13:
      return 13;
    case compression_write_d_mantissa_21:
      return 21;
    case compression_write_half_float:
      return 10;
    case compression_write_bfloat_16:
      return 7;
    default:
      return 0;
  }
}

/* Range of log10(|x|) to test for the mantissa-type filters. These match the
 * range each filter claims to support. */
static void mantissa_log_range(const enum lossy_compression_schemes c,
                               double *lo, double *hi) {
  switch (c) {
    case compression_write_d_mantissa_9:
    case compression_write_d_mantissa_13:
    case compression_write_d_mantissa_21:
      /* Same exponent range as a double */
      *lo = -300.;
      *hi = 300.;
      break;
    case compression_write_half_float:
      /* Stay within the normal half-float range [6.1e-5, 6.5e4] */
      *lo = -4.;
      *hi = 4.;
      break;
    default:
      /* Same exponent range as a float */
      *lo = -35.;
      *hi = 35.;
      break;
  }
}

/**
 * @brief Returns the list of field types a given filter should be tested
 * with. Returns the number of types.
 */
static int types_for_filter(const enum lossy_compression_schemes c,
                            enum IO_DATA_TYPE types[2]) {
  if (is_dscale(c)) {
    types[0] = FLOAT;
    types[1] = DOUBLE;
    return 2;
  }
  if (is_nbit(c)) {
    types[0] = LONGLONG;
    types[1] = ULONGLONG;
    return 2;
  }
  switch (c) {
    case compression_write_d_mantissa_9:
    case compression_write_d_mantissa_13:
    case compression_write_d_mantissa_21:
      types[0] = DOUBLE;
      return 1;
    case compression_write_f_mantissa_9:
    case compression_write_f_mantissa_13:
    case compression_write_half_float:
    case compression_write_bfloat_16:
      types[0] = FLOAT;
      return 1;
    default:
      return 0;
  }
}

static const char *type_name(const enum IO_DATA_TYPE t) {
  switch (t) {
    case FLOAT:
      return "FLOAT";
    case DOUBLE:
      return "DOUBLE";
    case LONGLONG:
      return "LONGLONG";
    case ULONGLONG:
      return "ULONGLONG";
    default:
      return "UNKNOWN";
  }
}

/**
 * @brief Fills the input array with test values appropriate for the given
 * filter and type. All values lie in the range the filter should support.
 */
static void make_data(const enum lossy_compression_schemes c,
                      const enum IO_DATA_TYPE t, void *data) {

  if (is_nbit(c)) {
    const int n = nbit_bits(c);
    if (t == ULONGLONG) {
      unsigned long long *d = (unsigned long long *)data;
      const unsigned long long max = (1ULL << n) - 1ULL;
      for (int i = 0; i < test_N; ++i) d[i] = rng_next() & max;
      /* Edge cases */
      d[0] = 0ULL;
      d[1] = 1ULL;
      d[2] = max;
      d[3] = max - 1ULL;
      /* Something typical of particle IDs: monotonic run */
      for (int i = 4; i < 100; ++i) d[i] = 574435569ULL + (unsigned)i;
    } else {
      long long *d = (long long *)data;
      const long long max = (long long)((1ULL << (n - 1)) - 1ULL);
      const long long min = -max - 1LL;
      const unsigned long long mask = (1ULL << n) - 1ULL;
      for (int i = 0; i < test_N; ++i)
        d[i] = (long long)(rng_next() & mask) + min;
      /* Edge cases */
      d[0] = 0LL;
      d[1] = 1LL;
      d[2] = -1LL;
      d[3] = max;
      d[4] = min;
      for (int i = 5; i < 100; ++i) d[i] = 574435569LL + i;
    }
    return;
  }

  if (is_dscale(c)) {
    /* Values in [-A, A] */
    const double A =
        (t == FLOAT) ? dscale_amplitude_float : dscale_amplitude_double;
    for (int i = 0; i < test_N; ++i) {
      const double x = A * (2. * rng_uniform() - 1.);
      if (t == FLOAT)
        ((float *)data)[i] = (float)x;
      else
        ((double *)data)[i] = x;
    }
    if (t == FLOAT)
      ((float *)data)[0] = 0.f;
    else
      ((double *)data)[0] = 0.;
    return;
  }

  /* Mantissa-type filters: random sign and log-uniform magnitude */
  double lo, hi;
  mantissa_log_range(c, &lo, &hi);
  for (int i = 0; i < test_N; ++i) {
    const double mag = pow(10., lo + (hi - lo) * rng_uniform());
    const double x = (rng_next() & 1ULL) ? mag : -mag;
    if (t == FLOAT)
      ((float *)data)[i] = (float)x;
    else
      ((double *)data)[i] = x;
  }
  if (t == FLOAT) {
    ((float *)data)[0] = 0.f;
    ((float *)data)[1] = 1.f;
    ((float *)data)[2] = -1.f;
  } else {
    ((double *)data)[0] = 0.;
    ((double *)data)[1] = 1.;
    ((double *)data)[2] = -1.;
  }
}

/**
 * @brief Compares the data read back with the original. Returns the number of
 * elements outside the tolerance (and prints the first few).
 */
static int compare_data(const enum lossy_compression_schemes c,
                        const enum IO_DATA_TYPE t, const void *in,
                        const void *out, const char *name) {

  int n_bad = 0;

  for (int i = 0; i < test_N; ++i) {

    int ok;
    double x_in = 0., x_out = 0., tol = 0.;

    if (t == LONGLONG) {
      ok = ((const long long *)in)[i] == ((const long long *)out)[i];
      if (!ok && n_bad < 5)
        message("%s: element %d: wrote %lld, read %lld", name, i,
                ((const long long *)in)[i], ((const long long *)out)[i]);
    } else if (t == ULONGLONG) {
      ok = ((const unsigned long long *)in)[i] ==
           ((const unsigned long long *)out)[i];
      if (!ok && n_bad < 5)
        message("%s: element %d: wrote %llu, read %llu", name, i,
                ((const unsigned long long *)in)[i],
                ((const unsigned long long *)out)[i]);
    } else {

      const double eps = (t == FLOAT) ? FLT_EPSILON : DBL_EPSILON;
      if (t == FLOAT) {
        x_in = ((const float *)in)[i];
        x_out = ((const float *)out)[i];
      } else {
        x_in = ((const double *)in)[i];
        x_out = ((const double *)out)[i];
      }

      if (is_dscale(c)) {
        /* Absolute accuracy of one unit in the last kept decimal digit
         * (HDF5's D-scale filter does not always round to nearest: e.g.
         * values with |x| < 10^-D are flushed to zero), plus a few ulps (at
         * the scale of the data range, since the filter subtracts the
         * minimum) of round-off in the native type. */
        const double A =
            (t == FLOAT) ? dscale_amplitude_float : dscale_amplitude_double;
        tol = pow(10., -dscale_digits(c)) + 4. * eps * A;
      } else {
        /* Relative accuracy set by the number of mantissa bits kept */
        tol = ldexp(1., -mantissa_bits(c)) * fabs(x_in);
      }
      ok = fabs(x_out - x_in) <= tol;
      if (!ok && n_bad < 5)
        message("%s: element %d: wrote %.17g, read %.17g (|diff|=%g, tol=%g)",
                name, i, x_in, x_out, fabs(x_out - x_in), tol);
    }

    if (!ok) n_bad++;
  }

  return n_bad;
}

/**
 * @brief Checks the type of the dataset on disk. Returns 0 on success.
 */
static int check_disk_type(const enum lossy_compression_schemes c,
                           const enum IO_DATA_TYPE t, hid_t h_data,
                           const char *name) {

  int n_bad = 0;
  const hid_t h_type = H5Dget_type(h_data);
  if (h_type < 0) error("Error getting type of dataset '%s'", name);
  const H5T_class_t cls = H5Tget_class(h_type);

  if (t == LONGLONG || t == ULONGLONG) {

    if (cls != H5T_INTEGER) {
      message("%s: on-disk type is not an integer", name);
      n_bad++;
    } else {
      /* The signedness of the original type must be preserved */
      const H5T_sign_t expected = (t == LONGLONG) ? H5T_SGN_2 : H5T_SGN_NONE;
      const H5T_sign_t sign = H5Tget_sign(h_type);
      if (sign != expected) {
        message("%s: on-disk type is %s but field is %s", name,
                sign == H5T_SGN_NONE ? "unsigned" : "signed",
                t == LONGLONG ? "signed" : "unsigned");
        n_bad++;
      }
      const size_t precision = H5Tget_precision(h_type);
      if (precision != (size_t)nbit_bits(c)) {
        message("%s: on-disk precision is %zu, expected %d", name, precision,
                nbit_bits(c));
        n_bad++;
      }
      if (H5Tget_size(h_type) != 8) {
        message("%s: on-disk size is %zu, expected 8", name,
                H5Tget_size(h_type));
        n_bad++;
      }
    }

  } else {
    if (cls != H5T_FLOAT) {
      message("%s: on-disk type is not a float", name);
      n_bad++;
    }
  }

  H5Tclose(h_type);
  return n_bad;
}

/**
 * @brief Writes, reads back and checks one (filter, type) combination.
 * Returns the number of failures.
 */
static int test_one(hid_t h_file, const enum lossy_compression_schemes c,
                    const enum IO_DATA_TYPE t) {

  char name[64];
  snprintf(name, sizeof(name), "%s_%s", lossy_compression_schemes_names[c],
           type_name(t));

  /* This will call error() if the combination is not allowed */
  io_check_field_compression(t, c, name);

  const size_t type_size = io_sizeof_type(t);
  void *in = malloc(test_N * type_size);
  void *out = malloc(test_N * type_size);
  if (in == NULL || out == NULL) error("Failed to allocate memory");
  memset(out, 0, test_N * type_size);
  make_data(c, t, in);

  /* Same set-up as in the snapshot writers */
  const hsize_t shape[1] = {test_N};
  const hsize_t chunk_shape[1] = {test_chunk};
  const hid_t h_space = H5Screate_simple(1, shape, NULL);
  if (h_space < 0) error("Error creating dataspace for '%s'", name);

  hid_t h_type = H5Tcopy(io_hdf5_type(t));
  hid_t h_prop = H5Pcreate(H5P_DATASET_CREATE);
  if (H5Pset_chunk(h_prop, 1, chunk_shape) < 0)
    error("Error setting chunk size for '%s'", name);

  char filter_name[32] = "None";
  set_hdf5_lossy_compression(&h_prop, &h_type, c, name, filter_name);

  if (H5Pset_shuffle(h_prop) < 0) error("Error setting shuffle for '%s'", name);
  if (H5Pset_deflate(h_prop, 4) < 0)
    error("Error setting deflate for '%s'", name);
  if (H5Pset_fletcher32(h_prop) < 0)
    error("Error setting fletcher32 for '%s'", name);

  hid_t h_data = H5Dcreate(h_file, name, h_type, h_space, H5P_DEFAULT, h_prop,
                           H5P_DEFAULT);
  if (h_data < 0) {
    /* Report and move on so that all the filters get tested */
    message("%s: error creating dataset (see HDF5 error stack above)", name);
    message("%-24s FAILED", name);
    H5Pclose(h_prop);
    H5Tclose(h_type);
    H5Sclose(h_space);
    free(in);
    free(out);
    return 1;
  }
  if (H5Dwrite(h_data, io_hdf5_type(t), H5S_ALL, H5S_ALL, H5P_DEFAULT, in) < 0)
    error("Error writing dataset '%s'", name);
  H5Dclose(h_data);
  H5Pclose(h_prop);
  H5Tclose(h_type);
  H5Sclose(h_space);

  /* Flush so that we really go through the filter pipeline on reading */
  H5Fflush(h_file, H5F_SCOPE_GLOBAL);

  /* Read back */
  h_data = H5Dopen(h_file, name, H5P_DEFAULT);
  if (h_data < 0) error("Error opening dataset '%s'", name);
  if (H5Dread(h_data, io_hdf5_type(t), H5S_ALL, H5S_ALL, H5P_DEFAULT, out) < 0)
    error("Error reading dataset '%s'", name);

  int n_fail = 0;

  /* Check the name of the filter reported */
  if (strcmp(filter_name, lossy_compression_schemes_names[c]) != 0) {
    message("%s: filter name is '%s', expected '%s'", name, filter_name,
            lossy_compression_schemes_names[c]);
    n_fail++;
  }

  n_fail += check_disk_type(c, t, h_data, name);
  H5Dclose(h_data);

  const int n_bad = compare_data(c, t, in, out, name);
  if (n_bad > 0) {
    message("%s: %d/%d elements outside tolerance", name, n_bad, test_N);
    n_fail++;
  }

  message("%-24s %s", name, n_fail ? "FAILED" : "passed");

  free(in);
  free(out);
  return n_fail;
}

int main(int argc, char *argv[]) {

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* Use the same file format version bounds as the snapshot writers. This
   * matters: with the default (earliest) bounds the object headers carry no
   * checksum, and HDF5 >= 1.14.4 then refuses to create datasets whose type
   * has an "unusual" number of unused bits (e.g. DMantissa9/13). SWIFT's
   * files use checksummed (v2) object headers, which skip that check. */
  const hid_t h_props = H5Pcreate(H5P_FILE_ACCESS);
  if (H5Pset_libver_bounds(h_props, HDF5_LOWEST_FILE_FORMAT_VERSION,
                           HDF5_HIGHEST_FILE_FORMAT_VERSION) < 0)
    error("Error setting the file format version bounds");

  const hid_t h_file =
      H5Fcreate(test_filename, H5F_ACC_TRUNC, H5P_DEFAULT, h_props);
  if (h_file < 0) error("Error creating file '%s'", test_filename);
  H5Pclose(h_props);

  int n_fail = 0;
  int n_tests = 0;

  /* Loop over all the lossy compression filters */
  for (int i = compression_write_lossless + 1; i < compression_level_count;
       ++i) {

    const enum lossy_compression_schemes c = (enum lossy_compression_schemes)i;

    enum IO_DATA_TYPE types[2];
    const int n_types = types_for_filter(c, types);

    /* Make sure any filter added in the future gets a test */
    if (n_types == 0) {
      message("No test defined for filter '%s'",
              lossy_compression_schemes_names[c]);
      n_fail++;
      continue;
    }

    for (int k = 0; k < n_types; ++k) {
      n_fail += test_one(h_file, c, types[k]);
      n_tests++;
    }
  }

  H5Fclose(h_file);

  if (n_fail > 0) {
    message("%d failure(s) in %d tests. Keeping '%s' for inspection.", n_fail,
            n_tests, test_filename);
    return 1;
  }

  remove(test_filename);
  message("All %d compression tests passed.", n_tests);
  return 0;
}

#else

int main(int argc, char *argv[]) {
  /* Nothing to test without HDF5 */
  return 0;
}

#endif /* HAVE_HDF5 */
