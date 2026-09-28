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
#include <config.h>

/* This whole file only makes sense against the single-rank HDF5 reader
 * (read_array_single, in src/single_io.c). The MPI readers (serial_io.c,
 * parallel_io.c) carry the identical fix but are not exercised here: the
 * canonical GEAR build and a bare ./configure are both non-MPI, so under
 * WITH_MPI this file compiles to a trivial no-op, matching the precedent
 * of testFeedback under a non-EAGLE feedback model. */
#if defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Some standard headers. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* glibc's own uninitialised-heap poisoning: the programmatic form of the
 * MALLOC_PERTURB_ environment variable. Using mallopt() instead of the env
 * var makes this test self-contained: it does not depend on how the test
 * runner invokes the binary. */
#ifdef __GLIBC__
#include <malloc.h>
#endif

/* Local headers. */
#include "io_properties.h"
#include "swift.h"

/* read_array_single() is deliberately not declared in single_io.h (only
 * read_ic_single()/write_output_single() are): it is an internal helper of
 * single_io.c that nonetheless has external linkage. This prototype must
 * match its definition there exactly. */
void read_array_single(hid_t h_grp, const struct io_props props, size_t N,
                       const struct unit_system *internal_units,
                       const struct unit_system *ic_units, int cleanup_h,
                       int cleanup_sqrt_a, double h, double a);

/* A stand-in particle array. `energy` mimics an OPTIONAL DOUBLE input field
 * (the ISRF PESpecificEnergy/LWSpecificEnergy/LWPhotonSpecificEnergy
 * fields that first exposed the bug this test guards). `stellar_type`
 * mimics the one OPTIONAL INT input field in the tree,
 * StellarParticleType, which reaches this exact default-fill branch on
 * every IC that omits it. */
struct fake_particle {
  double energy;
  int stellar_type;
};

/* Number of particles to fill. Large enough that an index-dependent
 * regression (e.g. only the first particle written correctly) would show
 * up, though the corruption mechanism itself is deterministic once the
 * allocator is poisoned: a single malloc() call supplies the default
 * value, and its bytes are then memcpy'd, unchanged, into every element. */
#define NUM_PARTICLES 1000

/**
 * @brief Read back every element of a field filled by read_array_single()'s
 * default-fill path and fail on the first non-zero or non-finite one.
 *
 * Compares raw bytes, not the floating-point value: under -ffast-math a
 * denormal or NaN-adjacent bit pattern can be flushed to a value that
 * reads as 0.0 in a floating comparison, which would silently pass on
 * exactly the corrupted input this test exists to catch. A direct bit
 * comparison against the all-zero pattern is exact under any codegen and
 * also catches a NaN encoding without ever calling isnan()/isinf() (which
 * do not compile under -Wnan-infinity-disabled).
 */
static void check_double_field_is_exactly_zero(
    const struct fake_particle *parts, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    uint64_t bits;
    memcpy(&bits, &parts[i].energy, sizeof(bits));
    if (bits != 0ULL) {
      error(
          "Particle %zu: OPTIONAL DOUBLE default-fill did not read back as "
          "exactly zero. Raw bits = 0x%016llx (value = %e). This is the "
          "float-width-fill-of-a-double-field regression fixed by "
          "c34605733.",
          i, (unsigned long long)bits, parts[i].energy);
    }
  }
}

static void check_int_field_is_exactly_zero(const struct fake_particle *parts,
                                            size_t n) {
  for (size_t i = 0; i < n; ++i) {
    if (parts[i].stellar_type != 0) {
      error(
          "Particle %zu: OPTIONAL INT default-fill did not read back as "
          "exactly zero (got %d).",
          i, parts[i].stellar_type);
    }
  }
}

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;

#ifdef __GLIBC__
  /* Force every subsequent malloc() in this process to return
   * uninitialised-looking memory (filled with ~byte) and every free() to
   * scrub with `byte`, exactly like MALLOC_PERTURB_=165 in the shell.
   * Without this, a fresh process's heap pages are usually already
   * zeroed by the kernel, so the pre-fix code often reads back 0.0 by
   * luck: that would make this a test that passes on broken code. Doing
   * the poisoning via mallopt() rather than requiring the test runner to
   * export MALLOC_PERTURB_ keeps the test self-contained; the tradeoff
   * is that this determinism is a glibc extension (guarded on __GLIBC__
   * above), so a non-glibc libc runs this test without the guaranteed
   * poisoning and could pass on broken code by chance of zeroed pages. */
  mallopt(M_PERTURB, 165);
#endif

  /* Build a real HDF5 group that OMITS both fields below, using HDF5's
   * in-memory ("core") driver so no file touches disk. */
  const hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  if (fapl < 0) error("Failed to create a file access property list.");
  if (H5Pset_fapl_core(fapl, 1 << 20, /*backing_store=*/0) < 0)
    error("Failed to select the HDF5 in-memory (core) driver.");
  const hid_t h_file = H5Fcreate("testIOOptionalFieldDefaultWidth_scratch.h5",
                                 H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  if (h_file < 0) error("Failed to create the in-memory HDF5 file.");
  const hid_t h_grp =
      H5Gcreate(h_file, "/PartType0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Failed to create the '/PartType0' group.");

  struct fake_particle *parts =
      (struct fake_particle *)malloc(NUM_PARTICLES * sizeof(struct fake_particle));
  if (parts == NULL) error("Failed to allocate the particle array.");

  /* An OPTIONAL DOUBLE field absent from the ICs: this is the exact
   * shape of PESpecificEnergy/LWSpecificEnergy/LWPhotonSpecificEnergy
   * that exposed the bug. Default value is 0 (io_make_input_field with no
   * explicit default). */
  struct io_props energy_props = io_make_input_field(
      "PESpecificEnergy", DOUBLE, /*dim=*/1, OPTIONAL, UNIT_CONV_NO_UNITS,
      parts, energy);

  /* An OPTIONAL INT field absent from the ICs, matching
   * StellarParticleType: this locks in that the fix's `props.type ==
   * DOUBLE` check, not io_is_double_precision() (which error()s on a
   * non-float/non-double type), guards the default-fill site. Reusing
   * io_is_double_precision() here would abort this test. */
  struct io_props stellar_type_props = io_make_input_field(
      "StellarParticleType", INT, /*dim=*/1, OPTIONAL, UNIT_CONV_NO_UNITS,
      parts, stellar_type);

  read_array_single(h_grp, energy_props, NUM_PARTICLES,
                    /*internal_units=*/NULL, /*ic_units=*/NULL,
                    /*cleanup_h=*/0, /*cleanup_sqrt_a=*/0, /*h=*/1.0,
                    /*a=*/1.0);
  check_double_field_is_exactly_zero(parts, NUM_PARTICLES);

  read_array_single(h_grp, stellar_type_props, NUM_PARTICLES,
                    /*internal_units=*/NULL, /*ic_units=*/NULL,
                    /*cleanup_h=*/0, /*cleanup_sqrt_a=*/0, /*h=*/1.0,
                    /*a=*/1.0);
  check_int_field_is_exactly_zero(parts, NUM_PARTICLES);

  free(parts);
  H5Gclose(h_grp);
  H5Fclose(h_file);
  H5Pclose(fapl);

  message("Optional-field default-fill widths (DOUBLE and INT) are correct.");
  return 0;
}

#else

int main(int argc, char *argv[]) {
  (void)argc;
  (void)argv;
  return 0;
}

#endif /* HAVE_HDF5 && !WITH_MPI */
