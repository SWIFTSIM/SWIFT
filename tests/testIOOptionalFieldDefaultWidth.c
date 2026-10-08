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

/* Exercises read_array_single() (src/single_io.c) only; under WITH_MPI
 * this file compiles to a trivial no-op instead. */
#if defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Some standard headers. */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* mallopt(M_PERTURB) is the programmatic form of MALLOC_PERTURB_. It only
 * governs glibc's own malloc(): if a replacement allocator was selected at
 * configure time, skip rather than trust an unpoisoned heap for a PASS. */
#if defined(__GLIBC__) && !defined(HAVE_TCMALLOC) && \
    !defined(HAVE_JEMALLOC) && !defined(HAVE_TBBMALLOC)
#include <malloc.h>
#define IO_TEST_CAN_POISON_HEAP 1
#else
#define IO_TEST_CAN_POISON_HEAP 0
#endif

/* Local headers. */
#include "io_properties.h"
#include "swift.h"

/* Not declared in single_io.h (only read_ic_single()/write_output_single()
 * are); must match its definition in single_io.c exactly. */
void read_array_single(hid_t h_grp, const struct io_props props, size_t N,
                       const struct unit_system *internal_units,
                       const struct unit_system *ic_units, int cleanup_h,
                       int cleanup_sqrt_a, double h, double a);

/* `energy` stands in for an OPTIONAL DOUBLE input field; `stellar_type`
 * for an OPTIONAL INT one. */
struct fake_particle {
  double energy;
  int stellar_type;
};

#define NUM_PARTICLES 1000

/**
 * @brief Fail on the first element whose raw bits are not exactly zero.
 *
 * Compares bits, not the floating-point value: under -ffast-math a
 * corrupted-but-subnormal or NaN-adjacent pattern can be flushed to
 * something that reads as 0.0 in a floating comparison. A raw bit
 * comparison is exact under any codegen and also catches a NaN encoding,
 * without calling isnan()/isinf() (banned under -Wnan-infinity-disabled).
 */
static void check_double_field_is_exactly_zero(
    const struct fake_particle *parts, size_t n) {
  for (size_t i = 0; i < n; ++i) {
    uint64_t bits;
    memcpy(&bits, &parts[i].energy, sizeof(bits));
    if (bits != 0ULL) {
      error(
          "Particle %zu: OPTIONAL DOUBLE default-fill did not read back as "
          "exactly zero. Raw bits = 0x%016llx (value = %e).",
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

#if IO_TEST_CAN_POISON_HEAP
  /* Without this, a fresh process's heap pages are usually already
   * zeroed by the kernel, and a broken default-fill would read back 0.0
   * by luck instead of by correctness. */
  mallopt(M_PERTURB, 165);
#else
  message(
      "No verified heap-poisoning allocator available; skipping (a PASS "
      "here would not be a reliable signal).");
  return 77;
#endif

  /* HDF5's in-memory ("core") driver, so no file touches disk. */
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

  struct fake_particle *parts = (struct fake_particle *)malloc(
      NUM_PARTICLES * sizeof(struct fake_particle));
  if (parts == NULL) error("Failed to allocate the particle array.");

  /* Default value is 0 (io_make_input_field with no explicit default). */
  struct io_props energy_props =
      io_make_input_field("PESpecificEnergy", DOUBLE, /*dim=*/1, OPTIONAL,
                          UNIT_CONV_NO_UNITS, parts, energy);

  /* io_is_double_precision() error()s on a non-float/non-double type, so
   * this locks in that the default-fill site does not call it. */
  struct io_props stellar_type_props =
      io_make_input_field("StellarParticleType", INT, /*dim=*/1, OPTIONAL,
                          UNIT_CONV_NO_UNITS, parts, stellar_type);

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
