/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

/* Local headers. */
#include "custom_mpi.h"
#include "swift.h"

#ifdef WITH_MPI

/* Byte written into every receive buffer before a call, and into the guard
 * region after it. Anything still carrying it after the call was not written
 * to; anything in the guard region *not* carrying it was overrun. */
#define POISON 0xAAu
#define GUARD_BYTES 4096

/* Number of back-to-back calls per case. Leaked requests or datatypes, or a
 * missed message, show up as a deadlock or an error here. */
#define REPEATS 50

/**
 * @brief Expected value of byte k of the payload sent by rank r.
 *
 * Defined on raw bytes so that the same checks work for any datatype.
 */
static unsigned char pattern(const int rank, const size_t k) {
  return (unsigned char)(131u * (unsigned)rank + 7u * k + 1u);
}

/**
 * @brief Gathers the per-rank counts and computes the displacements.
 *
 * @return The total number of elements across all ranks.
 */
static size_t gather_counts(const size_t mine, size_t *counts, size_t *displs,
                            MPI_Comm comm) {
  int size;
  MPI_Comm_size(comm, &size);
  MPI_Allgather(&mine, 1, MPI_UNSIGNED_LONG, counts, 1, MPI_UNSIGNED_LONG,
                comm);
  displs[0] = 0;
  for (int i = 1; i < size; ++i) displs[i] = displs[i - 1] + counts[i - 1];
  return displs[size - 1] + counts[size - 1];
}

/**
 * @brief Runs swift_mpi_allgatherv_sizet() and checks the result byte-for-byte
 * against the native MPI_Allgatherv(), against the known pattern, and against
 * a guard region past the end of the receive buffer.
 *
 * @param name Label used in the messages.
 * @param comm The communicator to run on.
 * @param mine Number of elements contributed by this rank.
 * @param type Datatype of the elements (same for send and receive).
 */
static void check_against_reference(const char *name, MPI_Comm comm,
                                    const size_t mine, MPI_Datatype type) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_Aint lb, extent;
  MPI_Type_get_extent(type, &lb, &extent);
  const size_t elem = (size_t)extent;

  size_t *counts = (size_t *)malloc(size * sizeof(size_t));
  size_t *displs = (size_t *)malloc(size * sizeof(size_t));
  int *icounts = (int *)malloc(size * sizeof(int));
  int *idispls = (int *)malloc(size * sizeof(int));
  if (counts == NULL || displs == NULL || icounts == NULL || idispls == NULL)
    error("[%s] Failed to allocate the count arrays", name);

  const size_t total = gather_counts(mine, counts, displs, comm);
  const size_t bytes = total * elem;

  /* The native reference takes int arguments: this check is only ever run on
   * cases that fit. The > INT_MAX cases are covered separately below. */
  for (int i = 0; i < size; ++i) {
    if (counts[i] > INT_MAX || displs[i] > INT_MAX)
      error("[%s] Reference case does not fit in an int", name);
    icounts[i] = (int)counts[i];
    idispls[i] = (int)displs[i];
  }

  /* Payload: a rank-dependent byte pattern. +1 so a rank with no data still
   * gets a valid (non-NULL) pointer. */
  unsigned char *send = (unsigned char *)malloc(mine * elem + 1);
  unsigned char *ref = (unsigned char *)malloc(bytes + GUARD_BYTES);
  unsigned char *got = (unsigned char *)malloc(bytes + GUARD_BYTES);
  if (send == NULL || ref == NULL || got == NULL)
    error("[%s] Failed to allocate the buffers", name);
  for (size_t k = 0; k < mine * elem; ++k) send[k] = pattern(rank, k);

  /* Reference result from the native collective. */
  memset(ref, POISON, bytes + GUARD_BYTES);
  MPI_Allgatherv(send, (int)mine, type, ref, icounts, idispls, type, comm);

  for (int rep = 0; rep < REPEATS; ++rep) {

    memset(got, POISON, bytes + GUARD_BYTES);
    const int err = swift_mpi_allgatherv_sizet(send, mine, type, got, counts,
                                               displs, type, comm);
    if (err != MPI_SUCCESS)
      error("[%s] Call returned %d on repeat %d", name, err, rep);

    /* Same as the native collective? */
    for (size_t k = 0; k < bytes; ++k) {
      if (got[k] != ref[k])
        error(
            "[%s] Byte %zu of %zu differs from MPI_Allgatherv (got 0x%02x, "
            "want 0x%02x) on repeat %d",
            name, k, bytes, got[k], ref[k], rep);
    }

    /* Nothing written past the end? */
    for (size_t k = 0; k < GUARD_BYTES; ++k) {
      if (got[bytes + k] != POISON)
        error("[%s] Wrote %zu bytes past the end of the receive buffer", name,
              k + 1);
    }
  }

  /* And, independently of the reference, is it the pattern we expect? */
  for (int r = 0; r < size; ++r) {
    for (size_t k = 0; k < counts[r] * elem; ++k) {
      const size_t off = displs[r] * elem + k;
      if (got[off] != pattern(r, k))
        error(
            "[%s] Byte %zu of rank %d's data (offset %zu) is 0x%02x, want "
            "0x%02x",
            name, k, r, off, got[off], pattern(r, k));
    }
  }

  if (rank == 0)
    message("%-22s %d ranks, %8zu elements, %10zu bytes: OK", name, size, total,
            bytes);

  free(send);
  free(ref);
  free(got);
  free(counts);
  free(displs);
  free(icounts);
  free(idispls);
}

/**
 * @brief Checks the layout built by create_large_count_type() without moving
 * any data: its size and extent must be exactly count * (size, extent) of the
 * base type, including for counts that do not fit in an int.
 */
static void check_large_count_type(const size_t count, MPI_Datatype oldtype,
                                   const char *tname) {

  int oldsize;
  MPI_Aint lb, oldextent;
  MPI_Type_size(oldtype, &oldsize);
  MPI_Type_get_extent(oldtype, &lb, &oldextent);

  MPI_Datatype t;
  if (create_large_count_type(count, oldtype, &t) != MPI_SUCCESS)
    error("[%s] create_large_count_type(%zu) failed", tname, count);

  /* A zero count must yield the null handle and nothing to free. */
  if (count == 0) {
    if (t != MPI_DATATYPE_NULL)
      error("[%s] count 0 did not return MPI_DATATYPE_NULL", tname);
    return;
  }

  MPI_Type_commit(&t);

  MPI_Count size_x, lb_x, extent_x;
  MPI_Type_size_x(t, &size_x);
  MPI_Type_get_extent_x(t, &lb_x, &extent_x);

  const MPI_Count want_size = (MPI_Count)count * (MPI_Count)oldsize;
  const MPI_Count want_extent = (MPI_Count)count * (MPI_Count)oldextent;

  if (size_x != want_size)
    error("[%s] count %zu: type size is %lld, want %lld", tname, count,
          (long long)size_x, (long long)want_size);
  if (lb_x != 0 || extent_x != want_extent)
    error("[%s] count %zu: lb/extent are %lld/%lld, want 0/%lld", tname, count,
          (long long)lb_x, (long long)extent_x, (long long)want_extent);

  MPI_Type_free(&t);
}

/**
 * @brief Exercises byte displacements beyond INT_MAX without needing GBs of
 * RAM -- the case the size_t interface exists for.
 *
 * The receive type is a single byte resized to an extent of 1 MiB, so a
 * displacement of d *elements* is a byte offset of d MiB. With >= 2 ranks and
 * 2100 elements per rank, everything from rank 1 onwards lands beyond 2 GiB.
 * The receive buffer is a MAP_NORESERVE mapping: only the pages actually
 * written to get physical memory, i.e. one page per element.
 */
static void check_huge_displacements(MPI_Comm comm) {
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (size < 2) {
    if (rank == 0)
      message("%-22s skipped (needs >= 2 ranks)", "huge displacements");
    return;
  }

  const size_t stride = (size_t)1 << 20; /* 1 MiB extent per element */
  const size_t mine = 2100;

  size_t *counts = (size_t *)malloc(size * sizeof(size_t));
  size_t *displs = (size_t *)malloc(size * sizeof(size_t));
  if (counts == NULL || displs == NULL)
    error("Failed to allocate the count arrays");
  const size_t total = gather_counts(mine, counts, displs, comm);
  const size_t bytes = total * stride;

  /* Make sure the test actually tests what it claims to. */
  if (displs[1] * stride <= (size_t)INT_MAX)
    error(
        "Test setup: rank 1's displacement (%zu bytes) does not exceed "
        "INT_MAX",
        displs[1] * stride);

  void *map = mmap(NULL, bytes + GUARD_BYTES, PROT_READ | PROT_WRITE,
                   MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  if (map == MAP_FAILED) {
    message("%-22s skipped: mmap of %zu GiB of address space failed (%s)",
            "huge displacements", (bytes + GUARD_BYTES) >> 30, strerror(errno));
    free(counts);
    free(displs);
    return;
  }
  unsigned char *recv = (unsigned char *)map;

  MPI_Datatype byte_1mib;
  MPI_Type_create_resized(MPI_BYTE, 0, (MPI_Aint)stride, &byte_1mib);
  MPI_Type_commit(&byte_1mib);

  unsigned char *send = (unsigned char *)malloc(mine);
  if (send == NULL) error("Failed to allocate the send buffer");
  for (size_t k = 0; k < mine; ++k) send[k] = pattern(rank, k);

  /* Sent as plain bytes, received into the strided layout: identical type
   * signatures, so this is a legal (and deliberately asymmetric) call. */
  const int err = swift_mpi_allgatherv_sizet(send, mine, MPI_BYTE, recv, counts,
                                             displs, byte_1mib, comm);
  if (err != MPI_SUCCESS) error("Call returned %d", err);

  for (int r = 0; r < size; ++r) {
    for (size_t k = 0; k < counts[r]; ++k) {
      const size_t off = (displs[r] + k) * stride;
      if (recv[off] != pattern(r, k))
        error(
            "Element %zu of rank %d at byte offset %zu (%.2f GiB) is "
            "0x%02x, want 0x%02x",
            k, r, off, (double)off / (double)((size_t)1 << 30), recv[off],
            pattern(r, k));

      /* The byte next to it was never written: still a zero page. */
      if (recv[off + 1] != 0)
        error("Byte %zu next to element %zu of rank %d was written to", off + 1,
              k, r);
    }
  }

  if (rank == 0)
    message("%-22s %d ranks, largest byte displacement %zu (INT_MAX is %d): OK",
            "huge displacements", size, displs[size - 1] * stride, INT_MAX);

  MPI_Type_free(&byte_1mib);
  free(send);
  free(counts);
  free(displs);
  munmap(map, bytes + GUARD_BYTES);
}

/* Per-rank element counts for the reference comparisons. */
static size_t law_uniform(const int rank, const int size) { return 5; }
static size_t law_varying(const int rank, const int size) {
  return (size_t)(37 * rank + 11) % 23;
}
static size_t law_one_empty(const int rank, const int size) {
  return (rank == 1) ? 0 : 3 * (size_t)(rank + 1);
}
static size_t law_all_empty(const int rank, const int size) { return 0; }
static size_t law_last_only(const int rank, const int size) {
  return (rank == size - 1) ? 17 : 0;
}
/* Large enough to leave the eager protocol behind (3.2 MB/rank for 32 B). */
static size_t law_large(const int rank, const int size) { return 100000; }

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  engine_rank = rank;

  /* Initialize CPU frequency, this also starts time. */
  unsigned long long cpufreq = 0;
  clocks_set_cpufreq(cpufreq);

  /* A 32-byte struct sent as raw bytes: exactly how fof_mpi_type is built. */
  MPI_Datatype struct32;
  MPI_Type_contiguous(32, MPI_BYTE, &struct32);
  MPI_Type_commit(&struct32);

  const struct {
    const char *name;
    size_t (*law)(int, int);
  } laws[] = {{"uniform", law_uniform},     {"varying", law_varying},
              {"one_empty", law_one_empty}, {"all_empty", law_all_empty},
              {"last_only", law_last_only}, {"large", law_large}};
  const int n_laws = (int)(sizeof(laws) / sizeof(laws[0]));

  const MPI_Datatype types[] = {struct32, MPI_DOUBLE};
  const char *tnames[] = {"struct32", "double"};
  const int n_types = (int)(sizeof(types) / sizeof(types[0]));

  /* Every count law with every datatype, against the native collective. */
  for (int t = 0; t < n_types; ++t) {
    for (int l = 0; l < n_laws; ++l) {
      char name[64];
      snprintf(name, sizeof(name), "%s/%s", tnames[t], laws[l].name);
      check_against_reference(name, MPI_COMM_WORLD, laws[l].law(rank, size),
                              types[t]);
    }
  }

  /* On other communicators: a duplicate, and one with the ranks reversed, to
   * make sure the ranks of 'comm' are used rather than those of COMM_WORLD. */
  MPI_Comm dup, rev;
  MPI_Comm_dup(MPI_COMM_WORLD, &dup);
  MPI_Comm_split(MPI_COMM_WORLD, 0, size - 1 - rank, &rev);
  int rev_rank;
  MPI_Comm_rank(rev, &rev_rank);
  check_against_reference("dup/varying", dup, law_varying(rank, size),
                          struct32);
  check_against_reference("reversed/varying", rev, law_varying(rev_rank, size),
                          struct32);
  MPI_Comm_free(&dup);
  MPI_Comm_free(&rev);

  /* Layouts for counts that do not fit in an int (structural, no data). */
  if (rank == 0) {
    const size_t counts[] = {0,
                             1,
                             (size_t)INT_MAX,
                             (size_t)INT_MAX + 1,
                             2 * (size_t)INT_MAX,
                             2 * (size_t)INT_MAX + 12345,
                             3 * (size_t)INT_MAX - 1};
    const int n_counts = (int)(sizeof(counts) / sizeof(counts[0]));
    for (int i = 0; i < n_counts; ++i) {
      check_large_count_type(counts[i], MPI_BYTE, "byte");
      check_large_count_type(counts[i], struct32, "struct32");
    }
    message("%-22s layouts for up to %zu elements (%zu x INT_MAX): OK",
            "large count types", counts[n_counts - 1],
            counts[n_counts - 1] / (size_t)INT_MAX);
  }

  /* Byte displacements beyond INT_MAX. */
  check_huge_displacements(MPI_COMM_WORLD);

  MPI_Type_free(&struct32);

  if (rank == 0) message("All tests passed.");

  MPI_Finalize();
  return 0;
}

#else

int main(int argc, char *argv[]) {
  printf("SWIFT was not compiled with MPI support, nothing to test.\n");
  return 0;
}

#endif /* WITH_MPI */
