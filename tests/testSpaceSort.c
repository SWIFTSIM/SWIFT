/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Will Roper (w.roper@sussex.ac.uk)
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local headers. */
#include "swift.h"

/**
 * @file testSpaceSort.c
 * @brief Exercise the displaced-only parallel bucket sort.
 *
 * The sort takes the parallel gather/scatter path when few particles are
 * out of place and falls back to the serial cycle walk when most are, so
 * the cases here are chosen to hit both paths and the D == 0 early exit:
 * nearly-sorted inputs exercise the parallel path (the production regime),
 * fully random and reversed inputs exercise the fallback, and the repeats
 * shake out races in the parallel phases.
 */

/**
 * @brief Fill a linked part/gpart set, each particle remembering its slot.
 */
static void make_parts(struct part *parts, struct xpart *xparts,
                       struct gpart *gparts, const size_t nr) {

  bzero(parts, nr * sizeof(struct part));
  bzero(xparts, nr * sizeof(struct xpart));
  bzero(gparts, nr * sizeof(struct gpart));
  for (size_t i = 0; i < nr; i++) {
    parts[i].id = (long long)i;
    parts[i].gpart = &gparts[i];
    gparts[i].type = swift_type_gas;
    gparts[i].id_or_neg_offset = -(long long)i;
  }
}

/**
 * @brief Check a sorted parts array against the bins it was sorted into.
 *
 * @param parts The sorted particles.
 * @param gparts The (unmoved) gravity counterparts.
 * @param ind The bin of each slot after the sort.
 * @param assign The bin each particle id belonged to.
 * @param counts Expected number of particles per bin.
 * @param nr Number of particles.
 * @param num_bins Number of bins.
 */
static void check_parts(const struct part *parts, const struct gpart *gparts,
                        const int *ind, const int *assign, const int *counts,
                        const size_t nr, const int num_bins) {

  int *seen = (int *)calloc(num_bins, sizeof(int));
  char *found = (char *)calloc(nr, sizeof(char));
  if (seen == NULL || found == NULL) error("Failed to allocate check arrays.");

  int prev_bin = 0;
  for (size_t k = 0; k < nr; k++) {
    const long long id = parts[k].id;
    if (id < 0 || (size_t)id >= nr)
      error("Corrupted id %lld at slot %zu (a lost or duplicated write).", id,
            k);
    if (found[id]) error("Particle %lld appears twice; slot %zu.", id, k);
    found[id] = 1;

    const int bin = assign[id];

    /* Bins must be contiguous and ascending, and ind must agree. */
    if (bin < prev_bin)
      error("Bin order violated at %zu: %d after %d.", k, bin, prev_bin);
    if (ind[k] != bin)
      error("ind[%zu]=%d but the particle there belongs to %d.", k, ind[k],
            bin);
    prev_bin = bin;
    seen[bin]++;

    /* The link must survive the shuffle in both directions. */
    if (parts[k].gpart != &gparts[id])
      error("Wrong gpart link at slot %zu (id %lld).", k, id);
    if (parts[k].gpart->id_or_neg_offset != -(long long)k)
      error("Stale back-link at slot %zu: %lld.", k,
            (long long)parts[k].gpart->id_or_neg_offset);
  }

  for (int b = 0; b < num_bins; b++)
    if (seen[b] != counts[b])
      error("Bin %d holds %d particles, expected %d.", b, seen[b], counts[b]);

  free(seen);
  free(found);
}

/**
 * @brief Run one parts sort for a given bin assignment.
 *
 * @param tp The #threadpool.
 * @param assign The bin of each particle id.
 * @param nr Number of particles.
 * @param num_bins Number of bins.
 * @param out_ids If non-NULL, filled with the resulting id order (for
 *        determinism checks).
 */
static void run_parts(struct threadpool *tp, const int *assign,
                      const size_t nr, const int num_bins,
                      long long *out_ids) {

  struct part *parts = NULL;
  struct xpart *xparts = NULL;
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&parts, part_align,
                     nr * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&xparts, xpart_align,
                     nr * sizeof(struct xpart)) != 0 ||
      posix_memalign((void **)&gparts, gpart_align,
                     nr * sizeof(struct gpart)) != 0)
    error("Failed to allocate particle arrays.");
  make_parts(parts, xparts, gparts, nr);

  int *ind = (int *)malloc(nr * sizeof(int));
  int *counts = (int *)calloc(num_bins, sizeof(int));
  int *counts_copy = (int *)malloc(num_bins * sizeof(int));
  if (ind == NULL || counts == NULL || counts_copy == NULL)
    error("Failed to allocate index arrays.");
  memcpy(ind, assign, nr * sizeof(int));
  for (size_t i = 0; i < nr; i++) counts[assign[i]]++;
  memcpy(counts_copy, counts, num_bins * sizeof(int));

  space_parts_sort(tp, parts, xparts, ind, counts, num_bins,
                   /*parts_offset=*/0, /*verbose=*/0);
  check_parts(parts, gparts, ind, assign, counts_copy, nr, num_bins);

  if (out_ids != NULL)
    for (size_t k = 0; k < nr; k++) out_ids[k] = parts[k].id;

  free(parts);
  free(xparts);
  free(gparts);
  free(ind);
  free(counts);
  free(counts_copy);
}

/**
 * @brief Run one gparts sort, checking the reverse links survive.
 */
static void run_gparts(struct threadpool *tp, const int *assign,
                       const size_t nr, const int num_bins) {

  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  if (posix_memalign((void **)&parts, part_align,
                     nr * sizeof(struct part)) != 0 ||
      posix_memalign((void **)&gparts, gpart_align,
                     nr * sizeof(struct gpart)) != 0)
    error("Failed to allocate particle arrays.");
  bzero(parts, nr * sizeof(struct part));
  bzero(gparts, nr * sizeof(struct gpart));

  /* Half gas, linked both ways; half dark matter, unlinked. */
  for (size_t i = 0; i < nr; i++) {
    if (i % 2 == 0) {
      gparts[i].type = swift_type_gas;
      gparts[i].id_or_neg_offset = -(long long)i;
      parts[i].id = (long long)i;
      parts[i].gpart = &gparts[i];
    } else {
      gparts[i].type = swift_type_dark_matter;
      gparts[i].id_or_neg_offset = (long long)i;
    }
  }

  int *ind = (int *)malloc(nr * sizeof(int));
  int *counts = (int *)calloc(num_bins, sizeof(int));
  if (ind == NULL || counts == NULL) error("Failed to allocate indices.");
  memcpy(ind, assign, nr * sizeof(int));
  for (size_t i = 0; i < nr; i++) counts[assign[i]]++;

  space_gparts_sort(tp, gparts, parts, /*sinks=*/NULL, /*sparts=*/NULL,
                    /*bparts=*/NULL, ind, counts, num_bins, /*verbose=*/0);

  /* Every gas part must point at a gpart that points back at it. */
  for (size_t i = 0; i < nr; i += 2) {
    if (parts[i].gpart == NULL) error("Lost the gpart link at %zu.", i);
    if (parts[i].gpart->type != swift_type_gas)
      error("Link at %zu points at a non-gas gpart.", i);
    if (parts[i].gpart->id_or_neg_offset != -(long long)i)
      error("Broken reverse link at %zu.", i);
  }
  for (size_t k = 1; k < nr; k++)
    if (ind[k] < ind[k - 1]) error("gpart bin order violated at %zu.", k);

  free(parts);
  free(gparts);
  free(ind);
  free(counts);
}

/**
 * @brief Fill a sorted assignment, then reassign a fraction of the
 *        particles to random bins - the production regime.
 */
static void make_nearly_sorted(int *assign, const size_t nr,
                               const int num_bins, const double frac) {

  for (size_t i = 0; i < nr; i++) assign[i] = (int)(i * num_bins / nr);
  const size_t nr_move = (size_t)(frac * nr);
  for (size_t m = 0; m < nr_move; m++)
    assign[rand() % nr] = rand() % num_bins;
}

int main(int argc, char *argv[]) {

  /* More threads than the machine will schedule at once is fine; what
   * matters is exercising the parallel phases at all. */
  struct threadpool tp;
  threadpool_init(&tp, 8);

  const size_t nr = 1 << 16;
  int *assign = (int *)malloc(nr * sizeof(int));
  if (assign == NULL) error("Failed to allocate assignment.");

  /* 1. Nearly sorted, ~1% movers: the parallel path, repeatedly, to shake
   * out races in the scan/gather/scatter phases. */
  for (int rep = 0; rep < 100; rep++) {
    srand(rep + 1);
    make_nearly_sorted(assign, nr, 64, 0.01);
    run_parts(&tp, assign, nr, 64, NULL);
  }
  message("Nearly sorted (1%% movers): OK");

  /* 2. Nearly sorted, ~10% movers: parallel path with substantial D. */
  for (int rep = 0; rep < 100; rep++) {
    srand(rep + 200);
    make_nearly_sorted(assign, nr, 64, 0.10);
    run_parts(&tp, assign, nr, 64, NULL);
  }
  message("Nearly sorted (10%% movers): OK");

  /* 3. Determinism: the same input must give the same output. */
  long long *ids_a = (long long *)malloc(nr * sizeof(long long));
  long long *ids_b = (long long *)malloc(nr * sizeof(long long));
  if (ids_a == NULL || ids_b == NULL) error("Failed to allocate id copies.");
  srand(4242);
  make_nearly_sorted(assign, nr, 64, 0.05);
  run_parts(&tp, assign, nr, 64, ids_a);
  run_parts(&tp, assign, nr, 64, ids_b);
  for (size_t k = 0; k < nr; k++)
    if (ids_a[k] != ids_b[k])
      error("Non-deterministic result at %zu: %lld vs %lld.", k, ids_a[k],
            ids_b[k]);
  free(ids_a);
  free(ids_b);
  message("Determinism: OK");

  /* 4. Already sorted: the D == 0 exit; nothing may move or break. */
  for (size_t i = 0; i < nr; i++) assign[i] = (int)(i * 64 / nr);
  run_parts(&tp, assign, nr, 64, NULL);
  message("Already sorted: OK");

  /* 5. Fully random: nearly everything displaced - the serial fallback. */
  for (int rep = 0; rep < 20; rep++) {
    srand(rep + 500);
    for (size_t i = 0; i < nr; i++) assign[i] = rand() % 64;
    run_parts(&tp, assign, nr, 64, NULL);
  }
  message("Fully random (fallback): OK");

  /* 6. Reversed bins: everything displaced - the fallback again. */
  for (size_t i = 0; i < nr; i++) assign[i] = 63 - (int)(i * 64 / nr);
  run_parts(&tp, assign, nr, 64, NULL);
  message("Reversed bins (fallback): OK");

  /* 7. One particle out of place at each end: D = 2, and the two moves are
   * each other's holes. (For a cycle-following sort this same input is one
   * giant cycle; here it is two copies.) */
  for (size_t i = 0; i < nr; i++) assign[i] = 1;
  assign[nr - 1] = 0;
  run_parts(&tp, assign, nr, 2, NULL);
  message("Two swapped ends: OK");

  /* 8. Empty bins between occupied ones, nearly sorted. */
  for (size_t i = 0; i < nr; i++) assign[i] = ((int)(i * 32 / nr)) * 2;
  for (size_t m = 0; m < nr / 100; m++)
    assign[rand() % nr] = (rand() % 32) * 2;
  run_parts(&tp, assign, nr, 64, NULL);
  message("Empty bins interleaved: OK");

  /* 9. The gravity sort on both paths. */
  for (int rep = 0; rep < 50; rep++) {
    srand(rep + 900);
    make_nearly_sorted(assign, nr, 64, 0.02);
    run_gparts(&tp, assign, nr, 64);
  }
  srand(1000);
  for (size_t i = 0; i < nr; i++) assign[i] = rand() % 64;
  run_gparts(&tp, assign, nr, 64);
  message("Gravity sort links (both paths): OK");

  free(assign);
  threadpool_clean(&tp);

  message("All space sort tests passed.");
  return 0;
}
