/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2026 Will Roper (w.roper@sussex.ac.uk)
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
#include <string.h>

/* This object's header. */
#include "error.h"
#include "memswap.h"
#include "memuse.h"
#include "space.h"
#include "threadpool.h"

/**
 * @file space_sort.c
 * @brief Bucket sort of the particle arrays, moving only the displaced
 *        particles.
 *
 * Between rebuilds almost every particle is still in the right top-level
 * cell, so the sort's job is to fix the few that are not. A slot whose
 * particle does not belong to the bin covering it is *displaced*, and every
 * displaced slot is two things at once: a hole (its bin is owed a particle)
 * and an orphan (the particle in it belongs elsewhere). Per bin, holes and
 * orphans are equinumerous, so a bijection between the displaced slots
 * exists and can be fixed before anything moves. That removes the need for
 * cycle-following entirely: gather the orphans into a buffer, scatter each
 * to its assigned hole. Every destination is written exactly once, by
 * exactly one thread, so the whole thing runs without a single atomic.
 *
 * The phases, separated by threadpool barriers:
 *
 *   1. Scan (parallel, over all N slots): find the displaced slots and
 *      compact them, in ascending slot order, into one list. Fixed-size
 *      chunks indexed by chunk id keep the result identical for any thread
 *      count and any scheduling. This is the only pass that touches all N.
 *   2. Match (serial, over the D displaced): ascending slot order means the
 *      list is already grouped by home bin, so the holes of each bin are a
 *      contiguous run. Give the i-th orphan targeting bin b the i-th hole
 *      of bin b.
 *   3. Gather (parallel, over D): copy the orphans into a buffer.
 *      Destinations are themselves displaced slots, so nothing may be
 *      scattered until everything is gathered - the barrier provides this.
 *   4. Scatter (parallel, over D): write each orphan to its hole, update
 *      the index array, and repair the particle links.
 *
 * If most of the array is displaced (the first-ever rebuild, say) the
 * gather buffer would grow towards the size of the particle arrays, so
 * above a threshold the sort falls back to the original serial
 * cycle-following walk, which moves everything and allocates nothing.
 */

/*! Slots per scan chunk. Fixed, so chunk boundaries (and therefore the
 * order of the compacted displaced list) do not depend on the number of
 * threads or their scheduling. */
#define space_sort_chunk_size ((size_t)131072)

/*! Fall back to the serial sort when more than 1/this of the particles are
 * displaced, bounding the gather buffer. */
#define space_sort_fallback_denominator ((size_t)4)

/**
 * @brief The move list of one sort: who is out of place, and where each
 *        one goes.
 */
struct space_sort_moves {

  /*! Displaced slots, in ascending order. */
  size_t *displaced;

  /*! Destination slot of each displaced particle. */
  size_t *dest;

  /*! Destination bin of each displaced particle. */
  int *target;

  /*! Number of displaced particles. */
  size_t nr_moves;

  /*! Total number of particles. */
  size_t nr;
};

/**
 * @brief Shared state for the scan mappers.
 */
struct space_sort_scan_data {

  /*! Bin index of each particle. */
  const int *ind;

  /*! First slot of each bin (num_bins + 1 entries). */
  const size_t *offsets;

  /*! Per-chunk displaced counts, then (after the prefix sum) write starts. */
  size_t *chunk_cursors;

  /*! The compacted displaced list (NULL during the counting pass). */
  size_t *displaced;

  /*! Total number of particles. */
  size_t nr;

  /*! Number of bins. */
  int num_bins;
};

/**
 * @brief Find the bin whose range contains a slot.
 *
 * @param offsets The bin offsets.
 * @param num_bins Number of bins.
 * @param k The slot.
 */
static int space_sort_home_bin(const size_t *offsets, const int num_bins,
                               const size_t k) {

  int lo = 0, hi = num_bins;
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (offsets[mid] <= k)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

/**
 * @brief #threadpool mapper counting the displaced slots in each chunk.
 *
 * @param map_data Pointer to a range of chunk ids.
 * @param num_chunks Number of chunk ids in this batch.
 * @param extra_data The #space_sort_scan_data.
 */
static void space_sort_count_mapper(void *map_data, int num_chunks,
                                    void *extra_data) {

  struct space_sort_scan_data *data = (struct space_sort_scan_data *)extra_data;
  const int *chunks = (const int *)map_data;
  const int *ind = data->ind;
  const size_t *offsets = data->offsets;

  for (int i = 0; i < num_chunks; i++) {
    const size_t c = (size_t)chunks[i];
    const size_t start = c * space_sort_chunk_size;
    size_t end = start + space_sort_chunk_size;
    if (end > data->nr) end = data->nr;

    /* Track the home bin monotonically from the chunk's first slot. */
    int b = space_sort_home_bin(offsets, data->num_bins, start);
    size_t count = 0;
    for (size_t k = start; k < end; k++) {
      while (offsets[b + 1] <= k) b++;
#ifdef SWIFT_DEBUG_CHECKS
      if (ind[k] < 0 || ind[k] >= data->num_bins)
        error("Invalid bin index %d (num_bins=%d) at slot %zu.", ind[k],
              data->num_bins, k);
#endif
      if (ind[k] != b) count++;
    }
    data->chunk_cursors[c] = count;
  }
}

/**
 * @brief #threadpool mapper compacting the displaced slots of each chunk.
 *
 * Each chunk writes into its own exclusive range of the list, handed to it
 * by the prefix sum, so there is no contention and the list comes out in
 * ascending slot order.
 *
 * @param map_data Pointer to a range of chunk ids.
 * @param num_chunks Number of chunk ids in this batch.
 * @param extra_data The #space_sort_scan_data.
 */
static void space_sort_fill_mapper(void *map_data, int num_chunks,
                                   void *extra_data) {

  struct space_sort_scan_data *data = (struct space_sort_scan_data *)extra_data;
  const int *chunks = (const int *)map_data;
  const int *ind = data->ind;
  const size_t *offsets = data->offsets;

  for (int i = 0; i < num_chunks; i++) {
    const size_t c = (size_t)chunks[i];
    const size_t start = c * space_sort_chunk_size;
    size_t end = start + space_sort_chunk_size;
    if (end > data->nr) end = data->nr;

    int b = space_sort_home_bin(offsets, data->num_bins, start);
    size_t cursor = data->chunk_cursors[c];
    for (size_t k = start; k < end; k++) {
      while (offsets[b + 1] <= k) b++;
      if (ind[k] != b) data->displaced[cursor++] = k;
    }
  }
}

/**
 * @brief Build the move list for one sort: find the displaced slots and
 *        assign each orphan a hole.
 *
 * @param tp The #threadpool.
 * @param ind Bin index of each particle.
 * @param counts Number of particles per bin (not modified).
 * @param num_bins Number of bins.
 * @param moves (return) The move list. Only valid when 1 is returned;
 *        nr_moves and nr are set either way, for reporting.
 * @return 1 if the move list was built, 0 if the caller should fall back to
 *         the serial sort.
 */
static int space_sort_plan_moves(struct threadpool *tp, const int *ind,
                                 const int *counts, const int num_bins,
                                 struct space_sort_moves *moves) {

  moves->displaced = NULL;
  moves->dest = NULL;
  moves->target = NULL;

  /* Bin offsets from the (untouched) counts. */
  size_t *offsets = NULL;
  if (swift_memalign("sort_offsets", (void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");
  size_t nr = 0;
  for (int k = 0; k < num_bins; k++) {
    offsets[k] = nr;
    nr += counts[k];
  }
  offsets[num_bins] = nr;
  moves->nr = nr;
  moves->nr_moves = 0;
  if (nr == 0) {
    swift_free("sort_offsets", offsets);
    return 1;
  }

  /* Cut the slots into fixed chunks and count the displaced in each. */
  const size_t nr_chunks =
      (nr + space_sort_chunk_size - 1) / space_sort_chunk_size;
  size_t *chunk_cursors = NULL;
  int *chunk_ids = NULL;
  if (swift_memalign("sort_chunks", (void **)&chunk_cursors,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(size_t) * nr_chunks) != 0 ||
      swift_memalign("sort_chunk_ids", (void **)&chunk_ids,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(int) * nr_chunks) != 0)
    error("Failed to allocate sort chunk arrays.");
  for (size_t c = 0; c < nr_chunks; c++) chunk_ids[c] = (int)c;

  struct space_sort_scan_data data;
  data.ind = ind;
  data.offsets = offsets;
  data.chunk_cursors = chunk_cursors;
  data.displaced = NULL;
  data.nr = nr;
  data.num_bins = num_bins;

  threadpool_map(tp, space_sort_count_mapper, chunk_ids, nr_chunks,
                 sizeof(int), threadpool_auto_chunk_size, &data);

  /* Prefix-sum the chunk counts into exclusive write starts, giving D. */
  size_t nr_moves = 0;
  for (size_t c = 0; c < nr_chunks; c++) {
    const size_t count = chunk_cursors[c];
    chunk_cursors[c] = nr_moves;
    nr_moves += count;
  }
  moves->nr_moves = nr_moves;

  /* Nothing out of place, or too much? */
  if (nr_moves == 0 || nr_moves > nr / space_sort_fallback_denominator) {
    swift_free("sort_offsets", offsets);
    swift_free("sort_chunks", chunk_cursors);
    swift_free("sort_chunk_ids", chunk_ids);
    return nr_moves == 0;
  }

  /* Compact the displaced slots, in ascending order. */
  if (swift_memalign("sort_displaced", (void **)&moves->displaced,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(size_t) * nr_moves) != 0 ||
      swift_memalign("sort_dest", (void **)&moves->dest, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * nr_moves) != 0 ||
      swift_memalign("sort_target", (void **)&moves->target,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(int) * nr_moves) != 0)
    error("Failed to allocate sort move list.");
  data.displaced = moves->displaced;

  threadpool_map(tp, space_sort_fill_mapper, chunk_ids, nr_chunks, sizeof(int),
                 threadpool_auto_chunk_size, &data);

  swift_free("sort_chunks", chunk_cursors);
  swift_free("sort_chunk_ids", chunk_ids);

  /* Match holes to orphans. The list is slot-ascending, so each bin's holes
   * are one contiguous run; find the runs and capture each orphan's target
   * while ind is still untouched. */
  size_t *hole_start = NULL;
  if (swift_memalign("sort_holes", (void **)&hole_start, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * 2 * num_bins) != 0)
    error("Failed to allocate sort hole runs.");
  size_t *hole_cursor = hole_start + num_bins;

  int b = 0;
  for (int k = 0; k < num_bins; k++) hole_start[k] = 0;
  for (int k = 0; k < num_bins; k++) hole_cursor[k] = 0;
  for (size_t i = 0; i < nr_moves; i++) {
    const size_t slot = moves->displaced[i];
    while (offsets[b + 1] <= slot) b++;
    if (hole_cursor[b] == 0) hole_start[b] = i;
    hole_cursor[b]++;
    moves->target[i] = ind[slot];
  }

  /* Hand the i-th orphan of bin t the i-th hole of bin t. */
#ifdef SWIFT_DEBUG_CHECKS
  size_t *hole_count = NULL;
  if (swift_memalign("sort_hole_count", (void **)&hole_count,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(size_t) * num_bins) != 0)
    error("Failed to allocate hole counts.");
  memcpy(hole_count, hole_cursor, sizeof(size_t) * num_bins);
#endif
  for (int k = 0; k < num_bins; k++) hole_cursor[k] = 0;
  for (size_t i = 0; i < nr_moves; i++) {
    const int t = moves->target[i];
    moves->dest[i] = moves->displaced[hole_start[t] + hole_cursor[t]++];
  }
#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (hole_cursor[k] != hole_count[k])
      error("Bin %d has %zu orphans but %zu holes.", k, hole_cursor[k],
            hole_count[k]);
  swift_free("sort_hole_count", hole_count);
#endif

  swift_free("sort_holes", hole_start);
  swift_free("sort_offsets", offsets);
  return 1;
}

/**
 * @brief Release a move list.
 */
static void space_sort_free_moves(struct space_sort_moves *moves) {

  if (moves->displaced != NULL) swift_free("sort_displaced", moves->displaced);
  if (moves->dest != NULL) swift_free("sort_dest", moves->dest);
  if (moves->target != NULL) swift_free("sort_target", moves->target);
}

/**
 * @brief Shared state for the gather/scatter mappers.
 *
 * The copies are family-agnostic (a memcpy of one or two element arrays);
 * only the link repair is per-family and lives in its own pass.
 */
struct space_sort_copy_data {

  /*! The particle array being sorted, as bytes. */
  char *arr;

  /*! The gather buffer for it. */
  char *buf;

  /*! Element size of the first array. */
  size_t size;

  /*! Optional second array moved in lock-step (the xparts), or NULL. */
  char *arr2;

  /*! The gather buffer for it. */
  char *buf2;

  /*! Element size of the second array. */
  size_t size2;

  /*! Bin index of each particle, updated by the scatter. */
  int *ind;

  /*! The move list. */
  const struct space_sort_moves *moves;
};

/**
 * @brief #threadpool mapper gathering the displaced particles.
 *
 * @param map_data Pointer to a range of the displaced list.
 * @param num Number of entries in this batch.
 * @param extra_data The #space_sort_copy_data.
 */
static void space_sort_gather_mapper(void *map_data, int num,
                                     void *extra_data) {

  struct space_sort_copy_data *data = (struct space_sort_copy_data *)extra_data;
  const size_t *displaced = (const size_t *)map_data;
  const size_t first = displaced - data->moves->displaced;

  for (int i = 0; i < num; i++) {
    const size_t src = displaced[i];
    const size_t out = first + i;
    memcpy(data->buf + out * data->size, data->arr + src * data->size,
           data->size);
    if (data->arr2 != NULL)
      memcpy(data->buf2 + out * data->size2, data->arr2 + src * data->size2,
             data->size2);
  }
}

/**
 * @brief #threadpool mapper scattering the gathered particles to their
 *        holes.
 *
 * Every destination is a distinct displaced slot, so no two writes collide
 * and every stale slot is overwritten. The full gather has completed before
 * this runs.
 *
 * @param map_data Pointer to a range of the displaced list.
 * @param num Number of entries in this batch.
 * @param extra_data The #space_sort_copy_data.
 */
static void space_sort_scatter_mapper(void *map_data, int num,
                                      void *extra_data) {

  struct space_sort_copy_data *data = (struct space_sort_copy_data *)extra_data;
  const size_t *displaced = (const size_t *)map_data;
  const size_t first = displaced - data->moves->displaced;

  for (int i = 0; i < num; i++) {
    const size_t out = first + i;
    const size_t d = data->moves->dest[out];
    memcpy(data->arr + d * data->size, data->buf + out * data->size,
           data->size);
    if (data->arr2 != NULL)
      memcpy(data->arr2 + d * data->size2, data->buf2 + out * data->size2,
             data->size2);
    data->ind[d] = data->moves->target[out];
  }
}

/**
 * @brief Run the gather and scatter for one family.
 *
 * @param tp The #threadpool.
 * @param data The #space_sort_copy_data, with the buffers allocated.
 */
static void space_sort_move(struct threadpool *tp,
                            struct space_sort_copy_data *data) {

  threadpool_map(tp, space_sort_gather_mapper, (void *)data->moves->displaced,
                 data->moves->nr_moves, sizeof(size_t),
                 threadpool_auto_chunk_size, data);
  threadpool_map(tp, space_sort_scatter_mapper, (void *)data->moves->displaced,
                 data->moves->nr_moves, sizeof(size_t),
                 threadpool_auto_chunk_size, data);
}

/**
 * @brief Serial cycle-following sort of the particles, used above the
 *        displaced-fraction threshold.
 *
 * Kept verbatim from the original implementation: it moves everything but
 * allocates nothing, which is the right trade when most of the array is out
 * of place (the first-ever rebuild, a fresh redistribution).
 */
static void space_parts_sort_serial(struct part *parts, struct xpart *xparts,
                                    int *restrict ind, int *restrict counts,
                                    int num_bins, ptrdiff_t parts_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (swift_memalign("parts_offsets", (void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct part temp_part = parts[k];
      struct xpart temp_xpart = xparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&parts[j], &temp_part, sizeof(struct part));
        memswap(&xparts[j], &temp_xpart, sizeof(struct xpart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (parts[j].gpart)
          parts[j].gpart->id_or_neg_offset = -(j + parts_offset);
      }
      parts[k] = temp_part;
      xparts[k] = temp_xpart;
      ind[k] = target_cid;
      if (parts[k].gpart)
        parts[k].gpart->id_or_neg_offset = -(k + parts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  swift_free("parts_offsets", offsets);
}

/**
 * @brief Shared state for the hydro link repair.
 */
struct space_parts_relink_data {
  struct part *parts;
  const struct space_sort_moves *moves;
  ptrdiff_t parts_offset;
};

/**
 * @brief #threadpool mapper repairing the gpart back-links of the moved
 *        hydro particles. Each move touches its own gpart, so no two
 *        entries write the same counterpart.
 */
static void space_parts_relink_mapper(void *map_data, int num,
                                      void *extra_data) {

  struct space_parts_relink_data *data =
      (struct space_parts_relink_data *)extra_data;
  const size_t *dest = (const size_t *)map_data;

  for (int i = 0; i < num; i++) {
    const size_t d = dest[i];
    if (data->parts[d].gpart)
      data->parts[d].gpart->id_or_neg_offset = -(d + data->parts_offset);
  }
}

/**
 * @brief Sort the particles and condensed particles according to the given
 * indices.
 *
 * @param tp The #threadpool to parallelise over.
 * @param parts The array of #part to sort.
 * @param xparts The corresponding #xpart array to sort as well.
 * @param ind The indices with respect to which the parts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of count).
 * @param parts_offset Offset of the #part array from the global #part array.
 * @param verbose Are we talkative?
 */
void space_parts_sort(struct threadpool *tp, struct part *parts,
                      struct xpart *xparts, int *restrict ind,
                      int *restrict counts, int num_bins,
                      ptrdiff_t parts_offset, int verbose) {

  struct space_sort_moves moves;
  if (!space_sort_plan_moves(tp, ind, counts, num_bins, &moves)) {

    /* Too much of the array is out of place: move everything, in place. */
    space_parts_sort_serial(parts, xparts, ind, counts, num_bins,
                            parts_offset);
  } else if (moves.nr_moves > 0) {

    struct space_sort_copy_data data;
    data.arr = (char *)parts;
    data.size = sizeof(struct part);
    data.arr2 = (char *)xparts;
    data.size2 = sizeof(struct xpart);
    data.ind = ind;
    data.moves = &moves;
    if (swift_memalign("sort_buff", (void **)&data.buf, part_align,
                       moves.nr_moves * sizeof(struct part)) != 0 ||
        swift_memalign("sort_buff_x", (void **)&data.buf2, xpart_align,
                       moves.nr_moves * sizeof(struct xpart)) != 0)
      error("Failed to allocate particle gather buffer.");

    space_sort_move(tp, &data);

    struct space_parts_relink_data rdata = {parts, &moves, parts_offset};
    threadpool_map(tp, space_parts_relink_mapper, moves.dest, moves.nr_moves,
                   sizeof(size_t), threadpool_auto_chunk_size, &rdata);

    swift_free("sort_buff", data.buf);
    swift_free("sort_buff_x", data.buf2);
    space_sort_free_moves(&moves);
  }

  if (verbose)
    message("%zu of %zu parts were displaced (%.2f %%).", moves.nr_moves,
            moves.nr, moves.nr > 0 ? 100. * moves.nr_moves / moves.nr : 0.);
}

/**
 * @brief Serial cycle-following sort of the s-particles (see
 *        space_parts_sort_serial()).
 */
static void space_sparts_sort_serial(struct spart *sparts, int *restrict ind,
                                     int *restrict counts, int num_bins,
                                     ptrdiff_t sparts_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (swift_memalign("sparts_offsets", (void **)&offsets,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct spart temp_spart = sparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&sparts[j], &temp_spart, sizeof(struct spart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (sparts[j].gpart)
          sparts[j].gpart->id_or_neg_offset = -(j + sparts_offset);
      }
      sparts[k] = temp_spart;
      ind[k] = target_cid;
      if (sparts[k].gpart)
        sparts[k].gpart->id_or_neg_offset = -(k + sparts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  swift_free("sparts_offsets", offsets);
}

/**
 * @brief Shared state for the star link repair.
 */
struct space_sparts_relink_data {
  struct spart *sparts;
  const struct space_sort_moves *moves;
  ptrdiff_t sparts_offset;
};

/**
 * @brief #threadpool mapper repairing the gpart back-links of the moved
 *        star particles.
 */
static void space_sparts_relink_mapper(void *map_data, int num,
                                       void *extra_data) {

  struct space_sparts_relink_data *data =
      (struct space_sparts_relink_data *)extra_data;
  const size_t *dest = (const size_t *)map_data;

  for (int i = 0; i < num; i++) {
    const size_t d = dest[i];
    if (data->sparts[d].gpart)
      data->sparts[d].gpart->id_or_neg_offset = -(d + data->sparts_offset);
  }
}

/**
 * @brief Sort the s-particles according to the given indices.
 *
 * @param tp The #threadpool to parallelise over.
 * @param sparts The array of #spart to sort.
 * @param ind The indices with respect to which the #spart are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param sparts_offset Offset of the #spart array from the global #spart.
 * array.
 * @param verbose Are we talkative?
 */
void space_sparts_sort(struct threadpool *tp, struct spart *sparts,
                       int *restrict ind, int *restrict counts, int num_bins,
                       ptrdiff_t sparts_offset, int verbose) {

  struct space_sort_moves moves;
  if (!space_sort_plan_moves(tp, ind, counts, num_bins, &moves)) {

    space_sparts_sort_serial(sparts, ind, counts, num_bins, sparts_offset);
  } else if (moves.nr_moves > 0) {

    struct space_sort_copy_data data;
    data.arr = (char *)sparts;
    data.size = sizeof(struct spart);
    data.arr2 = NULL;
    data.buf2 = NULL;
    data.size2 = 0;
    data.ind = ind;
    data.moves = &moves;
    if (swift_memalign("sort_buff", (void **)&data.buf, spart_align,
                       moves.nr_moves * sizeof(struct spart)) != 0)
      error("Failed to allocate particle gather buffer.");

    space_sort_move(tp, &data);

    struct space_sparts_relink_data rdata = {sparts, &moves, sparts_offset};
    threadpool_map(tp, space_sparts_relink_mapper, moves.dest, moves.nr_moves,
                   sizeof(size_t), threadpool_auto_chunk_size, &rdata);

    swift_free("sort_buff", data.buf);
    space_sort_free_moves(&moves);
  }

  if (verbose)
    message("%zu of %zu sparts were displaced (%.2f %%).", moves.nr_moves,
            moves.nr, moves.nr > 0 ? 100. * moves.nr_moves / moves.nr : 0.);
}

/**
 * @brief Serial cycle-following sort of the b-particles (see
 *        space_parts_sort_serial()).
 */
static void space_bparts_sort_serial(struct bpart *bparts, int *restrict ind,
                                     int *restrict counts, int num_bins,
                                     ptrdiff_t bparts_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (swift_memalign("bparts_offsets", (void **)&offsets,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct bpart temp_bpart = bparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&bparts[j], &temp_bpart, sizeof(struct bpart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (bparts[j].gpart)
          bparts[j].gpart->id_or_neg_offset = -(j + bparts_offset);
      }
      bparts[k] = temp_bpart;
      ind[k] = target_cid;
      if (bparts[k].gpart)
        bparts[k].gpart->id_or_neg_offset = -(k + bparts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  swift_free("bparts_offsets", offsets);
}

/**
 * @brief Shared state for the black hole link repair.
 */
struct space_bparts_relink_data {
  struct bpart *bparts;
  const struct space_sort_moves *moves;
  ptrdiff_t bparts_offset;
};

/**
 * @brief #threadpool mapper repairing the gpart back-links of the moved
 *        black holes.
 */
static void space_bparts_relink_mapper(void *map_data, int num,
                                       void *extra_data) {

  struct space_bparts_relink_data *data =
      (struct space_bparts_relink_data *)extra_data;
  const size_t *dest = (const size_t *)map_data;

  for (int i = 0; i < num; i++) {
    const size_t d = dest[i];
    if (data->bparts[d].gpart)
      data->bparts[d].gpart->id_or_neg_offset = -(d + data->bparts_offset);
  }
}

/**
 * @brief Sort the b-particles according to the given indices.
 *
 * @param tp The #threadpool to parallelise over.
 * @param bparts The array of #bpart to sort.
 * @param ind The indices with respect to which the #bpart are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param bparts_offset Offset of the #bpart array from the global #bpart.
 * array.
 * @param verbose Are we talkative?
 */
void space_bparts_sort(struct threadpool *tp, struct bpart *bparts,
                       int *restrict ind, int *restrict counts, int num_bins,
                       ptrdiff_t bparts_offset, int verbose) {

  struct space_sort_moves moves;
  if (!space_sort_plan_moves(tp, ind, counts, num_bins, &moves)) {

    space_bparts_sort_serial(bparts, ind, counts, num_bins, bparts_offset);
  } else if (moves.nr_moves > 0) {

    struct space_sort_copy_data data;
    data.arr = (char *)bparts;
    data.size = sizeof(struct bpart);
    data.arr2 = NULL;
    data.buf2 = NULL;
    data.size2 = 0;
    data.ind = ind;
    data.moves = &moves;
    if (swift_memalign("sort_buff", (void **)&data.buf, bpart_align,
                       moves.nr_moves * sizeof(struct bpart)) != 0)
      error("Failed to allocate particle gather buffer.");

    space_sort_move(tp, &data);

    struct space_bparts_relink_data rdata = {bparts, &moves, bparts_offset};
    threadpool_map(tp, space_bparts_relink_mapper, moves.dest, moves.nr_moves,
                   sizeof(size_t), threadpool_auto_chunk_size, &rdata);

    swift_free("sort_buff", data.buf);
    space_sort_free_moves(&moves);
  }

  if (verbose)
    message("%zu of %zu bparts were displaced (%.2f %%).", moves.nr_moves,
            moves.nr, moves.nr > 0 ? 100. * moves.nr_moves / moves.nr : 0.);
}

/**
 * @brief Serial cycle-following sort of the sink-particles (see
 *        space_parts_sort_serial()).
 */
static void space_sinks_sort_serial(struct sink *sinks, int *restrict ind,
                                    int *restrict counts, int num_bins,
                                    ptrdiff_t sinks_offset) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (swift_memalign("sinks_offsets", (void **)&offsets, SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct sink temp_sink = sinks[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap(&sinks[j], &temp_sink, sizeof(struct sink));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (sinks[j].gpart)
          sinks[j].gpart->id_or_neg_offset = -(j + sinks_offset);
      }
      sinks[k] = temp_sink;
      ind[k] = target_cid;
      if (sinks[k].gpart)
        sinks[k].gpart->id_or_neg_offset = -(k + sinks_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  swift_free("sinks_offsets", offsets);
}

/**
 * @brief Shared state for the sink link repair.
 */
struct space_sinks_relink_data {
  struct sink *sinks;
  const struct space_sort_moves *moves;
  ptrdiff_t sinks_offset;
};

/**
 * @brief #threadpool mapper repairing the gpart back-links of the moved
 *        sinks.
 */
static void space_sinks_relink_mapper(void *map_data, int num,
                                      void *extra_data) {

  struct space_sinks_relink_data *data =
      (struct space_sinks_relink_data *)extra_data;
  const size_t *dest = (const size_t *)map_data;

  for (int i = 0; i < num; i++) {
    const size_t d = dest[i];
    if (data->sinks[d].gpart)
      data->sinks[d].gpart->id_or_neg_offset = -(d + data->sinks_offset);
  }
}

/**
 * @brief Sort the sink-particles according to the given indices.
 *
 * @param tp The #threadpool to parallelise over.
 * @param sinks The array of #sink to sort.
 * @param ind The indices with respect to which the #sink are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param sinks_offset Offset of the #sink array from the global #sink.
 * array.
 * @param verbose Are we talkative?
 */
void space_sinks_sort(struct threadpool *tp, struct sink *sinks,
                      int *restrict ind, int *restrict counts, int num_bins,
                      ptrdiff_t sinks_offset, int verbose) {

  struct space_sort_moves moves;
  if (!space_sort_plan_moves(tp, ind, counts, num_bins, &moves)) {

    space_sinks_sort_serial(sinks, ind, counts, num_bins, sinks_offset);
  } else if (moves.nr_moves > 0) {

    struct space_sort_copy_data data;
    data.arr = (char *)sinks;
    data.size = sizeof(struct sink);
    data.arr2 = NULL;
    data.buf2 = NULL;
    data.size2 = 0;
    data.ind = ind;
    data.moves = &moves;
    if (swift_memalign("sort_buff", (void **)&data.buf, sink_align,
                       moves.nr_moves * sizeof(struct sink)) != 0)
      error("Failed to allocate particle gather buffer.");

    space_sort_move(tp, &data);

    struct space_sinks_relink_data rdata = {sinks, &moves, sinks_offset};
    threadpool_map(tp, space_sinks_relink_mapper, moves.dest, moves.nr_moves,
                   sizeof(size_t), threadpool_auto_chunk_size, &rdata);

    swift_free("sort_buff", data.buf);
    space_sort_free_moves(&moves);
  }

  if (verbose)
    message("%zu of %zu sinks were displaced (%.2f %%).", moves.nr_moves,
            moves.nr, moves.nr > 0 ? 100. * moves.nr_moves / moves.nr : 0.);
}

/**
 * @brief Serial cycle-following sort of the g-particles (see
 *        space_parts_sort_serial()).
 */
static void space_gparts_sort_serial(struct gpart *gparts, struct part *parts,
                                     struct sink *sinks, struct spart *sparts,
                                     struct bpart *bparts, int *restrict ind,
                                     int *restrict counts, int num_bins) {
  /* Create the offsets array. */
  size_t *offsets = NULL;
  if (swift_memalign("gparts_offsets", (void **)&offsets,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(size_t) * (num_bins + 1)) != 0)
    error("Failed to allocate temporary cell offsets array.");

  offsets[0] = 0;
  for (int k = 1; k <= num_bins; k++) {
    offsets[k] = offsets[k - 1] + counts[k - 1];
    counts[k - 1] = 0;
  }

  /* Loop over local cells. */
  for (int cid = 0; cid < num_bins; cid++) {
    for (size_t k = offsets[cid] + counts[cid]; k < offsets[cid + 1]; k++) {
      counts[cid]++;
      int target_cid = ind[k];
      if (target_cid == cid) {
        continue;
      }
      struct gpart temp_gpart = gparts[k];
      while (target_cid != cid) {
        size_t j = offsets[target_cid] + counts[target_cid]++;
        while (ind[j] == target_cid) {
          j = offsets[target_cid] + counts[target_cid]++;
        }
        memswap_unaligned(&gparts[j], &temp_gpart, sizeof(struct gpart));
        memswap(&ind[j], &target_cid, sizeof(int));
        if (gparts[j].type == swift_type_gas) {
          parts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        } else if (gparts[j].type == swift_type_stars) {
          sparts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        } else if (gparts[j].type == swift_type_black_hole) {
          bparts[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        } else if (gparts[j].type == swift_type_sink) {
          sinks[-gparts[j].id_or_neg_offset].gpart = &gparts[j];
        }
      }
      gparts[k] = temp_gpart;
      ind[k] = target_cid;
      if (gparts[k].type == swift_type_gas) {
        parts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_stars) {
        sparts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_black_hole) {
        bparts[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_sink) {
        sinks[-gparts[k].id_or_neg_offset].gpart = &gparts[k];
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < num_bins; k++)
    if (offsets[k + 1] != offsets[k] + counts[k])
      error("Bad offsets after shuffle.");
#endif /* SWIFT_DEBUG_CHECKS */

  swift_free("gparts_offsets", offsets);
}

/**
 * @brief Shared state for the gravity link repair.
 */
struct space_gparts_relink_data {
  struct gpart *gparts;
  struct part *parts;
  struct sink *sinks;
  struct spart *sparts;
  struct bpart *bparts;
  const struct space_sort_moves *moves;
};

/**
 * @brief #threadpool mapper pointing the counterparts of the moved gravity
 *        particles back at their new slots. Each move touches its own
 *        counterpart, so no two entries write the same one.
 */
static void space_gparts_relink_mapper(void *map_data, int num,
                                       void *extra_data) {

  struct space_gparts_relink_data *data =
      (struct space_gparts_relink_data *)extra_data;
  const size_t *dest = (const size_t *)map_data;

  for (int i = 0; i < num; i++) {
    struct gpart *gp = &data->gparts[dest[i]];
    if (gp->type == swift_type_gas) {
      data->parts[-gp->id_or_neg_offset].gpart = gp;
    } else if (gp->type == swift_type_stars) {
      data->sparts[-gp->id_or_neg_offset].gpart = gp;
    } else if (gp->type == swift_type_black_hole) {
      data->bparts[-gp->id_or_neg_offset].gpart = gp;
    } else if (gp->type == swift_type_sink) {
      data->sinks[-gp->id_or_neg_offset].gpart = gp;
    }
  }
}

/**
 * @brief Sort the g-particles according to the given indices.
 *
 * @param tp The #threadpool to parallelise over.
 * @param gparts The array of #gpart to sort.
 * @param parts Global #part array for re-linking.
 * @param sinks Global #sink array for re-linking.
 * @param sparts Global #spart array for re-linking.
 * @param bparts Global #bpart array for re-linking.
 * @param ind The indices with respect to which the gparts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param verbose Are we talkative?
 */
void space_gparts_sort(struct threadpool *tp, struct gpart *gparts,
                       struct part *parts, struct sink *sinks,
                       struct spart *sparts, struct bpart *bparts,
                       int *restrict ind, int *restrict counts, int num_bins,
                       int verbose) {

  struct space_sort_moves moves;
  if (!space_sort_plan_moves(tp, ind, counts, num_bins, &moves)) {

    space_gparts_sort_serial(gparts, parts, sinks, sparts, bparts, ind, counts,
                             num_bins);
  } else if (moves.nr_moves > 0) {

    struct space_sort_copy_data data;
    data.arr = (char *)gparts;
    data.size = sizeof(struct gpart);
    data.arr2 = NULL;
    data.buf2 = NULL;
    data.size2 = 0;
    data.ind = ind;
    data.moves = &moves;
    if (swift_memalign("sort_buff", (void **)&data.buf, gpart_align,
                       moves.nr_moves * sizeof(struct gpart)) != 0)
      error("Failed to allocate particle gather buffer.");

    space_sort_move(tp, &data);

    struct space_gparts_relink_data rdata = {gparts, parts,  sinks,
                                             sparts, bparts, &moves};
    threadpool_map(tp, space_gparts_relink_mapper, moves.dest, moves.nr_moves,
                   sizeof(size_t), threadpool_auto_chunk_size, &rdata);

    swift_free("sort_buff", data.buf);
    space_sort_free_moves(&moves);
  }

  if (verbose)
    message("%zu of %zu gparts were displaced (%.2f %%).", moves.nr_moves,
            moves.nr, moves.nr > 0 ? 100. * moves.nr_moves / moves.nr : 0.);
}
