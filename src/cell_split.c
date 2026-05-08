/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* This object's header. */
#include "cell.h"

/* Local headers. */
#include "lock.h"
#include "memswap.h"
#include "threadpool.h"

/**
 * @brief Per-thread scratch for the bucket ids used by #cell_split.
 *
 * The split now computes bucket ids directly from the particle arrays rather
 * than receiving pre-filled position buffers from the caller. We keep one
 * grow-only integer array per threadpool worker and reuse it across calls,
 * sized to the largest particle family count that worker has seen so far.
 */
struct cell_split_bid_scratch {
  int *bids;
  int capacity;
};

static struct cell_split_bid_scratch *cell_split_bid_scratch = NULL;
static int cell_split_bid_scratch_count = 0;
static swift_lock_type cell_split_bid_scratch_lock = lock_static_initializer;

/**
 * @brief Ensure the global per-thread bid scratch table exists and return the
 *        entry for the current worker, grown to at least @c needed entries.
 *
 * @param needed Required capacity for the returned worker-local @c bids array.
 * @return Pointer to the current worker's scratch entry.
 */
static struct cell_split_bid_scratch *cell_split_get_bid_scratch(
    const int needed) {

  const int tid = threadpool_gettid();

  if (tid >= cell_split_bid_scratch_count) {
    lock_lock(&cell_split_bid_scratch_lock);

    if (tid >= cell_split_bid_scratch_count) {
      const int new_count = tid + 1;
      struct cell_split_bid_scratch *new_scratch =
          (struct cell_split_bid_scratch *)swift_malloc(
              "cell_split_bid_scratch",
              sizeof(struct cell_split_bid_scratch) * new_count);
      if (new_scratch == NULL)
        error("Failed to allocate cell_split bid scratch table.");

      bzero(new_scratch, sizeof(struct cell_split_bid_scratch) * new_count);

      if (cell_split_bid_scratch != NULL) {
        memcpy(new_scratch, cell_split_bid_scratch,
               sizeof(struct cell_split_bid_scratch) *
                   cell_split_bid_scratch_count);
        swift_free("cell_split_bid_scratch", cell_split_bid_scratch);
      }

      cell_split_bid_scratch = new_scratch;
      cell_split_bid_scratch_count = new_count;
    }

    if (lock_unlock(&cell_split_bid_scratch_lock) != 0)
      error("Failed to unlock cell_split bid scratch lock.");
  }

  struct cell_split_bid_scratch *scratch = &cell_split_bid_scratch[tid];
  if (needed > scratch->capacity) {
    if (scratch->bids != NULL) swift_free("cell_split_bids", scratch->bids);
    scratch->bids =
        (int *)swift_malloc("cell_split_bids", sizeof(int) * needed);
    if (scratch->bids == NULL)
      error("Failed to allocate cell_split bid scratch of size %d.", needed);
    scratch->capacity = needed;
  }

  return scratch;
}

/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 * @param parts_offset Offset of the cell parts array relative to the
 *        space's parts array, i.e. c->hydro.parts - s->parts.
 * @param sparts_offset Offset of the cell sparts array relative to the
 *        space's sparts array, i.e. c->stars.parts - s->stars.parts.
 * @param bparts_offset Offset of the cell bparts array relative to the
 *        space's bparts array, i.e. c->black_holes.parts -
 * s->black_holes.parts.
 * @param sinks_offset Offset of the cell sink array relative to the
 *        space's sink array, i.e. c->sinks.parts - s->sinks.parts.
 */
void cell_split(struct cell *c, const ptrdiff_t parts_offset,
                const ptrdiff_t sparts_offset,
                const ptrdiff_t bparts_offset,
                const ptrdiff_t sinks_offset) {

  const int count = c->hydro.count, gcount = c->grav.count,
            scount = c->stars.count, bcount = c->black_holes.count,
            sink_count = c->sinks.count;
  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct sink *sinks = c->sinks.parts;
  const double pivot[3] = {c->loc[0] + c->width[0] / 2,
                           c->loc[1] + c->width[1] / 2,
                           c->loc[2] + c->width[2] / 2};
  int bucket_count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int bucket_offset[9];

  int max_count = count;
  max_count = max(max_count, gcount);
  max_count = max(max_count, scount);
  max_count = max(max_count, bcount);
  max_count = max(max_count, sink_count);
  struct cell_split_bid_scratch *scratch = cell_split_get_bid_scratch(max_count);
  int *restrict bids = scratch->bids;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < count; k++) {
    const int bid = (parts[k].x[0] >= pivot[0]) * 4 +
                    (parts[k].x[1] >= pivot[1]) * 2 + (parts[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    bids[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = bids[k];
      if (bid != bucket) {
        struct part part = parts[k];
        struct xpart xpart = xparts[k];
        int held_bid = bids[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (bids[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&parts[j], &part, sizeof(struct part));
          memswap(&xparts[j], &xpart, sizeof(struct xpart));
          const int tmp_bid = bids[j];
          bids[j] = held_bid;
          held_bid = tmp_bid;
          if (parts[j].gpart)
            parts[j].gpart->id_or_neg_offset = -(j + parts_offset);
          bid = held_bid;
        }
        parts[k] = part;
        xparts[k] = xpart;
        bids[k] = held_bid;
        if (parts[k].gpart)
          parts[k].gpart->id_or_neg_offset = -(k + parts_offset);
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->hydro.count = bucket_count[k];
    c->progeny[k]->hydro.count_total = c->progeny[k]->hydro.count;
    c->progeny[k]->hydro.parts = &c->hydro.parts[bucket_offset[k]];
    c->progeny[k]->hydro.xparts = &c->hydro.xparts[bucket_offset[k]];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the buffs are OK. */
  for (int k = 1; k < count; k++) {
    if (bids[k] < bids[k - 1]) error("Hydro bid array not sorted.");
  }

  /* Verify that _all_ the parts have been assigned to a cell. */
  for (int k = 1; k < 8; k++)
    if (&c->progeny[k - 1]->hydro.parts[c->progeny[k - 1]->hydro.count] !=
        c->progeny[k]->hydro.parts)
      error("Particle sorting failed (internal consistency).");
  if (c->progeny[0]->hydro.parts != c->hydro.parts)
    error("Particle sorting failed (left edge).");
  if (&c->progeny[7]->hydro.parts[c->progeny[7]->hydro.count] !=
      &c->hydro.parts[count])
    error("Particle sorting failed (right edge).");

  /* Verify a few sub-cells. */
  for (int k = 0; k < c->progeny[0]->hydro.count; k++)
    if (c->progeny[0]->hydro.parts[k].x[0] >= pivot[0] ||
        c->progeny[0]->hydro.parts[k].x[1] >= pivot[1] ||
        c->progeny[0]->hydro.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=0).");
  for (int k = 0; k < c->progeny[1]->hydro.count; k++)
    if (c->progeny[1]->hydro.parts[k].x[0] >= pivot[0] ||
        c->progeny[1]->hydro.parts[k].x[1] >= pivot[1] ||
        c->progeny[1]->hydro.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=1).");
  for (int k = 0; k < c->progeny[2]->hydro.count; k++)
    if (c->progeny[2]->hydro.parts[k].x[0] >= pivot[0] ||
        c->progeny[2]->hydro.parts[k].x[1] < pivot[1] ||
        c->progeny[2]->hydro.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=2).");
  for (int k = 0; k < c->progeny[3]->hydro.count; k++)
    if (c->progeny[3]->hydro.parts[k].x[0] >= pivot[0] ||
        c->progeny[3]->hydro.parts[k].x[1] < pivot[1] ||
        c->progeny[3]->hydro.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=3).");
  for (int k = 0; k < c->progeny[4]->hydro.count; k++)
    if (c->progeny[4]->hydro.parts[k].x[0] < pivot[0] ||
        c->progeny[4]->hydro.parts[k].x[1] >= pivot[1] ||
        c->progeny[4]->hydro.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=4).");
  for (int k = 0; k < c->progeny[5]->hydro.count; k++)
    if (c->progeny[5]->hydro.parts[k].x[0] < pivot[0] ||
        c->progeny[5]->hydro.parts[k].x[1] >= pivot[1] ||
        c->progeny[5]->hydro.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=5).");
  for (int k = 0; k < c->progeny[6]->hydro.count; k++)
    if (c->progeny[6]->hydro.parts[k].x[0] < pivot[0] ||
        c->progeny[6]->hydro.parts[k].x[1] < pivot[1] ||
        c->progeny[6]->hydro.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=6).");
  for (int k = 0; k < c->progeny[7]->hydro.count; k++)
    if (c->progeny[7]->hydro.parts[k].x[0] < pivot[0] ||
        c->progeny[7]->hydro.parts[k].x[1] < pivot[1] ||
        c->progeny[7]->hydro.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=7).");
#endif

  /* Now do the same song and dance for the sparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < scount; k++) {
    const int bid = (sparts[k].x[0] > pivot[0]) * 4 +
                    (sparts[k].x[1] > pivot[1]) * 2 + (sparts[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    bids[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = bids[k];
      if (bid != bucket) {
        struct spart spart = sparts[k];
        int held_bid = bids[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (bids[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&sparts[j], &spart, sizeof(struct spart));
          const int tmp_bid = bids[j];
          bids[j] = held_bid;
          held_bid = tmp_bid;
          if (sparts[j].gpart)
            sparts[j].gpart->id_or_neg_offset = -(j + sparts_offset);
          bid = held_bid;
        }
        sparts[k] = spart;
        bids[k] = held_bid;
        if (sparts[k].gpart)
          sparts[k].gpart->id_or_neg_offset = -(k + sparts_offset);
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->stars.count = bucket_count[k];
    c->progeny[k]->stars.count_total = c->progeny[k]->stars.count;
    c->progeny[k]->stars.parts = &c->stars.parts[bucket_offset[k]];
    c->progeny[k]->stars.parts_rebuild = c->progeny[k]->stars.parts;
  }

  /* Now do the same song and dance for the bparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < bcount; k++) {
    const int bid = (bparts[k].x[0] > pivot[0]) * 4 +
                    (bparts[k].x[1] > pivot[1]) * 2 + (bparts[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    bids[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = bids[k];
      if (bid != bucket) {
        struct bpart bpart = bparts[k];
        int held_bid = bids[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (bids[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&bparts[j], &bpart, sizeof(struct bpart));
          const int tmp_bid = bids[j];
          bids[j] = held_bid;
          held_bid = tmp_bid;
          if (bparts[j].gpart)
            bparts[j].gpart->id_or_neg_offset = -(j + bparts_offset);
          bid = held_bid;
        }
        bparts[k] = bpart;
        bids[k] = held_bid;
        if (bparts[k].gpart)
          bparts[k].gpart->id_or_neg_offset = -(k + bparts_offset);
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->black_holes.count = bucket_count[k];
    c->progeny[k]->black_holes.count_total = c->progeny[k]->black_holes.count;
    c->progeny[k]->black_holes.parts = &c->black_holes.parts[bucket_offset[k]];
  }

  /* Now do the same song and dance for the sinks. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < sink_count; k++) {
    const int bid = (sinks[k].x[0] > pivot[0]) * 4 +
                    (sinks[k].x[1] > pivot[1]) * 2 +
                    (sinks[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    bids[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = bids[k];
      if (bid != bucket) {
        struct sink sink = sinks[k];
        int held_bid = bids[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (bids[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&sinks[j], &sink, sizeof(struct sink));
          const int tmp_bid = bids[j];
          bids[j] = held_bid;
          held_bid = tmp_bid;
          if (sinks[j].gpart)
            sinks[j].gpart->id_or_neg_offset = -(j + sinks_offset);
          bid = held_bid;
        }
        sinks[k] = sink;
        bids[k] = held_bid;
        if (sinks[k].gpart)
          sinks[k].gpart->id_or_neg_offset = -(k + sinks_offset);
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->sinks.count = bucket_count[k];
    c->progeny[k]->sinks.count_total = c->progeny[k]->sinks.count;
    c->progeny[k]->sinks.parts = &c->sinks.parts[bucket_offset[k]];
    c->progeny[k]->sinks.parts_rebuild = c->progeny[k]->sinks.parts;
  }

  /* Finally, do the same song and dance for the gparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < gcount; k++) {
    const int bid = (gparts[k].x[0] > pivot[0]) * 4 +
                    (gparts[k].x[1] > pivot[1]) * 2 + (gparts[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    bids[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap particles to their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = bids[k];
      if (bid != bucket) {
        struct gpart gpart = gparts[k];
        int held_bid = bids[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (bids[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap_unaligned(&gparts[j], &gpart, sizeof(struct gpart));
          const int tmp_bid = bids[j];
          bids[j] = held_bid;
          held_bid = tmp_bid;
          if (gparts[j].type == swift_type_gas) {
            parts[-gparts[j].id_or_neg_offset - parts_offset].gpart =
                &gparts[j];
          } else if (gparts[j].type == swift_type_stars) {
            sparts[-gparts[j].id_or_neg_offset - sparts_offset].gpart =
                &gparts[j];
          } else if (gparts[j].type == swift_type_sink) {
            sinks[-gparts[j].id_or_neg_offset - sinks_offset].gpart =
                &gparts[j];
          } else if (gparts[j].type == swift_type_black_hole) {
            bparts[-gparts[j].id_or_neg_offset - bparts_offset].gpart =
                &gparts[j];
          }
          bid = held_bid;
        }
        gparts[k] = gpart;
        bids[k] = held_bid;
        if (gparts[k].type == swift_type_gas) {
          parts[-gparts[k].id_or_neg_offset - parts_offset].gpart = &gparts[k];
        } else if (gparts[k].type == swift_type_stars) {
          sparts[-gparts[k].id_or_neg_offset - sparts_offset].gpart =
              &gparts[k];
        } else if (gparts[k].type == swift_type_sink) {
          sinks[-gparts[k].id_or_neg_offset - sinks_offset].gpart = &gparts[k];
        } else if (gparts[k].type == swift_type_black_hole) {
          bparts[-gparts[k].id_or_neg_offset - bparts_offset].gpart =
              &gparts[k];
        }
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->grav.count = bucket_count[k];
    c->progeny[k]->grav.count_total = c->progeny[k]->grav.count;
    c->progeny[k]->grav.parts = &c->grav.parts[bucket_offset[k]];
    c->progeny[k]->grav.parts_rebuild = c->progeny[k]->grav.parts;
  }
}

/**
 * @brief Re-arrange the #part in a top-level cell such that all the extra
 * ones for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param parts_offset The offset between the first #part in the array and the
 * first #part in the global array in the space structure (for re-linking).
 */
void cell_reorder_extra_parts(struct cell *c, const ptrdiff_t parts_offset) {
  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  const int count_real = c->hydro.count;

  if (c->depth != 0 || c->nodeID != engine_rank)
    error("This function should only be called on local top-level cells!");

  int first_not_extra = count_real;

  /* Find extra particles */
  for (int i = 0; i < count_real; ++i) {
    if (parts[i].time_bin == time_bin_not_created) {
      /* Find the first non-extra particle after the end of the
         real particles */
      while (parts[first_not_extra].time_bin == time_bin_not_created) {
        ++first_not_extra;
      }

#ifdef SWIFT_DEBUG_CHECKS
      if (first_not_extra >= count_real + space_extra_parts)
        error("Looking for extra particles beyond this cell's range!");
#endif

      /* Swap everything, including g-part pointer */
      memswap(&parts[i], &parts[first_not_extra], sizeof(struct part));
      memswap(&xparts[i], &xparts[first_not_extra], sizeof(struct xpart));
      if (parts[i].gpart)
        parts[i].gpart->id_or_neg_offset = -(i + parts_offset);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < c->hydro.count_total; ++i) {
    if (parts[i].time_bin == time_bin_not_created && i < c->hydro.count) {
      error("Extra particle before the end of the regular array");
    }
    if (parts[i].time_bin != time_bin_not_created && i >= c->hydro.count) {
      error("Regular particle after the end of the regular array");
    }
  }
#endif
}

/**
 * @brief Re-arrange the #spart in a top-level cell such that all the extra
 * ones for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param sparts_offset The offset between the first #spart in the array and
 * the first #spart in the global array in the space structure (for
 * re-linking).
 */
void cell_reorder_extra_sparts(struct cell *c, const ptrdiff_t sparts_offset) {
  struct spart *sparts = c->stars.parts;
  const int count_real = c->stars.count;

  if (c->depth != 0 || c->nodeID != engine_rank)
    error("This function should only be called on local top-level cells!");

  int first_not_extra = count_real;

  /* Find extra particles */
  for (int i = 0; i < count_real; ++i) {
    if (sparts[i].time_bin == time_bin_not_created) {
      /* Find the first non-extra particle after the end of the
         real particles */
      while (sparts[first_not_extra].time_bin == time_bin_not_created) {
        ++first_not_extra;
      }

#ifdef SWIFT_DEBUG_CHECKS
      if (first_not_extra >= count_real + space_extra_sparts)
        error("Looking for extra particles beyond this cell's range!");
#endif

      /* Swap everything, including g-part pointer */
      memswap(&sparts[i], &sparts[first_not_extra], sizeof(struct spart));
      if (sparts[i].gpart)
        sparts[i].gpart->id_or_neg_offset = -(i + sparts_offset);
      sparts[first_not_extra].gpart = NULL;
#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[first_not_extra].time_bin != time_bin_not_created)
        error("Incorrect swap occured!");
#endif
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < c->stars.count_total; ++i) {
    if (sparts[i].time_bin == time_bin_not_created && i < c->stars.count) {
      error("Extra particle before the end of the regular array");
    }
    if (sparts[i].time_bin != time_bin_not_created && i >= c->stars.count) {
      error("Regular particle after the end of the regular array");
    }
  }
#endif
}

/**
 * @brief Re-arrange the #sink in a top-level cell such that all the extra
 * ones for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param sinks_offset The offset between the first #sink in the array and
 * the first #sink in the global array in the space structure (for
 * re-linking).
 */
void cell_reorder_extra_sinks(struct cell *c, const ptrdiff_t sinks_offset) {
  struct sink *sinks = c->sinks.parts;
  const int count_real = c->sinks.count;

  if (c->depth != 0 || c->nodeID != engine_rank)
    error("This function should only be called on local top-level cells!");

  int first_not_extra = count_real;

  /* Find extra particles */
  for (int i = 0; i < count_real; ++i) {
    if (sinks[i].time_bin == time_bin_not_created) {
      /* Find the first non-extra particle after the end of the
         real particles */
      while (sinks[first_not_extra].time_bin == time_bin_not_created) {
        ++first_not_extra;
      }

#ifdef SWIFT_DEBUG_CHECKS
      if (first_not_extra >= count_real + space_extra_sinks)
        error("Looking for extra particles beyond this cell's range!");
#endif

      /* Swap everything, including g-part pointer */
      memswap(&sinks[i], &sinks[first_not_extra], sizeof(struct sink));
      if (sinks[i].gpart)
        sinks[i].gpart->id_or_neg_offset = -(i + sinks_offset);
      sinks[first_not_extra].gpart = NULL;
#ifdef SWIFT_DEBUG_CHECKS
      if (sinks[first_not_extra].time_bin != time_bin_not_created)
        error("Incorrect swap occured!");
#endif
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < c->sinks.count_total; ++i) {
    if (sinks[i].time_bin == time_bin_not_created && i < c->sinks.count) {
      error("Extra particle before the end of the regular array");
    }
    if (sinks[i].time_bin != time_bin_not_created && i >= c->sinks.count) {
      error("Regular particle after the end of the regular array");
    }
  }
#endif
}

/**
 * @brief Re-arrange the #gpart in a top-level cell such that all the extra
 * ones for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param parts The global array of #part (for re-linking).
 * @param sparts The global array of #spart (for re-linking).
 * @param sinks The global array of #sink (for re-linking).
 */
void cell_reorder_extra_gparts(struct cell *c, struct part *parts,
                               struct spart *sparts, struct sink *sinks,
                               struct bpart *bparts) {
  struct gpart *gparts = c->grav.parts;
  const int count_real = c->grav.count;

  if (c->depth != 0 || c->nodeID != engine_rank)
    error("This function should only be called on local top-level cells!");

  int first_not_extra = count_real;

  /* Find extra particles */
  for (int i = 0; i < count_real; ++i) {
    if (gparts[i].time_bin == time_bin_not_created) {
      /* Find the first non-extra particle after the end of the
         real particles */
      while (gparts[first_not_extra].time_bin == time_bin_not_created) {
        ++first_not_extra;
      }

#ifdef SWIFT_DEBUG_CHECKS
      if (first_not_extra >= count_real + space_extra_gparts)
        error("Looking for extra particles beyond this cell's range!");
#endif

      /* Swap everything (including pointers) */
      memswap_unaligned(&gparts[i], &gparts[first_not_extra],
                        sizeof(struct gpart));
      if (gparts[i].type == swift_type_gas) {
        parts[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      } else if (gparts[i].type == swift_type_sink) {
        sinks[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      } else if (gparts[i].type == swift_type_stars) {
        sparts[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      } else if (gparts[i].type == swift_type_black_hole) {
        bparts[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < c->grav.count_total; ++i) {
    if (gparts[i].time_bin == time_bin_not_created && i < c->grav.count) {
      error("Extra particle before the end of the regular array");
    }
    if (gparts[i].time_bin != time_bin_not_created && i >= c->grav.count) {
      error("Regular particle after the end of the regular array");
    }
  }
#endif
}
