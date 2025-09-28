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
#include "memswap.h"
#include "part_type.h"

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
 * @param buff A buffer with at least max(c->hydro.count, c->grav.count)
 * entries, used for sorting indices.
 * @param sbuff A buffer with at least max(c->stars.count, c->grav.count)
 * entries, used for sorting indices for the sparts.
 * @param bbuff A buffer with at least max(c->black_holes.count, c->grav.count)
 * entries, used for sorting indices for the sparts.
 * @param gbuff A buffer with at least max(c->hydro.count, c->grav.count)
 * entries, used for sorting indices for the gparts.
 * @param sinkbuff A buffer with at least max(c->sinks.count, c->grav.count)
 * entries, used for sorting indices for the sinks.
 */
void cell_split(struct cell *c, const ptrdiff_t parts_offset,
                const ptrdiff_t sparts_offset, const ptrdiff_t bparts_offset,
                const ptrdiff_t sinks_offset, struct cell_buff *restrict buff,
                struct cell_buff *restrict sbuff,
                struct cell_buff *restrict bbuff,
                struct cell_buff *restrict gbuff,
                struct cell_buff *restrict sinkbuff) {

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

#define HYDRO_BUCKET(x)                                        \
  ((((x)[0] >= pivot[0]) << 2) | (((x)[1] >= pivot[1]) << 1) | \
   ((x)[2] >= pivot[2]))
#define STRICT_BUCKET(x)                                     \
  ((((x)[0] > pivot[0]) << 2) | (((x)[1] > pivot[1]) << 1) | \
   ((x)[2] > pivot[2]))

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k < count; k++) {
    if (buff[k].x[0] != parts[k].x[0] || buff[k].x[1] != parts[k].x[1] ||
        buff[k].x[2] != parts[k].x[2])
      error("Inconsistent buff contents.");
  }
  for (int k = 0; k < gcount; k++) {
    if (gbuff[k].x[0] != gparts[k].x[0] || gbuff[k].x[1] != gparts[k].x[1] ||
        gbuff[k].x[2] != gparts[k].x[2])
      error("Inconsistent gbuff contents.");
  }
  for (int k = 0; k < scount; k++) {
    if (sbuff[k].x[0] != sparts[k].x[0] || sbuff[k].x[1] != sparts[k].x[1] ||
        sbuff[k].x[2] != sparts[k].x[2])
      error("Inconsistent sbuff contents.");
  }
  for (int k = 0; k < bcount; k++) {
    if (bbuff[k].x[0] != bparts[k].x[0] || bbuff[k].x[1] != bparts[k].x[1] ||
        bbuff[k].x[2] != bparts[k].x[2])
      error("Inconsistent bbuff contents.");
  }
  for (int k = 0; k < sink_count; k++) {
    if (sinkbuff[k].x[0] != sinks[k].x[0] ||
        sinkbuff[k].x[1] != sinks[k].x[1] || sinkbuff[k].x[2] != sinks[k].x[2])
      error("Inconsistent sinkbuff contents.");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Hydro */
  for (int k = 0; k < count; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
    const int bid = HYDRO_BUCKET(buff[k].x);
    buff[k].ind = bid;
#else
    const int bid = HYDRO_BUCKET(parts[k].x);
#endif
    bucket_count[bid]++;
  }

  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {

#if defined(SWIFT_DEBUG_CHECKS)
      int bid = buff[k].ind;
#else
      int bid = HYDRO_BUCKET(parts[k].x);
#endif
      if (bid != bucket) {
        struct part part = parts[k];
        struct xpart xpart = xparts[k];
#if defined(SWIFT_DEBUG_CHECKS)
        struct cell_buff temp_buff = buff[k];
#endif
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (1) {
#if defined(SWIFT_DEBUG_CHECKS)
            if (buff[j].ind != bid) break;
#else
            if (HYDRO_BUCKET(parts[j].x) != bid) break;
#endif
            j++;
            bucket_count[bid]++;
          }
          memswap(&parts[j], &part, sizeof(struct part));
          memswap(&xparts[j], &xpart, sizeof(struct xpart));
#if defined(SWIFT_DEBUG_CHECKS)
          memswap(&buff[j], &temp_buff, sizeof(struct cell_buff));
#endif
          if (parts[j].gpart)
            parts[j].gpart->id_or_neg_offset = -(j + parts_offset);
#if defined(SWIFT_DEBUG_CHECKS)
          bid = temp_buff.ind;
#else
          bid = HYDRO_BUCKET(part.x);
#endif
        }
        parts[k] = part;
        xparts[k] = xpart;
#if defined(SWIFT_DEBUG_CHECKS)
        buff[k] = temp_buff;
#endif
        if (parts[k].gpart)
          parts[k].gpart->id_or_neg_offset = -(k + parts_offset);
      }
      bucket_count[bid]++;
    }
  }

  for (int k = 0; k < 8; k++) {
    c->progeny[k]->hydro.count = bucket_count[k];
    c->progeny[k]->hydro.count_total = c->progeny[k]->hydro.count;
    c->progeny[k]->hydro.parts = &c->hydro.parts[bucket_offset[k]];
    c->progeny[k]->hydro.xparts = &c->hydro.xparts[bucket_offset[k]];
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 1; k < count; k++) {
    if (buff[k].ind < buff[k - 1].ind) error("Buff not sorted.");
    if (buff[k].x[0] != parts[k].x[0] || buff[k].x[1] != parts[k].x[1] ||
        buff[k].x[2] != parts[k].x[2])
      error("Inconsistent buff contents (k=%i).", k);
  }
  for (int k = 1; k < 8; k++)
    if (&c->progeny[k - 1]->hydro.parts[c->progeny[k - 1]->hydro.count] !=
        c->progeny[k]->hydro.parts)
      error("Particle sorting failed (internal consistency).");
  if (c->progeny[0]->hydro.parts != c->hydro.parts)
    error("Particle sorting failed (left edge).");
  if (&c->progeny[7]->hydro.parts[c->progeny[7]->hydro.count] !=
      &c->hydro.parts[count])
    error("Particle sorting failed (right edge).");

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

  /* Stars */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;
  for (int k = 0; k < scount; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
    const int bid = STRICT_BUCKET(sbuff[k].x);
    sbuff[k].ind = bid;
#else
    const int bid = STRICT_BUCKET(sparts[k].x);
#endif
    bucket_count[bid]++;
  }

  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
      int bid = sbuff[k].ind;
#else
      int bid = STRICT_BUCKET(sparts[k].x);
#endif
      if (bid != bucket) {
        struct spart spart = sparts[k];
#if defined(SWIFT_DEBUG_CHECKS)
        struct cell_buff temp_buff = sbuff[k];
#endif
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (1) {
#if defined(SWIFT_DEBUG_CHECKS)
            if (sbuff[j].ind != bid) break;
#else
            if (STRICT_BUCKET(sparts[j].x) != bid) break;
#endif
            j++;
            bucket_count[bid]++;
          }
          memswap(&sparts[j], &spart, sizeof(struct spart));
#if defined(SWIFT_DEBUG_CHECKS)
          memswap(&sbuff[j], &temp_buff, sizeof(struct cell_buff));
#endif
          if (sparts[j].gpart)
            sparts[j].gpart->id_or_neg_offset = -(j + sparts_offset);
#if defined(SWIFT_DEBUG_CHECKS)
          bid = temp_buff.ind;
#else
          bid = STRICT_BUCKET(spart.x);
#endif
        }
        sparts[k] = spart;
#if defined(SWIFT_DEBUG_CHECKS)
        sbuff[k] = temp_buff;
#endif
        if (sparts[k].gpart)
          sparts[k].gpart->id_or_neg_offset = -(k + sparts_offset);
      }
      bucket_count[bid]++;
    }
  }

  for (int k = 0; k < 8; k++) {
    c->progeny[k]->stars.count = bucket_count[k];
    c->progeny[k]->stars.count_total = c->progeny[k]->stars.count;
    c->progeny[k]->stars.parts = &c->stars.parts[bucket_offset[k]];
    c->progeny[k]->stars.parts_rebuild = c->progeny[k]->stars.parts;
  }

  /* Black holes */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;
  for (int k = 0; k < bcount; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
    const int bid = STRICT_BUCKET(bbuff[k].x);
    bbuff[k].ind = bid;
#else
    const int bid = STRICT_BUCKET(bparts[k].x);
#endif
    bucket_count[bid]++;
  }

  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
      int bid = bbuff[k].ind;
#else
      int bid = STRICT_BUCKET(bparts[k].x);
#endif
      if (bid != bucket) {
        struct bpart bpart = bparts[k];
#if defined(SWIFT_DEBUG_CHECKS)
        struct cell_buff temp_buff = bbuff[k];
#endif
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (1) {
#if defined(SWIFT_DEBUG_CHECKS)
            if (bbuff[j].ind != bid) break;
#else
            if (STRICT_BUCKET(bparts[j].x) != bid) break;
#endif
            j++;
            bucket_count[bid]++;
          }
          memswap(&bparts[j], &bpart, sizeof(struct bpart));
#if defined(SWIFT_DEBUG_CHECKS)
          memswap(&bbuff[j], &temp_buff, sizeof(struct cell_buff));
#endif
          if (bparts[j].gpart)
            bparts[j].gpart->id_or_neg_offset = -(j + bparts_offset);
#if defined(SWIFT_DEBUG_CHECKS)
          bid = temp_buff.ind;
#else
          bid = STRICT_BUCKET(bpart.x);
#endif
        }
        bparts[k] = bpart;
#if defined(SWIFT_DEBUG_CHECKS)
        bbuff[k] = temp_buff;
#endif
        if (bparts[k].gpart)
          bparts[k].gpart->id_or_neg_offset = -(k + bparts_offset);
      }
      bucket_count[bid]++;
    }
  }

  for (int k = 0; k < 8; k++) {
    c->progeny[k]->black_holes.count = bucket_count[k];
    c->progeny[k]->black_holes.count_total = c->progeny[k]->black_holes.count;
    c->progeny[k]->black_holes.parts = &c->black_holes.parts[bucket_offset[k]];
  }

  /* Sinks */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;
  for (int k = 0; k < sink_count; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
    const int bid = STRICT_BUCKET(sinkbuff[k].x);
    sinkbuff[k].ind = bid;
#else
    const int bid = STRICT_BUCKET(sinks[k].x);
#endif
    bucket_count[bid]++;
  }

  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
      int bid = sinkbuff[k].ind;
#else
      int bid = STRICT_BUCKET(sinks[k].x);
#endif
      if (bid != bucket) {
        struct sink sink = sinks[k];
#if defined(SWIFT_DEBUG_CHECKS)
        struct cell_buff temp_buff = sinkbuff[k];
#endif
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (1) {
#if defined(SWIFT_DEBUG_CHECKS)
            if (sinkbuff[j].ind != bid) break;
#else
            if (STRICT_BUCKET(sinks[j].x) != bid) break;
#endif
            j++;
            bucket_count[bid]++;
          }
          memswap(&sinks[j], &sink, sizeof(struct sink));
#if defined(SWIFT_DEBUG_CHECKS)
          memswap(&sinkbuff[j], &temp_buff, sizeof(struct cell_buff));
#endif
          if (sinks[j].gpart)
            sinks[j].gpart->id_or_neg_offset = -(j + sinks_offset);
#if defined(SWIFT_DEBUG_CHECKS)
          bid = temp_buff.ind;
#else
          bid = STRICT_BUCKET(sink.x);
#endif
        }
        sinks[k] = sink;
#if defined(SWIFT_DEBUG_CHECKS)
        sinkbuff[k] = temp_buff;
#endif
        if (sinks[k].gpart)
          sinks[k].gpart->id_or_neg_offset = -(k + sinks_offset);
      }
      bucket_count[bid]++;
    }
  }

  for (int k = 0; k < 8; k++) {
    c->progeny[k]->sinks.count = bucket_count[k];
    c->progeny[k]->sinks.count_total = c->progeny[k]->sinks.count;
    c->progeny[k]->sinks.parts = &c->sinks.parts[bucket_offset[k]];
    c->progeny[k]->sinks.parts_rebuild = c->progeny[k]->sinks.parts;
  }

  /* Gravity */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;
  for (int k = 0; k < gcount; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
    const int bid = STRICT_BUCKET(gbuff[k].x);
    gbuff[k].ind = bid;
#else
    const int bid = STRICT_BUCKET(gparts[k].x);
#endif
    bucket_count[bid]++;
  }

  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
#if defined(SWIFT_DEBUG_CHECKS)
      int bid = gbuff[k].ind;
#else
      int bid = STRICT_BUCKET(gparts[k].x);
#endif
      if (bid != bucket) {
        struct gpart gpart = gparts[k];
#if defined(SWIFT_DEBUG_CHECKS)
        struct cell_buff temp_buff = gbuff[k];
#endif
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (1) {
#if defined(SWIFT_DEBUG_CHECKS)
            if (gbuff[j].ind != bid) break;
#else
            if (STRICT_BUCKET(gparts[j].x) != bid) break;
#endif
            j++;
            bucket_count[bid]++;
          }
          memswap_unaligned(&gparts[j], &gpart, sizeof(struct gpart));
#if defined(SWIFT_DEBUG_CHECKS)
          memswap(&gbuff[j], &temp_buff, sizeof(struct cell_buff));
#endif
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
#if defined(SWIFT_DEBUG_CHECKS)
          bid = temp_buff.ind;
#else
          bid = STRICT_BUCKET(gpart.x);
#endif
        }
        gparts[k] = gpart;
#if defined(SWIFT_DEBUG_CHECKS)
        gbuff[k] = temp_buff;
#endif
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

  for (int k = 0; k < 8; k++) {
    c->progeny[k]->grav.count = bucket_count[k];
    c->progeny[k]->grav.count_total = c->progeny[k]->grav.count;
    c->progeny[k]->grav.parts = &c->grav.parts[bucket_offset[k]];
    c->progeny[k]->grav.parts_rebuild = c->progeny[k]->grav.parts;
  }

#undef HYDRO_BUCKET
#undef STRICT_BUCKET
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
        error(
            "Looking for extra particles beyond this cell's range! "
            "(cell_type/cell_subtype=%s/%s, depth=%d, "
            "c->grav.count=%d, first_not_extra=%d, space_extra_gparts=%d, "
            "part_type=%s, i=%d, gparts[i].type=%s)",
            cellID_names[c->type], subcellID_names[c->subtype], c->depth,
            count_real, first_not_extra, space_extra_gparts,
            part_type_names[gparts[first_not_extra].type], i,
            part_type_names[gparts[i].type]);
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
