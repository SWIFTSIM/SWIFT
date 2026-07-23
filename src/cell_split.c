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

/**
 * @brief @brief Sort the parts into eight bins along the given pivots.
 *
 * This only reorders the buffers (buff, sbuff, bbuff, gbuff and sinkbuff), the
 * particle arrays themselves are sorted later on in the process based on the
 * location of the buffers (which contain the indices of the particles they
 * represent).
 *
 * @param c The #cell array to be sorted.
 * @param ind Scratch space with at least
 * max(c->hydro.count, c->grav.count, c->stars.count, c->black_holes.count,
 * c->sinks.count) entries, reused for each family in turn.
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
int cell_split(struct cell *c, int *restrict ind,
               struct cell_buff *restrict buff,
               struct cell_buff *restrict sbuff,
               struct cell_buff *restrict bbuff,
               struct cell_buff *restrict gbuff,
               struct cell_buff *restrict sinkbuff) {

  const int count = c->hydro.count, gcount = c->grav.count,
            scount = c->stars.count, bcount = c->black_holes.count,
            sink_count = c->sinks.count;
  const double pivot[3] = {c->loc[0] + c->width[0] / 2,
                           c->loc[1] + c->width[1] / 2,
                           c->loc[2] + c->width[2] / 2};
  int bucket_count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int bucket_offset[9];

  /* Track whether any particle was actually moved out of its incoming slot.
   * If nothing moves anywhere in a top-level cell's whole subtree, its leaves
   * are already in final order and the gather/write-back can be skipped. */
  int moved = 0;

  /* Fill ind with the bucket indices. */
  for (int k = 0; k < count; k++) {
    const int bid = (buff[k].x[0] >= pivot[0]) * 4 +
                    (buff[k].x[1] >= pivot[1]) * 2 + (buff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    ind[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap buffer entries (and their ind) to
   * their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = ind[k];
      if (bid != bucket) {
        moved = 1;
        struct cell_buff temp_buff = buff[k];
        int temp_ind = ind[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (ind[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&buff[j], &temp_buff, sizeof(struct cell_buff));
          const int swap_ind = ind[j];
          ind[j] = temp_ind;
          temp_ind = swap_ind;
          bid = temp_ind;
        }
        buff[k] = temp_buff;
        ind[k] = temp_ind;
      }
      bucket_count[bid]++;
    }
  }

  /* Store the counts and offsets. The #part/#xpart arrays have not moved
   * yet, so these pointers describe where the particles will end up. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->hydro.count = bucket_count[k];
    c->progeny[k]->hydro.count_total = c->progeny[k]->hydro.count;
    c->progeny[k]->hydro.parts = &c->hydro.parts[bucket_offset[k]];
    c->progeny[k]->hydro.xparts = &c->hydro.xparts[bucket_offset[k]];
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* The buffer should now be grouped into non-decreasing bucket order. We
   * cannot check this against the #part array itself here, since the
   * particles have not been physically moved yet -- only the buffer has. */
  for (int k = 1; k < count; k++) {
    if (ind[k] < ind[k - 1]) error("Buff not sorted.");
  }
  for (int k = 0; k < count; k++) {
    const int bid = (buff[k].x[0] >= pivot[0]) * 4 +
                    (buff[k].x[1] >= pivot[1]) * 2 + (buff[k].x[2] >= pivot[2]);
    if (bid != ind[k]) error("Buff ind inconsistent with position.");
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
#endif

  /* Now do the same song and dance for the sparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill ind with the bucket indices. */
  for (int k = 0; k < scount; k++) {
    const int bid = (sbuff[k].x[0] >= pivot[0]) * 4 +
                    (sbuff[k].x[1] >= pivot[1]) * 2 +
                    (sbuff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    ind[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap buffer entries (and their ind) to
   * their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = ind[k];
      if (bid != bucket) {
        moved = 1;
        struct cell_buff temp_buff = sbuff[k];
        int temp_ind = ind[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (ind[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&sbuff[j], &temp_buff, sizeof(struct cell_buff));
          const int swap_ind = ind[j];
          ind[j] = temp_ind;
          temp_ind = swap_ind;
          bid = temp_ind;
        }
        sbuff[k] = temp_buff;
        ind[k] = temp_ind;
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

  /* Fill ind with the bucket indices. */
  for (int k = 0; k < bcount; k++) {
    const int bid = (bbuff[k].x[0] >= pivot[0]) * 4 +
                    (bbuff[k].x[1] >= pivot[1]) * 2 +
                    (bbuff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    ind[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap buffer entries (and their ind) to
   * their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = ind[k];
      if (bid != bucket) {
        moved = 1;
        struct cell_buff temp_buff = bbuff[k];
        int temp_ind = ind[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (ind[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&bbuff[j], &temp_buff, sizeof(struct cell_buff));
          const int swap_ind = ind[j];
          ind[j] = temp_ind;
          temp_ind = swap_ind;
          bid = temp_ind;
        }
        bbuff[k] = temp_buff;
        ind[k] = temp_ind;
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

  /* Fill ind with the bucket indices. */
  for (int k = 0; k < sink_count; k++) {
    const int bid = (sinkbuff[k].x[0] >= pivot[0]) * 4 +
                    (sinkbuff[k].x[1] >= pivot[1]) * 2 +
                    (sinkbuff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    ind[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap buffer entries (and their ind) to
   * their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = ind[k];
      if (bid != bucket) {
        moved = 1;
        struct cell_buff temp_buff = sinkbuff[k];
        int temp_ind = ind[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (ind[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&sinkbuff[j], &temp_buff, sizeof(struct cell_buff));
          const int swap_ind = ind[j];
          ind[j] = temp_ind;
          temp_ind = swap_ind;
          bid = temp_ind;
        }
        sinkbuff[k] = temp_buff;
        ind[k] = temp_ind;
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

  /* Fill ind with the bucket indices. */
  for (int k = 0; k < gcount; k++) {
    const int bid = (gbuff[k].x[0] >= pivot[0]) * 4 +
                    (gbuff[k].x[1] >= pivot[1]) * 2 +
                    (gbuff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    ind[k] = bid;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
    bucket_count[k - 1] = 0;
  }

  /* Run through the buckets, and swap buffer entries (and their ind) to
   * their correct spot. */
  for (int bucket = 0; bucket < 8; bucket++) {
    for (int k = bucket_offset[bucket] + bucket_count[bucket];
         k < bucket_offset[bucket + 1]; k++) {
      int bid = ind[k];
      if (bid != bucket) {
        moved = 1;
        struct cell_buff temp_buff = gbuff[k];
        int temp_ind = ind[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (ind[j] == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&gbuff[j], &temp_buff, sizeof(struct cell_buff));
          const int swap_ind = ind[j];
          ind[j] = temp_ind;
          temp_ind = swap_ind;
          bid = temp_ind;
        }
        gbuff[k] = temp_buff;
        ind[k] = temp_ind;
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

  return moved;
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
