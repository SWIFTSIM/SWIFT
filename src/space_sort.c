/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* This object's header. */
#include "error.h"
#include "memswap.h"
#include "memuse.h"
#include "space.h"

/**
 * @brief Sort the particles and condensed particles according to the given
 * indices.
 *
 * @param parts The array of #part to sort.
 * @param xparts The corresponding #xpart array to sort as well.
 * @param ind The indices with respect to which the parts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of count).
 * @param parts_offset Offset of the #part array from the global #part array.
 */
void space_parts_sort(struct part *parts, struct xpart *xparts,
                      int *restrict ind, int *restrict counts, int num_bins,
                      ptrdiff_t parts_offset) {
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
 * @brief Sort the s-particles according to the given indices.
 *
 * @param sparts The array of #spart to sort.
 * @param ind The indices with respect to which the #spart are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param sparts_offset Offset of the #spart array from the global #spart.
 * array.
 */
void space_sparts_sort(struct spart *sparts, int *restrict ind,
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
 * @brief Sort the b-particles according to the given indices.
 *
 * @param bparts The array of #bpart to sort.
 * @param ind The indices with respect to which the #bpart are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param bparts_offset Offset of the #bpart array from the global #bpart.
 * array.
 */
void space_bparts_sort(struct bpart *bparts, int *restrict ind,
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
 * @brief Sort the sink-particles according to the given indices.
 *
 * @param sinks The array of #sink to sort.
 * @param ind The indices with respect to which the #sink are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 * @param sinks_offset Offset of the #sink array from the global #sink.
 * array.
 */
void space_sinks_sort(struct sink *sinks, int *restrict ind,
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
 * @brief Sort the g-particles according to the given indices.
 *
 * @param gparts The array of #gpart to sort.
 * @param parts Global #part array for re-linking.
 * @param sinks Global #sink array for re-linking.
 * @param sparts Global #spart array for re-linking.
 * @param bparts Global #bpart array for re-linking.
 * @param ind The indices with respect to which the gparts are sorted.
 * @param counts Number of particles per index.
 * @param num_bins Total number of bins (length of counts).
 */
void space_gparts_sort(struct gpart *gparts, struct part *parts,
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
