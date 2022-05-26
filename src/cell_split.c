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
#include "../config.h"

/* This object's header. */
#include "active.h"
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "hilbert.h"
#include "multipole.h"
#include "star_formation_logger.h"

/* Local headers. */
#include "memswap.h"

/**
 * @brief Compute hilbert indices for particles and sort them.
 *        Then split cells into octants recursively.
 *
 * @param c The #cell array to be sorted.
 */
void cell_split(struct cell *c) {

  const int count = c->hydro.count, gcount = c->grav.count,
    scount = c->stars.count, bcount = c->black_holes.count,
    sink_count = c->sinks.count;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct sink *sinks = c->sinks.parts;
  const int depth = c->depth;

  /* Set up buckets for the progeny */
  int bucket_count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int bucket_offset[9];

  /* Lets split the parts */
  /* Fill the buckets using the hilbert keys */
  for (int k = 0; k < count; k++) {
    unsigned long key = parts[k].hilb_key;

    /* Shift bits to the correct depth and mask to get the final 3
       bits which are the bin at this depth */
    int cell_ind = (key << (depth * 3)) & 7;
    bucket_count[cell_ind]++;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->hydro.count = bucket_count[k];
    c->progeny[k]->hydro.count_total = c->progeny[k]->hydro.count;
    c->progeny[k]->hydro.parts = &c->hydro.parts[bucket_offset[k]];
    c->progeny[k]->hydro.xparts = &c->hydro.xparts[bucket_offset[k]];
  }

  /* Now do the same song and dance for the sparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buckets using the hilbert keys */
  for (int k = 0; k < scount; k++) {
    unsigned long key = sparts[k].hilb_key;

    /* Shift bits to the correct depth and mask to get the final 3
       bits which are the bin at this depth */
    int cell_ind = (key << (depth * 3)) & 7;
    bucket_count[cell_ind]++;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
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

   /* Fill the buckets using the hilbert keys */
  for (int k = 0; k < bcount; k++) {
    unsigned long key = bparts[k].hilb_key;

    /* Shift bits to the correct depth and mask to get the final 3
       bits which are the bin at this depth */
    int cell_ind = (key << (depth * 3)) & 7;
    bucket_count[cell_ind]++;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->black_holes.count = bucket_count[k];
    c->progeny[k]->black_holes.count_total = c->progeny[k]->black_holes.count;
    c->progeny[k]->black_holes.parts = &c->black_holes.parts[bucket_offset[k]];
  }

  /* Now do the same song and dance for the sinks. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buckets using the hilbert keys */
  for (int k = 0; k < sink_count; k++) {
    unsigned long key = sinks[k].hilb_key;

    /* Shift bits to the correct depth and mask to get the final 3
       bits which are the bin at this depth */
    int cell_ind = (key << (depth * 3)) & 7;
    bucket_count[cell_ind]++;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->sinks.count = bucket_count[k];
    c->progeny[k]->sinks.count_total = c->progeny[k]->sinks.count;
    c->progeny[k]->sinks.parts = &c->sinks.parts[bucket_offset[k]];
  }

  /* Finally, do the same song and dance for the gparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buckets using the hilbert keys */
  for (int k = 0; k < gcount; k++) {
    unsigned long key = gparts[k].hilb_key;

    /* Shift bits to the correct depth and mask to get the final 3
       bits which are the bin at this depth */
    int cell_ind = key & 7;
    for (int order = 1; order < depth; order ++) {
      cell_ind += (key << (depth * 3)) & 7;
    }
    bucket_count[cell_ind % 8]++;
  }

  /* Set the buffer offsets. */
  bucket_offset[0] = 0;
  for (int k = 1; k <= 8; k++) {
    bucket_offset[k] = bucket_offset[k - 1] + bucket_count[k - 1];
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->grav.count = bucket_count[k];
    c->progeny[k]->grav.count_total = c->progeny[k]->grav.count;
    c->progeny[k]->grav.parts = &c->grav.parts[bucket_offset[k]];
    c->progeny[k]->grav.parts_rebuild = c->progeny[k]->grav.parts;
  }

#ifdef SWIFT_DEBUG_CHECKS

  /* Get the pivot to validate we have sorted everything correctly */
  const double pivot[3] = {c->loc[0] + c->width[0] / 2,
                           c->loc[1] + c->width[1] / 2,
                           c->loc[2] + c->width[2] / 2};

  /* Verify that _all_ the parts have been assigned to a cell. */
  for (int k = 1; k < 8; k++)
    if (&c->progeny[k - 1]->grav.parts[c->progeny[k - 1]->grav.count] !=
        c->progeny[k]->grav.parts)
      error("Particle sorting failed (internal consistency).");
  if (c->progeny[0]->grav.parts != c->grav.parts)
    error("Particle sorting failed (left edge).");
  if (&c->progeny[7]->grav.parts[c->progeny[7]->grav.count] !=
      &c->grav.parts[gcount])
    error("Particle sorting failed (right edge).");

  /* Verify a few sub-cells. */
  for (int k = 0; k < c->progeny[0]->grav.count; k++)
    if (c->progeny[0]->grav.parts[k].x[0] >= pivot[0] ||
        c->progeny[0]->grav.parts[k].x[1] >= pivot[1] ||
        c->progeny[0]->grav.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=0: depth=%d, (grav.parts[%d].x - pivot)=[%e, %e, %e]).",
            c->depth, k,
            c->progeny[0]->grav.parts[k].x[0] - pivot[0],
            c->progeny[0]->grav.parts[k].x[1] - pivot[1],
            c->progeny[0]->grav.parts[k].x[2] - pivot[2]);
  for (int k = 0; k < c->progeny[1]->grav.count; k++)
    if (c->progeny[1]->grav.parts[k].x[0] >= pivot[0] ||
        c->progeny[1]->grav.parts[k].x[1] < pivot[1] ||
        c->progeny[1]->grav.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=1).");
  for (int k = 0; k < c->progeny[2]->grav.count; k++)
    if (c->progeny[2]->grav.parts[k].x[0] >= pivot[0] ||
        c->progeny[2]->grav.parts[k].x[1] < pivot[1] ||
        c->progeny[2]->grav.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=2).");
  for (int k = 0; k < c->progeny[3]->grav.count; k++)
    if (c->progeny[3]->grav.parts[k].x[0] >= pivot[0] ||
        c->progeny[3]->grav.parts[k].x[1] >= pivot[1] ||
        c->progeny[3]->grav.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=3).");
  for (int k = 0; k < c->progeny[4]->grav.count; k++)
    if (c->progeny[4]->grav.parts[k].x[0] < pivot[0] ||
        c->progeny[4]->grav.parts[k].x[1] >= pivot[1] ||
        c->progeny[4]->grav.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=4).");
  for (int k = 0; k < c->progeny[5]->grav.count; k++)
    if (c->progeny[5]->grav.parts[k].x[0] < pivot[0] ||
        c->progeny[5]->grav.parts[k].x[1] < pivot[1] ||
        c->progeny[5]->grav.parts[k].x[2] < pivot[2])
      error("Sorting failed (progeny=5).");
  for (int k = 0; k < c->progeny[6]->grav.count; k++)
    if (c->progeny[6]->grav.parts[k].x[0] < pivot[0] ||
        c->progeny[6]->grav.parts[k].x[1] < pivot[1] ||
        c->progeny[6]->grav.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=6).");
  for (int k = 0; k < c->progeny[7]->grav.count; k++)
    if (c->progeny[7]->grav.parts[k].x[0] < pivot[0] ||
        c->progeny[7]->grav.parts[k].x[1] >= pivot[1] ||
        c->progeny[7]->grav.parts[k].x[2] >= pivot[2])
      error("Sorting failed (progeny=7).");
#endif
  
}


/**
 * @brief Compute hilbert indices for particles and sort them.
 *        Then split cells into octants recursively.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell array to be sorted.
 */
void cell_split_recursive(struct space *s, struct cell *c) {

  /* Extract cell data */
  const int count = c->hydro.count, gcount = c->grav.count,
            scount = c->stars.count;
  const int with_self_gravity = s->with_self_gravity;

  /* Split or let it be? */
  /* Note this is currently limited to a depth of 21! */
  if (((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize)))
      && (c->depth + 1 < 21)){

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->hydro.count = 0;
      cp->grav.count = 0;
      cp->stars.count = 0;
      cp->sinks.count = 0;
      cp->black_holes.count = 0;
      cp->hydro.count_total = 0;
      cp->grav.count_total = 0;
      cp->sinks.count_total = 0;
      cp->stars.count_total = 0;
      cp->black_holes.count_total = 0;
      cp->hydro.ti_old_part = c->hydro.ti_old_part;
      cp->grav.ti_old_part = c->grav.ti_old_part;
      cp->grav.ti_old_multipole = c->grav.ti_old_multipole;
      cp->stars.ti_old_part = c->stars.ti_old_part;
      cp->sinks.ti_old_part = c->sinks.ti_old_part;
      cp->black_holes.ti_old_part = c->black_holes.ti_old_part;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->width[0] = c->width[0] / 2;
      cp->width[1] = c->width[1] / 2;
      cp->width[2] = c->width[2] / 2;
      cp->dmin = c->dmin / 2;
      if (k & 4) cp->loc[0] += cp->width[0];
      if (k == 1 || k == 2 || k == 5 || k == 6) cp->loc[1] += cp->width[1];
      if (k > 1 && k < 6) cp->loc[2] += cp->width[2];
      cp->depth = c->depth + 1;
      cp->split = 0;
      cp->hydro.h_max = 0.f;
      cp->hydro.h_max_active = 0.f;
      cp->hydro.dx_max_part = 0.f;
      cp->hydro.dx_max_sort = 0.f;
      cp->stars.h_max = 0.f;
      cp->stars.h_max_active = 0.f;
      cp->stars.dx_max_part = 0.f;
      cp->stars.dx_max_sort = 0.f;
      cp->sinks.r_cut_max = 0.f;
      cp->sinks.r_cut_max_active = 0.f;
      cp->sinks.dx_max_part = 0.f;
      cp->black_holes.h_max = 0.f;
      cp->black_holes.h_max_active = 0.f;
      cp->black_holes.dx_max_part = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->top = c->top;
      cp->super = NULL;
      cp->hydro.super = NULL;
      cp->grav.super = NULL;
      cp->flags = 0;
      star_formation_logger_init(&cp->stars.sfh);
#ifdef WITH_MPI
      cp->mpi.tag = -1;
#endif  // WITH_MPI
      cp->tl_cell_type = c->tl_cell_type;
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
      cell_assign_cell_index(cp, c);
#endif
    }

    /* Split the cell's particle data. */
    cell_split(c);

    for (int k = 0; k < 8; k++) {

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->hydro.count == 0 && cp->grav.count == 0
          && cp->stars.count == 0 && cp->black_holes.count == 0
          && cp->sinks.count == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {  /* Split again */

        /* Recurse */
        cell_split_recursive(s, cp);

      }
    }
  }   /* Split or let it be? */

  /* Otherwise, clean up progeny. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
  }
}

/**
 * @brief Compute hilbert indices for particles and sort them.
 *        Then split cells into octants recursively.
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
void cell_sort_and_split(struct space *s, struct cell *c,
                         const ptrdiff_t parts_offset,
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
  const double cell_loc[3] = {c->loc[0],
                              c->loc[1],
                              c->loc[2]};
  const double cell_width[3] = {c->width[0],
                                c->width[1],
                                c->width[2]};
  const int nbits = 21;

  /* Get hilbert keys for parts and sort them. */
  if (count > 0) {

     /* Set up keys and indices to sort the keys */
    unsigned long *part_keys =
      (unsigned long *)malloc(count * sizeof(unsigned long));
    int *part_sinds = (int *)malloc(count * sizeof(int));

    /* Loop over particles */
    for (int k = 0; k < count; k++) {

      /* Convert position to 21 bit integer coordinates */
      unsigned long bits[3];
      bits[0] = (1ul << (nbits - 1))
        * (parts[k].x[0] - cell_loc[0]) / cell_width[0];
      bits[1] = (1ul << (nbits - 1))
        * (parts[k].x[1] - cell_loc[1]) / cell_width[1];
      bits[2] = (1ul << (nbits - 1))
        * (parts[k].x[2] - cell_loc[2]) / cell_width[2];

      /* Get hilbert key */
      part_keys[k] = hilbert_get_key_3d(bits, nbits);
      parts[k].hilb_key = part_keys[k];

      /* Set index */
      part_sinds[k] = k;
    }

    /* Now we can sort */
    qsort_r(part_sinds, count, sizeof(int), sort_h_comp,
            part_keys);

    /* Finally, loop over the particles swapping particles to the
       correct place */
    for (int k = 0; k < count; k++) {
      struct part part = parts[k];
      struct xpart xpart = xparts[k];
      memswap(&parts[part_sinds[k]], &part, sizeof(struct part));
      memswap(&xparts[part_sinds[k]], &xpart, sizeof(struct xpart));
      if (parts[k].gpart)
        parts[k].gpart->id_or_neg_offset = -(k + parts_offset);
    }
    
    /* Set the memory free */
    free(part_sinds);
    free(part_keys);
  }

  /* Get hilbert keys for sparts and sort them. */
  if (scount > 0) {

     /* Set up keys and indices to sort the keys */
    unsigned long *spart_keys =
      (unsigned long *)malloc(scount * sizeof(unsigned long));
    int *spart_sinds = (int *)malloc(scount * sizeof(int));

    /* Loop over particles */
    for (int k = 0; k < scount; k++) {

      /* Convert position to 21 bit integer coordinates */
      unsigned long bits[3];
      bits[0] = (1ul << (nbits - 1))
        * (sparts[k].x[0] - cell_loc[0]) / cell_width[0];
      bits[1] = (1ul << (nbits - 1))
        * (sparts[k].x[1] - cell_loc[1]) / cell_width[1];
      bits[2] = (1ul << (nbits - 1))
        * (sparts[k].x[2] - cell_loc[2]) / cell_width[2];

      /* Get hilbert key */
      spart_keys[k] = hilbert_get_key_3d(bits, nbits);
      sparts[k].hilb_key = spart_keys[k];

      /* Set index */
      spart_sinds[k] = k;
    }

    /* Now we can sort */
    qsort_r(spart_sinds, scount, sizeof(int), sort_h_comp,
            spart_keys);

    /* Finally, loop over the particles swapping particles to the
       correct place */
    for (int k = 0; k < scount; k++) {
      struct spart spart = sparts[k];
      memswap(&sparts[spart_sinds[k]], &spart, sizeof(struct spart));
      if (sparts[k].gpart)
          sparts[k].gpart->id_or_neg_offset = -(k + sparts_offset);

    }

    /* Set the memory free */
    free(spart_sinds);
    free(spart_keys);
  }

  /* Get hilbert keys for bparts and sort them. */
  if (bcount > 0) {

     /* Set up keys and indices to sort the keys */
    unsigned long *bpart_keys =
      (unsigned long *)malloc(bcount * sizeof(unsigned long));
    int *bpart_sinds = (int *)malloc(bcount * sizeof(int));

    /* Loop over particles */
    for (int k = 0; k < bcount; k++) {

      /* Convert position to 21 bit integer coordinates */
      unsigned long bits[3];
      bits[0] = (1ul << (nbits - 1))
        * (bparts[k].x[0] - cell_loc[0]) / cell_width[0];
      bits[1] = (1ul << (nbits - 1))
        * (bparts[k].x[1] - cell_loc[1]) / cell_width[1];
      bits[2] = (1ul << (nbits - 1))
        * (bparts[k].x[2] - cell_loc[2]) / cell_width[2];

      /* Get hilbert key */
      bpart_keys[k] = hilbert_get_key_3d(bits, nbits);
      bparts[k].hilb_key = bpart_keys[k];

      /* Set index */
      bpart_sinds[k] = k;
    }

    /* Now we can sort */
    qsort_r(bpart_sinds, bcount, sizeof(int), sort_h_comp,
            bpart_keys);

    /* Finally, loop over the particles swapping particles to the
       correct place */
    for (int k = 0; k < bcount; k++) {
      struct bpart bpart = bparts[k];
      memswap(&bparts[bpart_sinds[k]], &bpart, sizeof(struct bpart));
      if (bparts[k].gpart)
          bparts[k].gpart->id_or_neg_offset = -(k + bparts_offset);

    }

    /* Set the memory free */
    free(bpart_sinds);
    free(bpart_keys);
  }

  /* Get hilbert keys for sinks and sort them. */
  if (sink_count > 0) {

     /* Set up keys and indices to sort the keys */
    unsigned long *sink_keys =
      (unsigned long *)malloc(sink_count * sizeof(unsigned long));
    int *sink_sinds = (int *)malloc(sink_count * sizeof(int));

    /* Loop over particles */
    for (int k = 0; k < sink_count; k++) {

      /* Convert position to 21 bit integer coordinates */
      unsigned long bits[3];
      bits[0] = (1ul << (nbits - 1))
        * (sinks[k].x[0] - cell_loc[0]) / cell_width[0];
      bits[1] = (1ul << (nbits - 1))
        * (sinks[k].x[1] - cell_loc[1]) / cell_width[1];
      bits[2] = (1ul << (nbits - 1))
        * (sinks[k].x[2] - cell_loc[2]) / cell_width[2];

      /* Get hilbert key */
      sink_keys[k] = hilbert_get_key_3d(bits, nbits);
      sinks[k].hilb_key = sink_keys[k];

      /* Set index */
      sink_sinds[k] = k;
    }

    /* Now we can sort */
    qsort_r(sink_sinds, sink_count, sizeof(int), sort_h_comp,
            sink_keys);

    /* Finally, loop over the particles swapping particles to the
       correct place */
    for (int k = 0; k < sink_count; k++) {
      struct sink sink = sinks[k];
      memswap(&sinks[sink_sinds[k]], &sink, sizeof(struct part));
      if (sinks[k].gpart)
        sinks[k].gpart->id_or_neg_offset = -(k + sinks_offset);

    }

    /* Set the memory free */
    free(sink_sinds);
    free(sink_keys);
  }

  /* Get hilbert keys for gparts and sort them. */
  if (gcount > 0) {

     /* Set up keys and indices to sort the keys */
    unsigned long *gpart_keys =
      (unsigned long *)malloc(gcount * sizeof(unsigned long));
    int *gpart_sinds = (int *)malloc(gcount * sizeof(int));

    /* Loop over particles */
    for (int k = 0; k < gcount; k++) {

      /* Convert position to 21 bit integer coordinates */
      unsigned long bits[3];
      bits[0] = (1ul << (nbits - 1))
        * (gparts[k].x[0] - cell_loc[0]) / cell_width[0];
      bits[1] = (1ul << (nbits - 1))
        * (gparts[k].x[1] - cell_loc[1]) / cell_width[1];
      bits[2] = (1ul << (nbits - 1))
        * (gparts[k].x[2] - cell_loc[2]) / cell_width[2];

      /* Get hilbert key */
      gpart_keys[k] = hilbert_get_key_3d(bits, nbits);
      gparts[k].hilb_key = gpart_keys[k];

      /* Set index */
      gpart_sinds[k] = k;
    }
    
    /* Now we can sort */
    qsort_r(gpart_sinds, gcount, sizeof(int), sort_h_comp,
            gpart_keys);
    
    /* Finally, loop over the particles swapping particles to the
       correct place */
    int j, k, sind;
    struct gpart temp_gpart;
    for (k = 0; k < gcount; k++) {

      /* Get the sorted index and swap particles if necessary. */
      if (k != gpart_sinds[k]) {

        /* Set up tempary variables */
        temp_gpart = gparts[k];
        j = k;

        /* Loop until particles are in the right place. */
        while (k != (sind = gpart_sinds[j])) {

          /* Swap particles in memory */
          memswap_unaligned(&gparts[j], &gparts[sind],
                            sizeof(struct gpart));

          /* Corrected the now sorted sind */
          gpart_sinds[j] = j;

          /* Move on to the next */
          j = sind;

          /* Swap sorting indices. */
          sind = gpart_sinds[sind];
        }

        /* Return the temporary particle and set index. */
        gparts[j] = temp_gpart;
        gpart_sinds[j] = j;
      }

      /* Make sure all hydro particles are pointing to the correct gpart. */
      if (gparts[k].type == swift_type_gas) {
        parts[-gparts[k].id_or_neg_offset - parts_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_stars) {
        sparts[-gparts[k].id_or_neg_offset - sparts_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_sink) {
        sinks[-gparts[k].id_or_neg_offset - sinks_offset].gpart = &gparts[k];
      } else if (gparts[k].type == swift_type_black_hole) {
        bparts[-gparts[k].id_or_neg_offset - bparts_offset].gpart = &gparts[k];
      }
    }

    /* Set the memory free */
    free(gpart_sinds);
    free(gpart_keys);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Lets make sure everything is sorted correctly with sensible keys */
  for (int k = 1; k < gcount; k++) {
    if (gparts[k].hilb_key == 0)
      error("Particle has improper key, it may not have been set!");
    if (gparts[k - 1].hilb_key > gparts[k].hilb_key)
      error("Sorting failed, keys are not in order! "
            "(hilb_key[k-1]=%lu, hilb_key[k]=%lu)",
            gparts[k - 1].hilb_key, gparts[k].hilb_key);
  }
#endif

  /* With all that done we are finally in a position to split the cells! */
  cell_split_recursive(s, c);

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
                               struct spart *sparts, struct sink *sinks) {
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

/**
 * @brief A function to calculate the properties in a cell recursively over
 *        a cell heirarchy.
 *
 * @param s The #space in which the cell lives.
 * @param c The cell to compute properties for.
 */
void cell_props_recursive(struct space *s, struct cell *c) {

  /* Lets extract some important information. */
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  const int with_self_gravity = s->with_self_gravity;
  int maxdepth = 0;
  float h_max = 0.0f;
  float h_max_active = 0.0f;
  float stars_h_max = 0.f;
  float stars_h_max_active = 0.f;
  float black_holes_h_max = 0.f;
  float black_holes_h_max_active = 0.f;
  float sinks_h_max = 0.f;
  float sinks_h_max_active = 0.f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_end_max = 0,
                ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_end_max = 0,
                ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_end_max = 0, ti_black_holes_beg_max = 0;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct sink *sinks = c->sinks.parts;
  struct engine *e = s->e;
  const integertime_t ti_current = e->ti_current;

  /* If the cell has been split recurse a level */
  if (c->split) {
    for (int k = 0; k < 8; k++) {

      /* Skip if this is an empty progeny cell. */
      if (c->progeny[k] == NULL) continue;

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Calculate properties for this progenitor cell */
      cell_props_recursive(s, cp);

      /* Update the cell-wide properties */
      h_max = max(h_max, cp->hydro.h_max);
      h_max_active = max(h_max_active, cp->hydro.h_max_active);
      stars_h_max = max(stars_h_max, cp->stars.h_max);
      stars_h_max_active = max(stars_h_max_active, cp->stars.h_max_active);
      black_holes_h_max = max(black_holes_h_max, cp->black_holes.h_max);
      black_holes_h_max_active =
        max(black_holes_h_max_active, cp->black_holes.h_max_active);
      sinks_h_max = max(sinks_h_max, cp->sinks.r_cut_max);
      sinks_h_max_active =
        max(sinks_h_max_active, cp->sinks.r_cut_max_active);

      /* Update time step information. */
      ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
      ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
      ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
      ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);
      ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);
      ti_stars_beg_max = max(ti_stars_beg_max, cp->stars.ti_beg_max);
      ti_sinks_end_min = min(ti_sinks_end_min, cp->sinks.ti_end_min);
      ti_sinks_beg_max = max(ti_sinks_beg_max, cp->sinks.ti_beg_max);
      ti_black_holes_end_min =
        min(ti_black_holes_end_min, cp->black_holes.ti_end_min);
      ti_black_holes_beg_max =
        max(ti_black_holes_beg_max, cp->black_holes.ti_beg_max);

      star_formation_logger_add(&c->stars.sfh, &cp->stars.sfh);

      /* Increase the depth */
      maxdepth = max(maxdepth, cp->maxdepth);
    }

    /* Deal with the gravity */
    if (with_self_gravity) {

      /* Reset everything */
      gravity_reset(c->grav.multipole);

      /* Compute CoM and bulk velocity from all progenies */
      double CoM[3] = {0., 0., 0.};
      double vel[3] = {0., 0., 0.};
      float max_delta_vel[3] = {0.f, 0.f, 0.f};
      float min_delta_vel[3] = {0.f, 0.f, 0.f};
      double mass = 0.;

      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct gravity_tensors *m = c->progeny[k]->grav.multipole;

          mass += m->m_pole.M_000;

          CoM[0] += m->CoM[0] * m->m_pole.M_000;
          CoM[1] += m->CoM[1] * m->m_pole.M_000;
          CoM[2] += m->CoM[2] * m->m_pole.M_000;

          vel[0] += m->m_pole.vel[0] * m->m_pole.M_000;
          vel[1] += m->m_pole.vel[1] * m->m_pole.M_000;
          vel[2] += m->m_pole.vel[2] * m->m_pole.M_000;

          max_delta_vel[0] = max(m->m_pole.max_delta_vel[0], max_delta_vel[0]);
          max_delta_vel[1] = max(m->m_pole.max_delta_vel[1], max_delta_vel[1]);
          max_delta_vel[2] = max(m->m_pole.max_delta_vel[2], max_delta_vel[2]);

          min_delta_vel[0] = min(m->m_pole.min_delta_vel[0], min_delta_vel[0]);
          min_delta_vel[1] = min(m->m_pole.min_delta_vel[1], min_delta_vel[1]);
          min_delta_vel[2] = min(m->m_pole.min_delta_vel[2], min_delta_vel[2]);
        }
      }

      /* Final operation on the CoM and bulk velocity */
      const double inv_mass = 1. / mass;
      c->grav.multipole->CoM[0] = CoM[0] * inv_mass;
      c->grav.multipole->CoM[1] = CoM[1] * inv_mass;
      c->grav.multipole->CoM[2] = CoM[2] * inv_mass;
      c->grav.multipole->m_pole.vel[0] = vel[0] * inv_mass;
      c->grav.multipole->m_pole.vel[1] = vel[1] * inv_mass;
      c->grav.multipole->m_pole.vel[2] = vel[2] * inv_mass;

      /* Min max velocity along each axis */
      c->grav.multipole->m_pole.max_delta_vel[0] = max_delta_vel[0];
      c->grav.multipole->m_pole.max_delta_vel[1] = max_delta_vel[1];
      c->grav.multipole->m_pole.max_delta_vel[2] = max_delta_vel[2];
      c->grav.multipole->m_pole.min_delta_vel[0] = min_delta_vel[0];
      c->grav.multipole->m_pole.min_delta_vel[1] = min_delta_vel[1];
      c->grav.multipole->m_pole.min_delta_vel[2] = min_delta_vel[2];

      /* Now shift progeny multipoles and add them up */
      struct multipole temp;
      double r_max = 0.;
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          const struct cell *cp = c->progeny[k];
          const struct multipole *m = &cp->grav.multipole->m_pole;

          /* Contribution to multipole */
          gravity_M2M(&temp, m, c->grav.multipole->CoM,
                      cp->grav.multipole->CoM);
          gravity_multipole_add(&c->grav.multipole->m_pole, &temp);

          /* Upper limit of max CoM<->gpart distance */
          const double dx =
              c->grav.multipole->CoM[0] - cp->grav.multipole->CoM[0];
          const double dy =
              c->grav.multipole->CoM[1] - cp->grav.multipole->CoM[1];
          const double dz =
              c->grav.multipole->CoM[2] - cp->grav.multipole->CoM[2];
          const double r2 = dx * dx + dy * dy + dz * dz;
          r_max = max(r_max, cp->grav.multipole->r_max + sqrt(r2));
        }
      }

      /* Alternative upper limit of max CoM<->gpart distance */
      const double dx =
          c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
              ? c->grav.multipole->CoM[0] - c->loc[0]
              : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
      const double dy =
          c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
              ? c->grav.multipole->CoM[1] - c->loc[1]
              : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
      const double dz =
          c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
              ? c->grav.multipole->CoM[2] - c->loc[2]
              : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];

      /* Take minimum of both limits */
      c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));

      /* Store the value at rebuild time */
      c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
      c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
      c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
      c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];

      /* Compute the multipole power */
      gravity_multipole_compute_power(&c->grav.multipole->m_pole);

    } /* Deal with gravity */
  } /* Cell is split */
  else {
    maxdepth = c->depth;

    ti_hydro_end_min = max_nr_timesteps;
    ti_hydro_end_max = 0;
    ti_hydro_beg_max = 0;

    ti_gravity_end_min = max_nr_timesteps;
    ti_gravity_end_max = 0;
    ti_gravity_beg_max = 0;

    ti_stars_end_min = max_nr_timesteps;
    ti_stars_end_max = 0;
    ti_stars_beg_max = 0;

    ti_black_holes_end_min = max_nr_timesteps;
    ti_black_holes_end_max = 0;
    ti_black_holes_beg_max = 0;

    /* parts: Get dt_min/dt_max and h_max. */
    for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (parts[k].time_bin == time_bin_not_created)
        error("Extra particle present in space_split()");
      if (parts[k].time_bin == time_bin_inhibited)
        error("Inhibited particle present in space_split()");
#endif

      /* When does this particle's time-step start and end? */
      const timebin_t time_bin = parts[k].time_bin;
      const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
      const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

      ti_hydro_end_min = min(ti_hydro_end_min, ti_end);
      ti_hydro_end_max = max(ti_hydro_end_max, ti_end);
      ti_hydro_beg_max = max(ti_hydro_beg_max, ti_beg);

      h_max = max(h_max, parts[k].h);

      if (part_is_active(&parts[k], e))
        h_max_active = max(h_max_active, parts[k].h);

      /* Collect SFR from the particles after rebuilt */
      star_formation_logger_log_inactive_part(&parts[k], &xparts[k],
                                              &c->stars.sfh);
    }

    /* xparts: Reset x_diff */
    for (int k = 0; k < count; k++) {
      xparts[k].x_diff[0] = 0.f;
      xparts[k].x_diff[1] = 0.f;
      xparts[k].x_diff[2] = 0.f;
    }

    /* gparts: Get dt_min/dt_max. */
    for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (gparts[k].time_bin == time_bin_not_created)
        error("Extra g-particle present in space_split()");
      if (gparts[k].time_bin == time_bin_inhibited)
        error("Inhibited g-particle present in space_split()");
#endif

      /* When does this particle's time-step start and end? */
      const timebin_t time_bin = gparts[k].time_bin;
      const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
      const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

      ti_gravity_end_min = min(ti_gravity_end_min, ti_end);
      ti_gravity_end_max = max(ti_gravity_end_max, ti_end);
      ti_gravity_beg_max = max(ti_gravity_beg_max, ti_beg);
    }

    /* sparts: Get dt_min/dt_max */
    for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[k].time_bin == time_bin_not_created)
        error("Extra s-particle present in space_split()");
      if (sparts[k].time_bin == time_bin_inhibited)
        error("Inhibited s-particle present in space_split()");
#endif

      /* When does this particle's time-step start and end? */
      const timebin_t time_bin = sparts[k].time_bin;
      const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
      const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

      ti_stars_end_min = min(ti_stars_end_min, ti_end);
      ti_stars_end_max = max(ti_stars_end_max, ti_end);
      ti_stars_beg_max = max(ti_stars_beg_max, ti_beg);

      stars_h_max = max(stars_h_max, sparts[k].h);

      if (spart_is_active(&sparts[k], e))
        stars_h_max_active = max(stars_h_max_active, sparts[k].h);

      /* Reset x_diff */
      sparts[k].x_diff[0] = 0.f;
      sparts[k].x_diff[1] = 0.f;
      sparts[k].x_diff[2] = 0.f;
    }

    /* sinks: Get dt_min/dt_max */
    for (int k = 0; k < sink_count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (sinks[k].time_bin == time_bin_not_created)
        error("Extra sink-particle present in space_split()");
      if (sinks[k].time_bin == time_bin_inhibited)
        error("Inhibited sink-particle present in space_split()");
#endif

      /* When does this particle's time-step start and end? */
      const timebin_t time_bin = sinks[k].time_bin;
      const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
      const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

      ti_sinks_end_min = min(ti_sinks_end_min, ti_end);
      ti_sinks_end_max = max(ti_sinks_end_max, ti_end);
      ti_sinks_beg_max = max(ti_sinks_beg_max, ti_beg);

      sinks_h_max = max(sinks_h_max, sinks[k].r_cut);

      if (sink_is_active(&sinks[k], e))
        sinks_h_max_active = max(sinks_h_max_active, sinks[k].r_cut);

      /* Reset x_diff */
      sinks[k].x_diff[0] = 0.f;
      sinks[k].x_diff[1] = 0.f;
      sinks[k].x_diff[2] = 0.f;
    }

    /* bparts: Get dt_min/dt_max */
    for (int k = 0; k < bcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
      if (bparts[k].time_bin == time_bin_not_created)
        error("Extra b-particle present in space_split()");
      if (bparts[k].time_bin == time_bin_inhibited)
        error("Inhibited b-particle present in space_split()");
#endif

      /* When does this particle's time-step start and end? */
      const timebin_t time_bin = bparts[k].time_bin;
      const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
      const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

      ti_black_holes_end_min = min(ti_black_holes_end_min, ti_end);
      ti_black_holes_end_max = max(ti_black_holes_end_max, ti_end);
      ti_black_holes_beg_max = max(ti_black_holes_beg_max, ti_beg);

      black_holes_h_max = max(black_holes_h_max, bparts[k].h);

      if (bpart_is_active(&bparts[k], e))
        black_holes_h_max_active = max(black_holes_h_max_active, bparts[k].h);

      /* Reset x_diff */
      bparts[k].x_diff[0] = 0.f;
      bparts[k].x_diff[1] = 0.f;
      bparts[k].x_diff[2] = 0.f;
    }

    /* Construct the multipole and the centre of mass*/
    if (with_self_gravity) {
      if (gcount > 0) {

        gravity_P2M(c->grav.multipole, c->grav.parts, c->grav.count,
                    e->gravity_properties);

        /* Compute the multipole power */
        gravity_multipole_compute_power(&c->grav.multipole->m_pole);

      } else {

        /* No gparts in that leaf cell */

        /* Set the values to something sensible */
        gravity_multipole_init(&c->grav.multipole->m_pole);
        if (c->nodeID == engine_rank) {
          c->grav.multipole->CoM[0] = c->loc[0] + c->width[0] / 2.;
          c->grav.multipole->CoM[1] = c->loc[1] + c->width[1] / 2.;
          c->grav.multipole->CoM[2] = c->loc[2] + c->width[2] / 2.;
          c->grav.multipole->r_max = 0.;
        }
      }

      /* Store the value at rebuild time */
      c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
      c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
      c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
      c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];
    }
  }

  /* Set the values for this cell. */
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max_active;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->stars.h_max = stars_h_max;
  c->stars.h_max_active = stars_h_max_active;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->sinks.r_cut_max = sinks_h_max;
  c->sinks.r_cut_max_active = sinks_h_max_active;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.h_max = black_holes_h_max;
  c->black_holes.h_max_active = black_holes_h_max_active;
  c->maxdepth = maxdepth;

  /* Set ownership according to the start of the parts array. */
  if (s->nr_parts > 0)
    c->owner = ((c->hydro.parts - s->parts) % s->nr_parts) * s->nr_queues /
               s->nr_parts;
  else if (s->nr_sinks > 0)
    c->owner = ((c->sinks.parts - s->sinks) % s->nr_sinks) * s->nr_queues /
               s->nr_sinks;
  else if (s->nr_sparts > 0)
    c->owner = ((c->stars.parts - s->sparts) % s->nr_sparts) * s->nr_queues /
               s->nr_sparts;
  else if (s->nr_bparts > 0)
    c->owner = ((c->black_holes.parts - s->bparts) % s->nr_bparts) *
               s->nr_queues / s->nr_bparts;
  else if (s->nr_gparts > 0)
    c->owner = ((c->grav.parts - s->gparts) % s->nr_gparts) * s->nr_queues /
               s->nr_gparts;
  else
    c->owner = 0; /* Ok, there is really nothing on this rank... */

  /* Store the global max depth */
  if (c->depth == 0) atomic_max(&s->maxdepth, maxdepth);
}

/**
 * @brief A mapper function to loop over cells calculating properties.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void cell_props_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    cell_props_recursive(s, c);
  }
  
#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
  }
#endif
}

#ifdef WITH_ZOOM_REGION
/**
 * @brief A wrapper for #threadpool mapper function to calcualte background
 * cell properties.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void bkg_cell_props_mapper(void *map_data, int num_cells, void *extra_data) {
  cell_props_mapper(map_data, num_cells, extra_data);
}

/**
 * @brief A wrapper for #threadpool mapper function to calcualte zoom
 * cell properties.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void zoom_cell_props_mapper(void *map_data, int num_cells, void *extra_data) {
  cell_props_mapper(map_data, num_cells, extra_data);
}

#endif
