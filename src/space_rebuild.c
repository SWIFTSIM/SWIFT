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
#include "space.h"

/* Local headers. */
#include "cell.h"
#include "engine.h"
#include "memswap.h"

/*! Expected maximal number of strays received at a rebuild */
extern int space_expected_max_nr_strays;

/*! Counter for cell IDs (when debugging) */
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
extern long long last_cell_id;
#endif

/**
 * @brief Re-build the top-level cells as well as the whole hierarchy.
 *
 * @param s The #space in which to update the cells.
 * @param repartitioned Did we just repartition?
 * @param verbose Print messages to stdout or not
 */
void space_rebuild(struct space *s, int repartitioned, int verbose) {

  const ticks tic = getticks();

/* Be verbose about this. */
#ifdef SWIFT_DEBUG_CHECKS
  if (s->e->nodeID == 0 || verbose) message("(re)building space");
  fflush(stdout);
#endif
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  /* Reset the cell counter */
  last_cell_id = 1;
#endif

  /* Re-grid if necessary, or just re-set the cell data. */
  space_regrid(s, verbose);

  /* Allocate extra space for particles that will be created */
  if (s->with_star_formation || s->e->policy & engine_policy_sinks)
    space_allocate_extras(s, verbose);

  struct cell *cells_top = s->cells_top;
  const integertime_t ti_current = (s->e != NULL) ? s->e->ti_current : 0;
  const int local_nodeID = s->e->nodeID;

  /* The current number of particles */
  size_t nr_parts = s->nr_parts;
  size_t nr_gparts = s->nr_gparts;
  size_t nr_sparts = s->nr_sparts;
  size_t nr_bparts = s->nr_bparts;
  size_t nr_sinks = s->nr_sinks;

  /* The number of particles we allocated memory for */
  size_t size_parts = s->size_parts;
  size_t size_gparts = s->size_gparts;
  size_t size_sparts = s->size_sparts;
  size_t size_bparts = s->size_bparts;
  size_t size_sinks = s->size_sinks;

  /* Counter for the number of inhibited particles found on the node */
  size_t count_inhibited_parts = 0;
  size_t count_inhibited_gparts = 0;
  size_t count_inhibited_sparts = 0;
  size_t count_inhibited_bparts = 0;
  size_t count_inhibited_sinks = 0;

  /* Counter for the number of extra particles found on the node */
  size_t count_extra_parts = 0;
  size_t count_extra_gparts = 0;
  size_t count_extra_sparts = 0;
  size_t count_extra_bparts = 0;
  size_t count_extra_sinks = 0;

  /* Number of particles we expect to have after strays exchange */
  const size_t h_index_size = size_parts + space_expected_max_nr_strays;
  const size_t g_index_size = size_gparts + space_expected_max_nr_strays;
  const size_t s_index_size = size_sparts + space_expected_max_nr_strays;
  const size_t b_index_size = size_bparts + space_expected_max_nr_strays;
  const size_t sink_index_size = size_sinks + space_expected_max_nr_strays;

  /* Allocate arrays to store the indices of the cells where particles
     belong. We allocate extra space to allow for particles we may
     receive from other nodes */
  int *h_index = (int *)swift_malloc("h_index", sizeof(int) * h_index_size);
  int *g_index = (int *)swift_malloc("g_index", sizeof(int) * g_index_size);
  int *s_index = (int *)swift_malloc("s_index", sizeof(int) * s_index_size);
  int *b_index = (int *)swift_malloc("b_index", sizeof(int) * b_index_size);
  int *sink_index =
      (int *)swift_malloc("sink_index", sizeof(int) * sink_index_size);
  if (h_index == NULL || g_index == NULL || s_index == NULL ||
      b_index == NULL || sink_index == NULL)
    error("Failed to allocate temporary particle indices.");

  /* Allocate counters of particles that will land in each cell */
  int *cell_part_counts =
      (int *)swift_malloc("cell_part_counts", sizeof(int) * s->nr_cells);
  int *cell_gpart_counts =
      (int *)swift_malloc("cell_gpart_counts", sizeof(int) * s->nr_cells);
  int *cell_spart_counts =
      (int *)swift_malloc("cell_spart_counts", sizeof(int) * s->nr_cells);
  int *cell_bpart_counts =
      (int *)swift_malloc("cell_bpart_counts", sizeof(int) * s->nr_cells);
  int *cell_sink_counts =
      (int *)swift_malloc("cell_sink_counts", sizeof(int) * s->nr_cells);

  if (cell_part_counts == NULL || cell_gpart_counts == NULL ||
      cell_spart_counts == NULL || cell_bpart_counts == NULL ||
      cell_sink_counts == NULL)
    error("Failed to allocate cell particle count buffer.");

  /* Initialise the counters, including buffer space for future particles */
  for (int i = 0; i < s->nr_cells; ++i) {
    cell_part_counts[i] = 0;
    cell_gpart_counts[i] = 0;
    cell_spart_counts[i] = 0;
    cell_bpart_counts[i] = 0;
    cell_sink_counts[i] = 0;
  }

  /* Run through the particles and get their cell index. */
  if (nr_parts > 0)
    space_parts_get_cell_index(s, h_index, cell_part_counts,
                               &count_inhibited_parts, &count_extra_parts,
                               verbose);
  if (nr_gparts > 0)
    space_gparts_get_cell_index(s, g_index, cell_gpart_counts,
                                &count_inhibited_gparts, &count_extra_gparts,
                                verbose);
  if (nr_sparts > 0)
    space_sparts_get_cell_index(s, s_index, cell_spart_counts,
                                &count_inhibited_sparts, &count_extra_sparts,
                                verbose);
  if (nr_bparts > 0)
    space_bparts_get_cell_index(s, b_index, cell_bpart_counts,
                                &count_inhibited_bparts, &count_extra_bparts,
                                verbose);
  if (nr_sinks > 0)
    space_sinks_get_cell_index(s, sink_index, cell_sink_counts,
                               &count_inhibited_sinks, &count_extra_sinks,
                               verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Some safety checks */
  if (repartitioned && count_inhibited_parts)
    error("We just repartitioned but still found inhibited parts.");
  if (repartitioned && count_inhibited_sparts)
    error("We just repartitioned but still found inhibited sparts.");
  if (repartitioned && count_inhibited_gparts)
    error("We just repartitioned but still found inhibited gparts.");
  if (repartitioned && count_inhibited_bparts)
    error("We just repartitioned but still found inhibited bparts.");
  if (repartitioned && count_inhibited_sinks)
    error("We just repartitioned but still found inhibited sinks.");

  if (count_extra_parts != s->nr_extra_parts)
    error(
        "Number of extra parts in the part array not matching the space "
        "counter.");
  if (count_extra_gparts != s->nr_extra_gparts)
    error(
        "Number of extra gparts in the gpart array not matching the space "
        "counter.");
  if (count_extra_sparts != s->nr_extra_sparts)
    error(
        "Number of extra sparts in the spart array not matching the space "
        "counter.");
  if (count_extra_bparts != s->nr_extra_bparts)
    error(
        "Number of extra bparts in the bpart array not matching the space "
        "counter.");
  if (count_extra_sinks != s->nr_extra_sinks)
    error(
        "Number of extra sinks in the sink array not matching the space "
        "counter.");
#endif

  const ticks tic2 = getticks();

  /* Move non-local parts and inhibited parts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_parts > 0)) {

    for (size_t k = 0; k < nr_parts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (h_index[k] == -1 || cells_top[h_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_parts -= 1;

        /* Swap the particle */
        memswap(&s->parts[k], &s->parts[nr_parts], sizeof(struct part));

        /* Swap the link with the gpart */
        if (s->parts[k].gpart != NULL) {
          s->parts[k].gpart->id_or_neg_offset = -k;
        }
        if (s->parts[nr_parts].gpart != NULL) {
          s->parts[nr_parts].gpart->id_or_neg_offset = -nr_parts;
        }

        /* Swap the xpart */
        memswap(&s->xparts[k], &s->xparts[nr_parts], sizeof(struct xpart));
        /* Swap the index */
        memswap(&h_index[k], &h_index[nr_parts], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all parts are in the correct places. */
  size_t check_count_inhibited_part = 0;
  for (size_t k = 0; k < nr_parts; k++) {
    if (h_index[k] == -1 || cells_top[h_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local parts to send list");
    }
  }
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    if (h_index[k] != -1 && cells_top[h_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local parts from send list");
    }
    if (h_index[k] == -1) ++check_count_inhibited_part;
  }
  if (check_count_inhibited_part != count_inhibited_parts)
    error("Counts of inhibited particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local sparts and inhibited sparts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_sparts > 0)) {

    for (size_t k = 0; k < nr_sparts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (s_index[k] == -1 || cells_top[s_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_sparts -= 1;

        /* Swap the particle */
        memswap(&s->sparts[k], &s->sparts[nr_sparts], sizeof(struct spart));

        /* Swap the link with the gpart */
        if (s->sparts[k].gpart != NULL) {
          s->sparts[k].gpart->id_or_neg_offset = -k;
        }
        if (s->sparts[nr_sparts].gpart != NULL) {
          s->sparts[nr_sparts].gpart->id_or_neg_offset = -nr_sparts;
        }

        /* Swap the index */
        memswap(&s_index[k], &s_index[nr_sparts], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all sparts are in the correct place. */
  size_t check_count_inhibited_spart = 0;
  for (size_t k = 0; k < nr_sparts; k++) {
    if (s_index[k] == -1 || cells_top[s_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local sparts to send list");
    }
  }
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    if (s_index[k] != -1 && cells_top[s_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local sparts from send list");
    }
    if (s_index[k] == -1) ++check_count_inhibited_spart;
  }
  if (check_count_inhibited_spart != count_inhibited_sparts)
    error("Counts of inhibited s-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local bparts and inhibited bparts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_bparts > 0)) {

    for (size_t k = 0; k < nr_bparts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (b_index[k] == -1 || cells_top[b_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_bparts -= 1;

        /* Swap the particle */
        memswap(&s->bparts[k], &s->bparts[nr_bparts], sizeof(struct bpart));

        /* Swap the link with the gpart */
        if (s->bparts[k].gpart != NULL) {
          s->bparts[k].gpart->id_or_neg_offset = -k;
        }
        if (s->bparts[nr_bparts].gpart != NULL) {
          s->bparts[nr_bparts].gpart->id_or_neg_offset = -nr_bparts;
        }

        /* Swap the index */
        memswap(&b_index[k], &b_index[nr_bparts], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all bparts are in the correct place. */
  size_t check_count_inhibited_bpart = 0;
  for (size_t k = 0; k < nr_bparts; k++) {
    if (b_index[k] == -1 || cells_top[b_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local bparts to send list");
    }
  }
  for (size_t k = nr_bparts; k < s->nr_bparts; k++) {
    if (b_index[k] != -1 && cells_top[b_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local bparts from send list");
    }
    if (b_index[k] == -1) ++check_count_inhibited_bpart;
  }
  if (check_count_inhibited_bpart != count_inhibited_bparts)
    error("Counts of inhibited b-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local sinks and inhibited sinks to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_sinks > 0)) {

    for (size_t k = 0; k < nr_sinks; /* void */) {

      /* Inhibited particle or foreign particle */
      if (sink_index[k] == -1 ||
          cells_top[sink_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_sinks -= 1;

        /* Swap the particle */
        memswap(&s->sinks[k], &s->sinks[nr_sinks], sizeof(struct sink));

        /* Swap the link with the gpart */
        if (s->sinks[k].gpart != NULL) {
          s->sinks[k].gpart->id_or_neg_offset = -k;
        }
        if (s->sinks[nr_sinks].gpart != NULL) {
          s->sinks[nr_sinks].gpart->id_or_neg_offset = -nr_sinks;
        }

        /* Swap the index */
        memswap(&sink_index[k], &sink_index[nr_sinks], sizeof(int));

      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all sinks are in the correct place. */
  size_t check_count_inhibited_sinks = 0;
  for (size_t k = 0; k < nr_sinks; k++) {
    if (sink_index[k] == -1 ||
        cells_top[sink_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local sinks to send list");
    }
  }
  for (size_t k = nr_sinks; k < s->nr_sinks; k++) {
    if (sink_index[k] != -1 &&
        cells_top[sink_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local sinks from send list");
    }
    if (sink_index[k] == -1) ++check_count_inhibited_sinks;
  }
  if (check_count_inhibited_sinks != count_inhibited_sinks)
    error("Counts of inhibited sink-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

  /* Move non-local gparts and inhibited parts to the end of the list. */
  if (!repartitioned && (s->e->nr_nodes > 1 || count_inhibited_gparts > 0)) {

    for (size_t k = 0; k < nr_gparts; /* void */) {

      /* Inhibited particle or foreign particle */
      if (g_index[k] == -1 || cells_top[g_index[k]].nodeID != local_nodeID) {

        /* One fewer particle */
        nr_gparts -= 1;

        /* Swap the particle */
        memswap_unaligned(&s->gparts[k], &s->gparts[nr_gparts],
                          sizeof(struct gpart));

        /* Swap the link with part/spart */
        if (s->gparts[k].type == swift_type_gas) {
          s->parts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        } else if (s->gparts[k].type == swift_type_stars) {
          s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        } else if (s->gparts[k].type == swift_type_sink) {
          s->sparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        } else if (s->gparts[k].type == swift_type_black_hole) {
          s->bparts[-s->gparts[k].id_or_neg_offset].gpart = &s->gparts[k];
        }

        if (s->gparts[nr_gparts].type == swift_type_gas) {
          s->parts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        } else if (s->gparts[nr_gparts].type == swift_type_stars) {
          s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        } else if (s->gparts[nr_gparts].type == swift_type_sink) {
          s->sparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        } else if (s->gparts[nr_gparts].type == swift_type_black_hole) {
          s->bparts[-s->gparts[nr_gparts].id_or_neg_offset].gpart =
              &s->gparts[nr_gparts];
        }

        /* Swap the index */
        memswap(&g_index[k], &g_index[nr_gparts], sizeof(int));
      } else {
        /* Increment when not exchanging otherwise we need to retest "k".*/
        k++;
      }
    }
  }

  if (verbose)
    message("Moving non-local particles took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all gparts are in the correct place. */
  size_t check_count_inhibited_gpart = 0;
  for (size_t k = 0; k < nr_gparts; k++) {
    if (g_index[k] == -1 || cells_top[g_index[k]].nodeID != local_nodeID) {
      error("Failed to move all non-local gparts to send list");
    }
  }
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    if (g_index[k] != -1 && cells_top[g_index[k]].nodeID == local_nodeID) {
      error("Failed to remove local gparts from send list");
    }
    if (g_index[k] == -1) ++check_count_inhibited_gpart;
  }
  if (check_count_inhibited_gpart != count_inhibited_gparts)
    error("Counts of inhibited g-particles do not match!");
#endif /* SWIFT_DEBUG_CHECKS */

#ifdef WITH_MPI

  /* Exchange the strays, note that this potentially re-allocates
     the parts arrays. This can be skipped if we just repartitioned space
     as there should be no strays in that case */
  if (!repartitioned) {

    size_t nr_parts_exchanged = s->nr_parts - nr_parts;
    size_t nr_gparts_exchanged = s->nr_gparts - nr_gparts;
    size_t nr_sparts_exchanged = s->nr_sparts - nr_sparts;
    size_t nr_bparts_exchanged = s->nr_bparts - nr_bparts;
    engine_exchange_strays(s->e, nr_parts, &h_index[nr_parts],
                           &nr_parts_exchanged, nr_gparts, &g_index[nr_gparts],
                           &nr_gparts_exchanged, nr_sparts, &s_index[nr_sparts],
                           &nr_sparts_exchanged, nr_bparts, &b_index[nr_bparts],
                           &nr_bparts_exchanged);

    /* Set the new particle counts. */
    s->nr_parts = nr_parts + nr_parts_exchanged;
    s->nr_gparts = nr_gparts + nr_gparts_exchanged;
    s->nr_sparts = nr_sparts + nr_sparts_exchanged;
    s->nr_bparts = nr_bparts + nr_bparts_exchanged;

  } else {
#ifdef SWIFT_DEBUG_CHECKS
    if (s->nr_parts != nr_parts)
      error("Number of parts changing after repartition");
    if (s->nr_sparts != nr_sparts)
      error("Number of sparts changing after repartition");
    if (s->nr_gparts != nr_gparts)
      error("Number of gparts changing after repartition");
#endif
  }

  /* Clear non-local cell counts. */
  for (int k = 0; k < s->nr_cells; k++) {
    if (s->cells_top[k].nodeID != local_nodeID) {
      cell_part_counts[k] = 0;
      cell_spart_counts[k] = 0;
      cell_gpart_counts[k] = 0;
      cell_bpart_counts[k] = 0;
    }
  }

  /* Re-allocate the index array for the parts if needed.. */
  if (s->nr_parts + 1 > h_index_size) {
    int *ind_new;
    if ((ind_new = (int *)swift_malloc(
             "h_index", sizeof(int) * (s->nr_parts + 1))) == NULL)
      error("Failed to allocate temporary particle indices.");
    memcpy(ind_new, h_index, sizeof(int) * nr_parts);
    swift_free("h_index", h_index);
    h_index = ind_new;
  }

  /* Re-allocate the index array for the sparts if needed.. */
  if (s->nr_sparts + 1 > s_index_size) {
    int *sind_new;
    if ((sind_new = (int *)swift_malloc(
             "s_index", sizeof(int) * (s->nr_sparts + 1))) == NULL)
      error("Failed to allocate temporary s-particle indices.");
    memcpy(sind_new, s_index, sizeof(int) * nr_sparts);
    swift_free("s_index", s_index);
    s_index = sind_new;
  }

  /* Re-allocate the index array for the bparts if needed.. */
  if (s->nr_bparts + 1 > s_index_size) {
    int *bind_new;
    if ((bind_new = (int *)swift_malloc(
             "b_index", sizeof(int) * (s->nr_bparts + 1))) == NULL)
      error("Failed to allocate temporary s-particle indices.");
    memcpy(bind_new, b_index, sizeof(int) * nr_bparts);
    swift_free("b_index", b_index);
    b_index = bind_new;
  }

  const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
  const double ih[3] = {s->iwidth[0], s->iwidth[1], s->iwidth[2]};

  /* Assign each received part to its cell. */
  for (size_t k = nr_parts; k < s->nr_parts; k++) {
    const struct part *const p = &s->parts[k];
    h_index[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cell_part_counts[h_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[h_index[k]].nodeID != local_nodeID)
      error("Received part that does not belong to me (nodeID=%i).",
            cells_top[h_index[k]].nodeID);
#endif
  }
  nr_parts = s->nr_parts;

  /* Assign each received spart to its cell. */
  for (size_t k = nr_sparts; k < s->nr_sparts; k++) {
    const struct spart *const sp = &s->sparts[k];
    s_index[k] =
        cell_getid(cdim, sp->x[0] * ih[0], sp->x[1] * ih[1], sp->x[2] * ih[2]);
    cell_spart_counts[s_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[s_index[k]].nodeID != local_nodeID)
      error("Received s-part that does not belong to me (nodeID=%i).",
            cells_top[s_index[k]].nodeID);
#endif
  }
  nr_sparts = s->nr_sparts;

  /* Assign each received bpart to its cell. */
  for (size_t k = nr_bparts; k < s->nr_bparts; k++) {
    const struct bpart *const bp = &s->bparts[k];
    b_index[k] =
        cell_getid(cdim, bp->x[0] * ih[0], bp->x[1] * ih[1], bp->x[2] * ih[2]);
    cell_bpart_counts[b_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[b_index[k]].nodeID != local_nodeID)
      error("Received s-part that does not belong to me (nodeID=%i).",
            cells_top[b_index[k]].nodeID);
#endif
  }
  nr_bparts = s->nr_bparts;

#else /* WITH_MPI */

  /* Update the part, spart and bpart counters */
  s->nr_parts = nr_parts;
  s->nr_sparts = nr_sparts;
  s->nr_bparts = nr_bparts;
  s->nr_sinks = nr_sinks;

#endif /* WITH_MPI */

  /* Sort the parts according to their cells. */
  if (nr_parts > 0)
    space_parts_sort(s->parts, s->xparts, h_index, cell_part_counts,
                     s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the part have been sorted correctly. */
  for (size_t k = 0; k < nr_parts; k++) {
    const struct part *p = &s->parts[k];

    if (p->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_ind =
        cell_getid(s->cdim, p->x[0] * s->iwidth[0], p->x[1] * s->iwidth[1],
                   p->x[2] * s->iwidth[2]);

    /* New cell of this part */
    const struct cell *c = &s->cells_top[new_ind];

    if (h_index[k] != new_ind)
      error("part's new cell index not matching sorted index.");

    if (p->x[0] < c->loc[0] || p->x[0] > c->loc[0] + c->width[0] ||
        p->x[1] < c->loc[1] || p->x[1] > c->loc[1] + c->width[1] ||
        p->x[2] < c->loc[2] || p->x[2] > c->loc[2] + c->width[2])
      error("part not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Sort the sparts according to their cells. */
  if (nr_sparts > 0)
    space_sparts_sort(s->sparts, s_index, cell_spart_counts, s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the spart have been sorted correctly. */
  for (size_t k = 0; k < nr_sparts; k++) {
    const struct spart *sp = &s->sparts[k];

    if (sp->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_sind =
        cell_getid(s->cdim, sp->x[0] * s->iwidth[0], sp->x[1] * s->iwidth[1],
                   sp->x[2] * s->iwidth[2]);

    /* New cell of this spart */
    const struct cell *c = &s->cells_top[new_sind];

    if (s_index[k] != new_sind)
      error("spart's new cell index not matching sorted index.");

    if (sp->x[0] < c->loc[0] || sp->x[0] > c->loc[0] + c->width[0] ||
        sp->x[1] < c->loc[1] || sp->x[1] > c->loc[1] + c->width[1] ||
        sp->x[2] < c->loc[2] || sp->x[2] > c->loc[2] + c->width[2])
      error("spart not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Sort the bparts according to their cells. */
  if (nr_bparts > 0)
    space_bparts_sort(s->bparts, b_index, cell_bpart_counts, s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the bpart have been sorted correctly. */
  for (size_t k = 0; k < nr_bparts; k++) {
    const struct bpart *bp = &s->bparts[k];

    if (bp->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_bind =
        cell_getid(s->cdim, bp->x[0] * s->iwidth[0], bp->x[1] * s->iwidth[1],
                   bp->x[2] * s->iwidth[2]);

    /* New cell of this bpart */
    const struct cell *c = &s->cells_top[new_bind];

    if (b_index[k] != new_bind)
      error("bpart's new cell index not matching sorted index.");

    if (bp->x[0] < c->loc[0] || bp->x[0] > c->loc[0] + c->width[0] ||
        bp->x[1] < c->loc[1] || bp->x[1] > c->loc[1] + c->width[1] ||
        bp->x[2] < c->loc[2] || bp->x[2] > c->loc[2] + c->width[2])
      error("bpart not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Sort the sink according to their cells. */
  if (nr_sinks > 0)
    space_sinks_sort(s->sinks, sink_index, cell_sink_counts, s->nr_cells, 0);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the sink have been sorted correctly. */
  for (size_t k = 0; k < nr_sinks; k++) {
    const struct sink *sink = &s->sinks[k];

    if (sink->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_bind =
        cell_getid(s->cdim, sink->x[0] * s->iwidth[0],
                   sink->x[1] * s->iwidth[1], sink->x[2] * s->iwidth[2]);

    /* New cell of this sink */
    const struct cell *c = &s->cells_top[new_bind];

    if (sink_index[k] != new_bind)
      error("sink's new cell index not matching sorted index.");

    if (sink->x[0] < c->loc[0] || sink->x[0] > c->loc[0] + c->width[0] ||
        sink->x[1] < c->loc[1] || sink->x[1] > c->loc[1] + c->width[1] ||
        sink->x[2] < c->loc[2] || sink->x[2] > c->loc[2] + c->width[2])
      error("sink not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_index = 0;
  h_index[nr_parts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_parts; k++) {
    if (h_index[k] < h_index[k + 1]) {
      cells_top[h_index[k]].hydro.count =
          k - last_index + 1 - space_extra_parts;
      last_index = k + 1;
    }
  }

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_sindex = 0;
  s_index[nr_sparts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_sparts; k++) {
    if (s_index[k] < s_index[k + 1]) {
      cells_top[s_index[k]].stars.count =
          k - last_sindex + 1 - space_extra_sparts;
      last_sindex = k + 1;
    }
  }

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_bindex = 0;
  b_index[nr_bparts] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_bparts; k++) {
    if (b_index[k] < b_index[k + 1]) {
      cells_top[b_index[k]].black_holes.count =
          k - last_bindex + 1 - space_extra_bparts;
      last_bindex = k + 1;
    }
  }

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_sink_index = 0;
  sink_index[nr_sinks] = s->nr_cells;  // sentinel.
  for (size_t k = 0; k < nr_sinks; k++) {
    if (sink_index[k] < sink_index[k + 1]) {
      cells_top[sink_index[k]].sinks.count =
          k - last_sink_index + 1 - space_extra_sinks;
      last_sink_index = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  swift_free("h_index", h_index);
  swift_free("cell_part_counts", cell_part_counts);
  swift_free("s_index", s_index);
  swift_free("cell_spart_counts", cell_spart_counts);
  swift_free("b_index", b_index);
  swift_free("cell_bpart_counts", cell_bpart_counts);
  swift_free("sink_index", sink_index);
  swift_free("cell_sink_counts", cell_sink_counts);

  /* Update the slice of unique IDs. */
  space_update_unique_id(s);

#ifdef WITH_MPI

  /* Re-allocate the index array for the gparts if needed.. */
  if (s->nr_gparts + 1 > g_index_size) {
    int *gind_new;
    if ((gind_new = (int *)swift_malloc(
             "g_index", sizeof(int) * (s->nr_gparts + 1))) == NULL)
      error("Failed to allocate temporary g-particle indices.");
    memcpy(gind_new, g_index, sizeof(int) * nr_gparts);
    swift_free("g_index", g_index);
    g_index = gind_new;
  }

  /* Assign each received gpart to its cell. */
  for (size_t k = nr_gparts; k < s->nr_gparts; k++) {
    const struct gpart *const p = &s->gparts[k];
    g_index[k] =
        cell_getid(cdim, p->x[0] * ih[0], p->x[1] * ih[1], p->x[2] * ih[2]);
    cell_gpart_counts[g_index[k]]++;
#ifdef SWIFT_DEBUG_CHECKS
    if (cells_top[g_index[k]].nodeID != s->e->nodeID)
      error("Received g-part that does not belong to me (nodeID=%i).",
            cells_top[g_index[k]].nodeID);
#endif
  }
  nr_gparts = s->nr_gparts;

#else /* WITH_MPI */

  /* Update the gpart counter */
  s->nr_gparts = nr_gparts;

#endif /* WITH_MPI */

  /* Mark that there are no inhibited particles left */
  s->nr_inhibited_parts = 0;
  s->nr_inhibited_gparts = 0;
  s->nr_inhibited_sparts = 0;
  s->nr_inhibited_bparts = 0;
  s->nr_inhibited_sinks = 0;

  /* Sort the gparts according to their cells. */
  if (nr_gparts > 0)
    space_gparts_sort(s->gparts, s->parts, s->sinks, s->sparts, s->bparts,
                      g_index, cell_gpart_counts, s->nr_cells);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the gpart have been sorted correctly. */
  for (size_t k = 0; k < nr_gparts; k++) {
    const struct gpart *gp = &s->gparts[k];

    if (gp->time_bin == time_bin_inhibited)
      error("Inhibited particle sorted into a cell!");

    /* New cell index */
    const int new_gind =
        cell_getid(s->cdim, gp->x[0] * s->iwidth[0], gp->x[1] * s->iwidth[1],
                   gp->x[2] * s->iwidth[2]);

    /* New cell of this gpart */
    const struct cell *c = &s->cells_top[new_gind];

    if (g_index[k] != new_gind)
      error("gpart's new cell index not matching sorted index.");

    if (gp->x[0] < c->loc[0] || gp->x[0] > c->loc[0] + c->width[0] ||
        gp->x[1] < c->loc[1] || gp->x[1] > c->loc[1] + c->width[1] ||
        gp->x[2] < c->loc[2] || gp->x[2] > c->loc[2] + c->width[2])
      error("gpart not sorted into the right top-level cell!");
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Extract the cell counts from the sorted indices. Deduct the extra
   * particles. */
  size_t last_gindex = 0;
  g_index[nr_gparts] = s->nr_cells;
  for (size_t k = 0; k < nr_gparts; k++) {
    if (g_index[k] < g_index[k + 1]) {
      cells_top[g_index[k]].grav.count =
          k - last_gindex + 1 - space_extra_gparts;
      last_gindex = k + 1;
    }
  }

  /* We no longer need the indices as of here. */
  swift_free("g_index", g_index);
  swift_free("cell_gpart_counts", cell_gpart_counts);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that the links are correct */
  if ((nr_gparts > 0 && nr_parts > 0) || (nr_gparts > 0 && nr_sparts > 0) ||
      (nr_gparts > 0 && nr_bparts > 0) || (nr_gparts > 0 && nr_sinks > 0))
    part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                      nr_parts, nr_gparts, nr_sinks, nr_sparts, nr_bparts,
                      verbose);
#endif

  /* Hook the cells up to the parts. Make list of local and non-empty cells */
  const ticks tic3 = getticks();
  struct part *finger = s->parts;
  struct xpart *xfinger = s->xparts;
  struct gpart *gfinger = s->gparts;
  struct spart *sfinger = s->sparts;
  struct bpart *bfinger = s->bparts;
  struct sink *sink_finger = s->sinks;
  s->nr_cells_with_particles = 0;
  s->nr_local_cells_with_particles = 0;
  s->nr_local_cells = 0;
  for (int k = 0; k < s->nr_cells; k++) {
    struct cell *restrict c = &cells_top[k];
    c->hydro.ti_old_part = ti_current;
    c->grav.ti_old_part = ti_current;
    c->grav.ti_old_multipole = ti_current;
    c->stars.ti_old_part = ti_current;
    c->sinks.ti_old_part = ti_current;
    c->black_holes.ti_old_part = ti_current;

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
    cell_assign_top_level_cell_index(c, s->cdim, s->dim, s->width);
#endif

    const int is_local = (c->nodeID == engine_rank);
    const int has_particles =
        (c->hydro.count > 0) || (c->grav.count > 0) || (c->stars.count > 0) ||
        (c->black_holes.count > 0) || (c->sinks.count > 0);

    if (is_local) {
      c->hydro.parts = finger;
      c->hydro.xparts = xfinger;
      c->grav.parts = gfinger;
      c->stars.parts = sfinger;
      c->black_holes.parts = bfinger;
      c->sinks.parts = sink_finger;

      /* Store the state at rebuild time */
      c->stars.parts_rebuild = c->stars.parts;
      c->grav.parts_rebuild = c->grav.parts;

      c->hydro.count_total = c->hydro.count + space_extra_parts;
      c->grav.count_total = c->grav.count + space_extra_gparts;
      c->stars.count_total = c->stars.count + space_extra_sparts;
      c->sinks.count_total = c->sinks.count + space_extra_sinks;
      c->black_holes.count_total = c->black_holes.count + space_extra_bparts;

      finger = &finger[c->hydro.count_total];
      xfinger = &xfinger[c->hydro.count_total];
      gfinger = &gfinger[c->grav.count_total];
      sfinger = &sfinger[c->stars.count_total];
      bfinger = &bfinger[c->black_holes.count_total];
      sink_finger = &sink_finger[c->sinks.count_total];

      /* Add this cell to the list of local cells */
      s->local_cells_top[s->nr_local_cells] = k;
      s->nr_local_cells++;
    }

    if (is_local && has_particles) {

      /* Add this cell to the list of non-empty cells */
      s->local_cells_with_particles_top[s->nr_local_cells_with_particles] = k;
      s->nr_local_cells_with_particles++;
    }
  }
  if (verbose) {
    message("Have %d local top-level cells with particles (total=%d)",
            s->nr_local_cells_with_particles, s->nr_cells);
    message("Have %d local top-level cells (total=%d)", s->nr_local_cells,
            s->nr_cells);
    message("hooking up cells took %.3f %s.",
            clocks_from_ticks(getticks() - tic3), clocks_getunit());
  }

  /* Re-order the extra particles such that they are at the end of their cell's
     memory pool. */
  if (s->with_star_formation || s->e->policy & engine_policy_sinks)
    space_reorder_extras(s, verbose);

  /* At this point, we have the upper-level cells. Now recursively split each
     cell to get the full AMR grid. */
  space_split(s, verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the multipole construction went OK */
  if (s->with_self_gravity)
    for (int k = 0; k < s->nr_cells; k++)
      cell_check_multipole(&s->cells_top[k], s->e->gravity_properties);
#endif

  /* Clean up any stray sort indices in the cell buffer. */
  space_free_buff_sort_indices(s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
