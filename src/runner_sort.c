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
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"

/*! The size of the sorting stack used at the leaf level */
const int sort_stack_size = 10;

/**
 * @brief Sorts again all the stars in a given cell hierarchy.
 *
 * This is intended to be used after the star formation task has been run
 * to get the cells back into a state where self/pair star tasks can be run.
 *
 * @param r The thread #runner.
 * @param c The top-level cell to run on.
 * @param timer Are we timing this?
 */
void runner_do_stars_resort(struct runner *r, struct cell *c, const int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != r->e->nodeID) error("Task must be run locally!");
#endif

  TIMER_TIC;

  /* Did we demand a recalculation of the stars'sorts? */
  if (cell_get_flag(c, cell_flag_do_stars_resort)) {
    runner_do_all_stars_sort(r, c);
    cell_clear_flag(c, cell_flag_do_stars_resort);
  }

  if (timer) TIMER_TOC(timer_do_stars_resort);
}

/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */
void runner_do_sort_ascending(struct sort_entry *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[sort_stack_size];
  int qpos, i, j, lo, hi, imin;
  struct sort_entry temp;
  float pivot;

  if (N >= (1LL << sort_stack_size)) {
    error(
        "The stack size for sorting is too small."
        "Either increase it or reduce the number of parts per cell.");
  }

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    /* Do we have a low number of element to sort? */
    if (hi - lo < 15) {
      /* Sort the last elements. */
      for (i = lo; i < hi; i++) {
        imin = i;
        /* Find the minimal value. */
        for (j = i + 1; j <= hi; j++) {
          if (sort[j].d < sort[imin].d) {
            imin = j;
          }
        }
        /* Swap the elements if a smaller element exists. */
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      /* Select a pivot */
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      /* Ensure that the elements before/after the pivot
         are smaller/larger than the pivot. */
      while (i <= j) {
        /* Find the first elements that do not respect
           the order. */
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        /* Did we get two different elements */
        if (i <= j) {
          if (i < j) {
            /* Swap the elements */
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      /* Add the next operations to the stack.
       * The order is important in order to decrease the stack size.
       */
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Recursively checks that the flags are consistent in a cell hierarchy.
 *
 * Debugging function. Exists in two flavours: hydro & stars.
 */
#define RUNNER_CHECK_SORTS(TYPE)                                               \
  void runner_check_sorts_##TYPE(struct cell *c, int flags) {                  \
                                                                               \
    if (flags & ~c->TYPE.sorted) error("Inconsistent sort flags (downward)!"); \
    if (c->split)                                                              \
      for (int k = 0; k < 8; k++)                                              \
        if (c->progeny[k] != NULL && c->progeny[k]->TYPE.count > 0)            \
          runner_check_sorts_##TYPE(c->progeny[k], c->TYPE.sorted);            \
  }
#else
#define RUNNER_CHECK_SORTS(TYPE)                                       \
  void runner_check_sorts_##TYPE(struct cell *c, int flags) {          \
    error("Calling debugging code without debugging flag activated."); \
  }
#endif

RUNNER_CHECK_SORTS(hydro)
RUNNER_CHECK_SORTS(stars)

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param cleanup If true, re-build the sorts for the selected flags instead
 *        of just adding them.
 * @param rt_requests_sort whether this sort was requested for RT. If true,
 *        this cell is allowed to be undrifted.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
void runner_do_hydro_sort(struct runner *r, struct cell *c, int flags,
                          int cleanup, int rt_requests_sort, int clock) {

  struct sort_entry *fingers[8];
  const int count = c->hydro.count;
  const struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  float buff[8];

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.super == NULL) error("Task called above the super level!!!");
#endif

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->hydro.do_sort;
  if (cleanup) {
    c->hydro.sorted = 0;
  } else {
    flags &= ~c->hydro.sorted;
  }
  if (flags == 0 && !cell_get_flag(c, cell_flag_do_hydro_sub_sort) &&
      !cell_get_flag(c, cell_flag_do_rt_sub_sort))
    return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_part_drifted(c, r->e)) {
    /* If the sort was requested by RT, cell may be intentionally
     * undrifted. */
    if (!rt_requests_sort) error("Sorting un-drifted cell");
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_hydro(c, c->hydro.sorted);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->hydro.sorted & ~c->hydro.sorted)
      error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->hydro.sorted == 0) c->hydro.ti_sort = r->e->ti_current;
#endif

  /* Allocate memory for sorting. */
  cell_malloc_hydro_sorts(c, flags);

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        if (c->progeny[k]->hydro.count > 0) {

          /* Only propagate cleanup if the progeny is stale. */
          runner_do_hydro_sort(
              r, c->progeny[k], flags,
              cleanup && (c->progeny[k]->hydro.dx_max_sort_old >
                          space_maxreldx * c->progeny[k]->dmin),
              rt_requests_sort, 0);
          dx_max_sort = max(dx_max_sort, c->progeny[k]->hydro.dx_max_sort);
          dx_max_sort_old =
              max(dx_max_sort_old, c->progeny[k]->hydro.dx_max_sort_old);
        } else {

          /* We need to clean up the unused flags that were in case the
             number of particles in the cell would change */
          cell_clear_hydro_sort_flags(c->progeny[k], /*clear_unused_flags=*/1);
        }
      }
    }
    c->hydro.dx_max_sort = dx_max_sort;
    c->hydro.dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->hydro.count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->hydro.count > 0) {
          fingers[k] = cell_get_hydro_sorts(c->progeny[k], j);
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (int i = 0; i < 7; i++)
        for (int k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            int temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      struct sort_entry *finger = cell_get_hydro_sorts(c, j);
      for (int ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (int k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          int temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */

      struct sort_entry *entries = cell_get_hydro_sorts(c, j);
      entries[count].d = FLT_MAX;
      entries[count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->hydro.sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->hydro.sorted == 0) {
#ifdef SWIFT_DEBUG_CHECKS
      if (xparts != NULL && c->nodeID != engine_rank)
        error("Have non-NULL xparts in foreign cell");
#endif

      /* And the individual sort distances if we are a local cell */
      if (xparts != NULL) {
        for (int k = 0; k < count; k++) {
          xparts[k].x_diff_sort[0] = 0.0f;
          xparts[k].x_diff_sort[1] = 0.0f;
          xparts[k].x_diff_sort[2] = 0.0f;
        }
      }
      c->hydro.dx_max_sort_old = 0.f;
      c->hydro.dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {parts[k].x[0], parts[k].x[1], parts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          struct sort_entry *entries = cell_get_hydro_sorts(c, j);
          entries[k].i = k;
          entries[k].d = px[0] * runner_shift[j][0] +
                         px[1] * runner_shift[j][1] +
                         px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        struct sort_entry *entries = cell_get_hydro_sorts(c, j);
        entries[count].d = FLT_MAX;
        entries[count].i = 0;
        runner_do_sort_ascending(entries, count);
        atomic_or(&c->hydro.sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct sort_entry *finger = cell_get_hydro_sorts(c, j);
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_hydro(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->hydro.sorted & ~c->hydro.sorted)
      error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->hydro.do_sort = 0;
  cell_clear_flag(c, cell_flag_do_hydro_sub_sort);
  cell_clear_flag(c, cell_flag_do_rt_sub_sort);
  cell_clear_flag(c, cell_flag_rt_requests_sort);
  c->hydro.requires_sorts = 0;

  if (clock) TIMER_TOC(timer_dosort);
}

/**
 * @brief Sort the stars particles in the given cell along all cardinal
 * directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param cleanup If true, re-build the sorts for the selected flags instead
 *        of just adding them.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
void runner_do_stars_sort(struct runner *r, struct cell *c, int flags,
                          int cleanup, int clock) {

  struct sort_entry *fingers[8];
  const int count = c->stars.count;
  struct spart *sparts = c->stars.parts;
  float buff[8];

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.super == NULL) error("Task called above the super level!!!");
#endif

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->stars.do_sort;
  if (cleanup) {
    c->stars.sorted = 0;
  } else {
    flags &= ~c->stars.sorted;
  }
  if (flags == 0 && !cell_get_flag(c, cell_flag_do_stars_sub_sort)) return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_spart_drifted(c, r->e)) {
    error("Sorting un-drifted cell c->nodeID=%d", c->nodeID);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_stars(c, c->stars.sorted);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->stars.sorted == 0) c->stars.ti_sort = r->e->ti_current;
#endif

  /* start by allocating the entry arrays in the requested dimensions. */
  cell_malloc_stars_sorts(c, flags);

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {

        if (c->progeny[k]->stars.count > 0) {

          /* Only propagate cleanup if the progeny is stale. */
          const int cleanup_prog =
              cleanup && (c->progeny[k]->stars.dx_max_sort_old >
                          space_maxreldx * c->progeny[k]->dmin);
          runner_do_stars_sort(r, c->progeny[k], flags, cleanup_prog, 0);
          dx_max_sort = max(dx_max_sort, c->progeny[k]->stars.dx_max_sort);
          dx_max_sort_old =
              max(dx_max_sort_old, c->progeny[k]->stars.dx_max_sort_old);
        } else {

          /* We need to clean up the unused flags that were in case the
             number of particles in the cell would change */
          cell_clear_stars_sort_flags(c->progeny[k], /*clear_unused_flags=*/1);
        }
      }
    }
    c->stars.dx_max_sort = dx_max_sort;
    c->stars.dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->stars.count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
          fingers[k] = cell_get_stars_sorts(c->progeny[k], j);
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (int i = 0; i < 7; i++)
        for (int k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            int temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      struct sort_entry *finger = cell_get_stars_sorts(c, j);
      for (int ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (int k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          int temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      struct sort_entry *entries = cell_get_stars_sorts(c, j);
      entries[count].d = FLT_MAX;
      entries[count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->stars.sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->stars.sorted == 0) {

      /* And the individual sort distances if we are a local cell */
      for (int k = 0; k < count; k++) {
        sparts[k].x_diff_sort[0] = 0.0f;
        sparts[k].x_diff_sort[1] = 0.0f;
        sparts[k].x_diff_sort[2] = 0.0f;
      }
      c->stars.dx_max_sort_old = 0.f;
      c->stars.dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {sparts[k].x[0], sparts[k].x[1], sparts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          struct sort_entry *entries = cell_get_stars_sorts(c, j);
          entries[k].i = k;
          entries[k].d = px[0] * runner_shift[j][0] +
                         px[1] * runner_shift[j][1] +
                         px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        struct sort_entry *entries = cell_get_stars_sorts(c, j);
        entries[count].d = FLT_MAX;
        entries[count].i = 0;
        runner_do_sort_ascending(entries, count);
        atomic_or(&c->stars.sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct sort_entry *finger = cell_get_stars_sorts(c, j);
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_stars(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->stars.do_sort = 0;
  cell_clear_flag(c, cell_flag_do_stars_sub_sort);
  c->stars.requires_sorts = 0;

  if (clock) TIMER_TOC(timer_do_stars_sort);
}

/**
 * @brief Recurse into a cell until reaching the super level and call
 * the hydro sorting function there.
 *
 * This function must be called at or above the super level!
 *
 * This function will sort the particles in all 13 directions.
 *
 * @param r the #runner.
 * @param c the #cell.
 */
void runner_do_all_hydro_sort(struct runner *r, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Function called on a foreign cell!");
#endif

  if (!cell_is_active_hydro(c, r->e)) return;

  /* Shall we sort at this level? */
  if (c->hydro.super == c) {

    /* Sort everything */
    runner_do_hydro_sort(r, c, 0x1FFF, /*cleanup=*/0, /*rt_requests_sort=*/0,
                         /*timer=*/0);

  } else {

#ifdef SWIFT_DEBUG_CHECKS
    if (c->hydro.super != NULL) error("Function called below the super level!");
#endif

    /* Ok, then, let's try lower */
    if (c->split) {
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) runner_do_all_hydro_sort(r, c->progeny[k]);
      }
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Reached a leaf without encountering a hydro super cell!");
#endif
    }
  }
}

/**
 * @brief Recurse into a cell until reaching the super level and call
 * the star sorting function there.
 *
 * This function must be called at or above the super level!
 *
 * This function will sort the particles in all 13 directions.
 *
 * @param r the #runner.
 * @param c the #cell.
 */
void runner_do_all_stars_sort(struct runner *r, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Function called on a foreign cell!");
#endif

  if (!cell_is_active_stars(c, r->e) && !cell_is_active_hydro(c, r->e)) return;

  /* Shall we sort at this level? */
  if (c->hydro.super == c) {

    /* Sort everything */
    runner_do_stars_sort(r, c, 0x1FFF, /*cleanup=*/0, /*timer=*/0);

  } else {

#ifdef SWIFT_DEBUG_CHECKS
    if (c->hydro.super != NULL) error("Function called below the super level!");
#endif

    /* Ok, then, let's try lower */
    if (c->split) {
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) runner_do_all_stars_sort(r, c->progeny[k]);
      }
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Reached a leaf without encountering a hydro super cell!");
#endif
    }
  }
}
