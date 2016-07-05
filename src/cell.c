/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Switch off timers. */
#ifdef TIMER
#undef TIMER
#endif

/* This object's header. */
#include "cell.h"

/* Local headers. */
#include "atomic.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "space.h"
#include "timers.h"

/* Global variables. */
int cell_next_tag = 0;

/**
 * @brief Get the size of the cell subtree.
 *
 * @param c The #cell.
 */

int cell_getsize(struct cell *c) {

  /* Number of cells in this subtree. */
  int count = 1;

  /* Sum up the progeny if split. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) count += cell_getsize(c->progeny[k]);

  /* Return the final count. */
  return count;
}

/**
 * @brief Unpack the data of a given cell and its sub-cells.
 *
 * @param pc An array of packed #pcell.
 * @param c The #cell in which to unpack the #pcell.
 * @param s The #space in which the cells are created.
 *
 * @return The number of cells created.
 */

int cell_unpack(struct pcell *pc, struct cell *c, struct space *s) {

  /* Unpack the current pcell. */
  c->h_max = pc->h_max;
  c->ti_end_min = pc->ti_end_min;
  c->ti_end_max = pc->ti_end_max;
  c->count = pc->count;
  c->gcount = pc->gcount;
  c->tag = pc->tag;

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  for (int k = 0; k < 8; k++)
    if (pc->progeny[k] >= 0) {
      struct cell *temp = space_getcell(s);
      temp->count = 0;
      temp->gcount = 0;
      temp->loc[0] = c->loc[0];
      temp->loc[1] = c->loc[1];
      temp->loc[2] = c->loc[2];
      temp->h[0] = c->h[0] / 2;
      temp->h[1] = c->h[1] / 2;
      temp->h[2] = c->h[2] / 2;
      temp->dmin = c->dmin / 2;
      if (k & 4) temp->loc[0] += temp->h[0];
      if (k & 2) temp->loc[1] += temp->h[1];
      if (k & 1) temp->loc[2] += temp->h[2];
      temp->depth = c->depth + 1;
      temp->split = 0;
      temp->dx_max = 0.f;
      temp->nodeID = c->nodeID;
      temp->parent = c;
      c->progeny[k] = temp;
      c->split = 1;
      count += cell_unpack(&pc[pc->progeny[k]], temp, s);
    }

  /* Return the total number of unpacked cells. */
  c->pcell_size = count;
  return count;
}

/**
 * @brief Link the cells recursively to the given #part array.
 *
 * @param c The #cell.
 * @param parts The #part array.
 *
 * @return The number of particles linked.
 */

int cell_link_parts(struct cell *c, struct part *parts) {

  c->parts = parts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_parts(c->progeny[k], &parts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->count;
}

/**
 * @brief Link the cells recursively to the given #gpart array.
 *
 * @param c The #cell.
 * @param gparts The #gpart array.
 *
 * @return The number of particles linked.
 */

int cell_link_gparts(struct cell *c, struct gpart *gparts) {

  c->gparts = gparts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_gparts(c->progeny[k], &gparts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->gcount;
}

/**
 * @brief Pack the data of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param pc Pointer to an array of packed cells in which the
 *      cells will be packed.
 *
 * @return The number of packed cells.
 */

int cell_pack(struct cell *c, struct pcell *pc) {

  /* Start by packing the data of the current cell. */
  pc->h_max = c->h_max;
  pc->ti_end_min = c->ti_end_min;
  pc->ti_end_max = c->ti_end_max;
  pc->count = c->count;
  pc->gcount = c->gcount;
  c->tag = pc->tag = atomic_inc(&cell_next_tag) % cell_max_tag;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      pc->progeny[k] = count;
      count += cell_pack(c->progeny[k], &pc[count]);
    } else
      pc->progeny[k] = -1;

  /* Return the number of packed cells used. */
  c->pcell_size = count;
  return count;
}

int cell_pack_ti_ends(struct cell *c, int *ti_ends) {

  /* Pack this cell's data. */
  ti_ends[0] = c->ti_end_min;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_ti_ends(c->progeny[k], &ti_ends[count]);
    }

  /* Return the number of packed values. */
  return count;
}

int cell_unpack_ti_ends(struct cell *c, int *ti_ends) {

  /* Unpack this cell's data. */
  c->ti_end_min = ti_ends[0];

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_ti_ends(c->progeny[k], &ti_ends[count]);
    }

  /* Return the number of packed values. */
  return count;
}

/**
 * @brief Lock a cell and hold its parents.
 *
 * @param c The #cell.
 */

int cell_locktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to lock this cell. */
  if (c->hold || lock_trylock(&c->lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->hold) {

    /* Unlock this cell. */
    if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {

    /* Lock this cell. */
    if (lock_trylock(&finger->lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->lock) != 0) error("Failed to unlock cell.");
  }

  /* If we reached the top of the tree, we're done. */
  if (finger == NULL) {
    TIMER_TOC(timer_locktree);
    return 0;
  }

  /* Otherwise, we hit a snag. */
  else {

    /* Undo the holds up to finger. */
    for (struct cell *finger2 = c->parent; finger2 != finger;
         finger2 = finger2->parent)
      atomic_dec(&finger2->hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

int cell_glocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to lock this cell. */
  if (c->ghold || lock_trylock(&c->glock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->ghold) {

    /* Unlock this cell. */
    if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {

    /* Lock this cell. */
    if (lock_trylock(&finger->glock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->ghold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->glock) != 0) error("Failed to unlock cell.");
  }

  /* If we reached the top of the tree, we're done. */
  if (finger == NULL) {
    TIMER_TOC(timer_locktree);
    return 0;
  }

  /* Otherwise, we hit a snag. */
  else {

    /* Undo the holds up to finger. */
    for (struct cell *finger2 = c->parent; finger2 != finger;
         finger2 = finger2->parent)
      atomic_dec(&finger2->ghold);

    /* Unlock this cell. */
    if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Unlock a cell's parents.
 *
 * @param c The #cell.
 */

void cell_unlocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->hold);

  TIMER_TOC(timer_locktree);
}

void cell_gunlocktree(struct cell *c) {

  TIMER_TIC

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->glock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->ghold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 * @param parts_offset Offset of the cell parts array relative to the
 *        space's parts array, i.e. c->parts - s->parts.
 */

void cell_split(struct cell *c, ptrdiff_t parts_offset) {

  int i, j;
  const int count = c->count, gcount = c->gcount;
  struct part *parts = c->parts;
  struct xpart *xparts = c->xparts;
  struct gpart *gparts = c->gparts;
  int left[8], right[8];
  double pivot[3];

  /* Init the pivots. */
  for (int k = 0; k < 3; k++) pivot[k] = c->loc[k] + c->h[k] / 2;

  /* Split along the x-axis. */
  i = 0;
  j = count - 1;
  while (i <= j) {
    while (i <= count - 1 && parts[i].x[0] <= pivot[0]) i += 1;
    while (j >= 0 && parts[j].x[0] > pivot[0]) j -= 1;
    if (i < j) {
      struct part temp = parts[i];
      parts[i] = parts[j];
      parts[j] = temp;
      struct xpart xtemp = xparts[i];
      xparts[i] = xparts[j];
      xparts[j] = xtemp;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int k = 0; k <= j; k++)
    if (parts[k].x[0] > pivot[0]) error("cell_split: sorting failed.");
  for (int k = i; k < count; k++)
    if (parts[k].x[0] < pivot[0]) error("cell_split: sorting failed.");
#endif

  left[1] = i;
  right[1] = count - 1;
  left[0] = 0;
  right[0] = j;

  /* Split along the y axis, twice. */
  for (int k = 1; k >= 0; k--) {
    i = left[k];
    j = right[k];
    while (i <= j) {
      while (i <= right[k] && parts[i].x[1] <= pivot[1]) i += 1;
      while (j >= left[k] && parts[j].x[1] > pivot[1]) j -= 1;
      if (i < j) {
        struct part temp = parts[i];
        parts[i] = parts[j];
        parts[j] = temp;
        struct xpart xtemp = xparts[i];
        xparts[i] = xparts[j];
        xparts[j] = xtemp;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    for (int kk = left[k]; kk <= j; kk++)
      if (parts[kk].x[1] > pivot[1]) {
        message("ival=[%i,%i], i=%i, j=%i.", left[k], right[k], i, j);
        error("sorting failed (left).");
      }
    for (int kk = i; kk <= right[k]; kk++)
      if (parts[kk].x[1] < pivot[1]) error("sorting failed (right).");
#endif

    left[2 * k + 1] = i;
    right[2 * k + 1] = right[k];
    left[2 * k] = left[k];
    right[2 * k] = j;
  }

  /* Split along the z axis, four times. */
  for (int k = 3; k >= 0; k--) {
    i = left[k];
    j = right[k];
    while (i <= j) {
      while (i <= right[k] && parts[i].x[2] <= pivot[2]) i += 1;
      while (j >= left[k] && parts[j].x[2] > pivot[2]) j -= 1;
      if (i < j) {
        struct part temp = parts[i];
        parts[i] = parts[j];
        parts[j] = temp;
        struct xpart xtemp = xparts[i];
        xparts[i] = xparts[j];
        xparts[j] = xtemp;
      }
    }

#ifdef SWIFT_DEBUG_CHECKS
    for (int kk = left[k]; kk <= j; kk++)
      if (parts[kk].x[2] > pivot[2]) {
        message("ival=[%i,%i], i=%i, j=%i.", left[k], right[k], i, j);
        error("sorting failed (left).");
      }
    for (int kk = i; kk <= right[k]; kk++)
      if (parts[kk].x[2] < pivot[2]) {
        message("ival=[%i,%i], i=%i, j=%i.", left[k], right[k], i, j);
        error("sorting failed (right).");
      }
#endif

    left[2 * k + 1] = i;
    right[2 * k + 1] = right[k];
    left[2 * k] = left[k];
    right[2 * k] = j;
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->count = right[k] - left[k] + 1;
    c->progeny[k]->parts = &c->parts[left[k]];
    c->progeny[k]->xparts = &c->xparts[left[k]];
  }

  /* Re-link the gparts. */
  part_relink_gparts(parts, count, parts_offset);

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that _all_ the parts have been assigned to a cell. */
  for (int k = 1; k < 8; k++)
    if (&c->progeny[k - 1]->parts[c->progeny[k - 1]->count] !=
        c->progeny[k]->parts)
      error("Particle sorting failed (internal consistency).");
  if (c->progeny[0]->parts != c->parts)
    error("Particle sorting failed (left edge).");
  if (&c->progeny[7]->parts[c->progeny[7]->count] != &c->parts[count])
    error("Particle sorting failed (right edge).");

  /* Verify a few sub-cells. */
  for (int k = 0; k < c->progeny[0]->count; k++)
    if (c->progeny[0]->parts[k].x[0] > pivot[0] ||
        c->progeny[0]->parts[k].x[1] > pivot[1] ||
        c->progeny[0]->parts[k].x[2] > pivot[2])
      error("Sorting failed (progeny=0).");
  for (int k = 0; k < c->progeny[1]->count; k++)
    if (c->progeny[1]->parts[k].x[0] > pivot[0] ||
        c->progeny[1]->parts[k].x[1] > pivot[1] ||
        c->progeny[1]->parts[k].x[2] <= pivot[2])
      error("Sorting failed (progeny=1).");
  for (int k = 0; k < c->progeny[2]->count; k++)
    if (c->progeny[2]->parts[k].x[0] > pivot[0] ||
        c->progeny[2]->parts[k].x[1] <= pivot[1] ||
        c->progeny[2]->parts[k].x[2] > pivot[2])
      error("Sorting failed (progeny=2).");
#endif

  /* Now do the same song and dance for the gparts. */

  /* Split along the x-axis. */
  i = 0;
  j = gcount - 1;
  while (i <= j) {
    while (i <= gcount - 1 && gparts[i].x[0] <= pivot[0]) i += 1;
    while (j >= 0 && gparts[j].x[0] > pivot[0]) j -= 1;
    if (i < j) {
      struct gpart gtemp = gparts[i];
      gparts[i] = gparts[j];
      gparts[j] = gtemp;
    }
  }
  left[1] = i;
  right[1] = gcount - 1;
  left[0] = 0;
  right[0] = j;

  /* Split along the y axis, twice. */
  for (int k = 1; k >= 0; k--) {
    i = left[k];
    j = right[k];
    while (i <= j) {
      while (i <= right[k] && gparts[i].x[1] <= pivot[1]) i += 1;
      while (j >= left[k] && gparts[j].x[1] > pivot[1]) j -= 1;
      if (i < j) {
        struct gpart gtemp = gparts[i];
        gparts[i] = gparts[j];
        gparts[j] = gtemp;
      }
    }
    left[2 * k + 1] = i;
    right[2 * k + 1] = right[k];
    left[2 * k] = left[k];
    right[2 * k] = j;
  }

  /* Split along the z axis, four times. */
  for (int k = 3; k >= 0; k--) {
    i = left[k];
    j = right[k];
    while (i <= j) {
      while (i <= right[k] && gparts[i].x[2] <= pivot[2]) i += 1;
      while (j >= left[k] && gparts[j].x[2] > pivot[2]) j -= 1;
      if (i < j) {
        struct gpart gtemp = gparts[i];
        gparts[i] = gparts[j];
        gparts[j] = gtemp;
      }
    }
    left[2 * k + 1] = i;
    right[2 * k + 1] = right[k];
    left[2 * k] = left[k];
    right[2 * k] = j;
  }

  /* Store the counts and offsets. */
  for (int k = 0; k < 8; k++) {
    c->progeny[k]->gcount = right[k] - left[k] + 1;
    c->progeny[k]->gparts = &c->gparts[left[k]];
  }

  /* Re-link the parts. */
  part_relink_parts(gparts, gcount, parts - parts_offset);
}

/**
 * @brief Initialises all particles to a valid state even if the ICs were stupid
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_init_parts(struct cell *c, void *data) {

  struct part *p = c->parts;
  struct xpart *xp = c->xparts;
  const int count = c->count;

  for (int i = 0; i < count; ++i) {
    p[i].ti_begin = 0;
    p[i].ti_end = 0;
    xp[i].v_full[0] = p[i].v[0];
    xp[i].v_full[1] = p[i].v[1];
    xp[i].v_full[2] = p[i].v[2];
    hydro_first_init_part(&p[i], &xp[i]);
    hydro_init_part(&p[i]);
    hydro_reset_acceleration(&p[i]);
  }
  c->ti_end_min = 0;
  c->ti_end_max = 0;
}

/**
 * @brief Initialises all g-particles to a valid state even if the ICs were
 *stupid
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_init_gparts(struct cell *c, void *data) {

  struct gpart *gp = c->gparts;
  const int gcount = c->gcount;

  for (int i = 0; i < gcount; ++i) {
    gp[i].ti_begin = 0;
    gp[i].ti_end = 0;
    gravity_first_init_gpart(&gp[i]);
  }
  c->ti_end_min = 0;
  c->ti_end_max = 0;
}

/**
 * @brief Converts hydro quantities to a valid state after the initial density
 *calculation
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_convert_hydro(struct cell *c, void *data) {

  struct part *p = c->parts;

  for (int i = 0; i < c->count; ++i) {
    hydro_convert_quantities(&p[i]);
  }
}

/**
 * @brief Cleans the links in a given cell.
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_clean_links(struct cell *c, void *data) {
  c->density = NULL;
  c->nr_density = 0;

  c->force = NULL;
  c->nr_force = 0;
}
