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
#include "cell.h"

/* Local headers. */
#include "timers.h"

/**
 * @brief Lock a cell for access to its array of #part and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_locktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->hydro.hold || lock_trylock(&c->hydro.lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->hydro.hold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->hydro.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->hydro.lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->hydro.hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->hydro.lock) != 0) error("Failed to unlock cell.");
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
      atomic_dec(&finger2->hydro.hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->hydro.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its array of #gpart and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_glocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->grav.phold || lock_trylock(&c->grav.plock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->grav.phold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->grav.plock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->grav.plock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->grav.phold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->grav.plock) != 0) error("Failed to unlock cell.");
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
      atomic_dec(&finger2->grav.phold);

    /* Unlock this cell. */
    if (lock_unlock(&c->grav.plock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its #multipole and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_mlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->grav.mhold || lock_trylock(&c->grav.mlock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->grav.mhold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->grav.mlock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->grav.mlock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->grav.mhold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->grav.mlock) != 0) error("Failed to unlock cell.");
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
      atomic_dec(&finger2->grav.mhold);

    /* Unlock this cell. */
    if (lock_unlock(&c->grav.mlock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its array of #spart and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_slocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->stars.hold || lock_trylock(&c->stars.lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->stars.hold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->stars.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->stars.lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->stars.hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->stars.lock) != 0) error("Failed to unlock cell.");
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
      atomic_dec(&finger2->stars.hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->stars.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its array of #sink and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_sink_locktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->sinks.hold || lock_trylock(&c->sinks.lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->sinks.hold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->sinks.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->sinks.lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->sinks.hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->sinks.lock) != 0) error("Failed to unlock cell.");
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
      atomic_dec(&finger2->sinks.hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->sinks.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Lock a cell for access to its array of #bpart and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_blocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to lock this cell. */
  if (c->black_holes.hold || lock_trylock(&c->black_holes.lock) != 0) {
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Did somebody hold this cell in the meantime? */
  if (c->black_holes.hold) {
    /* Unlock this cell. */
    if (lock_unlock(&c->black_holes.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }

  /* Climb up the tree and lock/hold/unlock. */
  struct cell *finger;
  for (finger = c->parent; finger != NULL; finger = finger->parent) {
    /* Lock this cell. */
    if (lock_trylock(&finger->black_holes.lock) != 0) break;

    /* Increment the hold. */
    atomic_inc(&finger->black_holes.hold);

    /* Unlock the cell. */
    if (lock_unlock(&finger->black_holes.lock) != 0)
      error("Failed to unlock cell.");
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
      atomic_dec(&finger2->black_holes.hold);

    /* Unlock this cell. */
    if (lock_unlock(&c->black_holes.lock) != 0) error("Failed to unlock cell.");

    /* Admit defeat. */
    TIMER_TOC(timer_locktree);
    return 1;
  }
}

/**
 * @brief Unlock a cell's parents for access to #part array.
 *
 * @param c The #cell.
 */
void cell_unlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->hydro.lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->hydro.hold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to #gpart array.
 *
 * @param c The #cell.
 */
void cell_gunlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->grav.plock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->grav.phold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to its #multipole.
 *
 * @param c The #cell.
 */
void cell_munlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->grav.mlock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->grav.mhold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to #spart array.
 *
 * @param c The #cell.
 */
void cell_sunlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->stars.lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->stars.hold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to #sink array.
 *
 * @param c The #cell.
 */
void cell_sink_unlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->sinks.lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->sinks.hold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Unlock a cell's parents for access to #bpart array.
 *
 * @param c The #cell.
 */
void cell_bunlocktree(struct cell *c) {
  TIMER_TIC;

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->black_holes.lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->black_holes.hold);

  TIMER_TOC(timer_locktree);
}
